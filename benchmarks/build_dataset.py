"""Build the ProteinLens benchmark set from CATH (one command, cached, resumable).

Procedure (all parameters are in CONFIG below and saved to benchmarks/data/selection.json):

1. Download CATH release files: the domain classification list, the names file and
   the S40 non-redundant list (no two domains share >= 40 % sequence identity).
2. Keep S40 domains that are X-ray structures at <= 2.5 Å resolution, 50-300
   residues long, and in CATH classes 1 (mainly alpha), 2 (mainly beta) or 3 (alpha/beta).
3. Fold = CATH topology (C.A.T). A fold is *eligible* if it has >= 10 such domains
   and >= 2 superfamilies (C.A.T.H) with >= 2 domains each, so that both the fold
   level and the superfamily level can be tested.
4. Pick 20 folds: 5 well-known "anchor" folds chosen in advance for interpretability
   (globin-like, ubiquitin-like roll, TIM barrel, immunoglobulin-like sandwich,
   Rossmann fold) plus 15 drawn at random (seed 42) from the remaining eligible
   folds, 5 from each class.
5. In each fold pick up to 4 superfamilies at random among those with >= 2 domains,
   then take domains round-robin across them (random order) until 15 per fold.
6. Download each domain's coordinates from CATH and parse them; a domain that fails
   to download/parse is replaced by the next candidate in the same order.

Outputs (committed): benchmarks/data/dataset.csv, benchmarks/data/reference_set.npz,
benchmarks/data/selection.json. Coordinates are cached in data/cache/cath/ (gitignored).
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from proteinlens.data import (  # noqa: E402
    ProteinLensError,
    _download,
    cache_dir,
    fetch_cath_domain,
    parse_structure,
)

CATH_RELEASE = "https://download.cathdb.info/cath/releases/latest-release"
MAIN_CONFIG = {
    "name": "main",
    "cath_release_url": CATH_RELEASE,
    "cath_domain_version": "v4_4_0",
    "seed": 42,
    "max_resolution": 2.5,
    "min_length": 50,
    "max_length": 300,
    "classes": [1, 2, 3],
    "min_domains_per_fold": 10,
    "min_superfamilies_with_2": 2,
    "anchor_folds": ["1.10.490", "3.10.20", "3.20.20", "2.60.40", "3.40.50"],
    "random_folds_per_class": 5,
    "superfamilies_per_fold": 4,
    "domains_per_fold": 15,
}
# A harder, pre-specified second benchmark: many folds, few examples each. Same quality filters.
HARD_CONFIG = {
    **MAIN_CONFIG,
    "name": "hard",
    "seed": 2026,
    "min_domains_per_fold": 5,
    "min_superfamilies_with_2": 0,
    "anchor_folds": [],
    "n_random_folds": 120,  # drawn uniformly from all eligible folds (any class 1-3)
    "superfamilies_per_fold": None,  # None = domains drawn uniformly at random within the fold
    "domains_per_fold": 5,
}
# Held-out TEST set for the score selected on main + hard (the "development" sets). Same design as
# the hard set, new seed, and no domain that appears in the main or hard sets.
TEST_CONFIG = {
    **HARD_CONFIG,
    "name": "test",
    "seed": 7777,
    "exclude_domains_in": ["data", "data_hard"],
}
CONFIGS = {"main": MAIN_CONFIG, "hard": HARD_CONFIG, "test": TEST_CONFIG}
CONFIG = MAIN_CONFIG  # replaced in main() according to --set
DATA_DIRS = {"main": "data", "hard": "data_hard", "test": "data_test"}


def download_cath_files() -> dict[str, Path]:
    folder = cache_dir() / "cath_meta"
    files = {
        "domains": "cath-classification-data/cath-domain-list.txt",
        "names": "cath-classification-data/cath-names.txt",
        "s40": "non-redundant-data-sets/cath-dataset-nonredundant-S40.list",
    }
    out = {}
    for key, rel in files.items():
        print(f"  {rel}", flush=True)
        out[key] = _download(f"{CATH_RELEASE}/{rel}", folder / Path(rel).name, timeout=300)
    return out


def read_names(path: Path) -> dict[str, str]:
    names = {}
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        code = line.split()[0]
        names[code] = line.split(":", 1)[1].strip() if ":" in line else ""
    return names


def eligible_domains(files: dict[str, Path]) -> pd.DataFrame:
    cols = ["domain", "c", "a", "t", "h", "s35", "s60", "s95", "s100", "s100n", "length", "res"]
    d = pd.read_csv(files["domains"], sep=r"\s+", comment="#", names=cols)
    s40 = set(files["s40"].read_text().split())
    d = d[d["domain"].isin(s40)]
    d = d[(d["res"] <= CONFIG["max_resolution"])]  # NMR = 999, obsolete = 1000
    d = d[d["length"].between(CONFIG["min_length"], CONFIG["max_length"])]
    d = d[d["c"].isin(CONFIG["classes"])].copy()
    d["fold"] = d["c"].astype(str) + "." + d["a"].astype(str) + "." + d["t"].astype(str)
    d["superfamily"] = d["fold"] + "." + d["h"].astype(str)
    for other in CONFIG.get("exclude_domains_in", []):
        used = set(pd.read_csv(REPO / "benchmarks" / other / "dataset.csv")["domain"])
        d = d[~d["domain"].isin(used)]
    return d.sort_values("domain").reset_index(drop=True)


def choose_folds(d: pd.DataFrame, rng: np.random.Generator) -> list[str]:
    sf_sizes = d.groupby(["fold", "superfamily"]).size()
    n_sf2 = (sf_sizes >= 2).groupby(level=0).sum()
    n_dom = d.groupby("fold").size()
    eligible = sorted(
        f
        for f in n_dom.index
        if n_dom[f] >= CONFIG["min_domains_per_fold"]
        and n_sf2.get(f, 0) >= CONFIG["min_superfamilies_with_2"]
    )
    for f in CONFIG["anchor_folds"]:
        assert f in eligible, f"anchor fold {f} is not eligible"
    if CONFIG.get("n_random_folds"):
        return sorted(rng.choice(eligible, CONFIG["n_random_folds"], replace=False))
    folds = list(CONFIG["anchor_folds"])
    for cls in CONFIG["classes"]:
        pool = [f for f in eligible if f.startswith(f"{cls}.") and f not in folds]
        folds += sorted(rng.choice(pool, CONFIG["random_folds_per_class"], replace=False))
    return folds


def candidate_order(d: pd.DataFrame, fold: str, rng: np.random.Generator) -> list[str]:
    """Round-robin order of candidate domains across up to N random superfamilies."""
    sub = d[d["fold"] == fold]
    if CONFIG["superfamilies_per_fold"] is None:
        return [str(x) for x in rng.permutation(sub["domain"].to_numpy())]
    sizes = sub.groupby("superfamily").size()
    sfs = sorted(sizes[sizes >= 2].index)
    k = min(CONFIG["superfamilies_per_fold"], len(sfs))
    chosen = sorted(rng.choice(sfs, k, replace=False))
    queues = [
        list(rng.permutation(sub[sub["superfamily"] == sf]["domain"].to_numpy())) for sf in chosen
    ]
    order = []
    while any(queues):
        for q in queues:
            if q:
                order.append(str(q.pop(0)))
    return order


def fetch_fold(fold: str, candidates: list[str]) -> list[dict]:
    """Download/parse candidates in order until enough domains are accepted."""
    accepted = []
    for dom in candidates:
        if len(accepted) == CONFIG["domains_per_fold"]:
            break
        try:
            s = parse_structure(
                fetch_cath_domain(dom, CONFIG["cath_domain_version"]), name=dom, origin="cath"
            )
        except (ProteinLensError, ValueError) as err:
            print(f"    skip {dom}: {err}", flush=True)
            continue
        if len(s) < 40:
            print(f"    skip {dom}: only {len(s)} residues with C-alpha", flush=True)
            continue
        accepted.append({"domain": dom, "sequence": s.sequence, "ca": s.ca_coords})
    print(f"  fold {fold}: {len(accepted)} domains", flush=True)
    return accepted


def main() -> None:
    global CONFIG
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=sorted(CONFIGS), default="main", help="which benchmark set")
    CONFIG = CONFIGS[ap.parse_args().set]
    OUT = REPO / "benchmarks" / DATA_DIRS[CONFIG["name"]]
    OUT.mkdir(parents=True, exist_ok=True)
    print("Downloading CATH classification files (cached after the first run)...", flush=True)
    files = download_cath_files()
    names = read_names(files["names"])
    d = eligible_domains(files)
    print(f"Eligible S40 domains after filters: {len(d)}", flush=True)

    rng = np.random.default_rng(CONFIG["seed"])
    folds = choose_folds(d, rng)
    orders = {f: candidate_order(d, f, rng) for f in folds}  # all randomness happens here

    print("Downloading domain coordinates from CATH (4 parallel workers)...", flush=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = dict(zip(folds, pool.map(lambda f: fetch_fold(f, orders[f]), folds), strict=True))

    rows, seqs, coords = [], [], []
    info = d.set_index("domain")
    for fold in folds:
        for item in results[fold]:
            r = info.loc[item["domain"]]
            rows.append(
                {
                    "domain": item["domain"],
                    "cath_class": int(r["c"]),
                    "fold": fold,
                    "superfamily": r["superfamily"],
                    "s35_cluster": f"{r['superfamily']}.{int(r['s35'])}",
                    "fold_name": names.get(fold, ""),
                    "superfamily_name": names.get(r["superfamily"], ""),
                    "cath_length": int(r["length"]),
                    "n_residues": len(item["sequence"]),
                    "resolution": float(r["res"]),
                    "sequence": item["sequence"],
                }
            )
            seqs.append(item["sequence"])
            coords.append(item["ca"])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "dataset.csv", index=False)
    offsets = np.cumsum([0] + [len(s) for s in seqs])
    np.savez_compressed(
        OUT / "reference_set.npz",
        domain=np.array(df["domain"], dtype=str),
        sequence=np.array(seqs, dtype=str),
        offsets=offsets,
        ca=np.concatenate(coords).astype(np.float32),
    )
    (OUT / "selection.json").write_text(
        json.dumps({**{k: v for k, v in CONFIG.items() if k != "name"}, "folds": folds}, indent=2)
    )
    print(
        f"Wrote {len(df)} domains, {df['fold'].nunique()} folds, "
        f"{df['superfamily'].nunique()} superfamilies to {OUT}",
        flush=True,
    )


if __name__ == "__main__":
    main()
