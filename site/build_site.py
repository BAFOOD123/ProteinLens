"""Build the ProteinLens results website from the committed result files.

    python site/build_site.py            # writes site/index.html (self-contained)

Every number on the page comes from this script: it reads benchmarks/results*/summary.json and
the other result files, recomputes the example comparisons with ProteinLens itself (from the
bundled offline fixtures), and runs the real faithfulness checker on a deliberately corrupted
report. Nothing is typed in by hand. The data are injected into site/template.html.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

SETS = {
    "main": ("results", "data", "Main set"),
    "hard": ("results_hard", "data_hard", "Hard set"),
    "test": ("results_test", "data_test", "Held-out test set"),
}
EXAMPLES = [
    ("Hemoglobin α vs β", "4HHB:A", "4HHB:B"),
    ("Ubiquitin vs SUMO", "1UBQ", "1WM3"),
    ("Ubiquitin vs TIM barrel", "1UBQ", "8TIM"),
    ("AlphaFold vs crystal (Hb α)", "P69905", "4HHB:A"),
]
BACKBONE = {"N", "CA", "C", "O"}  # enough for 3Dmol.js to draw a cartoon, at a fraction of the size


def r4(x: float) -> float:
    return round(float(x), 4)


# --------------------------------------------------------------------------------------
# Benchmark sets
# --------------------------------------------------------------------------------------


def tm_histograms(res: Path, data: Path, bins: int = 40) -> dict:
    """Density of pairwise TM-scores (normalized by the shorter protein) per relationship."""
    pw = np.load(res / "pairwise.npz")
    names = [str(x) for x in pw["domain"]]
    t = pd.read_csv(data / "dataset.csv").set_index("domain").loc[names]
    fold, sf = t["fold"].to_numpy(), t["superfamily"].to_numpy()
    tm = pw["tm"]
    iu = np.triu_indices(len(names), k=1)
    tm_max = np.maximum(tm, tm.T)[iu]
    same_sf = (sf[:, None] == sf[None, :])[iu]
    same_fold = (fold[:, None] == fold[None, :])[iu]
    groups = {
        "same superfamily": same_sf,
        "same fold, different superfamily": same_fold & ~same_sf,
        "different fold": ~same_fold,
    }
    edges = np.linspace(0, 1, bins + 1)
    out = {"edges": [r4(e) for e in edges]}
    for name, mask in groups.items():
        counts, _ = np.histogram(tm_max[mask], edges)
        out[name] = {"n": int(mask.sum()), "density": [r4(c) for c in counts / max(1, mask.sum())]}
    return out


def load_set(key: str) -> dict | None:
    res_name, data_name, label = SETS[key]
    res, data = ROOT / "benchmarks" / res_name, ROOT / "benchmarks" / data_name
    if not (res / "summary.json").exists():
        return None
    s = json.loads((res / "summary.json").read_text())
    out = {"key": key, "label": label, "summary": s, "tm_hist": tm_histograms(res, data)}
    fails = res / "failures.csv"
    if fails.exists():
        f = pd.read_csv(fails).head(5)
        out["top_failures"] = f[
            [
                "domain",
                "length",
                "fold",
                "struct_nn",
                "nn_length",
                "struct_pred_fold",
                "struct_nn_tm",
                "best_same_fold_tm",
            ]
        ].to_dict(orient="records")
    return out


# --------------------------------------------------------------------------------------
# Example comparisons (recomputed with ProteinLens from the offline fixtures)
# --------------------------------------------------------------------------------------


def backbone_pdb(s, rotation=None, translation=None) -> str:
    """PDB text with only backbone atoms, optionally moved by x' = R x + t."""
    from proteinlens.data import Structure

    keep = [i for i, n in enumerate(s.atom_names) if n in BACKBONE]
    light = Structure(
        name=s.name,
        source=s.source,
        chain=s.chain,
        sequence=s.sequence,
        ca_coords=s.ca_coords,
        residue_ids=s.residue_ids,
        bfactors=s.bfactors,
        atom_names=[s.atom_names[i] for i in keep],
        atom_elements=[s.atom_elements[i] for i in keep],
        atom_coords=s.atom_coords[keep],
        atom_residue_index=s.atom_residue_index[keep],
        residue_names=s.residue_names,
    )
    return light.to_pdb(rotation, translation)


def examples() -> list[dict]:
    from proteinlens.align import compare
    from proteinlens.data import load_structure
    from proteinlens.report import build_facts, template_report
    from proteinlens.verify import check_report

    out = []
    for label, ia, ib in EXAMPLES:
        a, b = load_structure(ia), load_structure(ib)
        r = compare(a, b)
        facts = build_facts(a, b, r)
        text = template_report(facts)
        check = check_report(text, facts)
        dist = r.per_residue_distance_a()
        out.append(
            {
                "label": label,
                "a": {
                    "name": a.name,
                    "length": len(a),
                    "description": a.metadata.get("description") or a.metadata.get("title", ""),
                },
                "b": {
                    "name": b.name,
                    "length": len(b),
                    "description": b.metadata.get("description") or b.metadata.get("title", ""),
                },
                "tm_a": r4(r.tm_score_a),
                "tm_b": r4(r.tm_score_b),
                "rmsd": round(r.rmsd, 2),
                "aligned": r.aligned_length,
                "identity_pct": round(100 * r.seq_identity, 1),
                "same_fold": r.likely_same_fold,
                "regions": [g.__dict__ for g in r.top_regions],
                "residue_ids": a.residue_ids,
                "distance": [None if np.isnan(d) else round(float(d), 2) for d in dist],
                "report": text,
                "verified": check.passed,
                "numbers_checked": check.numbers_checked,
                "pdb_a": backbone_pdb(a, np.asarray(r.rotation), np.asarray(r.translation)),
                "pdb_b": backbone_pdb(b),
            }
        )
    return out


def checker_demo() -> dict:
    """Run the real checker on a faithful report and on a corrupted copy of it."""
    from proteinlens.align import compare
    from proteinlens.data import load_structure
    from proteinlens.report import build_facts
    from proteinlens.verify import check_report

    a, b = load_structure("4HHB:A"), load_structure("4HHB:B")
    r = compare(a, b)
    facts = build_facts(a, b, r)
    al = facts["alignment"]
    good = (
        f"Hemoglobin subunit alpha (4HHB:A, {len(a)} residues) and hemoglobin subunit beta "
        f"(4HHB:B, {len(b)} residues) superimpose closely: {al['aligned_residue_pairs']} residue "
        f"pairs align with an RMSD of {al['rmsd_angstrom']:.2f} Å. The TM-score is "
        f"{al['tm_score_normalized_by_a']:.3f}, above 0.5, so they likely share the same fold."
    )
    bad = (
        good.replace(f"{al['rmsd_angstrom']:.2f} Å", "0.84 Å").replace(
            "likely share the same fold", "likely have different folds"
        )
        + " Both chains bind heme to transport oxygen, like 2DHB."
    )
    res_good, res_bad = check_report(good, facts), check_report(bad, facts)
    return {
        "good": good,
        "good_passed": res_good.passed,
        "good_numbers": res_good.numbers_checked,
        "bad": bad,
        "bad_passed": res_bad.passed,
        "issues": [{"kind": i.kind, "text": i.text, "detail": i.detail} for i in res_bad.issues],
    }


def encode_coords(ca: np.ndarray) -> str:
    """C-alpha coordinates at 0.1 Å, delta-encoded (first point absolute), comma-separated."""
    q = np.round(np.asarray(ca) * 10).astype(int)
    d = np.vstack([q[:1], np.diff(q, axis=0)])
    return ",".join(map(str, d.ravel()))


def _superpose(args: tuple[str, str, str]) -> tuple[list[float], list[float]]:
    """Rotation/translation moving domain a onto domain b (runs in a worker process)."""
    from proteinlens.align import compare
    from proteinlens.classify import load_reference

    data_dir, a, b = args
    by_name = {st.name: st for st in load_reference(data_dir).structures}
    r = compare(by_name[a], by_name[b])
    return (
        [round(float(x), 4) for x in np.asarray(r.rotation).ravel()],
        [round(float(x), 2) for x in r.translation],
    )


def explorer() -> dict:
    """Every benchmark domain, its top-5 neighbours, and the superposition onto its nearest one.

    Neighbours are ranked by the current default (coverage-aware) score, with ties broken by the
    lowest index exactly as in the benchmark. The original-score neighbour is kept for comparison.
    """
    from concurrent.futures import ProcessPoolExecutor

    from proteinlens.classify import load_reference

    domains: dict[str, int] = {}
    records: list[list] = []  # [id, fold, superfamily, sequence, coords]
    fold_names: dict[str, str] = {}
    sets_out = {}
    for key, (res_name, data_name, label) in SETS.items():
        res, data = ROOT / "benchmarks" / res_name, ROOT / "benchmarks" / data_name
        if not (res / "pairwise.npz").exists():
            continue
        pw = np.load(res / "pairwise.npz")
        names = [str(x) for x in pw["domain"]]
        ref = load_reference(str(data))
        by_name = {st.name: st for st in ref.structures}
        table = ref.table.set_index("domain")
        for name in names:
            if name not in domains:
                row = table.loc[name]
                domains[name] = len(records)
                st = by_name[name]
                records.append(
                    [
                        name,
                        row["fold"],
                        row["superfamily"],
                        st.sequence,
                        encode_coords(st.ca_coords),
                    ]
                )
                fold_names[row["fold"]] = (
                    row["fold_name"] if isinstance(row["fold_name"], str) else ""
                )
        tm = pw["tm"].astype(float)
        cov = np.minimum(tm, tm.T)
        np.fill_diagonal(cov, -1.0)
        tmq = tm.copy()
        np.fill_diagonal(tmq, -1.0)
        orders = [np.argsort(-cov[i], kind="stable")[:5] for i in range(len(names))]
        jobs = [(str(data), name, names[int(orders[i][0])]) for i, name in enumerate(names)]
        with ProcessPoolExecutor(4) as pool:
            transforms = list(pool.map(_superpose, jobs, chunksize=16))
        rows = []
        n_ok = 0
        for i, name in enumerate(names):
            order = orders[i]
            nn = int(order[0])
            nn_orig = int(np.argmax(tmq[i]))
            n_ok += table.loc[names[nn], "fold"] == table.loc[name, "fold"]
            rot, trans = transforms[i]  # moves this domain onto its closest match
            rows.append(
                [
                    domains[name],
                    [[domains[names[j]], round(float(cov[i, j]), 3)] for j in order],
                    domains[names[nn_orig]],
                    round(float(tmq[i, nn_orig]), 3),
                    rot,
                    trans,
                ]
            )
        # The explorer must agree with the reported benchmark numbers.
        summary = json.loads((res / "summary.json").read_text())
        expected = summary["coverage_aware_score"]["fold"]["accuracy"]
        assert round(n_ok / len(names), 4) == expected, (key, n_ok / len(names), expected)
        sets_out[key] = {"label": label, "rows": rows}
    return {"domains": records, "folds": fold_names, "sets": sets_out}


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "site"), help="output folder")
    out_dir = Path(ap.parse_args().out)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Use a temporary copy of the bundled fixtures so this works offline and never writes
    # into the repository's test data.
    tmp = Path(tempfile.mkdtemp())
    shutil.copytree(ROOT / "tests" / "fixtures" / "cache", tmp / "cache")
    os.environ["PROTEINLENS_CACHE"] = str(tmp / "cache")

    data = {
        "sets": [s for s in (load_set(k) for k in SETS) if s],
        "examples": examples(),
        "checker": checker_demo(),
        "explorer": explorer(),
        "llm": json.loads((ROOT / "benchmarks/results/llm_experiment/summary.json").read_text()),
    }
    sel = ROOT / "benchmarks" / "score_selection.json"
    if sel.exists():
        data["score_selection"] = json.loads(sel.read_text())
    template = (ROOT / "site" / "template.html").read_text()
    blob = json.dumps(data, separators=(",", ":"), ensure_ascii=False)
    blob = blob.replace("</", "<\\/")  # never close the <script> tag early
    page = template.replace("/*__DATA__*/null", blob)
    (out_dir / "fragment.html").write_text(page)  # body content (for hosted previews)
    m = re.search(r"<title>(.*?)</title>", page)
    title = m.group(1) if m else "ProteinLens"
    # The standalone page carries title/meta in <head>; drop them from the body copy.
    body = re.sub(r"<title>.*?</title>\s*", "", page, count=1)
    body = re.sub(r'<meta name="claude:viewport-fit-sides"[^>]*>\s*', "", body, count=1)
    full = (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
        f"<title>{title}</title>\n</head>\n<body>\n{body}\n</body>\n</html>\n"
    )
    (out_dir / "index.html").write_text(full)
    shutil.rmtree(tmp, ignore_errors=True)
    kb = len(full.encode()) / 1024
    print(f"wrote {out_dir / 'index.html'} ({kb:.0f} KB) with {len(data['sets'])} benchmark sets")


if __name__ == "__main__":
    main()
