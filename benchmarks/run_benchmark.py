"""All-vs-all comparison of the benchmark set (structure with TM-align + sequence baseline).

Resumable: work is split into chunks of query rows; each finished chunk is saved
under data/cache/benchmark_chunks/ and skipped on the next run. If the time limit
(--max-minutes) is reached, the script exits with code 3 and you simply run it
again. When all chunks exist the matrices are assembled into
benchmarks/results/pairwise.npz (or benchmarks/smoke/ with --smoke).

Matrices (n x n, row = query i, column = reference j):
  tm        TM-score normalized by the length of the *query* (row) protein
  rmsd      TM-align RMSD (Å) over aligned pairs
  aligned   TM-align aligned length (pairs < 5 Å)
  sid_struct sequence identity over the structurally aligned pairs
  seq_score Smith-Waterman score (BLOSUM62, gaps 11/1)
  seq_id    identical residues in the local sequence alignment / shorter length
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from proteinlens.align import compare, sequence_compare  # noqa: E402
from proteinlens.classify import load_reference  # noqa: E402
from proteinlens.data import cache_dir  # noqa: E402

FIELDS = ["tm_ij", "tm_ji", "rmsd", "aligned", "sid_struct", "seq_score", "seq_id"]


def smoke_indices(table) -> list[int]:
    """First 2 domains of every fold (40 domains): a fast end-to-end check."""
    return [int(i) for i in table.groupby("fold", sort=False).head(2).index]


DATA_DIRS = {"main": "data", "hard": "data_hard"}
RESULT_DIRS = {"main": "results", "hard": "results_hard"}


def row_job(args: tuple[int, list[int], str]) -> tuple[int, np.ndarray]:
    """Compare query i against every j in ``cols`` (all j > i). Runs in a worker process."""
    i, cols, data_dir = args
    ref = load_reference(data_dir)
    a = ref.structures[i]
    out = np.zeros((len(cols), len(FIELDS)))
    for k, j in enumerate(cols):
        b = ref.structures[j]
        r = compare(a, b)
        s = sequence_compare(a.sequence, b.sequence)
        out[k] = [
            r.tm_score_a,
            r.tm_score_b,
            r.rmsd,
            r.aligned_length,
            r.seq_identity,
            s.score,
            s.identity,
        ]
    return i, out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=sorted(DATA_DIRS), default="main", help="benchmark set")
    ap.add_argument("--smoke", action="store_true", help="40-domain subset")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max-minutes", type=float, default=25.0)
    ap.add_argument("--chunk-rows", type=int, default=20)
    args = ap.parse_args()

    data_dir = str(REPO / "benchmarks" / DATA_DIRS[args.set])
    ref = load_reference(data_dir)
    idx = smoke_indices(ref.table) if args.smoke else list(range(len(ref)))
    names = [ref.structures[i].name for i in idx]
    n = len(idx)
    tag = hashlib.sha1(",".join(names).encode()).hexdigest()[:10]
    chunk_dir = cache_dir() / "benchmark_chunks" / f"{'smoke' if args.smoke else 'full'}_{tag}"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    out_dir = REPO / "benchmarks" / ("smoke" if args.smoke else RESULT_DIRS[args.set])
    out_dir.mkdir(parents=True, exist_ok=True)

    chunks = [list(range(s, min(s + args.chunk_rows, n))) for s in range(0, n, args.chunk_rows)]
    n_pairs = n * (n - 1) // 2
    print(f"{n} domains, {n_pairs} pairs, {len(chunks)} chunks -> {chunk_dir}", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(args.workers) as pool:
        for c, rows in enumerate(chunks):
            path = chunk_dir / f"chunk_{c:03d}.npz"
            if path.exists():
                continue
            if (time.time() - t0) / 60 > args.max_minutes:
                print("Time limit reached; run again to resume.", flush=True)
                return 3
            jobs = [(idx[r], [idx[j] for j in range(r + 1, n)], data_dir) for r in rows]
            results = dict(pool.map(row_job, jobs))
            data = {f"row_{r}": results[idx[r]] for r in rows}
            tmp = path.with_suffix(".tmp.npz")
            np.savez(tmp, **data)
            tmp.replace(path)  # atomic, so an interrupted run never leaves a broken chunk
            print(f"  chunk {c + 1}/{len(chunks)} done ({time.time() - t0:.0f}s)", flush=True)

    # Assemble the full matrices.
    m = {
        f: np.zeros((n, n)) for f in ("tm", "rmsd", "aligned", "sid_struct", "seq_score", "seq_id")
    }
    for c, rows in enumerate(chunks):
        data = np.load(chunk_dir / f"chunk_{c:03d}.npz")
        for r in rows:
            block = data[f"row_{r}"]
            cols = np.arange(r + 1, n)
            m["tm"][r, cols], m["tm"][cols, r] = block[:, 0], block[:, 1]
            for f, k in (
                ("rmsd", 2),
                ("aligned", 3),
                ("sid_struct", 4),
                ("seq_score", 5),
                ("seq_id", 6),
            ):
                m[f][r, cols] = m[f][cols, r] = block[:, k]
    np.fill_diagonal(m["tm"], 1.0)
    np.fill_diagonal(m["sid_struct"], 1.0)
    np.fill_diagonal(m["seq_id"], 1.0)
    np.savez_compressed(
        out_dir / "pairwise.npz",
        domain=np.array(names),
        **{k: v.astype(np.float32) for k, v in m.items()},
    )
    print(f"Saved {out_dir / 'pairwise.npz'} in {time.time() - t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
