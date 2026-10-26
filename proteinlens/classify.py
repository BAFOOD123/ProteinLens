"""Fold classification by structural nearest neighbour against the CATH reference set.

The reference set is the benchmark set built by ``benchmarks/build_dataset.py``
(about 300 CATH domains from 20 folds). Its C-alpha coordinates and sequences
are stored in ``benchmarks/data/reference_set.npz`` so classification works
offline.

Method (1-nearest-neighbour, "1-NN"): align the query against every reference
domain with TM-align and predict the fold and superfamily of the reference
domain with the highest TM-score, normalized by the *query's* length (so the
score answers "how much of my protein is covered by this reference?").
"""

from __future__ import annotations

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from proteinlens.align import compare
from proteinlens.data import Structure

DATA_DIR = Path(__file__).resolve().parent.parent / "benchmarks" / "data"


def light_structure(name: str, sequence: str, ca: np.ndarray) -> Structure:
    """A Structure with only C-alpha coordinates (all the alignment needs)."""
    return Structure(
        name=name,
        source="cath",
        chain="A",
        sequence=sequence,
        ca_coords=np.asarray(ca, dtype=float),
        residue_ids=[str(i + 1) for i in range(len(sequence))],
        bfactors=np.zeros(len(sequence)),
        atom_names=[],
        atom_elements=[],
        atom_coords=np.zeros((0, 3)),
        atom_residue_index=np.zeros(0, dtype=int),
        residue_names=[],
    )


@dataclass
class ReferenceSet:
    table: pd.DataFrame  # one row per domain (see benchmarks/data/dataset.csv)
    structures: list[Structure]

    def __len__(self) -> int:
        return len(self.structures)


@lru_cache(maxsize=2)
def load_reference(data_dir: str | None = None) -> ReferenceSet:
    folder = Path(data_dir) if data_dir else DATA_DIR
    table = pd.read_csv(folder / "dataset.csv")
    npz = np.load(folder / "reference_set.npz")
    off = npz["offsets"]
    structs = [
        light_structure(str(d), str(s), npz["ca"][off[i] : off[i + 1]])
        for i, (d, s) in enumerate(zip(npz["domain"], npz["sequence"], strict=True))
    ]
    if list(table["domain"]) != [s.name for s in structs]:
        raise ValueError("dataset.csv and reference_set.npz are out of sync; rebuild the dataset.")
    return ReferenceSet(table=table, structures=structs)


@dataclass
class Hit:
    domain: str
    fold: str
    fold_name: str
    superfamily: str
    superfamily_name: str
    tm_query: float  # TM-score normalized by the query length (used for ranking)
    tm_reference: float  # TM-score normalized by the reference length
    rmsd: float
    aligned_length: int
    seq_identity: float


def _score_one(args: tuple[Structure, Structure]) -> tuple[float, float, float, int, float]:
    q, ref = args
    r = compare(q, ref)
    return r.tm_score_a, r.tm_score_b, r.rmsd, r.aligned_length, r.seq_identity


def classify(
    query: Structure,
    reference: ReferenceSet | None = None,
    top_k: int = 5,
    exclude: set[str] | None = None,
    workers: int = 1,
) -> list[Hit]:
    """Return the ``top_k`` most similar reference domains, best first."""
    ref = reference or load_reference()
    keep = [i for i, s in enumerate(ref.structures) if not exclude or s.name not in exclude]
    jobs = [(query, ref.structures[i]) for i in keep]
    if workers > 1:
        with ProcessPoolExecutor(workers) as pool:
            scores = list(pool.map(_score_one, jobs, chunksize=16))
    else:
        scores = [_score_one(j) for j in jobs]
    hits = []
    for i, (tq, tr, rmsd, n, sid) in zip(keep, scores, strict=True):
        row = ref.table.iloc[i]
        hits.append(
            Hit(
                domain=row["domain"],
                fold=row["fold"],
                fold_name=row["fold_name"] if isinstance(row["fold_name"], str) else "",
                superfamily=row["superfamily"],
                superfamily_name=(
                    row["superfamily_name"] if isinstance(row["superfamily_name"], str) else ""
                ),
                tm_query=round(tq, 4),
                tm_reference=round(tr, 4),
                rmsd=round(rmsd, 3),
                aligned_length=n,
                seq_identity=round(sid, 4),
            )
        )
    hits.sort(key=lambda h: h.tm_query, reverse=True)
    return hits[:top_k]


# --------------------------------------------------------------------------------------
# Nearest-neighbour prediction from a precomputed score matrix (used by the benchmark)
# --------------------------------------------------------------------------------------


def nearest_neighbour(scores: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Leave-one-out 1-NN: for each row i, the label of the column j != i with the top score.

    ``scores[i, j]`` = similarity of query i to reference j (higher = more similar).
    Ties are broken by the lowest index, which is deterministic.
    """
    s = scores.astype(float).copy()
    np.fill_diagonal(s, -np.inf)  # leave one out: a protein may not vote for itself
    return labels[np.argmax(s, axis=1)]


def majority_baseline(labels: np.ndarray) -> np.ndarray:
    """Leave-one-out majority class: predict the most common label among the *other* items.

    Ties are broken alphabetically so the result is deterministic.
    """
    counts = Counter(labels.tolist())
    preds = []
    for lab in labels:
        c = counts.copy()
        c[lab] -= 1
        best = max(c.values())
        preds.append(min(k for k, v in c.items() if v == best))
    return np.array(preds)
