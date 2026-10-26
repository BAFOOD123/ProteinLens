"""Comparison engine: structural alignment (TM-align) and a sequence-alignment baseline.

Key numbers (all computed on C-alpha atoms):

* **TM-score** (0-1). Measures how much of a protein can be superimposed on the
  other. It is normalized by the length of one protein, so there are two
  versions: ``tm_score_a`` (divided by the length of A) and ``tm_score_b``.
  Rule of thumb (Xu & Zhang, 2010): TM-score > 0.5 usually means the same fold;
  around 0.17 is what two random, unrelated proteins give.
* **RMSD** (Å). Root-mean-square distance between aligned C-alpha atoms after the
  best superposition, over the residue pairs TM-align counts as aligned
  (pairs closer than 5 Å). Small RMSD over few residues can be misleading,
  which is why TM-score is the main metric.
* **Aligned length**: number of residue pairs closer than 5 Å after superposition
  (TM-align's own definition).
* **Sequence identity**: fraction of those aligned pairs with the same amino acid.

The structural alignment itself is done by TM-align (Zhang & Skolnick, 2005)
through the ``tmtools`` Python bindings.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from tmtools import tm_align

from proteinlens.data import Structure

SAME_FOLD_TM = 0.5  # TM-score above this: very likely the same fold (Xu & Zhang 2010)
RANDOM_TM = 0.17  # typical TM-score of unrelated proteins
DEVIATION_CUTOFF = 3.0  # Å; aligned pairs further apart than this count as "deviating"


@dataclass
class Region:
    """A stretch of consecutive aligned residue pairs that deviate after superposition."""

    residues_a: str  # author residue numbers in A, e.g. "45-52"
    residues_b: str
    length: int  # number of aligned pairs in the region
    mean_distance: float  # Å
    max_distance: float  # Å


@dataclass
class AlignmentResult:
    """Everything ProteinLens computes when comparing two structures."""

    name_a: str
    name_b: str
    length_a: int
    length_b: int
    tm_score_a: float  # TM-score normalized by the length of A
    tm_score_b: float  # TM-score normalized by the length of B
    rmsd: float  # Å, over the aligned pairs (TM-align definition)
    aligned_length: int  # pairs closer than 5 Å
    seq_identity: float  # fraction (0-1) of aligned pairs with identical residues
    rotation: list[list[float]]  # 3x3 matrix U; x_on_b = U @ x_a + t
    translation: list[float]
    aligned_seq_a: str  # alignment strings ('-' = gap)
    aligned_seq_b: str
    match_line: str  # ':' pair < 5 Å, '.' other aligned pair, ' ' gap
    pair_index_a: list[int] = field(repr=False)  # aligned pairs (indices into each structure)
    pair_index_b: list[int] = field(repr=False)
    pair_distance: list[float] = field(repr=False)  # Å, after superposition
    top_regions: list[Region] = field(default_factory=list)

    @property
    def tm_score_max(self) -> float:
        return max(self.tm_score_a, self.tm_score_b)

    @property
    def tm_score_min(self) -> float:
        return min(self.tm_score_a, self.tm_score_b)

    @property
    def likely_same_fold(self) -> bool:
        """Simple rule: the larger of the two TM-scores is above 0.5."""
        return self.tm_score_max > SAME_FOLD_TM

    def per_residue_distance_a(self) -> np.ndarray:
        """Distance (Å) for every residue of A to its partner in B; NaN if unaligned."""
        d = np.full(self.length_a, np.nan)
        d[self.pair_index_a] = self.pair_distance
        return d

    def per_residue_distance_b(self) -> np.ndarray:
        d = np.full(self.length_b, np.nan)
        d[self.pair_index_b] = self.pair_distance
        return d

    def to_dict(self, include_pairs: bool = False) -> dict[str, Any]:
        out = asdict(self)
        if not include_pairs:
            for k in ("pair_index_a", "pair_index_b", "pair_distance"):
                out.pop(k)
        out["tm_score_max"] = self.tm_score_max
        out["likely_same_fold"] = self.likely_same_fold
        return out


def _swap(r: AlignmentResult) -> AlignmentResult:
    """The same alignment seen from the other side (B onto A)."""
    u = np.asarray(r.rotation)
    t = np.asarray(r.translation)
    u_inv = u.T  # inverse of a rotation matrix is its transpose
    t_inv = -u_inv @ t
    swapped = AlignmentResult(
        name_a=r.name_b,
        name_b=r.name_a,
        length_a=r.length_b,
        length_b=r.length_a,
        tm_score_a=r.tm_score_b,
        tm_score_b=r.tm_score_a,
        rmsd=r.rmsd,
        aligned_length=r.aligned_length,
        seq_identity=r.seq_identity,
        rotation=u_inv.tolist(),
        translation=t_inv.tolist(),
        aligned_seq_a=r.aligned_seq_b,
        aligned_seq_b=r.aligned_seq_a,
        match_line=r.match_line,
        pair_index_a=r.pair_index_b,
        pair_index_b=r.pair_index_a,
        pair_distance=r.pair_distance,
        top_regions=[
            Region(g.residues_b, g.residues_a, g.length, g.mean_distance, g.max_distance)
            for g in r.top_regions
        ],
    )
    return swapped


def _range_label(ids: list[str]) -> str:
    return ids[0] if len(ids) == 1 else f"{ids[0]}-{ids[-1]}"


def find_deviating_regions(
    a: Structure,
    b: Structure,
    idx_a: list[int],
    idx_b: list[int],
    dist: list[float],
    cutoff: float = DEVIATION_CUTOFF,
    min_length: int = 3,
    top: int = 3,
) -> list[Region]:
    """Find runs of >= ``min_length`` consecutive aligned pairs all further apart than ``cutoff``.

    "Consecutive" means neighbouring positions in the alignment with no gap in
    either protein. Regions are ranked by mean distance (largest first).
    """
    regions: list[Region] = []
    run: list[int] = []

    def close_run() -> None:
        if len(run) >= min_length:
            d = [dist[k] for k in run]
            regions.append(
                Region(
                    residues_a=_range_label([a.residue_ids[idx_a[k]] for k in run]),
                    residues_b=_range_label([b.residue_ids[idx_b[k]] for k in run]),
                    length=len(run),
                    mean_distance=round(float(np.mean(d)), 2),
                    max_distance=round(float(np.max(d)), 2),
                )
            )
        run.clear()

    for k in range(len(dist)):
        contiguous = k > 0 and idx_a[k] == idx_a[k - 1] + 1 and idx_b[k] == idx_b[k - 1] + 1
        if run and not contiguous:
            close_run()
        if dist[k] > cutoff:
            run.append(k)
        else:
            close_run()
    close_run()
    regions.sort(key=lambda g: g.mean_distance, reverse=True)
    return regions[:top]


def _align_raw(a: Structure, b: Structure) -> AlignmentResult:
    """Run TM-align with A as chain 1 and B as chain 2."""
    res = tm_align(a.ca_coords, b.ca_coords, a.sequence, b.sequence)
    u = np.asarray(res.u)
    t = np.asarray(res.t)
    # Walk the alignment strings to recover which residue is paired with which.
    idx_a: list[int] = []
    idx_b: list[int] = []
    close: list[bool] = []  # pair is closer than 5 Å (':' in TM-align's match line)
    ia = ib = 0
    for xa, yb, m in zip(res.seqxA, res.seqyA, res.seqM, strict=True):
        if xa != "-" and yb != "-":
            idx_a.append(ia)
            idx_b.append(ib)
            close.append(m == ":")
        ia += xa != "-"
        ib += yb != "-"
    moved = a.ca_coords[idx_a] @ u.T + t
    dist = np.linalg.norm(moved - b.ca_coords[idx_b], axis=1)
    n_close = int(sum(close))
    identical = sum(
        1
        for i, j, c in zip(idx_a, idx_b, close, strict=True)
        if c and a.sequence[i] == b.sequence[j]
    )
    return AlignmentResult(
        name_a=a.name,
        name_b=b.name,
        length_a=len(a),
        length_b=len(b),
        tm_score_a=float(res.tm_norm_chain1),
        tm_score_b=float(res.tm_norm_chain2),
        rmsd=float(res.rmsd),
        aligned_length=n_close,
        seq_identity=identical / n_close if n_close else 0.0,
        rotation=u.tolist(),
        translation=t.tolist(),
        aligned_seq_a=res.seqxA,
        aligned_seq_b=res.seqyA,
        match_line=res.seqM,
        pair_index_a=idx_a,
        pair_index_b=idx_b,
        pair_distance=[float(x) for x in dist],
        top_regions=find_deviating_regions(a, b, idx_a, idx_b, [float(x) for x in dist]),
    )


def compare(a: Structure, b: Structure) -> AlignmentResult:
    """Structurally align A and B.

    TM-align's search is very slightly order-dependent (A vs B and B vs A can
    differ in the third decimal). To make ``compare(a, b)`` and ``compare(b, a)``
    give exactly the same numbers, we always run TM-align in one canonical
    order and then flip the result if needed.
    """
    if len(a) < 5 or len(b) < 5:
        raise ValueError("Both structures need at least 5 residues with C-alpha atoms.")
    key_a = (len(a), a.sequence, a.name)
    key_b = (len(b), b.sequence, b.name)
    if key_a <= key_b:
        return _align_raw(a, b)
    return _swap(_align_raw(b, a))


# --------------------------------------------------------------------------------------
# Sequence-only baseline
# --------------------------------------------------------------------------------------


def make_sequence_aligner() -> PairwiseAligner:
    """Smith-Waterman local alignment, BLOSUM62, gap open 11 / extend 1 (BLAST defaults)."""
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


_ALIGNER: PairwiseAligner | None = None


@dataclass
class SequenceResult:
    score: float  # raw local alignment score (BLOSUM62 units)
    identity: float  # identical aligned residues / length of the shorter sequence
    aligned_length: int  # aligned (non-gap) columns in the local alignment


def _clean_seq(s: str) -> str:
    # BLOSUM62 has no entry for unusual letters; map them to X (unknown).
    return "".join(c if c in "ARNDCQEGHILKMFPSTWYVBZX" else "X" for c in s.upper())


def sequence_compare(seq_a: str, seq_b: str, with_identity: bool = True) -> SequenceResult:
    """Sequence-only comparison used as the baseline the structure method must beat."""
    global _ALIGNER
    if _ALIGNER is None:
        _ALIGNER = make_sequence_aligner()
    sa, sb = _clean_seq(seq_a), _clean_seq(seq_b)
    if not with_identity:
        return SequenceResult(float(_ALIGNER.score(sa, sb)), float("nan"), 0)
    aln = _ALIGNER.align(sa, sb)[0]
    counts = aln.counts()
    aligned = counts.identities + counts.mismatches
    return SequenceResult(
        score=float(aln.score),
        identity=counts.identities / min(len(sa), len(sb)),
        aligned_length=int(aligned),
    )
