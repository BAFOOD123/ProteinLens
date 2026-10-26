"""Known-answer tests for the comparison engine."""

import numpy as np
import pytest

from proteinlens.align import SAME_FOLD_TM, compare, sequence_compare
from proteinlens.data import load_structure


@pytest.fixture
def S(offline_cache):
    names = ["1UBQ", "1WM3", "4HHB:A", "4HHB:B", "8TIM"]
    return {n: load_structure(n, with_metadata=False) for n in names}


def test_self_alignment_is_perfect(S):
    r = compare(S["1UBQ"], S["1UBQ"])
    assert r.tm_score_a == pytest.approx(1.0) and r.tm_score_b == pytest.approx(1.0)
    assert r.rmsd == pytest.approx(0.0, abs=1e-3)
    assert r.aligned_length == 76 and r.seq_identity == 1.0
    assert r.top_regions == []


def test_globins_score_high(S):
    r = compare(S["4HHB:A"], S["4HHB:B"])  # hemoglobin alpha vs beta
    assert r.tm_score_a > 0.8 and r.tm_score_b > 0.8
    assert r.rmsd < 2.0
    assert r.likely_same_fold


def test_unrelated_pair_scores_low(S):
    r = compare(S["1UBQ"], S["8TIM"])  # ubiquitin vs a TIM barrel
    assert r.tm_score_max < SAME_FOLD_TM
    assert not r.likely_same_fold


def test_structure_sees_what_sequence_misses(S):
    # Ubiquitin and SUMO share a fold but only ~14% of their residues.
    r = compare(S["1UBQ"], S["1WM3"])
    assert r.tm_score_max > 0.7
    assert r.seq_identity < 0.25
    assert sequence_compare(S["1UBQ"].sequence, S["1WM3"].sequence).identity < 0.25


def test_symmetric_under_swap(S):
    ab = compare(S["4HHB:A"], S["8TIM"])
    ba = compare(S["8TIM"], S["4HHB:A"])
    assert ab.tm_score_a == ba.tm_score_b and ab.tm_score_b == ba.tm_score_a
    assert ab.rmsd == ba.rmsd and ab.aligned_length == ba.aligned_length
    assert ab.pair_index_a == ba.pair_index_b
    # The rotations are inverses of each other.
    assert np.allclose(np.asarray(ab.rotation) @ np.asarray(ba.rotation), np.eye(3), atol=1e-6)


def test_superposition_reproduces_distances(S):
    a, b = S["4HHB:A"], S["4HHB:B"]
    r = compare(a, b)
    moved = a.ca_coords[r.pair_index_a] @ np.asarray(r.rotation).T + np.asarray(r.translation)
    d = np.linalg.norm(moved - b.ca_coords[r.pair_index_b], axis=1)
    assert np.allclose(d, r.pair_distance, atol=1e-4)
    # TM-align's RMSD is from a separate least-squares fit, so allow a small difference.
    assert np.sqrt(np.mean(d[d < 5] ** 2)) == pytest.approx(r.rmsd, abs=0.05)
    per_res = r.per_residue_distance_a()
    assert per_res.shape == (len(a),) and np.isnan(per_res).sum() == len(a) - len(r.pair_index_a)


def test_rigid_motion_does_not_change_scores(S):
    a = S["1UBQ"]
    moved = load_structure("1UBQ", with_metadata=False)
    theta = 0.7
    rot = np.array(
        [[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]]
    )
    moved.ca_coords = a.ca_coords @ rot.T + np.array([10.0, -5.0, 3.0])
    moved.name = "moved"
    r = compare(a, moved)
    assert r.tm_score_a == pytest.approx(1.0) and r.rmsd == pytest.approx(0.0, abs=1e-3)


def test_deviating_region_is_found(S):
    a = S["1UBQ"]
    bent = load_structure("1UBQ", with_metadata=False)
    bent.name = "bent"
    idx = [i for i, rid in enumerate(bent.residue_ids) if 20 <= int(rid) <= 25]
    bent.ca_coords = bent.ca_coords.copy()
    bent.ca_coords[idx] += np.array([6.0, 0.0, 0.0])  # displace one loop by 6 Å
    r = compare(a, bent)
    assert r.top_regions, "displaced loop should be reported"
    top = r.top_regions[0]
    start, end = (int(x) for x in top.residues_a.split("-"))
    # TM-align may leave the loop ends unaligned, so the region can be a bit shorter.
    assert 20 <= start and end <= 25 and top.length >= 3
    assert top.mean_distance > 3.0


def test_sequence_baseline():
    assert sequence_compare("MKTAYIAKQR", "MKTAYIAKQR").identity == 1.0
    hb = sequence_compare(
        "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF",
        "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFG",
    )
    assert 0.3 < hb.identity < 0.6 and hb.score > 0
    assert np.isnan(sequence_compare("MKT", "MKT", with_identity=False).identity)
