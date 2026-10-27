"""Tests for nearest-neighbour classification and the bundled CATH reference set."""

import numpy as np

from proteinlens.classify import classify, load_reference, majority_baseline, nearest_neighbour
from proteinlens.data import load_structure


def test_nearest_neighbour_leaves_self_out():
    scores = np.array([[1.0, 0.9, 0.1], [0.9, 1.0, 0.2], [0.1, 0.3, 1.0]])
    labels = np.array(["x", "x", "y"])
    # Item 2 is most similar to itself, but may not vote for itself.
    assert nearest_neighbour(scores, labels).tolist() == ["x", "x", "x"]


def test_majority_baseline_is_leave_one_out():
    labels = np.array(["a", "a", "b", "b", "b"])
    # For a "b" query the other items are a,a,b,b -> tie, broken alphabetically -> "a".
    assert majority_baseline(labels).tolist() == ["b", "b", "a", "a", "a"]


def test_reference_set_is_consistent():
    ref = load_reference()
    assert len(ref) == len(ref.table) >= 250
    assert ref.table["fold"].nunique() >= 15
    assert all(len(s.sequence) == len(s.ca_coords) for s in ref.structures)


def test_classify_ubiquitin_and_globin(offline_cache):
    hits = classify(load_structure("1UBQ", with_metadata=False), top_k=3)
    assert hits[0].fold == "3.10.20"  # ubiquitin-like (UB roll)
    assert hits[0].tm_query >= hits[1].tm_query >= hits[2].tm_query
    hb = classify(load_structure("4HHB:B", with_metadata=False), top_k=1)
    assert hb[0].fold == "1.10.490"  # globin-like


def test_classify_can_exclude_domains():
    ref = load_reference()
    q = ref.structures[0]
    hits = classify(q, ref, top_k=1, exclude={q.name})
    assert hits[0].domain != q.name


def test_coverage_score_ranks_by_smaller_tm(offline_cache):
    from proteinlens.classify import hit_score

    q = load_structure("1UBQ", with_metadata=False)
    hits = classify(q, top_k=5, score="coverage")
    scores = [hit_score(h, "coverage") for h in hits]
    assert scores == sorted(scores, reverse=True)
    assert all(hit_score(h, "coverage") == min(h.tm_query, h.tm_reference) for h in hits)
    assert hits[0].fold == "3.10.20"


def test_unknown_score_is_rejected():
    import pytest

    with pytest.raises(ValueError):
        classify(load_reference().structures[0], score="bogus")
