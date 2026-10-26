"""The LLM experiment pipeline, run end to end with a mock client (no API calls)."""

import json

import pytest

from proteinlens import experiment
from proteinlens.classify import load_reference
from proteinlens.llm import MockClient


@pytest.fixture(scope="module")
def small_pairs():
    return experiment.build_pairs(load_reference(), per_category=2, seed=7)


def test_relation():
    assert experiment.relation("1.10.490.10", "1.10.490.10") == "same superfamily"
    assert experiment.relation("1.10.490.10", "1.10.490.20") == "same fold, different superfamily"
    assert experiment.relation("1.10.490.10", "3.20.20.70") == "different fold"


def test_pairs_are_balanced_and_labelled(small_pairs):
    assert small_pairs["category"].value_counts().to_dict() == {c: 2 for c in experiment.CATEGORIES}
    for _, r in small_pairs.iterrows():
        assert r.same_fold_truth == (r.fold_a == r.fold_b)
        assert experiment.relation(r.superfamily_a, r.superfamily_b) == r.category


def test_pipeline_with_mock(small_pairs, monkeypatch):
    monkeypatch.setattr(experiment, "entry_title", lambda d: "")
    answer = json.dumps({"same_fold": True, "explanation": "The proteins bind DNA."})
    client = MockClient([answer])
    res = experiment.run_condition_rows(small_pairs, client)
    assert len(res) == 2 * len(small_pairs) and set(res["condition"]) == {"A", "B"}
    # "Always same fold" is right exactly for the two same-fold categories.
    assert res["correct"].sum() == 2 * 4
    # The explanation contains outside knowledge, so it must never pass the checker.
    assert not res["faithful"].any()
    s = experiment.summarize(res)
    assert s["mock / B"]["accuracy"] == pytest.approx(4 / 6, abs=1e-3)
    # Grounded prompts contain the facts, ungrounded ones do not.
    b_prompts = [c["user"] for c in client.calls if "Computed facts" in c["user"]]
    assert len(b_prompts) == len(small_pairs)


def test_tm_rule_and_pending(small_pairs, tmp_path):
    rule = experiment.tm_rule_rows(small_pairs)
    assert set(rule["tm_rule_same_fold"]) <= {True, False}
    experiment.write_pending(tmp_path, rule, "make eval-llm")
    status = json.loads((tmp_path / "summary.json").read_text())
    assert status["status"] == "PENDING" and "make eval-llm" in status["how_to_run"]


def test_cost_estimate_is_positive():
    assert experiment.estimate_cost(["x" * 3500], 1.0, 5.0, out_tokens=100) > 0
