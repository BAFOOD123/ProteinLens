"""End-to-end test of the Streamlit app without a browser (streamlit.testing.AppTest)."""

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parent.parent / "app" / "streamlit_app.py")


@pytest.fixture
def app(offline_cache, monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    return AppTest.from_file(APP, default_timeout=120)


def test_compare_example_without_api_key(app):
    app.run()
    assert not app.exception
    next(b for b in app.button if b.label == "Hemoglobin α vs β").click().run()
    assert not app.exception
    assert any("Verified against computed data" in s.value for s in app.success)
    assert any("Likely the same fold" in i.value for i in app.info)
    assert [m.value for m in app.metric][:2] == ["0.904", "0.876"]


def test_classify_page(app):
    app.run()
    app.sidebar.radio[0].set_value("Classify").run()
    next(b for b in app.button if b.label == "Example: 1UBQ").click().run()
    assert not app.exception
    assert any("3.10.20" in s.value for s in app.success)


def test_benchmarks_page(app):
    app.run()
    app.sidebar.radio[0].set_value("Benchmarks").run()
    assert not app.exception
    assert any("PENDING" in w.value for w in app.warning)


GOOD_NARRATION = (
    "Hemoglobin subunit alpha (4HHB:A, 141 residues) and hemoglobin subunit beta (4HHB:B, "
    "146 residues) superimpose closely: 139 residue pairs align with an RMSD of 1.41 Å. "
    "The TM-score is 0.904 normalized by chain A and 0.876 normalized by chain B, above 0.5, "
    "so they likely share the same fold."
)


def _with_mock_llm(monkeypatch, responses):
    import json

    import proteinlens.llm as llm
    from proteinlens.llm import MockClient

    texts = [json.dumps({"summary": t, "same_fold": True}) for t in responses]
    monkeypatch.setattr(llm, "get_client", lambda *a, **k: MockClient(texts, model="mock-llm"))


def test_compare_with_llm_narration(app, monkeypatch):
    _with_mock_llm(monkeypatch, [GOOD_NARRATION])
    app.run()
    next(b for b in app.button if b.label == "Hemoglobin α vs β").click().run()
    assert not app.exception
    assert any("LLM narration: on (mock-llm)" in s.value for s in app.sidebar.success)
    assert any("Verified against computed data" in s.value for s in app.success)
    assert any("Source: LLM narration" in c.value for c in app.caption)


def test_compare_llm_failure_falls_back(app, monkeypatch):
    _with_mock_llm(monkeypatch, [GOOD_NARRATION + " Both bind oxygen."])
    app.run()
    next(b for b in app.button if b.label == "Hemoglobin α vs β").click().run()
    assert not app.exception
    assert any("Source: template" in c.value for c in app.caption)
    assert any("Verified against computed data" in s.value for s in app.success)
