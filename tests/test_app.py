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
