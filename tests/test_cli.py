"""CLI smoke tests (offline)."""

import json

from typer.testing import CliRunner

from proteinlens.cli import app

runner = CliRunner()


def test_compare_text(offline_cache):
    res = runner.invoke(app, ["compare", "4HHB:A", "4HHB:B"])
    assert res.exit_code == 0, res.output
    assert "TM-score (normalized by A): 0.904" in res.output
    assert "likely the SAME fold" in res.output


def test_compare_json(offline_cache):
    res = runner.invoke(app, ["compare", "1UBQ", "8TIM", "--json"])
    assert res.exit_code == 0, res.output
    data = json.loads(res.output)
    assert data["likely_same_fold"] is False and data["length_a"] == 76


def test_bad_chain_exit_code(offline_cache):
    res = runner.invoke(app, ["compare", "1UBQ", "4HHB:Q"])
    assert res.exit_code == 1


def test_info(offline_cache):
    res = runner.invoke(app, ["info", "P69905"])
    assert res.exit_code == 0 and "AlphaFold" in res.output
