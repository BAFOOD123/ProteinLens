"""The headline numbers in README.md and docs/WRITEUP.md must match the saved result files."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SUMMARY = json.loads((ROOT / "benchmarks/results/summary.json").read_text())
TEST = json.loads((ROOT / "benchmarks/results_test/summary.json").read_text())
HARD = json.loads((ROOT / "benchmarks/results_hard/summary.json").read_text())
LLM = json.loads((ROOT / "benchmarks/results/llm_experiment/summary.json").read_text())


def pct(x: float) -> str:
    return f"{100 * x:.1f}"


def acc_string(e: dict) -> str:
    """e.g. '95.8 % [93.4, 97.9]' exactly as written in the docs."""
    return f"{pct(e['accuracy'])} % [{pct(e['ci95'][0])}, {pct(e['ci95'][1])}]"


def expected_strings() -> list[str]:
    out = []
    for block in ("all_queries", "remote_homologs"):
        for method in ("structure", "sequence"):
            for level in ("fold", "superfamily"):
                out.append(acc_string(SUMMARY[block][method][level]))
    for level in ("superfamily",):
        for subset in ("all", "remote"):
            d = SUMMARY["structure_minus_sequence"][level][subset]
            out.append(f"{pct(d['diff'])} points [{pct(d['ci95'][0])}, {pct(d['ci95'][1])}]")
    b = LLM["tm_rule_baseline"]
    out.append(f"{pct(b['accuracy'])} % [{pct(b['ci95'][0])}, {pct(b['ci95'][1])}]")
    out.append(f"{SUMMARY['n_domains']} CATH S40 domains")
    # Hard set
    for method in ("structure", "sequence"):
        for level in ("fold", "superfamily"):
            out.append(acc_string(HARD["all_queries"][method][level]))
        out.append(acc_string(HARD["queries_with_superfamily_partner"][method]))
    for method in ("structure", "sequence"):
        e = HARD["remote_homologs"][method]["superfamily"]
        out.append(f"{pct(e['accuracy'])} % [{pct(e['ci95'][0])}, {pct(e['ci95'][1])}]")
    ab = HARD["abstain_below_tm_0_5"]
    out.append(
        f"answers {pct(ab['coverage'])} % of hard-set queries with "
        f"{acc_string(ab['fold_accuracy_answered'])}"
    )
    ab_main = SUMMARY["abstain_below_tm_0_5"]
    out.append(f"{pct(ab_main['coverage'])} % of main-set queries")
    # Held-out test and the score change
    out.append(acc_string(TEST["all_queries"]["structure"]["fold"]))
    out.append(acc_string(TEST["coverage_aware_score"]["fold"]))
    for name in ("results", "results_hard", "results_test"):
        d = json.loads((ROOT / f"benchmarks/{name}/summary.json").read_text())
        diff = d["coverage_aware_score"]["minus_primary"]["fold"]
        sign = "+" if diff["diff"] > 0 else ""
        out.append(f"{sign}{pct(diff['diff'])} [{pct(diff['ci95'][0])}, {pct(diff['ci95'][1])}]")
    errs = HARD["structure_fold_errors_vs_neighbour_size"]["errors"]
    out.append(f"{pct(errs['fraction_neighbour_longer'])} % of the {errs['n']} fold errors")
    return out


def test_readme_numbers_match_results():
    readme = (ROOT / "README.md").read_text().replace("**", "")  # ignore bold markup
    missing = [s for s in expected_strings() if s not in readme]
    assert not missing, f"README is out of sync with results: {missing}"


def test_writeup_numbers_match_results():
    text = (ROOT / "docs/WRITEUP.md").read_text()
    for block in ("all_queries", "remote_homologs"):
        for method in ("structure", "sequence"):
            for level in ("fold", "superfamily"):
                assert acc_string(SUMMARY[block][method][level]) in text


def test_llm_experiment_is_pending_or_complete():
    assert LLM["status"] in {"PENDING", "COMPLETE"}
    if LLM["status"] == "PENDING":
        assert "PENDING" in (ROOT / "README.md").read_text()
