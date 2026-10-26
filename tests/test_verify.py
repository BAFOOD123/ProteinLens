"""The faithfulness checker must pass faithful reports and catch corrupted ones."""

import json

import pytest

from proteinlens.align import compare
from proteinlens.data import load_structure
from proteinlens.llm import MockClient
from proteinlens.report import build_facts, generate_report, template_report
from proteinlens.verify import check_report, fold_verdict_claims


@pytest.fixture
def globin_facts(offline_cache):
    a, b = load_structure("4HHB:A"), load_structure("4HHB:B")
    return build_facts(a, b, compare(a, b))


@pytest.fixture
def unrelated_facts(offline_cache):
    a, b = load_structure("1UBQ"), load_structure("8TIM")
    return build_facts(a, b, compare(a, b))


GOOD = (
    "Hemoglobin subunit alpha (4HHB:A, 141 residues) and hemoglobin subunit beta (4HHB:B, "
    "146 residues) superimpose closely: 139 residue pairs align with an RMSD of 1.41 Å. "
    "The TM-score is 0.904 normalized by chain A and 0.876 normalized by chain B, well "
    "above 0.5, so they likely share the same fold. Only 43.2% of aligned residues are "
    "identical."
)


def test_faithful_report_passes(globin_facts):
    res = check_report(GOOD, globin_facts, claimed_same_fold=True)
    assert res.passed, res.summary()
    assert res.numbers_checked >= 8


def test_rounding_is_allowed(globin_facts):
    text = "The TM-score is about 0.90 (or 90%), the RMSD is 1.4 Å and identity is 43%."
    assert check_report(text, globin_facts).passed


def test_template_reports_always_pass(globin_facts, unrelated_facts):
    for facts in (globin_facts, unrelated_facts):
        res = check_report(template_report(facts), facts)
        assert res.passed, res.summary()


@pytest.mark.parametrize(
    "corrupt, kind",
    [
        (GOOD.replace("0.904", "0.954"), "number"),  # wrong TM-score
        (GOOD.replace("1.41 Å", "0.91 Å"), "number"),  # wrong RMSD
        (GOOD.replace("139 residue pairs", "131 residue pairs"), "number"),  # wrong count
        (GOOD.replace("43.2%", "53.2%"), "number"),  # wrong identity
        (GOOD + " The TM-score difference is 0.028.", "number"),  # invented derived number
        (GOOD.replace("likely share the same fold", "likely have different folds"), "verdict"),
        (GOOD.replace("so they likely share", "so they do not share"), "verdict"),
        (GOOD + " Both resemble 2DHB.", "identifier"),  # PDB ID from memory
        (GOOD + " Both belong to CATH 1.10.490.10.", "identifier"),  # CATH code not in facts
        (GOOD + " Both are TIM barrels.", "fold_name"),  # fold name from memory
        (GOOD + " Both chains bind heme to transport oxygen.", "outside_knowledge"),
        (GOOD + " They diverged from a common ancestor.", "outside_knowledge"),
    ],
)
def test_corrupted_reports_are_caught(globin_facts, corrupt, kind):
    res = check_report(corrupt, globin_facts)
    assert not res.passed
    assert kind in {i.kind for i in res.issues}, res.summary()


def test_wrong_structured_verdict_is_caught(unrelated_facts):
    res = check_report("They were compared.", unrelated_facts, claimed_same_fold=True)
    assert not res.passed and res.issues[0].kind == "verdict"


def test_negation_is_local():
    claims = fold_verdict_claims(
        "Because the larger TM-score is not above 0.5, they likely have different folds. "
        "They do not share the same fold. They share a common fold."
    )
    assert [c[1] for c in claims] == [False, False, True]


def _llm(summary: str, same: bool) -> str:
    return json.dumps({"summary": summary, "same_fold": same})


def test_generate_uses_llm_when_it_passes(globin_facts):
    client = MockClient([_llm(GOOD, True)])
    out = generate_report(globin_facts, client)
    assert out.source == "llm" and out.verified and out.attempts == 1


def test_generate_regenerates_with_feedback(globin_facts):
    bad = _llm(GOOD.replace("0.904", "0.95"), True)
    client = MockClient([bad, _llm(GOOD, True)])
    out = generate_report(globin_facts, client)
    assert out.source == "llm" and out.attempts == 2
    assert "0.95" in client.calls[1]["user"]  # the checker's complaint was fed back


def test_generate_falls_back_to_template(globin_facts):
    client = MockClient([_llm(GOOD + " It binds oxygen.", True), "not json"])
    out = generate_report(globin_facts, client, max_attempts=3)
    assert out.source == "template" and out.verified
    assert out.text == template_report(globin_facts)
    assert len(client.calls) == 3


def test_no_client_gives_template(globin_facts):
    out = generate_report(globin_facts, None)
    assert out.source == "template" and out.verified and out.attempts == 0
