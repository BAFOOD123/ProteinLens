"""Grounded comparison reports.

1. :func:`build_facts` collects every computed number into one JSON-friendly dict.
2. :func:`template_report` turns the facts into plain English with fixed sentences
   (always available, always faithful).
3. :func:`generate_report` asks an LLM to narrate *only* the facts, checks the
   result with :mod:`proteinlens.verify`, regenerates up to ``max_attempts``
   times with the checker's complaints as feedback, and falls back to the
   template if no draft passes.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from typing import Any

from proteinlens.align import (
    DEVIATION_CUTOFF,
    RANDOM_TM,
    SAME_FOLD_TM,
    AlignmentResult,
)
from proteinlens.classify import Hit
from proteinlens.data import Structure, cache_dir
from proteinlens.llm import LLMClient, LLMError
from proteinlens.verify import CheckResult, check_report

ALIGNED_PAIR_CUTOFF = 5.0  # Å, TM-align's definition of an aligned pair


def _protein_facts(s: Structure) -> dict[str, Any]:
    m = s.metadata
    out: dict[str, Any] = {
        "name": s.name,
        "source": s.source,
        "chain": s.chain,
        "length_residues": len(s),
    }
    for key, fact in (
        ("description", "description"),
        ("organism", "organism"),
        ("title", "title"),
        ("method", "method"),
        ("resolution", "resolution_angstrom"),
        ("mean_plddt", "mean_plddt"),
    ):
        if m.get(key) not in (None, "", []):
            out[fact] = m[key]
    if m.get("n_chain_breaks"):
        out["chain_breaks"] = m["n_chain_breaks"]
    return out


def _hit_facts(hits: list[Hit] | None) -> list[dict[str, Any]]:
    return [
        {
            "cath_fold": h.fold,
            "fold_name": h.fold_name,
            "reference_domain": h.domain,
            "tm_score": round(h.tm_coverage, 3),  # normalized by the longer protein
        }
        for h in hits or []
    ]


def build_facts(
    a: Structure,
    b: Structure,
    result: AlignmentResult,
    hits_a: list[Hit] | None = None,
    hits_b: list[Hit] | None = None,
) -> dict[str, Any]:
    """Everything a report is allowed to say, rounded the way it should be printed."""
    facts: dict[str, Any] = {
        "protein_a": _protein_facts(a),
        "protein_b": _protein_facts(b),
        "alignment": {
            "tm_score_normalized_by_a": round(result.tm_score_a, 3),
            "tm_score_normalized_by_b": round(result.tm_score_b, 3),
            "rmsd_angstrom": round(result.rmsd, 2),
            "aligned_residue_pairs": result.aligned_length,
            "sequence_identity_percent": round(100 * result.seq_identity, 1),
        },
        "interpretation": {
            "same_fold_threshold": SAME_FOLD_TM,
            "likely_same_fold": result.likely_same_fold,
            "typical_tm_score_of_unrelated_proteins": RANDOM_TM,
            "aligned_pair_cutoff_angstrom": ALIGNED_PAIR_CUTOFF,
            "deviation_cutoff_angstrom": DEVIATION_CUTOFF,
        },
        "deviating_regions": [
            {
                "residues_a": g.residues_a,
                "residues_b": g.residues_b,
                "aligned_pairs": g.length,
                "mean_distance_angstrom": g.mean_distance,
                "max_distance_angstrom": g.max_distance,
            }
            for g in result.top_regions
        ],
    }
    if hits_a:
        facts["nearest_known_folds_a"] = _hit_facts(hits_a)
    if hits_b:
        facts["nearest_known_folds_b"] = _hit_facts(hits_b)
    return facts


# --------------------------------------------------------------------------------------
# Template report
# --------------------------------------------------------------------------------------


def _label(p: dict[str, Any]) -> str:
    desc = p.get("description")
    return f"{p['name']} ({desc})" if desc else p["name"]


def template_report(facts: dict[str, Any]) -> str:
    a, b = facts["protein_a"], facts["protein_b"]
    al, it = facts["alignment"], facts["interpretation"]
    tm_max = max(al["tm_score_normalized_by_a"], al["tm_score_normalized_by_b"])
    lines = [
        f"ProteinLens compared {_label(a)}, {a['length_residues']} residues, with "
        f"{_label(b)}, {b['length_residues']} residues.",
        f"TM-align superimposed {al['aligned_residue_pairs']} residue pairs closer than "
        f"{it['aligned_pair_cutoff_angstrom']:g} Å, with an RMSD of {al['rmsd_angstrom']:.2f} Å.",
        f"The TM-score is {al['tm_score_normalized_by_a']:.3f} normalized by {a['name']} and "
        f"{al['tm_score_normalized_by_b']:.3f} normalized by {b['name']}.",
    ]
    if it["likely_same_fold"]:
        lines.append(
            f"Because the larger TM-score ({tm_max:.3f}) is above {it['same_fold_threshold']:g}, "
            "the two structures likely share the same fold."
        )
    else:
        lines.append(
            f"Because the larger TM-score ({tm_max:.3f}) is not above "
            f"{it['same_fold_threshold']:g}, the two structures likely have different folds "
            "(unrelated proteins typically score about "
            f"{it['typical_tm_score_of_unrelated_proteins']:g})."
        )
    lines.append(
        f"Sequence identity over the aligned pairs is {al['sequence_identity_percent']:.1f}%."
    )
    regions = facts["deviating_regions"]
    if regions:
        parts = [
            f"residues {g['residues_a']} of {a['name']} versus {g['residues_b']} of {b['name']} "
            f"(mean {g['mean_distance_angstrom']:.2f} Å, max {g['max_distance_angstrom']:.2f} Å)"
            for g in regions
        ]
        lines.append("The largest deviations after superposition are " + "; ".join(parts) + ".")
    else:
        lines.append(
            "No stretch of 3 or more aligned residues deviates by more than "
            f"{it['deviation_cutoff_angstrom']:g} Å after superposition."
        )
    for key, p in (("nearest_known_folds_a", a), ("nearest_known_folds_b", b)):
        hits = facts.get(key)
        if hits:
            h = hits[0]
            lines.append(
                f"The closest CATH reference fold for {p['name']} is {h['cath_fold']} "
                f'"{h["fold_name"]}" (TM-score {h["tm_score"]:.3f} to {h["reference_domain"]}).'
            )
    return " ".join(lines)


# --------------------------------------------------------------------------------------
# LLM narration with verification
# --------------------------------------------------------------------------------------

SYSTEM_PROMPT = """You write short, plain-English comparisons of two protein structures \
for a general science audience.

Rules (a strict automatic checker enforces them):
- Use ONLY the facts in the JSON you are given. Do not add anything from memory: no \
biological function, no evolutionary history, no fold or family names that are not in the \
facts, no PDB IDs or other identifiers that are not in the facts.
- Every number you write must appear in the facts (you may round, e.g. 0.904 -> 0.90). \
Do not compute new numbers such as differences or averages.
- The fold verdict must agree with "likely_same_fold". Explain it with the TM-score rule \
(TM-score above 0.5 usually means the same fold).
- Mention the deviating regions if there are any.
- 3 to 6 sentences in "summary". Put the verdict in "same_fold"."""

REPORT_SCHEMA: dict[str, Any] = {
    "type": "object",
    "properties": {
        "summary": {"type": "string"},
        "same_fold": {"type": "boolean"},
    },
    "required": ["summary", "same_fold"],
    "additionalProperties": False,
}


@dataclass
class ReportOutcome:
    text: str
    source: str  # "llm" or "template"
    verified: bool  # passed the faithfulness checker
    attempts: int  # LLM attempts made (0 = no LLM)
    checks: list[CheckResult] = field(default_factory=list)
    note: str = ""
    cost_usd: float = 0.0


def _user_prompt(facts: dict[str, Any], feedback: str = "") -> str:
    msg = "Facts (JSON):\n" + json.dumps(facts, indent=1)
    if feedback:
        msg += (
            "\n\nYour previous draft failed the checker for these reasons; fix them and "
            "write it again:\n" + feedback
        )
    return msg


def _log(event: dict[str, Any]) -> None:
    """Append one line per generated report so the pass rate can be computed later."""
    path = cache_dir() / "llm" / "report_log.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps({"time": time.time(), **event}) + "\n")


def generate_report(
    facts: dict[str, Any], client: LLMClient | None, max_attempts: int = 3
) -> ReportOutcome:
    """LLM report that passed the checker, or the (verified) template report."""
    template = template_report(facts)
    if client is None:
        check = check_report(template, facts)
        return ReportOutcome(template, "template", check.passed, 0, [check], "no API key")

    checks: list[CheckResult] = []
    feedback = ""
    cost = 0.0
    for attempt in range(1, max_attempts + 1):
        try:
            resp = client.generate(SYSTEM_PROMPT, _user_prompt(facts, feedback), REPORT_SCHEMA)
            cost += resp.cost_usd if not resp.cached else 0.0
            data = json.loads(resp.text)
            text, verdict = str(data["summary"]).strip(), bool(data["same_fold"])
        except (LLMError, json.JSONDecodeError, KeyError, TypeError) as err:
            _log({"model": client.model, "attempt": attempt, "passed": False, "error": str(err)})
            checks.append(CheckResult(False))
            feedback = f"- response was not usable: {err}"
            continue
        check = check_report(text, facts, claimed_same_fold=verdict)
        checks.append(check)
        _log({"model": client.model, "attempt": attempt, "passed": check.passed})
        if check.passed:
            return ReportOutcome(text, "llm", True, attempt, checks, cost_usd=cost)
        feedback = "\n".join(f"- {i.text!r}: {i.detail}" for i in check.issues)

    final = check_report(template, facts)
    return ReportOutcome(
        template,
        "template",
        final.passed,
        max_attempts,
        checks + [final],
        note=f"no LLM draft passed the checker in {max_attempts} attempts",
        cost_usd=cost,
    )


def report_pass_rate() -> dict[str, float]:
    """Pass rate of LLM drafts logged so far (per attempt and per report)."""
    path = cache_dir() / "llm" / "report_log.jsonl"
    if not path.exists():
        return {"attempts": 0, "attempt_pass_rate": float("nan")}
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    passed = sum(r["passed"] for r in rows)
    return {"attempts": len(rows), "attempt_pass_rate": passed / len(rows) if rows else 0.0}
