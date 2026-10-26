"""Grounded vs ungrounded LLM experiment.

Question asked for each pair of CATH domains: *do they share the same fold (same CATH
topology)?*

* Condition A (ungrounded): the model sees only the domain IDs, PDB IDs and the PDB
  entry titles, and answers from memory.
* Condition B (grounded): the model additionally sees ProteinLens's computed facts
  for the pair (TM-scores, RMSD, aligned length, identity, deviating regions) and
  is told to use only those facts.

Answers are scored against CATH. B's explanations (and A's, for reference) are run
through the faithfulness checker. A no-LLM baseline, "same fold if the larger
TM-score > 0.5", is scored on the same pairs.

The reference-set facts deliberately leave out nearest-fold hits: the domains are
themselves in the reference set, so a nearest-fold lookup would leak the answer.
Residue numbers in these facts are positions within each domain (1, 2, 3 ...).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from proteinlens.align import compare
from proteinlens.classify import ReferenceSet, load_reference
from proteinlens.data import ProteinLensError, fetch_pdb_metadata
from proteinlens.llm import LLMClient, LLMError
from proteinlens.report import build_facts
from proteinlens.verify import check_report

CATEGORIES = ["same superfamily", "same fold, different superfamily", "different fold"]

SYSTEM_UNGROUNDED = """You are an expert in protein structure classification. Answer \
from your own knowledge."""

SYSTEM_GROUNDED = """You are an expert in protein structure classification. You are given \
facts computed by a structural alignment program. Base your answer ONLY on these computed \
facts: do not use outside knowledge about the proteins, their function or their \
classification. Every number you write must appear in the facts (rounding is fine). A \
TM-score above 0.5 usually indicates the same fold."""

QUESTION = """Question: do domain A and domain B share the same fold, i.e. the same CATH \
topology (C.A.T level)? Reply with "same_fold" (true/false) and a 2-4 sentence \
"explanation"."""

ANSWER_SCHEMA: dict[str, Any] = {
    "type": "object",
    "properties": {"same_fold": {"type": "boolean"}, "explanation": {"type": "string"}},
    "required": ["same_fold", "explanation"],
    "additionalProperties": False,
}


def relation(sf_a: str, sf_b: str) -> str:
    if sf_a == sf_b:
        return CATEGORIES[0]
    if sf_a.rsplit(".", 1)[0] == sf_b.rsplit(".", 1)[0]:
        return CATEGORIES[1]
    return CATEGORIES[2]


def build_pairs(ref: ReferenceSet, per_category: int = 20, seed: int = 7) -> pd.DataFrame:
    """Sample ``per_category`` distinct domain pairs for each relation (seeded)."""
    rng = np.random.default_rng(seed)
    t = ref.table
    n = len(t)
    by_cat: dict[str, list[tuple[int, int]]] = {c: [] for c in CATEGORIES}
    for i in range(n):
        for j in range(i + 1, n):
            by_cat[relation(t.superfamily[i], t.superfamily[j])].append((i, j))
    rows = []
    for cat in CATEGORIES:
        pool = by_cat[cat]
        used: set[int] = set()
        for k in rng.permutation(len(pool)):
            i, j = pool[k]
            if i in used or j in used:  # keep domains distinct within a category
                continue
            used.update((i, j))
            a, b = t.iloc[i], t.iloc[j]
            rows.append(
                {
                    "domain_a": a.domain,
                    "domain_b": b.domain,
                    "category": cat,
                    "same_fold_truth": a.fold == b.fold,
                    "fold_a": a.fold,
                    "fold_b": b.fold,
                    "superfamily_a": a.superfamily,
                    "superfamily_b": b.superfamily,
                }
            )
            if len(used) == 2 * per_category:
                break
    return pd.DataFrame(rows)


def entry_title(domain: str) -> str:
    try:
        return fetch_pdb_metadata(domain[:4]).get("title") or ""
    except ProteinLensError:
        return ""


def pair_facts(ref: ReferenceSet, a_id: str, b_id: str) -> dict[str, Any]:
    idx = {s.name: k for k, s in enumerate(ref.structures)}
    a, b = ref.structures[idx[a_id]], ref.structures[idx[b_id]]
    return build_facts(a, b, compare(a, b))


def identity_block(row: pd.Series, titles: dict[str, str]) -> str:
    lines = []
    for tag, dom in (("A", row.domain_a), ("B", row.domain_b)):
        lines.append(
            f"Domain {tag}: CATH domain {dom} (PDB entry {dom[:4].upper()}, chain {dom[4]}); "
            f'PDB entry title: "{titles.get(dom, "")}"'
        )
    return "\n".join(lines)


def prompts(row: pd.Series, titles: dict[str, str], facts: dict[str, Any]) -> dict[str, tuple]:
    ident = identity_block(row, titles)
    return {
        "A": (SYSTEM_UNGROUNDED, f"{ident}\n\n{QUESTION}"),
        "B": (
            SYSTEM_GROUNDED,
            f"{ident}\n\nComputed facts (JSON):\n{json.dumps(facts, indent=1)}\n\n{QUESTION}",
        ),
    }


def estimate_cost(
    prompt_texts: list[str], price_in: float, price_out: float, out_tokens: int = 600
) -> float:
    """Rough upper estimate: ~3.5 characters per token, plus a 30 % margin; output tokens
    include an allowance for the model's thinking."""
    tokens_in = sum(len(p) for p in prompt_texts) / 3.5 * 1.3
    tokens_out = out_tokens * len(prompt_texts)
    return (tokens_in * price_in + tokens_out * price_out) / 1e6


@dataclass
class Answer:
    same_fold: bool | None
    explanation: str
    error: str = ""
    cost_usd: float = 0.0


def ask(client: LLMClient, system: str, user: str) -> Answer:
    try:
        r = client.generate(system, user, ANSWER_SCHEMA)
        d = json.loads(r.text)
        return Answer(bool(d["same_fold"]), str(d["explanation"]), cost_usd=r.cost_usd)
    except (LLMError, json.JSONDecodeError, KeyError, TypeError) as err:
        return Answer(None, "", error=str(err))


def run_condition_rows(
    pairs: pd.DataFrame, client: LLMClient, ref: ReferenceSet | None = None
) -> pd.DataFrame:
    """Ask both conditions for every pair; score answers and faithfulness."""
    ref = ref or load_reference()
    titles = {d: entry_title(d) for d in set(pairs.domain_a) | set(pairs.domain_b)}
    out = []
    for _, row in pairs.iterrows():
        facts = pair_facts(ref, row.domain_a, row.domain_b)
        for cond, (system, user) in prompts(row, titles, facts).items():
            ans = ask(client, system, user)
            check = check_report(ans.explanation, facts) if ans.explanation else None
            out.append(
                {
                    "domain_a": row.domain_a,
                    "domain_b": row.domain_b,
                    "category": row.category,
                    "condition": cond,
                    "model": client.model,
                    "same_fold_truth": bool(row.same_fold_truth),
                    "answer": ans.same_fold,
                    "correct": (ans.same_fold == bool(row.same_fold_truth))
                    if ans.same_fold is not None
                    else False,
                    "faithful": check.passed if check else False,
                    "issue_kinds": ",".join(sorted({i.kind for i in check.issues}))
                    if check
                    else "",
                    "numbers_checked": check.numbers_checked if check else 0,
                    "explanation": ans.explanation,
                    "error": ans.error,
                    "cost_usd": ans.cost_usd,
                }
            )
    return pd.DataFrame(out)


def tm_rule_rows(pairs: pd.DataFrame, ref: ReferenceSet | None = None) -> pd.DataFrame:
    """The no-LLM baseline: same fold if the larger TM-score is above 0.5."""
    ref = ref or load_reference()
    rows = []
    for _, row in pairs.iterrows():
        f = pair_facts(ref, row.domain_a, row.domain_b)
        al = f["alignment"]
        pred = f["interpretation"]["likely_same_fold"]
        rows.append(
            {
                **row.to_dict(),
                "tm_score_a": al["tm_score_normalized_by_a"],
                "tm_score_b": al["tm_score_normalized_by_b"],
                "rmsd": al["rmsd_angstrom"],
                "seq_identity_percent": al["sequence_identity_percent"],
                "tm_rule_same_fold": pred,
                "tm_rule_correct": pred == bool(row.same_fold_truth),
            }
        )
    return pd.DataFrame(rows)


def bootstrap(values: np.ndarray, seed: int = 0, n_boot: int = 2000) -> list[float]:
    if len(values) == 0:
        return [float("nan"), float("nan")]
    rng = np.random.default_rng(seed)
    means = values[rng.integers(0, len(values), (n_boot, len(values)))].mean(axis=1)
    return [round(float(np.percentile(means, 2.5)), 4), round(float(np.percentile(means, 97.5)), 4)]


def summarize(results: pd.DataFrame) -> dict[str, Any]:
    """Accuracy and faithfulness per model x condition (x category), with bootstrap CIs."""
    out: dict[str, Any] = {}
    for (model, cond), g in results.groupby(["model", "condition"]):
        acc = g["correct"].to_numpy(dtype=float)
        faith = g["faithful"].to_numpy(dtype=float)
        out[f"{model} / {cond}"] = {
            "n": int(len(g)),
            "accuracy": round(float(acc.mean()), 4),
            "accuracy_ci95": bootstrap(acc),
            "faithfulness_pass_rate": round(float(faith.mean()), 4),
            "faithfulness_ci95": bootstrap(faith),
            "unparseable_answers": int((g["answer"].isna()).sum()),
            "by_category": {
                c: round(float(gc["correct"].mean()), 4) for c, gc in g.groupby("category")
            },
            "issue_kinds": g["issue_kinds"]
            .str.split(",")
            .explode()
            .value_counts()
            .drop("", errors="ignore")
            .to_dict(),
            "cost_usd": round(float(g["cost_usd"].sum()), 4),
        }
    return out


def write_pending(out_dir: Path, rule: pd.DataFrame, command: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    rule.to_csv(out_dir / "pairs.csv", index=False)
    acc = rule["tm_rule_correct"].to_numpy(dtype=float)
    by_cat = rule.groupby("category")["tm_rule_correct"].mean().round(4).to_dict()
    status = {
        "status": "PENDING",
        "reason": "no ANTHROPIC_API_KEY was available when this was run",
        "how_to_run": command,
        "n_pairs": int(len(rule)),
        "tm_rule_baseline": {
            "accuracy": round(float(acc.mean()), 4),
            "ci95": bootstrap(acc),
            "by_category": by_cat,
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(status, indent=2))
