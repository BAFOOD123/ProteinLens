"""Run the grounded vs ungrounded LLM experiment (see proteinlens/experiment.py).

    make eval-llm                 # needs ANTHROPIC_API_KEY
    python benchmarks/llm_experiment.py --max-cost 3 --yes

Plan: the cheap model (claude-haiku-5-5) answers all 60 pairs in both conditions
(120 calls); the strong model (claude-opus-5-5) answers a balanced 18-pair subset in
both conditions (36 calls). Every response is cached in data/cache/llm/, so a re-run
only pays for what is missing. A cost estimate is printed before any call, and the
run stops if it exceeds --max-cost or --max-calls.

Without an API key the script still builds the pairs, scores the no-LLM TM-score
rule, and writes results marked PENDING.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from proteinlens.classify import load_reference  # noqa: E402
from proteinlens.experiment import (  # noqa: E402
    build_pairs,
    entry_title,
    estimate_cost,
    pair_facts,
    prompts,
    run_condition_rows,
    summarize,
    tm_rule_rows,
    write_pending,
)
from proteinlens.llm import PRICES, AnthropicClient, CachedClient, api_key_available  # noqa: E402

OUT = REPO / "benchmarks" / "results" / "llm_experiment"
BULK_MODEL, STRONG_MODEL = "claude-haiku-5-5", "claude-opus-5-5"
COMMAND = "ANTHROPIC_API_KEY=... make eval-llm"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-category", type=int, default=20)
    ap.add_argument("--strong-per-category", type=int, default=6)
    ap.add_argument("--max-cost", type=float, default=3.0, help="USD; abort if estimate exceeds")
    ap.add_argument("--max-calls", type=int, default=300)
    ap.add_argument("--yes", action="store_true", help="run without the confirmation pause")
    args = ap.parse_args()

    ref = load_reference()
    pairs = build_pairs(ref, per_category=args.per_category)
    rule = tm_rule_rows(pairs, ref)
    strong = pairs.groupby("category", group_keys=False).head(args.strong_per_category)
    print(f"{len(pairs)} pairs ({pairs['category'].value_counts().to_dict()})")
    print(f"TM-score > 0.5 rule (no LLM): accuracy {rule['tm_rule_correct'].mean():.3f}")

    if not api_key_available():
        write_pending(OUT, rule, COMMAND)
        print(f"No ANTHROPIC_API_KEY: wrote PENDING results to {OUT}")
        return 0

    # ---- cost estimate before spending anything ----
    titles = {d: entry_title(d) for d in set(pairs.domain_a) | set(pairs.domain_b)}

    def texts(df: pd.DataFrame) -> list[str]:
        out = []
        for _, row in df.iterrows():
            for system, user in prompts(
                row, titles, pair_facts(ref, row.domain_a, row.domain_b)
            ).values():
                out.append(system + user)
        return out

    bulk_texts, strong_texts = texts(pairs), texts(strong)
    n_calls = len(bulk_texts) + len(strong_texts)
    est = estimate_cost(bulk_texts, *PRICES[BULK_MODEL]) + estimate_cost(
        strong_texts, *PRICES[STRONG_MODEL]
    )
    print(f"Planned calls: {len(bulk_texts)} x {BULK_MODEL} + {len(strong_texts)} x {STRONG_MODEL}")
    print(f"Estimated cost (upper, uncached): ${est:.2f}  (limit ${args.max_cost:.2f})")
    if n_calls > args.max_calls or est > args.max_cost:
        print("Aborting: over the call or cost limit.")
        return 2
    if not args.yes:
        input("Press Enter to start (Ctrl-C to cancel)... ")

    frames = []
    for model, subset in ((BULK_MODEL, pairs), (STRONG_MODEL, strong)):
        client = CachedClient(AnthropicClient(model=model, effort="low", max_tokens=4000))
        print(f"Running {model} on {len(subset)} pairs ...", flush=True)
        frames.append(run_condition_rows(subset, client, ref))
    results = pd.concat(frames, ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)
    rule.to_csv(OUT / "pairs.csv", index=False)
    results.to_csv(OUT / "responses.csv", index=False)
    summary = {
        "status": "COMPLETE",
        "n_pairs": len(pairs),
        "tm_rule_baseline_accuracy": round(float(rule["tm_rule_correct"].mean()), 4),
        "conditions": summarize(results),
        "total_cost_usd_this_run": round(float(results["cost_usd"].sum()), 4),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["conditions"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
