# Grounded vs ungrounded LLM experiment

**Status: PENDING — run `make eval-llm` with an API key** (`ANTHROPIC_API_KEY` was not available
when this repository was built, so no LLM answers exist yet and no LLM accuracy is reported).

```bash
export ANTHROPIC_API_KEY=sk-ant-...      # never commit this
make eval-llm                            # prints a cost estimate, then asks to continue
```

What is already real (computed without any LLM):

* `pairs.csv` — the 60 CATH domain pairs (seed 7): 20 from the same superfamily, 20 from the same
  fold but different superfamilies, 20 from different folds, with ProteinLens's computed TM-scores.
* `summary.json` → `tm_rule_baseline` — the no-LLM rule "same fold if the larger TM-score > 0.5"
  scores **83.3 % [73.3, 93.3]**: 100 % on same-superfamily pairs, 100 % on different-fold pairs and
  50 % on same-fold/different-superfamily pairs. Any grounded LLM should be compared against this
  rule, not only against the ungrounded LLM.

What `make eval-llm` will add (`responses.csv`, an updated `summary.json`):

| | Condition A — ungrounded | Condition B — grounded |
|---|---|---|
| Model sees | domain IDs, PDB IDs, PDB entry titles | the same + ProteinLens facts JSON |
| Instructions | answer from your own knowledge | use only the computed facts |
| Scored | same-fold accuracy vs CATH; checker run on explanations (for reference) | same-fold accuracy vs CATH; faithfulness checker pass rate |

* Models: `claude-haiku-5-5` on all 60 pairs × 2 conditions (120 calls); `claude-opus-5-5` on an
  18-pair balanced subset × 2 (36 calls). 156 calls in total (cap: 300).
* Every response is cached in `data/cache/llm/`; re-running costs nothing for cached calls.
* The script stops before calling the API if the estimated cost exceeds `--max-cost` (default $3).
* The grounded facts deliberately omit nearest-fold hits (the domains are in the reference set, so
  that would leak the answer).
* For condition A the checker is not a rule-violation test (A was allowed to use memory); its pass
  rate shows how much of an ungrounded explanation is supported by computed numbers.
