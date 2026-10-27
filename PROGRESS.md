# PROGRESS

Start: **2026-10-26 23:32 UTC**. Budget: 12 h working time (target ~11 h), ~$100 credits (target < $90).
Cost readout: none available in this environment → elapsed time and work done are used as a proxy.

## Network check (M0, 2026-10-26 23:33 UTC)
| Host | Result |
|---|---|
| files.rcsb.org (PDB files) | ✅ 200 |
| data.rcsb.org (metadata API) | ✅ 200 |
| alphafold.ebi.ac.uk (API) | ✅ 200 |
| rest.uniprot.org | ✅ 200 |
| scop.berkeley.edu (SCOPe) | ❌ 403 — blocked. Not needed: CATH used instead (see DECISIONS.md) |
| download.cathdb.info (CATH v4.4 release files) | ✅ 200 |
| www.cathdb.info (per-domain PDB files) | ✅ 200 |
| pypi.org | ✅ 200 (all packages installed cleanly, incl. tmtools 0.3.0) |
| api.anthropic.com | reachable (401 without key). **No `ANTHROPIC_API_KEY` in the environment.** |

Optional ACTION NEEDED (only if you want SCOPe instead of CATH): allow `scop.berkeley.edu`.

## Milestone log
| Milestone | Started (UTC) | Finished | Elapsed total | Notes |
|---|---|---|---|---|
| M0 Setup | 23:32 | 23:36 | 0:04 | scaffold, CI, network check, all packages install |
| M1 Data layer | 23:36 | 23:39 | 0:07 | RCSB/AlphaFold/CATH/uploads, cache, cleaning; 12 tests (offline + network) |
| M2 Comparison engine | 23:39 | 23:41 | 0:09 | tmtools TM-align, symmetric results, deviating regions, CLI compare |
| M3 Benchmark set | 23:41 | 23:44 | 0:12 | 288 CATH S40 domains / 20 folds / 60 superfamilies; deterministic rebuild |
| M4 Benchmark results | 23:44 | 23:54 | 0:22 | all-vs-all 41,328 pairs in 357 s; results + 4 figures + failure analysis |
| M5 Reports + checker | 23:48 | 23:56 | 0:24 | facts, template, LLM narration, checker; corruption tests; CLI classify/report |
| M6 LLM experiment | 23:56 | 23:57 | 0:25 | pipeline + mock tests; PENDING (no key); TM-rule baseline 83.3 % on 60 pairs |
| M7 Streamlit app | 23:57 | 00:02 | 0:30 | Compare/Classify/Benchmarks/About; AppTest e2e; Playwright screenshots |
| M8 Documentation | 00:02 | 00:06 | 0:34 | README, WRITEUP, UNDERSTANDING (+15 Q&A), CITATIONS, DEPLOY |

## Cost tracking
No usage/cost readout is exposed in this environment. Proxy: elapsed time and work done.
No live LLM API calls have been made (no `ANTHROPIC_API_KEY`). Spend so far is development only.

## Headline results so far (from `benchmarks/results/summary.json`)
- Structure 1-NN fold accuracy 99.7 % [98.96, 100], superfamily 95.8 % [93.4, 97.9] (n = 288)
- Sequence 1-NN fold 80.9 % [76.0, 85.4], superfamily 78.1 % [73.3, 83.0]
- Remote homologs (n = 133): superfamily structure 93.2 % vs sequence 54.9 %

## M9 — Final verification (2026-10-27 00:06–00:15 UTC)
- Fresh clone of the pushed branch → `make setup` (39 s) → `make lint` ✅ → `make test` ✅ (63 passed at that commit).
- `make dataset` from an empty cache in the fresh clone rebuilt `benchmarks/data/dataset.csv` **byte-identical**
  (md5 00e6750f8398edcac57f0a61bb4cc20b) and `reference_set.npz` unchanged (`git status` clean).
- `make benchmark-smoke` ✅ in the fresh clone.
- **Full all-vs-all benchmark re-run in the fresh clone (empty cache, 41,328 alignments): every result file
  (`summary.json`, `per_query.csv`, `per_fold.csv`, `failures.csv`, `superfamily_failures.csv`, `pair_stats.csv`,
  `accuracy_vs_identity.csv`) is byte-identical to the committed one, and all `pairwise.npz` matrices are identical.**
- Network tests (`make test-all`) ✅: RCSB, AlphaFold DB and CATH downloads.
- CI (GitHub Actions, Python 3.12 + 3.13) green on every milestone commit M3–M8.
- App: headless AppTest e2e (with no key, with a mock LLM, and with a failing LLM → template fallback) ✅;
  real-browser run with Playwright/Chromium for the screenshots in `docs/figures/` ✅.
- README / WRITEUP headline numbers are checked against `summary.json` by `tests/test_docs_numbers.py` ✅.

## Final summary
All milestones M0–M9 done. Elapsed ≈ 0:45 of the 12 h budget. No live LLM calls were made, so API spend
is $0; development cost is not visible from inside the environment (no readout).
**PENDING:** the live grounded-vs-ungrounded LLM run (`ANTHROPIC_API_KEY=... make eval-llm`; the script's own upper estimate is $0.54 for
120 Haiku 5.5 + 36 Opus 5.5 calls; it aborts above $3).

## Extension E1 — harder second benchmark + abstention (2026-10-27 02:24–03:05 UTC)
User asked to continue with the best-fit improvement. Chose depth on the core (M3/M4) over new scope.
- Design pre-specified and committed (commit "Extension: pre-specify harder second benchmark …") before any
  hard-set result existed: 600 domains / 120 random folds / 5 per fold, seed 2026.
- Main-set dataset byte-identical after refactoring the builder (md5 checked); main-set result files unchanged
  except two added keys in `summary.json` (and three more descriptive keys added afterwards).
- Hard run: 179,700 pairs; hit the 20-min limit once and resumed automatically (chunks on disk).
- Results (`benchmarks/results_hard/summary.json`): structure fold 88.0 % [85.5, 90.3] vs sequence 65.5 %;
  superfamily (queries with a partner, n = 493) 95.3 % vs 78.7 %; remote homologs 91.5 % vs 56.0 %.
- Failure analysis: 95.8 % of fold errors have a longer nearest neighbour (size bias of query-normalized TM).
- Abstention at TM 0.5: hard set 90.0 % coverage / 93.2 % accuracy; main set 98.3 % / 100.0 %.
- README, WRITEUP, UNDERSTANDING, app (Benchmarks tabs) updated; docs-number test covers hard-set numbers.

## Extension E2 — better score (dev/test) + results website (2026-10-27 03:10 → 19:30 UTC, VM paused in between)
User asked to "make it better" and create a simple website for the results.
- Main + hard became development sets; `select_score.py` (fixed rule) picked the coverage-aware score.
- Held-out test set (seed 7777, 600 domains, none shared) pre-registered and committed before running.
  The VM paused mid-run; the resumable runner continued from the saved chunks.
- Test result: fold accuracy 84.8 % → 89.5 % (+4.7 [2.7, 6.7]); adopted as default in classify/CLI/report/app.
- Website: `site/build_site.py` + `site/template.html` → `site/index.html`; checked at 1300 px and 400 px,
  light and dark (no horizontal overflow); freshness test in CI; GitHub Pages workflow (needs Pages enabled).
