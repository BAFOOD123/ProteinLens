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
| M8 Documentation | 00:02 | 00:10 | 0:38 | README, WRITEUP, UNDERSTANDING (+15 Q&A), CITATIONS, DEPLOY |

## Cost tracking
No usage/cost readout is exposed in this environment. Proxy: elapsed time and work done.
No live LLM API calls have been made (no `ANTHROPIC_API_KEY`). Spend so far is development only.

## Headline results so far (from `benchmarks/results/summary.json`)
- Structure 1-NN fold accuracy 99.7 % [98.96, 100], superfamily 95.8 % [93.4, 97.9] (n = 288)
- Sequence 1-NN fold 80.9 % [76.0, 85.4], superfamily 78.1 % [73.3, 83.0]
- Remote homologs (n = 133): superfamily structure 93.2 % vs sequence 54.9 %
