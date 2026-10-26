# CLAUDE.md — ProteinLens working memory

Re-read this file at every milestone boundary.

## Mission
Build **ProteinLens**: an open-source Python project that
1. compares protein 3D structures (TM-align-style alignment, computed — not recalled),
2. classifies a structure's fold against a ground-truth answer key (SCOPe/CATH) and *measures* accuracy,
3. writes plain-English comparison reports where an LLM only narrates computed facts and an
   automatic faithfulness checker verifies every number,
4. runs a grounded-vs-ungrounded LLM experiment scored against ground truth.

Audience: a high-school computational-biology student's flagship portfolio project.

## Hard rules
- **Everything real.** Every number in README / write-up / app comes from code that was run. No
  invented, estimated, rounded-up or cherry-picked results. Mediocre results are reported honestly.
- **Explainable.** Clear, simple, commented code. `docs/UNDERSTANDING.md` explains it all.
- **Reproducible.** Fixed seeds, pinned deps, one-command setup, cached downloads.
- No hidden test-set tuning (tune only on a validation split, say so).
- Don't claim the LLM "understands" proteins. Claim: structure is computed, narration is verified.
- No docking / MD / training neural nets. Depth over breadth.
- No secrets in repo; `.env.example` provided; large data gitignored (`data/cache/`).

## Operating rules
- Track elapsed time in `PROGRESS.md` (start: 2026-10-26 23:32 UTC). Milestone over by >50% → cut
  scope, note in `DECISIONS.md`, move on. Last 30 min = final verification.
- Commit after every meaningful step; **push after every milestone and at least hourly**
  (`git push -u origin claude/charming-lovelace-s8dipt`). Tests pass at every milestone boundary.
- Don't ask the user unless truly blocked. Record decisions in `DECISIONS.md` (one-line reason).
- Long jobs: resumable chunks ≤ ~20 min, incremental writes to disk, run in background with logs,
  read only log tails. Never rely on a process surviving.
- Verify URLs / packages / model names before relying on them.
- Python 3.11+, `pyproject.toml`, `ruff`, `pytest`, type hints. Modest dependencies.

## Budget
~$100 total. ~$85 build, ~$5 live LLM eval (M6), rest reserve. At ~70% before M7: stop adding
scope. Hard stop ~$90. Cost proxy = elapsed time/work (no cost readout). Cheapest capable model for
bulk LLM calls; print cost estimate before any live run; cap ~300 live calls; cache all responses.

Cut order (cut from the bottom): M1–M4 never cut → M8 docs → M5 report/checker → M7 app → M6 live.

## Milestones (time-boxes)
- M0 Setup (0.5h): scaffold, pyproject, Makefile (setup/test/lint/benchmark/eval-llm/app), CI,
  CLAUDE/PROGRESS/DECISIONS, network check.
- M1 Data layer (1h): fetch/cache by PDB ID, UniProt (AlphaFold), upload; chains, ranges, missing
  residues, altlocs, HETATM strip, metadata. Tests offline (fixtures) + online.
- M2 Comparison engine (1.5h): typed result (TM both norms, RMSD, aligned len, seq id, R/t,
  per-residue distances, top deviating regions). Known-answer tests; CLI `compare`.
- M3 Benchmark set + baselines (1.5h): SCOPe ASTRAL ≤40% or CATH S40; ~300 domains, ≥15 folds;
  documented selection + seed; sequence baseline + majority baseline.
- M4 Structural classification benchmark (1.5h): LOO 1-NN by TM-score, all-vs-all cached &
  parallel; fold + superfamily accuracy w/ 95% bootstrap CI; remote-homolog subset; confusion
  matrix; accuracy-vs-seq-id plot; honest failure analysis. `make benchmark`.
- M5 Grounded report + faithfulness checker (1.5h): facts object, LLM narration (structured), checker
  verifying every number/claim, regenerate up to N then template fallback; corruption tests.
- M6 Grounded vs ungrounded experiment (1h): ~60 pairs (same SF / same fold diff SF / diff fold),
  conditions A (IDs only) vs B (facts). Score vs ground truth; faithfulness on B. PENDING if no key.
- M7 Streamlit app (1.5h): Compare / Classify / Benchmarks pages; works without key; examples;
  deploy instructions; headless e2e test; screenshots.
- M8 Docs (1h): README (quickstart, screenshot, real numbers + CIs, limitations, "how built"),
  docs/WRITEUP.md, docs/UNDERSTANDING.md (+15 Q&A), CITATIONS.md.
- M9 Final verification (0.5h): fresh clone venv, setup/test/benchmark(smoke), lint, CI, app,
  README numbers match result files, final PROGRESS.

## Definition of done
- `make setup`, `make test`, `make lint` pass from fresh clone.
- Full benchmark run; results + figures committed and match README exactly.
- Faithfulness checker tested and catches corrupted reports.
- Grounded-vs-ungrounded run with real results, or clearly PENDING with the exact command.
- Streamlit app runs with and without API key; tested end to end.
- README, WRITEUP, UNDERSTANDING, CITATIONS, LICENSE, CI complete.
- No secrets; `.env.example`; large data gitignored. PROGRESS.md and DECISIONS.md up to date.

## Final hand-off message
(1) what was built, (2) headline benchmark numbers exactly as in README, (3) what is PENDING/limited,
(4) three decisions most worth review, (5) exact commands to run the app and reproduce the benchmark.
