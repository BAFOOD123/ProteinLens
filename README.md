# ProteinLens

**Compare protein 3D structures, classify their fold against CATH, and get plain-English
reports in which every number is checked against what was actually computed.**

![ProteinLens compare page: hemoglobin α vs β superimposed](docs/figures/compare.png)

A language model can talk about proteins from memory, and it can be confidently wrong. ProteinLens
keeps the model away from the science:

1. **Structural similarity is computed**, with TM-align on C-α atoms, not recalled.
2. **Fold classification is measured** against a ground-truth answer key (CATH v4.4) on a
   non-redundant benchmark, with confidence intervals and a failure analysis.
3. **The LLM only narrates a facts object**, and an automatic **faithfulness checker** verifies
   every number, identifier, fold name and same/different-fold claim in its output. Reports that
   fail are regenerated, then replaced by a template report. Without an API key the template is
   used, so everything works offline.
4. A **grounded vs ungrounded experiment** asks whether giving an LLM the computed facts changes
   its accuracy (pipeline built and tested; live run PENDING, see below).

## Headline results

Leave-one-out 1-nearest-neighbour classification of **288 CATH S40 domains** (20 folds,
60 superfamilies; no two domains share ≥ 40 % sequence identity by CATH's clustering).
95 % percentile-bootstrap confidence intervals (2,000 resamples, seed 0). Source:
[`benchmarks/results/summary.json`](benchmarks/results/summary.json).

| Method | Fold (CATH topology) | Superfamily |
|---|---|---|
| **Structure** (TM-score 1-NN) | **99.7 %** [99.0, 100.0] (287/288) | **95.8 %** [93.4, 97.9] (276/288) |
| Sequence (Smith-Waterman 1-NN) | 80.9 % [76.0, 85.4] (233/288) | 78.1 % [73.3, 83.0] (225/288) |
| Majority class (leave-one-out) | 0.0 % (largest class = 5.2 % of the set) | 0.0 % (largest class = 4.5 %) |

**Remote homologs** (n = 133 queries whose same-superfamily relatives all share < 25 % sequence
identity with them):

| Method | Fold | Superfamily |
|---|---|---|
| **Structure** | **99.2 %** [97.7, 100.0] (132/133) | **93.2 %** [88.7, 97.0] (124/133) |
| Sequence | 60.2 % [51.9, 69.2] (80/133) | 54.9 % [45.9, 63.2] (73/133) |

Structure minus sequence (paired bootstrap): +17.7 points [12.8, 22.6] at superfamily level on
all queries, and **+38.4 points [29.3, 47.4] on remote homologs**. As sequence identity falls below
20 %, the sequence baseline collapses (30.0 % superfamily accuracy below 15 % identity, n = 50) while
structure stays at 92.0 %:

<p float="left">
  <img src="benchmarks/figures/accuracy.png" width="58%" />
  <img src="benchmarks/figures/accuracy_vs_identity.png" width="40%" />
</p>

**Read these numbers with care** (details in [Limitations](#limitations) and the
[failure analysis](benchmarks/results/FAILURE_ANALYSIS.md)): the benchmark is clean and small
(20 well-populated folds, X-ray domains ≤ 2.5 Å), and CATH itself is built partly from structure
comparison, so agreement with CATH is partly expected. The one fold error and the 12 superfamily
errors are all between neighbouring topologies or superfamilies of the same fold.

**Grounded vs ungrounded LLM experiment: PENDING — run `make eval-llm` with an API key.** No LLM
accuracy is reported because no API key was available when this was built. Already computed on the
experiment's 60 pairs: the no-LLM rule "same fold if TM-score > 0.5" scores 83.3 % [73.3, 93.3]
(100 % on same-superfamily pairs, 50 % on same-fold/different-superfamily pairs, 100 % on
different-fold pairs). See [`benchmarks/results/llm_experiment/`](benchmarks/results/llm_experiment/README.md).

## Quickstart (60 seconds)

```bash
git clone https://github.com/bafood123/proteinlens.git && cd proteinlens
make setup                 # creates .venv and installs pinned dependencies (Python ≥ 3.12)
make test                  # offline tests (bundled fixtures; no network needed)
make app                   # web app at http://localhost:8501
```

Command line:

```bash
.venv/bin/proteinlens compare 4HHB:A 4HHB:B        # hemoglobin α vs β
.venv/bin/proteinlens compare 1UBQ 1WM3 --json      # ubiquitin vs SUMO, full result as JSON
.venv/bin/proteinlens classify 8TIM                 # nearest CATH folds
.venv/bin/proteinlens report 1UBQ 1WM3              # verified plain-English report
.venv/bin/proteinlens info P69905                   # AlphaFold model by UniProt ID
```

Inputs can be PDB IDs (`4HHB`, `4HHB:B`, `4HHB:B:1-100`), UniProt accessions (AlphaFold DB
models), CATH domain IDs (`1cukA01`) or local `.pdb`/`.cif` files. Downloads are cached in
`data/cache/`.

Optional LLM narration: copy `.env.example`, set `ANTHROPIC_API_KEY` (and optionally
`PROTEINLENS_MODEL`, default `claude-opus-5-5`). Reports are still checked; failures fall back to
the template.

Reproduce the benchmark: `make benchmark` (≈ 8 minutes on 4 CPUs; resumable) or
`make benchmark-smoke` (40 domains, ≈ 10 s).

## How it works

| Module | What it does |
|---|---|
| `proteinlens/data.py` | Fetches and caches RCSB PDB entries, AlphaFold DB models and CATH domains; selects a chain/range; strips water and ligands; keeps one copy of alternate-location atoms; detects chain breaks; pulls metadata (title, organism, resolution, ligands, pLDDT). |
| `proteinlens/align.py` | TM-align (via `tmtools`): TM-score normalized by each protein, RMSD, aligned length, identity, rotation/translation, per-residue distances and the most deviating regions. Results are exactly symmetric under swapping. Also the Smith-Waterman sequence baseline. |
| `proteinlens/classify.py` | 1-NN fold classification against the bundled 288-domain CATH reference set. |
| `proteinlens/report.py` | Builds the facts object, the template report and the generate → check → regenerate → fallback loop. |
| `proteinlens/verify.py` | The faithfulness checker. |
| `proteinlens/llm.py` | Claude client (official `anthropic` SDK, structured JSON output), mock client and response cache. |
| `proteinlens/experiment.py` | Grounded vs ungrounded experiment. |
| `app/streamlit_app.py` | Compare / Classify / Benchmarks pages, 3D superposition with 3Dmol.js. |
| `benchmarks/` | Dataset builder, all-vs-all runner, analysis, results and figures ([details](benchmarks/README.md)). |

The faithfulness checker is tested on deliberately corrupted reports: wrong TM-score, wrong RMSD,
wrong counts, an invented derived number, a flipped or negated fold verdict, a PDB ID or CATH code
from memory, a fold name from memory, and functional or evolutionary claims (`tests/test_verify.py`).

More: [plain-language explainer](docs/UNDERSTANDING.md) · [paper-style write-up](docs/WRITEUP.md) ·
[deploying the app](docs/DEPLOY.md) · [design decisions](DECISIONS.md) · [progress log](PROGRESS.md).

| Classify | Benchmarks |
|---|---|
| ![classify page](docs/figures/classify.png) | ![benchmarks page](docs/figures/benchmarks.png) |

## Limitations

* **Small, clean benchmark.** 20 folds with ≥ 10 eligible members each; X-ray domains at ≤ 2.5 Å,
  50–300 residues, already cut into domains by CATH. Real inputs (multi-domain chains, NMR models,
  low-confidence predictions, folds absent from the reference set) are harder. 1-NN cannot predict
  a fold that is not in the reference set; the app warns when the best TM-score is ≤ 0.5.
* **Circularity with the answer key.** CATH uses structural comparison when classifying domains,
  so high agreement between a structure method and CATH is partly built in. The meaningful
  comparison is structure vs sequence on the same queries.
* **Residual redundancy.** 9 of 41,328 pairs exceed 40 % identity by our measure (identical
  residues ÷ shorter length), because CATH's S40 clustering uses different identity/overlap
  rules. They were kept and disclosed rather than removed after seeing results.
* **One pre-specified scoring choice.** Ranking by TM-score normalized by the query length was
  fixed in advance. Normalizing by the shorter protein would have scored higher at superfamily
  level (97.2 % vs 95.8 %); we report the pre-specified number.
* **The superfamily boundary is partly a curation decision.** Several "errors" are domains from
  the same kind of enzyme that CATH places in different superfamilies of one fold.
* **The faithfulness checker is narrow.** It verifies numbers, identifiers, fold names from a
  vocabulary, the fold verdict, and a list of functional and evolutionary keywords. It cannot prove
  that a fluent sentence without such tokens is correct or well chosen. It does not make the LLM
  "understand" proteins; the honest claim is that the structure is computed and the narration is
  checked against it.
* **The LLM experiment has not been run** (no API key at build time).
* **TM-score > 0.5 is a rule of thumb**, not a law: 50.6 % of same-fold/different-superfamily
  pairs and 91.2 % of same-superfamily pairs exceed it, versus 0.44 % of different-fold pairs.

## How this was built

ProteinLens was built with AI coding assistance (Claude Code) from a written project brief. All
results in this README, the write-up and the app were produced by the scripts in this repository
(`make benchmark`, `make eval-llm`) and are read from the committed result files. None were typed
in by hand. The decision log is in [`DECISIONS.md`](DECISIONS.md).

## Data, software and licenses

Code: MIT ([LICENSE](LICENSE)). ProteinLens depends on (but does not include) **tmtools**, a Python
binding to **TM-align**, which is licensed **GPL-3.0**. If you redistribute a bundle that includes
tmtools, its license applies to that bundle. Data: **RCSB PDB** (structures and metadata), **AlphaFold
DB** (CC BY 4.0), **CATH** v4.4 (CC BY 4.0), **UniProt** accessions. Software: TM-align, Biopython,
3Dmol.js, Streamlit, Plotly. Full references are in [CITATIONS.md](CITATIONS.md).
