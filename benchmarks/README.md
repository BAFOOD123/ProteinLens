# Benchmarks

Everything here is produced by the scripts in this folder; nothing is typed in by hand.

| Command | What it does | Time (4 vCPU) |
|---|---|---|
| `make dataset` (`build_dataset.py`) | builds the CATH benchmark set | ~1 min first time, seconds when cached |
| `make benchmark` | dataset + all-vs-all alignment (`run_benchmark.py`) + analysis (`analyze.py`) | ~8 min |
| `make benchmark-smoke` | the same on a 40-domain subset → `benchmarks/smoke/` | ~10 s |
| `make benchmark-hard` | the harder second benchmark (below) → `results_hard/`, `figures_hard/` | ~30 min |
| `make eval-llm` (`llm_experiment.py`) | grounded vs ungrounded LLM experiment (needs an API key) | see `results/llm_experiment/` |

`run_benchmark.py` is resumable: finished chunks are saved under `data/cache/benchmark_chunks/`,
and if it hits its time limit (`--max-minutes`, default 25) it exits with code 3; run it again
to continue.

## Ground truth: CATH v4.4

SCOPe's download site (`scop.berkeley.edu`) was blocked from the build environment, so the answer
key is **CATH** (release `latest-release` = v4.4.0, files dated 14 Jan 2025, downloaded from
`download.cathdb.info`). CATH levels map onto the words used in this project as:

| CATH level | Code | Called here |
|---|---|---|
| Class | C | class (mainly α / mainly β / α-β) |
| Architecture | C.A | – |
| Topology | C.A.T | **fold** |
| Homologous superfamily | C.A.T.H | **superfamily** |

## Selection procedure (exact)

Parameters are in `build_dataset.py` (`CONFIG`) and saved with the chosen folds in
`data/selection.json`. Random seed: **42**.

1. Start from CATH's **S40 non-redundant domain list** (`cath-dataset-nonredundant-S40.list`,
   34,653 domains; no two share ≥ 40 % sequence identity).
2. Keep domains that are X-ray structures at **≤ 2.5 Å** resolution, **50–300 residues** long, in
   CATH classes **1, 2 or 3**.
3. A fold (topology) is eligible if it has **≥ 10** such domains **and ≥ 2 superfamilies with ≥ 2
   domains each** (so both levels can be tested). 146 folds qualify.
4. Choose **20 folds**: 5 **anchor folds** fixed in advance because they are textbook examples and
   make the demo interpretable — globin-like (1.10.490), ubiquitin-like roll (3.10.20), TIM barrel
   (3.20.20), immunoglobulin-like (2.60.40), Rossmann fold (3.40.50) — plus **5 random eligible folds
   from each of the three classes** (`numpy.random.default_rng(42)`).
5. In each fold choose up to **4 superfamilies** at random among those with ≥ 2 domains, then take
   domains **round-robin** across them (random order within each) until there are **15**.
6. Download each domain's coordinates from the CATH API
   (`https://www.cathdb.info/version/v4_4_0/api/rest/id/<domain>.pdb`) and parse its C-α atoms.

Result: **288 domains, 20 folds, 60 superfamilies** (`data/dataset.csv`). The total is below
20 × 15 = 300 because in some folds the randomly chosen superfamilies together had fewer than 15
eligible domains (e.g. Ig-like: 13, restriction endonuclease: 11); no domain failed to download. Every
superfamily in the set has at least 2 members. Rebuilding gives a byte-identical `dataset.csv`.

The C-α coordinates and sequences are committed in `data/reference_set.npz` (≈ 0.5 MB), so
classification and the app work without re-downloading anything.

## Methods

* **Structure**: TM-align (via `tmtools`) on C-α atoms for all 41,328 pairs. Classification is
  **leave-one-out 1-nearest-neighbour**: each domain is predicted to have the fold/superfamily of the
  *other* domain with the highest TM-score **normalized by the query's length** (chosen before any
  results were seen).
* **Sequence baseline**: Smith-Waterman local alignment (Biopython `PairwiseAligner`, BLOSUM62, gap
  open 11 / extend 1), 1-NN on the raw alignment score.
* **Majority baseline**: predict the most common label among the other domains.
* **Confidence intervals**: percentile bootstrap over queries, 2,000 resamples, seed 0. Differences
  between methods use a paired bootstrap (the same resampled queries for both).
* **Remote-homolog subset**: queries that have ≥ 1 other domain of their superfamily in the set but
  **< 25 % sequence identity to every one of them** (identity = identical residues in the local
  alignment ÷ length of the shorter sequence). CATH does not have SCOP-style families; because the set
  is S40 non-redundant, every same-superfamily pair is already in a different 40 % sequence cluster.
* No parameter was tuned on these results. There is nothing to tune in 1-NN with a fixed score;
  other TM-score normalizations are reported only as a secondary sensitivity analysis.

## Outputs

* `results/summary.json` — all headline numbers, CIs, subset sizes, sensitivity analysis
* `results/per_query.csv` — every query's nearest neighbours and predictions
* `results/per_fold.csv`, `results/pair_stats.csv`, `results/accuracy_vs_identity.csv`
* `results/failures.csv` + `results/FAILURE_ANALYSIS.md` — misclassifications, examined by hand
* `results/pairwise.npz` — the full 288 × 288 matrices (TM-score, RMSD, sequence scores)
* `figures/accuracy.png`, `figures/accuracy_vs_identity.png`, `figures/confusion_matrix.png`,
  `figures/tm_distributions.png`

## Second benchmark: the "hard" set (pre-specified before running)

The main set is easy in one important way: every query has 10–14 members of its own fold to match.
The hard set tests the realistic case of **many folds with few known examples each**. Its design
and analysis were written down and committed before any of its results were computed.

* Same quality filters as the main set (S40, X-ray ≤ 2.5 Å, 50–300 residues, classes 1–3).
* Eligible fold: ≥ 5 such domains (430 folds). **120 folds** drawn uniformly at random with seed
  **2026** (no anchor folds, no superfamily requirement).
* **5 domains per fold**, drawn uniformly at random within the fold (no superfamily balancing).
* Result: **600 domains, 120 folds (27 / 25 / 68 in classes 1 / 2 / 3), 236 superfamilies**, of which
  107 have a single member; 6 folds also occur in the main set. Each query has only 4 fold-mates.
* Same methods, same primary metric, same bootstrap, same remote-homolog definition as the main set.
* Superfamily accuracy is reported over all queries and over the queries that have at least one
  same-superfamily partner (a singleton superfamily can never be predicted correctly).
* Expected difficulty (stated in advance): fold accuracy below the main set's, because there are
  6× more classes and far fewer examples per class.

## Abstention analysis (both sets, pre-specified)

1-NN always answers, even when nothing in the reference set resembles the query. Using the
literature threshold (TM-score > 0.5 usually means the same fold; not tuned on these data), we also
report the **coverage** (fraction of queries whose best match has TM > 0.5) and the accuracy on the
answered and abstained queries (`summary.json` → `abstain_below_tm_0_5`).
