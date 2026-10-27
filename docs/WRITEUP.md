# ProteinLens: computed structural comparison, measured fold classification, and verified LLM narration

*All numbers below are read from `benchmarks/results/` and were produced by the scripts in this
repository.*

## 1. Introduction

Large language models (LLMs) can describe proteins fluently, but what they say comes from memory
and is not checked against the molecule in front of them. For structural questions such as "do
these two proteins share a fold?" a reliable answer requires *computing* a structural comparison.
ProteinLens separates the two jobs. Structure comparison and fold classification are done by
established algorithms and measured against a ground-truth classification. An LLM is used only to
narrate a structured set of computed facts, and an automatic checker verifies that narration. We
ask three questions: (1) how accurately does structural nearest-neighbour search recover CATH fold
and superfamily labels; (2) how much better is it than sequence similarity, especially for
remote homologs; and (3) can LLM-written reports be constrained to the computed facts?

## 2. Methods

**Data.** The ground truth is CATH v4.4 (SCOPe's server was unreachable from the build
environment). Starting from the CATH S40 non-redundant domain list, we kept X-ray domains at
≤ 2.5 Å resolution, 50–300 residues long, in classes 1–3. A fold (CATH topology) was eligible if it
had ≥ 10 such domains and ≥ 2 superfamilies with ≥ 2 domains each (146 folds). We chose 5 textbook
"anchor" folds in advance (globin-like, ubiquitin-like roll, TIM barrel, immunoglobulin-like,
Rossmann) plus 5 random eligible folds per class (seed 42). We then sampled up to 15 domains per
fold, round-robin across up to 4 random superfamilies. The result is **288 domains, 20 folds and
60 superfamilies**, and every superfamily has ≥ 2 members.

**Structure comparison.** TM-align (via `tmtools`) aligns C-α atoms and reports the TM-score
normalized by each protein's length, the RMSD and aligned length (pairs < 5 Å), and the
superposition. We also derive per-residue distances and the top deviating regions (≥ 3 consecutive
aligned pairs > 3 Å apart). Comparisons are made exactly symmetric by running TM-align in a
canonical order.

**Classification.** All 41,328 pairs were aligned (357 s on 4 CPUs). Each domain was classified by
leave-one-out 1-nearest-neighbour using the TM-score normalized by the query length, which was
fixed in advance. Baselines were 1-NN on the Smith-Waterman score (BLOSUM62, gap 11/1) and the
leave-one-out majority class. We report accuracy with 95 % percentile-bootstrap CIs (2,000
resamples over queries; paired for differences). The **remote-homolog subset** is the set of
queries that have at least one same-superfamily domain in the set but < 25 % sequence identity to
all of them (n = 133).

**Grounded reports.** A facts object holds every computed number and the metadata. The LLM
(Claude via the official SDK, JSON-schema output) is told to use only these facts. A checker
extracts every number (allowing for printed rounding), identifier, CATH code, known fold name,
fold-verdict sentence (with local negation handling) and functional/evolutionary keyword, and
verifies each against the facts. Failing drafts are regenerated up to 3 times, with the checker's
complaints as feedback, and then replaced by a deterministic template report, which always passes.

**LLM experiment.** We drew 60 domain pairs (seed 7): 20 from the same superfamily, 20 from the
same fold but different superfamilies, and 20 from different folds. Condition A gives the model
only domain/PDB IDs and entry titles; condition B adds the computed facts. Answers are scored
against CATH, and explanations are run through the checker. A no-LLM rule (same fold if the larger
TM-score > 0.5) is scored on the same pairs.

## 3. Results

| | Fold | Superfamily |
|---|---|---|
| Structure 1-NN, all (n = 288) | 99.7 % [99.0, 100.0] | 95.8 % [93.4, 97.9] |
| Sequence 1-NN, all | 80.9 % [76.0, 85.4] | 78.1 % [73.3, 83.0] |
| Structure 1-NN, remote homologs (n = 133) | 99.2 % [97.7, 100.0] | 93.2 % [88.7, 97.0] |
| Sequence 1-NN, remote homologs | 60.2 % [51.9, 69.2] | 54.9 % [45.9, 63.2] |

The leave-one-out majority baseline scored 0 % because the classes are nearly balanced (the largest
fold is 5.2 % of the set). Structure beat sequence by 17.7 points [12.8, 22.6] at superfamily level
overall and by 38.4 points [29.3, 47.4] on remote homologs. Binned by each query's highest
identity to a same-superfamily domain, sequence accuracy fell from 100 % (≥ 30 %, n = 98) to 30.0 %
(< 15 %, n = 50), while structure stayed between 90.9 % and 100 %. Pairwise TM-scores separate the
relations imperfectly: 91.2 % of same-superfamily pairs, 50.6 % of same-fold/different-superfamily
pairs and 0.44 % of different-fold pairs exceed 0.5.

There was one fold error. 1f06A02 (a diaminopimelate dehydrogenase domain) has a best TM-score of
only 0.385 to its own fold and was assigned to a neighbouring topology of the same architecture
with TM 0.423. All 12 superfamily errors chose a neighbour from a related superfamily, 11 of them in
the correct fold. Several involve domains whose PDB titles name the same kind of enzyme
(leucyl- or isoleucyl-tRNA synthetases) but which CATH places in different superfamilies.

The checker passed every template report tested. It rejected all 12 kinds of deliberately
corrupted report in the test suite: wrong values, an invented derived number, a flipped or negated
verdict, an identifier, CATH code or fold name from memory, and functional or evolutionary claims.
**The live LLM experiment is PENDING** (no API key at build time). On its 60 pairs the TM > 0.5 rule
alone scores 83.3 % [73.3, 93.3]: 100 % on same-superfamily and different-fold pairs, and 50 % on
same-fold/different-superfamily pairs. This rule is the bar a grounded LLM must beat to show it
adds anything.

## 4. Limitations

The benchmark is small and clean (20 well-populated folds, high-resolution single domains), so it
overstates performance on arbitrary inputs. CATH is partly derived from structural comparison,
so agreement with it is partly circular; the structure-vs-sequence contrast is the informative
result. Nine pairs exceed 40 % identity by our measure despite S40 filtering; we disclose them
rather than remove them post hoc. The pre-specified TM normalization was not the best in a
sensitivity analysis (shorter-length normalization: 97.2 % superfamily); we kept the
pre-specified choice. The checker verifies only tokens it can recognise. It cannot certify the
correctness of free-form wording, and it does not imply that the LLM understands protein
structure.

## 5. Future work

Run the PENDING LLM experiment. Extend the benchmark to more folds, multi-domain chains and
AlphaFold models stratified by pLDDT. Add a calibrated confidence score, e.g. abstain when the top
TM-score is ≤ 0.5. Compare against other structural search tools. Evaluate the checker against
human judgements of report faithfulness.
