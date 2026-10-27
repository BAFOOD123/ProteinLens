# Failure analysis — hard set (600 domains, 120 folds, 5 per fold)

Sources: `failures.csv`, `per_query.csv` and `summary.json` in this folder, all written by
`benchmarks/analyze.py --set hard`. Protein descriptions are RCSB entry titles.

Structure 1-NN made **72 fold-level errors** out of 600 queries (88.0 % accuracy). The errors are
not random. They follow one dominant pattern.

## Pattern 1: a small query "covered" by part of a much larger neighbour

| | Errors (n = 72) | Correct (n = 528) |
|---|---|---|
| Nearest neighbour longer than the query | **95.8 %** | 56.1 % |
| Median neighbour length ÷ query length | **1.9** | 1.02 |

(`summary.json` → `structure_fold_errors_vs_neighbour_size`.) The primary score is the TM-score
normalized by the *query's* length, so it asks "how much of my protein can be superimposed on the
reference?". A short domain, especially a small bundle of helices, can often be superimposed onto a
*piece* of a larger protein. It then scores well even though the two folds differ as wholes. The
most confident errors look exactly like this:

| Query (length) | Nearest neighbour (length) | TM (query-normalized) | Best TM to own fold |
|---|---|---|---|
| 2fcwA00 (106), 1.20.81 "Receptor-associated Protein"; entry: LDL-receptor modules with receptor-associated protein | 2q12A00 (242), 1.20.1270; entry: BAR domain of APPL1 | 0.767 | 0.478 |
| 5figA00 (100), 1.20.1270; entry: copper storage protein 3 | 3d19B00 (262), 1.20.1260 "Ferritin"; entry: conserved metalloprotein | 0.691 | 0.559 |
| 7s0rB01 (76), 1.20.81; entry: factor-H-binding fragment of a streptococcal protein | 4nqiD00 (232), 1.20.1270; entry: I-BAR domain | 0.656 | 0.532 |

Accuracy rises with query length (`structure_fold_accuracy_by_query_length_tercile`): **78.8 %**
for 51–108 residues, **89.9 %** for 109–165, **95.5 %** for 166–300. Mainly-α domains are hardest
(`structure_fold_accuracy_by_class`): **74.1 %** (class 1, n = 135), versus 92.8 % for mainly-β
(n = 125) and 91.8 % for α/β (n = 340). 35 of the 72 errors are in class 1.

### What this means for the choice of score

The sensitivity analysis (`sensitivity_tm_normalization`) shows that normalizing by the **longer**
of the two proteins, which forces the match to cover *both*, would have given **92.0 %** fold
accuracy here instead of 88.0 %. On the main set, a different normalization (by the shorter
protein) scored best at superfamily level. **No single normalization wins on both sets**, and we
keep the pre-specified one for the headline numbers. A principled follow-up, to be tested on new data
rather than tuned here, is to require the match to cover both proteins (for example, rank by the
smaller of the two TM-scores, or penalize large length differences).

## Pattern 2: the right answer is simply not similar enough

For 61 of the 72 errors, the query's best TM-score to *any* member of its own fold is ≤ 0.5. With
only 4 fold-mates per query, some folds are represented by members that are not globally
similar to each other. This is the situation the main set hid by providing 10–14 fold-mates.

## Abstention helps

Answering only when the best TM-score is above 0.5 (threshold from the literature, not tuned;
`abstain_below_tm_0_5`):

* coverage **90.0 %** (540 of 600 queries answered),
* fold accuracy on answered queries **93.2 %** [90.9, 95.2],
* fold accuracy on the 60 abstained queries **41.7 %** [30.0, 55.0].

35 of the 72 errors would have been withheld as "low confidence".

## Superfamily level

Superfamily accuracy is 78.3 % over all 600 queries, but 107 of the 236 superfamilies have only one
member here, so 107 queries *cannot* be right. Over the 493 queries that have a same-superfamily
partner, structure scores **95.3 %** [93.5, 97.0] vs sequence 78.7 % [75.0, 82.6]
(`queries_with_superfamily_partner`).

## Why the remote-homolog subset scores higher than "all queries"

The remote-homolog subset (n = 234) requires a same-superfamily partner in the set. Having such a
relative makes a query easier to place, so structure fold accuracy is 95.7 % there vs 88.0 % overall.
The informative comparison inside that subset is structure (91.5 % superfamily) vs sequence (56.0 %).
