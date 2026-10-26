# Failure analysis

Source files: `failures.csv` (fold level) and `superfamily_failures.csv` (superfamily level),
both written by `benchmarks/analyze.py`. Protein descriptions below are the RCSB entry titles
(fetched with `proteinlens.data.fetch_pdb_metadata`); superfamily names are from CATH's
`cath-names.txt` (several superfamilies have no name there).

Structure 1-NN (TM-score normalized by the query) made **1 fold-level error** and **12
superfamily-level errors** out of 288 queries.

## The one fold-level error

| Query | True fold | Nearest neighbour | Predicted fold | TM to neighbour | Best TM to own fold |
|---|---|---|---|---|---|
| 1f06A02 (138 res) — domain 2 of *C. glutamicum* diaminopimelate dehydrogenase | 3.30.360 "Dihydrodipicolinate Reductase; domain 2" | 6fgjA01 — *S. aureus* small alarmone synthetase 2 | 3.30.460 "Beta Polymerase; domain 2" | 0.423 | 0.385 |

This query does not look much like *anything* in the set: its best TM-score to any of its 14
fold-mates is 0.385, below the 0.5 "same fold" rule of thumb, and the wrong neighbour wins
with an equally weak 0.423. Both folds are in CATH architecture 3.30 (two-layer α/β sandwich),
so the confusion is between neighbouring topologies of the same architecture. A 1-NN classifier
can only be as good as the closest example it is given; here CATH's topology assignment is not
reflected by a high global TM-score to the other members we sampled. A practical safeguard
(not used in the benchmark numbers) is to flag predictions whose top TM-score is below 0.5
as low-confidence; this one would have been flagged.

## Superfamily-level errors (12)

**Every superfamily error picked a neighbour from the same fold** (11 of 12 are in the correct
fold; the 12th is the fold error above). This is exactly the hardest case: different
superfamilies that share a topology. The pair statistics show why: 50.6 % of
"same fold, different superfamily" pairs have TM-score > 0.5 (`pair_stats.csv`), so a global
shape score alone often cannot tell such superfamilies apart. Most of these wrong neighbours
are confident (TM 0.54–0.85), so this is not just noise near a threshold.

| Group (fold) | Query → wrong neighbour (TM) | What the entries are (RCSB titles) |
|---|---|---|
| 1.10.730 (two superfamilies, .10 and .20) | 1wkbA03 → 6lpfA01 (0.853); 6lpfA01 → 3zgzA04 (0.836); 1ileA03 → 1ffyA03 (0.809); 1ffyA03 → 1ileA03 (0.565) | domains of leucyl-tRNA synthetases (archaeal 1WKB, *E. coli* 3ZGZ, human 6LPF) and isoleucyl-tRNA synthetases (1ILE, 1FFY) |
| 3.20.20 TIM barrel | 2iqtA00 (class I aldolase, .70) → 3qxbA00 (putative xylose isomerase, .150) (0.581); 7lvlA01 (dihydrodipicolinate synthase, .70) → 3erpA01 (putative oxidoreductase, .100) (0.571) | TIM-barrel enzymes from different superfamilies |
| 3.40.91 restriction endonuclease-like | 2fokA03 (FokI, .30) → 2ewfA03 (nickase N.BspD6I, .50) (0.727); 1fzrA00 (T7 endonuclease I, .30) → 1knvB00 (Bse634I, .10) (0.459) | nucleases from different superfamilies |
| 2.40.40 Barwin-like | 1yleA02 (.20) → 3sumB00 (cerato-platanin, .10 "RlpA-like") (0.714); 1twfA03 (RNA polymerase II, .20) → 3sumB00 (0.636) | |
| 3.30.460 Beta polymerase domain 2 | 4at7B02 (NF90–NF45 complex, .10) → 7lt2A01 (cGAS-like receptor, .90) (0.544) | |
| 3.30.360 | 1f06A02 → 6fgjA01 (0.423) | the fold error above |

Observations:

* **The tRNA-synthetase cases are arguably label-boundary cases, not structural failures.**
  Domains whose PDB titles say they come from the same kind of enzyme (leucyl-tRNA synthetase:
  1WKB vs 6LPF; isoleucyl-tRNA synthetase: 1ILE vs 1FFY) sit in two different CATH superfamilies
  of fold 1.10.730 and superimpose with TM-scores of 0.81–0.85. We still count them as errors —
  the answer key is CATH — but they show that "superfamily" is a curated judgement that uses more
  evidence than global shape.
* **TIM barrels** are the classic example of one fold shared by many superfamilies; two of the
  15 TIM-barrel queries were assigned to the wrong TIM-barrel superfamily.
* The **sequence baseline got 3 of these 12 right** (1wkbA03, 1yleA02, 7lvlA01), so the two
  methods fail on partly different queries.
* Two wrong neighbours have TM < 0.5 (1fzrA00 at 0.459, 1f06A02 at 0.423): low-confidence
  predictions that a TM < 0.5 warning would flag.

## Why the headline numbers are high — and why that should be read with care

* **CATH itself is built partly from structure comparison.** Domains are assigned to topologies
  and superfamilies using structural similarity (plus sequence and curator judgement). A
  structure-based classifier agreeing with CATH is therefore partly expected; the benchmark
  measures agreement with CATH, not discovery of new biology.
* **Only 20 folds, each with ≥ 10 eligible members**, and every query has 10–14 fold-mates in
  the reference set. Real queries may come from folds with few or no known members, where 1-NN
  must fail.
* **Clean inputs**: X-ray domains at ≤ 2.5 Å, 50–300 residues, already cut into single domains by
  CATH. Multi-domain chains, NMR models or low-confidence predicted structures would be harder.
* **Residual sequence redundancy**: although the set is CATH S40, 9 of 41,328 pairs have ≥ 40 %
  identity by our measure (identical residues ÷ shorter length), because CATH clusters with its
  own identity and overlap rules. The extreme case is 1jkfA03 (60 residues) vs 1p1jA02 (115
  residues): 58 of the short domain's 60 residues are identical to the longer one (96.7 %) — it is essentially contained
  in the longer one. These pairs were kept (removing them after seeing results would be post-hoc
  selection); queries with such partners are excluded from the remote-homolog subset by its
  definition.
* **The pre-specified TM-score normalization was not the best one** at superfamily level:
  normalizing by the shorter protein gives 97.2 % vs 95.8 % (`summary.json`,
  `sensitivity_tm_normalization`). We report the pre-specified choice as the headline and do not
  switch after seeing the results.
