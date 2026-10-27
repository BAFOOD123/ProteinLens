# Understanding ProteinLens

A plain-language guide to every idea in this project, written for a curious high-school student.
Numbers quoted here come from `benchmarks/results/summary.json`.

---

## 1. What is a protein structure?

A protein is a chain of amino acids (20 kinds, written as letters such as M, K, L). The order of
letters is the **sequence**. In water the chain folds into a specific 3D shape, and that shape
determines what the protein can do.

Scientists measure these shapes with X-ray crystallography, cryo-electron microscopy or NMR, and
deposit them in the **Protein Data Bank (PDB)**. Each entry has a 4-character ID such as `4HHB`
(human hemoglobin). A PDB file is a list of atoms with x, y, z coordinates in ångströms
(1 Å = 10⁻¹⁰ m). Entries can contain several **chains** (4HHB has chains A–D) as well as
non-protein molecules such as water, heme or drugs (the "HETATM" records).

**AlphaFold DB** holds *predicted* structures for millions of proteins, looked up by UniProt
accession (e.g. `P69905`, hemoglobin α). Its files store a confidence score called **pLDDT**
(0–100) for every residue.

ProteinLens represents each amino acid by its **C-alpha (Cα) atom**, the backbone carbon at the
centre of each residue. That gives one point per residue: a protein becomes a 3D curve.

**Domains.** Large proteins are often built from several compact, independently folding units
called domains. CATH classifies domains, not whole proteins, which is why benchmark IDs look like
`1cukA01` (PDB 1CUK, chain A, domain 01).

**Cleaning** (in `data.py`): keep one model and one chain, remove water and ligands, keep the
most-occupied copy of atoms recorded at two positions ("alternate locations"), and flag
**chain breaks**. Neighbouring Cα atoms are ~3.8 Å apart, so a gap above 4.5 Å means residues are
missing from the experimental model.

## 2. Comparing two structures: superposition, RMSD and TM-score

To compare two shapes we (1) decide which residue of A matches which residue of B (the
*alignment*) and (2) rotate and shift A to lie on top of B (the *superposition*).

**RMSD** (root-mean-square deviation) is the average distance between matched Cα atoms after
superposition:  RMSD = √( (1/N) Σ dᵢ² ). Lower is better. RMSD has two problems:

* It is dominated by the worst parts. One floppy loop that sits 15 Å away can make two nearly
  identical proteins look different.
* It depends on how many residues you chose to match. Matching only the 20 best residues gives a
  tiny RMSD that says little about the whole protein.

**TM-score** fixes both. For a protein of length L:

  TM-score = (1/L) Σᵢ 1 / (1 + (dᵢ/d₀)²),  with d₀ = 1.24·∛(L − 15) − 1.8

Every matched pair contributes between 0 and 1: close pairs count almost fully and far pairs count
almost nothing, so one bad loop cannot ruin the score. Dividing by L (not by the number matched)
rewards covering the *whole* protein. d₀ grows with protein size, so scores are comparable
across sizes. Rules of thumb: about **0.17** for unrelated proteins, **> 0.5** usually the same fold
(Xu & Zhang, 2010), 1.0 identical.

**Two TM-scores.** Since L is one protein's length, there are two versions. If a small domain (70
residues) sits inside a big protein (250 residues), the score normalized by the small one can be
high ("all of me is covered") while the one normalized by the big one is low ("most of me is not").
For classification ProteinLens uses the score normalized by the *query*.

**TM-align** (Zhang & Skolnick, 2005) searches for the alignment and superposition that maximize
the TM-score, without needing the sequences to be similar. ProteinLens calls it through `tmtools`.
"Aligned length" counts matched pairs closer than 5 Å, and the reported RMSD is over those pairs.

**Deviating regions** are runs of at least 3 consecutive matched residues that are all more than
3 Å apart after superposition. They are where the two shapes differ most.

## 3. Comparing sequences

Two sequences are aligned by sliding them against each other and allowing gaps. **BLOSUM62** is a
table scoring how often each amino-acid swap is seen in related proteins: similar amino acids
score positive, dissimilar ones negative. **Smith-Waterman** finds the best-scoring *local*
alignment. **Sequence identity** is the fraction of aligned positions with the same letter (here:
identical residues ÷ length of the shorter sequence).

Above about 30–40 % identity, two proteins are almost always related. Below about 25 % (the
"twilight zone") sequence similarity becomes unreliable. Related proteins can drift to 10–15 %
identity while keeping the same shape.

## 4. How CATH classifies proteins

CATH sorts domains into a four-level hierarchy:

| Level | Meaning | Example |
|---|---|---|
| **C**lass | mostly α-helix (1), mostly β-strand (2), mixed α/β (3) | 3 |
| **A**rchitecture | overall arrangement of secondary structures | 3.20 (α-β barrel) |
| **T**opology (= "fold" here) | the specific connectivity of the chain | 3.20.20 TIM barrel |
| **H**omologous superfamily | evidence of a common evolutionary origin | 3.20.20.70 aldolase class I |

CATH assigns domains using structure-comparison algorithms, sequence evidence and human curators.
**S40** is CATH's non-redundant list in which domains are clustered at 40 % sequence identity, so
near-copies do not inflate results. SCOP/SCOPe is a similar classification. ProteinLens planned
to use SCOPe, but its website was blocked from the build machine, so it uses CATH.

## 5. Why structure beats sequence for remote homologs

During evolution mutations change the sequence quickly, but the 3D shape is under strong pressure
to stay the same because the protein must still fold and work. After a long time two related
proteins may share only 10–15 % of their letters, too little for sequence search to notice, yet
still superimpose well. These are **remote homologs**.

Ubiquitin vs SUMO is a classic example: 14.1 % identity over the structurally aligned residues,
but a TM-score of 0.84 (`proteinlens compare 1UBQ 1WM3`). In the benchmark, on queries whose
same-superfamily relatives all share < 25 % identity (n = 133), structure found the right
superfamily 93.2 % of the time and sequence 54.9 %.

## 6. Nearest-neighbour classification and leave-one-out

**1-nearest-neighbour (1-NN):** to classify a protein, find the most similar protein with a known
label and copy its label. "Most similar" means highest TM-score for the structure method and
highest alignment score for the sequence method.

**Leave-one-out:** each of the 288 domains takes a turn as the "unknown" query, and it is not
allowed to match itself. Every domain is tested once without being used to answer its own question.

**Majority baseline:** always predict the most common label. In leave-one-out on a balanced set
it scores 0 %. Removing the query makes its own fold one member smaller than the others, so the
"most common" fold is always a different one. The more intuitive chance level is the largest
fold's share of the set, 5.2 %. Both are reported.

**Why no training or tuning?** 1-NN with a fixed score has nothing to learn, and the one choice
(which TM-score normalization to use) was fixed before results were seen. That avoids "tuning on
the test set", which would make the numbers look better than they really are.

### A harder test, and knowing when not to answer

In the main benchmark every query has 10–14 members of its own fold to match, which makes 1-NN's
job easy. A second **hard set** (600 domains, 120 folds, only 5 per fold) was designed and committed
*before* it was run, so its result could not be steered. Fold accuracy dropped from 99.7 % to
**88.0 %**, which is still far above sequence (65.5 %).

Looking at the errors showed a **size bias**. The TM-score used for ranking is divided by the
*query's* length, so it measures "how much of my protein is covered". A small domain can be
covered by *part* of a much bigger protein, and in 95.8 % of the errors the wrong neighbour was
longer than the query. A score that requires covering *both* proteins would have done better on
this set, but not on the main set. Switching after seeing the results would be test-set tuning,
so the headline keeps the original choice.

**Abstention** means saying "I'm not sure" instead of guessing. If ProteinLens answers only when its
best match has TM-score > 0.5, it answers 90.0 % of hard-set queries and gets 93.2 % of those right.
On the 10 % it refuses, a forced guess would have been right only 41.7 % of the time. This trade of
coverage for accuracy is how a tool like this should behave on real, unknown proteins.

## 7. Bootstrap confidence intervals

288 queries is a sample, and a different sample would give a slightly different accuracy. The
**bootstrap** estimates how much: draw 288 queries *with replacement* from the 288 (some appear
twice, some not at all), recompute the accuracy, repeat 2,000 times, and take the middle 95 % of
those accuracies (2.5th to 97.5th percentile). That range is the **95 % confidence interval**.
Example: superfamily accuracy 95.8 % [93.4, 97.9].

To compare two methods we use a **paired** bootstrap: each resample uses the same queries for
both methods, and we record the *difference*. If the whole interval is above zero, the advantage
is unlikely to be luck. On remote homologs, structure minus sequence = +38.4 points [29.3, 47.4].

## 8. The faithfulness checker

**What it does.** ProteinLens builds a **facts object**, a JSON dictionary of everything computed
(TM-scores, RMSD, lengths, identity, deviating regions, metadata, nearest folds). An LLM is asked
to write a short report using only those facts. The checker then reads the report and:

1. finds every number and checks that it matches a fact, allowing for rounding (0.90 for 0.904 is
   fine; 0.95 is not);
2. checks every identifier (PDB IDs, domain IDs) and CATH code against the facts;
3. flags fold names ("TIM barrel", "Rossmann" …) that are not in the facts;
4. finds sentences claiming "same fold" or "different folds", handles nearby negation ("do not
   share"), and checks them against the computed verdict;
5. flags words about function or evolution ("binds", "oxygen", "ancestor" …) that the facts do
   not contain, because nothing computed supports them.

If a report fails, the LLM gets the list of problems and tries again (up to 3 times). After that,
ProteinLens uses a fixed-sentence **template report** that is correct by construction. The tests
feed the checker deliberately corrupted reports to prove it catches each kind of error.

**What it does NOT guarantee.**
* A sentence with no numbers, IDs or listed keywords is not checked for meaning. "The proteins look
  quite different" when TM = 0.9 would slip through unless it uses the fold-verdict phrases.
* It checks against the facts, not against reality. If the alignment or metadata were wrong, a
  "verified" report would faithfully repeat the error.
* Its keyword lists are finite, so outside knowledge phrased in unlisted words can pass.
* Passing it does not mean the LLM understands proteins. The honest claim is narrower: *the
  structure is computed; the narration is checked against the computation.*

## 9. The grounded vs ungrounded experiment

60 pairs of CATH domains: 20 from the same superfamily, 20 from the same fold but different
superfamilies (the hard case), and 20 from different folds. The LLM answers "same fold?" twice:
(A) seeing only IDs and names, so it must rely on memory; (B) also seeing ProteinLens's facts.
Accuracy is scored against CATH and the explanations go through the checker. The live run is
**PENDING** because no API key was available. Already computed: the simple rule "same fold if
TM > 0.5" gets 83.3 %. It is perfect on easy pairs and 50 % on same-fold/different-superfamily
pairs, so a grounded LLM must beat 83.3 % to show it adds anything beyond the rule.

## 10. The limitations, one by one

* **Clean benchmarks.** The main set has 20 folds with many members each. With few examples per
  fold (the hard set) accuracy falls to 88.0 %, and real-world inputs are messier still.
* **Size bias.** Ranking by the query-normalized TM-score lets small domains match pieces of big
  proteins, which causes most hard-set errors.
* **Circularity.** CATH uses structural comparison when building its classification, so a
  structure method agreeing with CATH is partly built in. The fair comparison is structure vs
  sequence on the same queries.
* **Unknown folds.** 1-NN can only answer with folds in its 288-domain reference set. A protein
  from another fold gets a wrong answer, which is why the app warns when the best TM-score ≤ 0.5.
* **Residual redundancy.** 9 pairs have ≥ 40 % identity by our measure (CATH measures identity
  differently). They were kept and disclosed, because removing them after seeing results would
  be cherry-picking.
* **A pre-registered choice that was not the best.** Normalizing TM-score by the shorter protein
  would have scored 97.2 % instead of 95.8 % at superfamily level. We report 95.8 %.
* **Superfamily labels are partly judgement calls.** Some "errors" are similar enzymes that CATH
  places in different superfamilies.
* **TM > 0.5 is a rule of thumb.** About half of same-fold/different-superfamily pairs fall below it.
* **The checker is narrow** (section 8), and **the LLM experiment has not been run yet**.

---

## 11. Practice interview questions (with answers)

1. **What does ProteinLens do, in one sentence?**
   It computes how similar two protein structures are, predicts a structure's fold by finding the
   most similar classified domain, measures how accurate that is against CATH, and writes reports
   whose numbers are automatically checked against the computation.

2. **Why not just ask ChatGPT or Claude whether two proteins have the same fold?**
   An LLM answers from memory: it may not know a protein, can mix up IDs, and cannot measure the
   actual coordinates. ProteinLens computes the answer and only lets the LLM describe computed
   facts, checking every number.

3. **What is the difference between RMSD and TM-score, and why prefer TM-score for folds?**
   RMSD averages distances, so it is dominated by the worst-matching parts and depends on how many
   residues are matched. TM-score weights each pair between 0 and 1, so far pairs barely count, and
   it divides by the full length, rewarding coverage. It is also size-adjusted, with known
   thresholds (≈ 0.17 random, > 0.5 same fold).

4. **Why are there two TM-scores for one comparison?**
   The score is divided by one protein's length. For proteins of different sizes, "how much of A is
   covered" and "how much of B is covered" differ.

5. **What are CATH's levels, and which ones did you predict?**
   Class, Architecture, Topology, Homologous superfamily. I predicted topology ("fold") and
   superfamily.

6. **How did you make sure the benchmark is not trivially easy?**
   I used CATH's S40 set (clustered at 40 % sequence identity), used leave-one-out so no protein
   matches itself, fixed the method before seeing results, compared against sequence and majority
   baselines, and analyzed a remote-homolog subset. I also reported the remaining caveats
   (circularity with CATH, 9 high-identity pairs).

7. **What is a remote homolog and why does structure help?**
   An evolutionary relative whose sequence has diverged so much (< 25 % identity here) that
   sequence comparison struggles. Shape changes much more slowly than sequence. On those 133
   queries, structure got 93.2 % of superfamilies right and sequence 54.9 %.

8. **Your majority baseline is 0 %. Isn't that a bug?**
   No. In leave-one-out on a nearly balanced set, removing the query always makes another class the
   majority, so it is always wrong. The fairer chance level is the largest class's share, 5.2 %.

9. **What is a bootstrap confidence interval?**
   Resample the queries with replacement many times (2,000), recompute accuracy each time, and take
   the middle 95 % of results. It shows how much the accuracy could vary with a different sample.

10. **Your fold accuracy is 99.7 %. Should I believe it?**
    It is real for that benchmark, but the benchmark is easy: 20 well-populated folds, clean single
    domains, and CATH itself uses structural comparison. On a harder, pre-specified set with 120
    folds and only 5 examples each, accuracy is 88.0 %, still well above sequence (65.5 %). The most
    meaningful result is the large gap over sequence on remote homologs in both sets.

11. **Describe a failure and what you learned from it.**
    1f06A02 was put in the wrong fold. Its best TM-score to its own fold was only 0.385, so it did not
    resemble its fold-mates. Several superfamily "errors" are domains from the same kind of enzyme
    (tRNA synthetases) that CATH places in different superfamilies. Global shape cannot always
    reproduce curated superfamily boundaries.

12. **How does the faithfulness checker work?**
    It extracts numbers, IDs, CATH codes, fold names, same/different-fold claims and function words
    from the text and checks each against the facts object, allowing for rounding. Failing reports
    are regenerated with the list of problems, then replaced by a template.

13. **What can the checker NOT catch?**
    Wrong claims made without numbers, IDs or listed keywords; outside knowledge in unlisted words;
    and errors already in the facts. It verifies consistency with the computation, not truth or
    understanding.

14. **Why is the LLM experiment marked PENDING?**
    No API key was available when the project was built, and I did not want to report invented
    numbers. The pipeline is implemented and tested with a mock model. Running `make eval-llm` with a
    key produces the real results, and a cost estimate is printed first.

15. **If you had more time, what would you do next?**
    Run the LLM experiment; test on harder data (more folds, multi-domain chains, AlphaFold models
    with low pLDDT); add a calibrated "I don't know" option when the best TM-score is low; and
    compare the checker's verdicts with human judgements of report accuracy.
