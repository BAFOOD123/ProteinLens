# Citations and attribution

ProteinLens computes nothing new about proteins by itself. It stands on these data resources and
tools; please cite them if you use results from ProteinLens.

## Data

* **RCSB Protein Data Bank** — experimental structures and entry metadata
  (`files.rcsb.org`, `data.rcsb.org`).
  Berman HM, Westbrook J, Feng Z, et al. *The Protein Data Bank.* Nucleic Acids Research
  28:235–242 (2000).
* **AlphaFold Protein Structure Database** — predicted structures (`alphafold.ebi.ac.uk`),
  licensed **CC BY 4.0**.
  Varadi M, et al. *AlphaFold Protein Structure Database: massively expanding the structural
  coverage of protein-sequence space with high-accuracy models.* Nucleic Acids Research
  50:D439–D444 (2022).
  Jumper J, et al. *Highly accurate protein structure prediction with AlphaFold.* Nature
  596:583–589 (2021).
* **CATH** v4.4 (release files dated 14 Jan 2025, `download.cathdb.info`; domain coordinates from
  `www.cathdb.info`), licensed **CC BY 4.0**. Used as the ground-truth answer key.
  Sillitoe I, et al. *CATH: increased structural coverage of functional space.* Nucleic Acids
  Research 49:D266–D273 (2021).
* **UniProt** accessions identify AlphaFold models.
  The UniProt Consortium. *UniProt: the Universal Protein Knowledgebase in 2023.* Nucleic Acids
  Research 51:D523–D531 (2023).

## Methods and software

* **TM-align** — structural alignment.
  Zhang Y, Skolnick J. *TM-align: a protein structure alignment algorithm based on the TM-score.*
  Nucleic Acids Research 33:2302–2309 (2005).
* **TM-score** — Zhang Y, Skolnick J. *Scoring function for automated assessment of protein
  structure template quality.* Proteins 57:702–710 (2004).
* **TM-score 0.5 threshold** — Xu J, Zhang Y. *How significant is a protein structure similarity
  with TM-score = 0.5?* Bioinformatics 26:889–895 (2010).
* **tmtools** — Python bindings to TM-align by Joris Vankerschaver,
  <https://github.com/jvkersch/tmtools> (**GPL-3.0**; ProteinLens imports it as a dependency and
  does not redistribute it).
* **Biopython** — structure parsing and pairwise sequence alignment.
  Cock PJA, et al. *Biopython: freely available Python tools for computational molecular biology
  and bioinformatics.* Bioinformatics 25:1422–1423 (2009).
* **Smith-Waterman local alignment** — Smith TF, Waterman MS. *Identification of common molecular
  subsequences.* Journal of Molecular Biology 147:195–197 (1981).
* **BLOSUM62** — Henikoff S, Henikoff JG. *Amino acid substitution matrices from protein blocks.*
  PNAS 89:10915–10919 (1992).
* **Bootstrap confidence intervals** — Efron B. *Bootstrap methods: another look at the
  jackknife.* Annals of Statistics 7:1–26 (1979).
* **3Dmol.js** — 3D viewer. Rego N, Koes D. *3Dmol.js: molecular visualization with WebGL.*
  Bioinformatics 31:1322–1324 (2015).
* **Streamlit**, **Plotly**, **NumPy**, **pandas**, **SciPy**, **Matplotlib**, **Typer** — app,
  charts and numerics.
* **Anthropic Claude** (official `anthropic` Python SDK) — optional report narration and the LLM
  experiment.
