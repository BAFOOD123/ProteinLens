# DECISIONS

One line each: decision — reason.

- **CATH v4.4 (latest-release) instead of SCOPe** — `scop.berkeley.edu` returns 403 from this environment; CATH files download fine.
- **"Fold" = CATH topology (C.A.T); "superfamily" = CATH homologous superfamily (C.A.T.H)** — standard CATH mapping to the SCOP-style levels in the brief.
- **Python ≥ 3.12** — pinned numpy 2.5 / scipy 1.18 require 3.12; CI tests 3.12 and 3.13.
- **Exact pins in `requirements.txt`, loose ranges in `pyproject.toml`** — reproducible installs via `make setup` while keeping the package installable elsewhere.
- **`tmtools` (Python bindings to TM-align) for structural alignment** — installed cleanly from PyPI; it is the reference TM-align algorithm, run on Cα atoms.
- **Chain-break threshold 4.5 Å (not 4.2 Å)** — 1984-era 4HHB chain B has 4.22 Å Cα–Cα gaps without missing residues; one missing residue gives > 5 Å.
- **`compare(a, b)` runs TM-align in a canonical order (by length, sequence, name) and flips the result** — makes results exactly symmetric under swapping; raw TM-align differs only in the 3rd decimal.
- **Aligned length / RMSD / identity follow TM-align's definitions** (pairs < 5 Å). Per-residue distances use the TM-score superposition, which differs slightly from TM-align's reported (Kabsch) RMSD.
- **Deviating region = ≥ 3 consecutive aligned pairs all > 3 Å apart**, top 3 by mean distance; none reported if nothing exceeds 3 Å (stated explicitly in reports).
- **5 "anchor" folds (globin, ubiquitin-like, TIM barrel, Ig-like, Rossmann) fixed before any results + 15 random folds (seed 42)** — makes the demo interpretable; anchors had to pass the same eligibility filters.
- **288 instead of 300 domains** — kept the seeded selection rather than re-rolling when the chosen superfamilies ran out of members.
- **Primary classifier: TM-score normalized by the query length, pre-specified.** Other normalizations reported only as a sensitivity analysis (one of them scores higher at superfamily level; we did not switch).
- **Sequence baseline ranks by raw Smith-Waterman score (BLOSUM62, 11/1)**; sequence identity = identities ÷ shorter length.
- **Remote-homolog subset = queries with a same-superfamily partner but < 25 % identity to all of them** (CATH has no SCOP-style family level; the S40 set already separates 40 % clusters).
- **Majority baseline is leave-one-out and scores 0 %** on this near-balanced set; the largest-class share (5.2 %) is reported as the chance reference.
- **Kept the 9 pairs with ≥ 40 % identity** (CATH S40 uses different identity/overlap rules); disclosed instead of removed post hoc.
- **No hyperparameters tuned on the benchmark**, so no validation split was needed; 1-NN with a fixed score has nothing to tune.
- **Default LLM: `claude-opus-5-5` (env `PROTEINLENS_MODEL`); bulk M6 calls: `claude-haiku-5-5`** — model IDs/prices from the Claude API reference (cached 2026-10-06), per the brief's "cheapest capable model for bulk".
- **No server-side refusal fallback** — a refusal or API error falls back to the verified template report instead, which keeps the code simple.
- **Checker is deliberately strict**: any functional/evolutionary word (bind, oxygen, ancestor, homolog …) or fold name not present in the facts fails a report, since nothing computed supports it.
- **3Dmol.js embedded directly instead of py3Dmol/stmol** — one fewer dependency; the CDN was blocked in the build VM, so `vendor_3dmol()` can inline a copy fetched from the npm registry (used for screenshots).
- **tmtools is GPL-3.0; ProteinLens code stays MIT** — tmtools is a pip dependency, not redistributed; this is stated in README and CITATIONS.
- **App seeds its cache from the bundled fixtures** — the four example pairs work offline and on first load.
- **LLM experiment pairs exclude nearest-fold hits from the grounded facts** — those domains are in the reference set, so a fold lookup would leak the answer.
- **Extension: a harder second benchmark (600 domains / 120 random folds / 5 per fold, seed 2026)** — the main set gives every query 10–14 fold-mates; this tests the few-examples case. Design committed before running.
- **Abstention at TM 0.5 reported for both sets** — threshold from the literature (Xu & Zhang 2010), not tuned; shows how often 1-NN can answer confidently.
- **Score improvement via a proper dev/test split** — main + hard are now "development" sets; a fixed rule (`select_score.py`, best mean fold accuracy) picked normalization by the longer protein. A new held-out test set (seed 7777, no shared domains) decides whether it is adopted (pre-registered rule: adopt if test fold accuracy ≥ original).
- **Coverage-aware score adopted as the default ranking** — met the pre-registered rule on the held-out test set (+4.7 fold-accuracy points [2.7, 6.7]); main/hard headline tables keep the pre-specified original score.
- **Results website generated, not hand-written** (`site/build_site.py` → self-contained `site/index.html`, Pages workflow, freshness test) — keeps every number on the site tied to the result files. Hand-drawn SVG charts (no chart library); 3Dmol.js from jsDelivr for the 3D view.
