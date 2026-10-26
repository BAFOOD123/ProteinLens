# DECISIONS

One line each: decision — reason.

- **CATH v4.4 (latest-release) instead of SCOPe** — `scop.berkeley.edu` returns 403 from this environment; CATH files download fine.
- **"Fold" = CATH topology (C.A.T); "superfamily" = CATH homologous superfamily (C.A.T.H)** — standard CATH mapping to the SCOP-style levels in the brief.
- **Python ≥ 3.12** — pinned numpy 2.5 / scipy 1.18 require 3.12; CI tests 3.12 and 3.13.
- **Exact pins in `requirements.txt`, loose ranges in `pyproject.toml`** — reproducible installs via `make setup` while keeping the package installable elsewhere.
- **`tmtools` (Python bindings to TM-align) for structural alignment** — installed cleanly from PyPI; it is the reference TM-align algorithm, run on Cα atoms.
