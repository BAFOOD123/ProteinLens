# ProteinLens — common tasks. Run `make setup` first.
PY ?= .venv/bin/python
PIP ?= .venv/bin/pip

.PHONY: setup test test-all lint format dataset benchmark benchmark-smoke benchmark-hard eval-llm app clean

setup:            ## Create .venv and install pinned dependencies + the package
	python3 -m venv .venv
	$(PIP) install -q --upgrade pip
	$(PIP) install -q -r requirements.txt
	$(PIP) install -q -e . --no-deps

test:             ## Offline unit tests (uses bundled fixtures, no network)
	$(PY) -m pytest -q -m "not network"

test-all:         ## All tests, including ones that download from RCSB / AlphaFold
	$(PY) -m pytest -q

lint:
	$(PY) -m ruff check .
	$(PY) -m ruff format --check .

format:
	$(PY) -m ruff format .
	$(PY) -m ruff check --fix .

dataset:          ## Build the CATH benchmark set (downloads, cached + resumable)
	$(PY) benchmarks/build_dataset.py

benchmark: dataset ## Full benchmark: all-vs-all TM-align + sequence baseline + analysis
	$(PY) benchmarks/run_benchmark.py
	$(PY) benchmarks/analyze.py

benchmark-smoke:  ## Quick end-to-end check on a 40-domain subset (writes to benchmarks/smoke/)
	$(PY) benchmarks/run_benchmark.py --smoke
	$(PY) benchmarks/analyze.py --smoke

benchmark-hard:   ## Harder second benchmark: 600 domains / 120 folds (resumable; re-run if it stops)
	$(PY) benchmarks/build_dataset.py --set hard
	$(PY) benchmarks/run_benchmark.py --set hard
	$(PY) benchmarks/analyze.py --set hard

eval-llm:         ## Grounded vs ungrounded LLM experiment (needs ANTHROPIC_API_KEY)
	$(PY) benchmarks/llm_experiment.py

app:
	.venv/bin/streamlit run app/streamlit_app.py

clean:
	rm -rf .pytest_cache .ruff_cache build *.egg-info
