# Deploying the ProteinLens app (free tiers)

The app needs no API key; with one, reports are narrated by Claude and still checked.
Both options below install `requirements.txt` (pinned versions).

## Streamlit Community Cloud
1. Push this repository to GitHub (public, or private with Streamlit access).
2. Go to <https://share.streamlit.io>, **Create app** → pick the repo and branch, main file
   `app/streamlit_app.py`, Python 3.12.
3. (Optional) **Advanced settings → Secrets**: add
   ```toml
   ANTHROPIC_API_KEY = "sk-ant-..."
   PROTEINLENS_MODEL = "claude-opus-5-5"
   ```
   The app copies `ANTHROPIC_API_KEY` from `st.secrets` into the environment.
4. Deploy. The first comparison of a new PDB ID downloads it from RCSB; the four example pairs
   work immediately from the bundled fixtures.

## Hugging Face Spaces
1. Create a new Space → SDK **Docker** (or **Streamlit** if offered) → link or push this repo.
2. For a Docker Space add a `Dockerfile`:
   ```dockerfile
   FROM python:3.12-slim
   WORKDIR /app
   COPY . .
   RUN pip install --no-cache-dir -r requirements.txt && pip install --no-cache-dir -e . --no-deps
   EXPOSE 7860
   CMD ["streamlit", "run", "app/streamlit_app.py", "--server.port=7860", "--server.address=0.0.0.0"]
   ```
3. (Optional) Settings → **Variables and secrets** → add `ANTHROPIC_API_KEY` as a secret.

## Notes
* The 3D viewer loads 3Dmol.js from the jsDelivr CDN in the visitor's browser. To serve it
  locally instead, run `python -c "from proteinlens.viewer import vendor_3dmol; vendor_3dmol()"`
  once (it downloads the npm package into `data/cache/3dmol/`).
* Downloads are cached in `data/cache/` (ephemeral on free tiers; that only costs re-downloads).
