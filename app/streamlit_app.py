"""ProteinLens web app.

Run locally with ``make app`` (or ``streamlit run app/streamlit_app.py``).
Works without an API key: reports then use the verified template text.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from proteinlens.align import DEVIATION_CUTOFF, SAME_FOLD_TM, compare  # noqa: E402
from proteinlens.classify import classify, load_reference  # noqa: E402
from proteinlens.data import (  # noqa: E402
    ProteinLensError,
    Structure,
    cache_dir,
    load_structure,
    parse_structure,
)
from proteinlens.llm import get_client  # noqa: E402
from proteinlens.report import build_facts, generate_report  # noqa: E402
from proteinlens.viewer import COLOR_A, COLOR_B, single_html, superposition_html  # noqa: E402

RESULTS = REPO / "benchmarks" / "results"
FIGURES = REPO / "benchmarks" / "figures"
EXAMPLES = {
    "Hemoglobin α vs β": ("4HHB:A", "4HHB:B"),
    "Ubiquitin vs SUMO": ("1UBQ", "1WM3"),
    "Ubiquitin vs TIM barrel": ("1UBQ", "8TIM"),
    "AlphaFold vs crystal (Hb α)": ("P69905", "4HHB:A"),
}

st.set_page_config(page_title="ProteinLens", page_icon="🔬", layout="wide")


# --------------------------------------------------------------------------------------
# Setup and cached helpers
# --------------------------------------------------------------------------------------


@st.cache_resource
def seed_cache() -> None:
    """Copy the bundled example structures into the cache so the examples work offline."""
    src = REPO / "tests" / "fixtures" / "cache"
    for f in src.rglob("*"):
        if f.is_file():
            dest = cache_dir() / f.relative_to(src)
            if not dest.exists():
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(f, dest)


@st.cache_resource(show_spinner=False)
def cached_load(identifier: str, chain: str | None, rng: str | None) -> Structure:
    return load_structure(identifier, chain=chain or None, residue_range=rng or None)


@st.cache_resource(show_spinner=False)
def cached_classify(key: str, _s: Structure, top_k: int = 5):
    return classify(_s, top_k=top_k)


def llm_client():
    # Streamlit Cloud: put ANTHROPIC_API_KEY in the app's secrets.
    try:
        if "ANTHROPIC_API_KEY" in st.secrets and not os.environ.get("ANTHROPIC_API_KEY"):
            os.environ["ANTHROPIC_API_KEY"] = st.secrets["ANTHROPIC_API_KEY"]
    except Exception:  # no secrets file: fine
        pass
    return get_client()


def structure_input(label: str, key: str, default: str) -> tuple[Structure | None, str]:
    """Text ID (+ chain, range) or an uploaded file. Returns (structure, cache key)."""
    st.markdown(f"**{label}**")
    ident = st.text_input(
        "PDB ID (e.g. 4HHB or 4HHB:B), UniProt ID or CATH domain",
        value=st.session_state.get(f"{key}_id", default),
        key=f"{key}_id_input",
    )
    c1, c2 = st.columns(2)
    chain = c1.text_input("Chain (optional)", key=f"{key}_chain")
    rng = c2.text_input("Residue range (optional, e.g. 1-100)", key=f"{key}_range")
    upload = st.file_uploader(
        "…or upload a .pdb / .cif file", type=["pdb", "cif", "ent", "gz"], key=f"{key}_upload"
    )
    try:
        if upload is not None:
            s = parse_structure(
                upload.getvalue(),
                chain=chain or None,
                name=upload.name.split(".")[0],
                filename=upload.name,
            )
            if rng:
                from proteinlens.data import parse_residue_range

                s = parse_structure(
                    upload.getvalue(),
                    chain=chain or None,
                    residue_range=parse_residue_range(rng),
                    name=upload.name.split(".")[0],
                    filename=upload.name,
                )
            return s, f"upload:{upload.name}:{upload.size}:{chain}:{rng}"
        if not ident.strip():
            return None, ""
        return cached_load(ident.strip(), chain, rng), f"{ident}:{chain}:{rng}"
    except ProteinLensError as err:
        st.error(str(err))
        return None, ""


def deviation_figure(a: Structure, b: Structure, result) -> go.Figure:
    d = result.per_residue_distance_a()
    x = np.arange(len(a))
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x,
            y=d,
            mode="lines",
            line=dict(color=COLOR_A, width=2),
            name="distance",
            customdata=np.array(a.residue_ids),
            hovertemplate="residue %{customdata} of A<br>%{y:.2f} Å<extra></extra>",
        )
    )
    fig.add_hline(
        y=DEVIATION_CUTOFF,
        line_dash="dash",
        line_color="#52514e",
        annotation_text=f"{DEVIATION_CUTOFF:g} Å",
        annotation_position="top left",
    )
    tick_step = max(1, len(a) // 12)
    fig.update_layout(
        height=280,
        margin=dict(l=10, r=10, t=30, b=10),
        showlegend=False,
        title=f"Distance after superposition, per residue of {a.name} (gaps = unaligned)",
        xaxis=dict(
            title="residue",
            tickmode="array",
            tickvals=x[::tick_step],
            ticktext=[a.residue_ids[i] for i in x[::tick_step]],
        ),
        yaxis=dict(title="Å", rangemode="tozero"),
        plot_bgcolor="white",
    )
    fig.update_yaxes(gridcolor="#e6e5e0")
    return fig


# --------------------------------------------------------------------------------------
# Pages
# --------------------------------------------------------------------------------------


def page_compare() -> None:
    st.header("Compare two structures")
    st.caption(
        "Structural alignment with TM-align; the report below is checked number by "
        "number against what was computed."
    )
    cols = st.columns(len(EXAMPLES))
    for col, (name, (x, y)) in zip(cols, EXAMPLES.items(), strict=True):
        if col.button(name, use_container_width=True):
            st.session_state["a_id"], st.session_state["b_id"] = x, y
            st.session_state["a_id_input"], st.session_state["b_id_input"] = x, y
            st.session_state["run_compare"] = True
    left, right = st.columns(2)
    with left:
        a, key_a = structure_input("Structure A", "a", "4HHB:A")
    with right:
        b, key_b = structure_input("Structure B", "b", "4HHB:B")
    with_folds = st.checkbox(
        "Also find the nearest CATH folds for both (adds a few seconds)", value=False
    )
    if st.button("Compare", type="primary"):
        st.session_state["run_compare"] = True
    if not st.session_state.get("run_compare") or a is None or b is None:
        return

    with st.spinner("Aligning with TM-align…"):
        r = compare(a, b)
        hits_a = cached_classify(key_a, a, 3) if with_folds else None
        hits_b = cached_classify(key_b, b, 3) if with_folds else None
    for s in (a, b):
        for w in s.warnings:
            st.warning(f"{s.name}: {w}")

    m = st.columns(5)
    m[0].metric(f"TM-score (norm. by A, {len(a)} res)", f"{r.tm_score_a:.3f}")
    m[1].metric(f"TM-score (norm. by B, {len(b)} res)", f"{r.tm_score_b:.3f}")
    m[2].metric("RMSD", f"{r.rmsd:.2f} Å")
    m[3].metric("Aligned pairs (< 5 Å)", f"{r.aligned_length}")
    m[4].metric("Sequence identity", f"{100 * r.seq_identity:.1f}%")
    if r.likely_same_fold:
        st.info(f"**Likely the same fold** — the larger TM-score is above {SAME_FOLD_TM}.")
    else:
        st.info(f"**Likely different folds** — neither TM-score is above {SAME_FOLD_TM}.")

    v1, v2 = st.columns([3, 2])
    with v1:
        st.markdown(
            f"**3D superposition** — <span style='color:{COLOR_A}'>■ A: {a.name}</span> "
            f"(moved onto B) &nbsp; <span style='color:{COLOR_B}'>■ B: {b.name}</span> "
            "&nbsp; <span style='color:#c0392b'>■ largest deviations</span>",
            unsafe_allow_html=True,
        )
        components.html(superposition_html(a, b, r), height=470)
    with v2:
        st.markdown("**Grounded report**")
        facts = build_facts(a, b, r, hits_a, hits_b)
        with st.spinner("Writing report…"):
            out = generate_report(facts, llm_client())
        st.write(out.text)
        if out.verified:
            st.success(
                "Verified against computed data ✓ "
                f"({out.checks[-1].numbers_checked} numbers checked)"
            )
        else:
            st.error("Report could not be verified ✗")
        src = "LLM narration" if out.source == "llm" else "template (no LLM)"
        st.caption(f"Source: {src}. {out.note}")
        with st.expander("Facts the report is allowed to use"):
            st.json(facts)

    st.plotly_chart(deviation_figure(a, b, r), use_container_width=True)
    if r.top_regions:
        st.markdown("**Largest deviations after superposition**")
        st.dataframe(pd.DataFrame([g.__dict__ for g in r.top_regions]), hide_index=True)
    with st.expander("Sequence alignment produced by TM-align (':' = pair closer than 5 Å)"):
        st.code(f"A {r.aligned_seq_a}\n  {r.match_line}\nB {r.aligned_seq_b}", language=None)
    if with_folds:
        st.markdown("**Nearest CATH reference folds**")
        rows = [
            {"structure": s.name, **h.__dict__}
            for s, hits in ((a, hits_a), (b, hits_b))
            for h in hits or []
        ]
        st.dataframe(pd.DataFrame(rows), hide_index=True)


def page_classify() -> None:
    st.header("Classify a structure")
    ref = load_reference()
    st.caption(
        f"Aligns your structure against {len(ref)} CATH reference domains from "
        f"{ref.table['fold'].nunique()} folds and reports the closest ones. "
        "Folds outside this reference set cannot be predicted."
    )
    ex = st.columns(4)
    for col, ident in zip(ex, ["1UBQ", "4HHB:B", "8TIM", "P69905"], strict=True):
        if col.button(f"Example: {ident}", use_container_width=True):
            st.session_state["c_id"] = ident
            st.session_state["c_id_input"] = ident
            st.session_state["run_classify"] = True
    s, key = structure_input("Structure", "c", "1UBQ")
    if st.button("Classify", type="primary"):
        st.session_state["run_classify"] = True
    if not st.session_state.get("run_classify") or s is None:
        return
    with st.spinner(f"Aligning against {len(ref)} reference domains…"):
        hits = cached_classify(key, s, 5)
    best = hits[0]
    if best.tm_query > SAME_FOLD_TM:
        st.success(
            f"Predicted fold: **{best.fold} — {best.fold_name}** (top TM-score {best.tm_query:.3f})"
        )
    else:
        st.warning(
            f"Low confidence: the best match {best.fold} — {best.fold_name} has "
            f"TM-score {best.tm_query:.3f} (not above {SAME_FOLD_TM}). The fold may not "
            "be in the reference set."
        )
    c1, c2 = st.columns([3, 2])
    with c1:
        cols = [
            "domain",
            "tm_query",
            "tm_reference",
            "rmsd",
            "fold",
            "fold_name",
            "superfamily",
            "superfamily_name",
            "aligned_length",
            "seq_identity",
        ]
        st.dataframe(pd.DataFrame([h.__dict__ for h in hits])[cols], hide_index=True)
        fig = go.Figure(
            go.Bar(
                x=[h.tm_query for h in hits][::-1],
                y=[f"{h.domain} ({h.fold})" for h in hits][::-1],
                orientation="h",
                marker_color=COLOR_A,
                hovertemplate="%{y}<br>TM-score %{x:.3f}<extra></extra>",
            )
        )
        fig.add_vline(x=SAME_FOLD_TM, line_dash="dash", line_color="#52514e")
        fig.update_layout(
            height=260,
            margin=dict(l=10, r=10, t=10, b=10),
            xaxis=dict(range=[0, 1], title="TM-score (normalized by your structure)"),
            plot_bgcolor="white",
        )
        st.plotly_chart(fig, use_container_width=True)
    with c2:
        components.html(single_html(s), height=390)


def page_benchmarks() -> None:
    st.header("Benchmark results")
    summary_path = RESULTS / "summary.json"
    if not summary_path.exists():
        st.warning("No results yet — run `make benchmark`.")
        return
    s = json.loads(summary_path.read_text())
    st.caption(
        f"Leave-one-out 1-nearest-neighbour classification of {s['n_domains']} CATH S40 domains "
        f"({s['n_folds']} folds, {s['n_superfamilies']} superfamilies). 95% bootstrap CIs. "
        "All numbers are read from benchmarks/results/summary.json."
    )
    rows = []
    for block, label in (
        (s["all_queries"], "all queries"),
        (s["remote_homologs"], "remote homologs"),
    ):
        for m in ("structure", "sequence", "majority"):
            for lv in ("fold", "superfamily"):
                e = block[m][lv]
                rows.append(
                    {
                        "subset": f"{label} (n={e['n']})",
                        "method": m,
                        "level": lv,
                        "accuracy": f"{100 * e['accuracy']:.1f}%",
                        "95% CI": f"{100 * e['ci95'][0]:.1f}–{100 * e['ci95'][1]:.1f}%",
                    }
                )
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    c1, c2 = st.columns(2)
    c1.image(str(FIGURES / "accuracy.png"))
    c2.image(str(FIGURES / "accuracy_vs_identity.png"))
    c1.image(str(FIGURES / "tm_distributions.png"))
    c2.image(str(FIGURES / "confusion_matrix.png"))
    with st.expander("Per-fold accuracy"):
        st.dataframe(pd.read_csv(RESULTS / "per_fold.csv"), hide_index=True)
    with st.expander("Failure analysis"):
        st.markdown((RESULTS / "FAILURE_ANALYSIS.md").read_text())
    st.subheader("Grounded vs ungrounded LLM experiment")
    exp = RESULTS / "llm_experiment" / "summary.json"
    if exp.exists():
        e = json.loads(exp.read_text())
        if e.get("status") == "PENDING":
            b = e["tm_rule_baseline"]
            st.warning(
                f"PENDING — run `{e['how_to_run']}`. Already computed: the no-LLM rule "
                f'"same fold if TM-score > 0.5" scores {100 * b["accuracy"]:.1f}% on the '
                f"{e['n_pairs']} experiment pairs."
            )
        else:
            st.json(e)


def page_about() -> None:
    st.header("About ProteinLens")
    intro = (REPO / "README.md").read_text().split("## Quickstart")[0]
    # Local image paths do not render inside Streamlit markdown; drop those lines.
    keep = [
        ln for ln in intro.splitlines() if not ln.lstrip().startswith(("![", "<p", "<img", "</p"))
    ]
    st.markdown("\n".join(keep))


# --------------------------------------------------------------------------------------

seed_cache()
st.sidebar.title("🔬 ProteinLens")
page = st.sidebar.radio("Page", ["Compare", "Classify", "Benchmarks", "About"])
client = llm_client()
if client is None:
    st.sidebar.info("LLM narration: off (no API key) — reports use the verified template.")
else:
    st.sidebar.success(f"LLM narration: on ({client.model}); every report is checked.")
st.sidebar.caption(
    "Structures: RCSB PDB, AlphaFold DB (CC BY 4.0). Folds: CATH v4.4. Alignment: TM-align."
)
{
    "Compare": page_compare,
    "Classify": page_classify,
    "Benchmarks": page_benchmarks,
    "About": page_about,
}[page]()
