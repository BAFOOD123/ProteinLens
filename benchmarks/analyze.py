"""Analyze the all-vs-all matrices: leave-one-out 1-NN fold/superfamily classification.

All analysis choices below were fixed before looking at results:

* Primary method: structure 1-NN using the TM-score normalized by the query length.
* Baselines: sequence 1-NN (Smith-Waterman score, BLOSUM62) and leave-one-out
  majority class.
* 95 % confidence intervals: percentile bootstrap over queries, 2000 resamples, seed 0.
* Remote-homolog subset: queries that have at least one other member of their CATH
  superfamily in the set, but < 25 % sequence identity to *every* such member
  (identity = identical residues in the local alignment / shorter length).
* Sensitivity analysis (secondary, reported for transparency only): other TM-score
  normalizations.

Outputs: benchmarks/results/{summary.json, per_query.csv, per_fold.csv, failures.csv,
pair_stats.csv} and figures in benchmarks/figures/ (or benchmarks/smoke/ with --smoke).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from proteinlens.classify import majority_baseline, nearest_neighbour  # noqa: E402

SEED = 0
N_BOOT = 2000
REMOTE_ID = 0.25
# Validated categorical palette slots (see dataviz reference palette): blue, orange, aqua.
C_STRUCT, C_SEQ, C_THIRD, C_MAJ = "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
INK, MUTED = "#0b0b0b", "#52514e"


def bootstrap_ci(correct: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """95 % percentile bootstrap CI of the mean of a 0/1 vector (resampling queries)."""
    n = len(correct)
    if n == 0:
        return (float("nan"), float("nan"))
    idx = rng.integers(0, n, size=(N_BOOT, n))
    means = correct[idx].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def paired_diff_ci(a: np.ndarray, b: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """95 % CI of mean(a) - mean(b), resampling the same queries for both (paired)."""
    n = len(a)
    idx = rng.integers(0, n, size=(N_BOOT, n))
    diffs = a[idx].mean(axis=1) - b[idx].mean(axis=1)
    return float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def acc_entry(correct: np.ndarray, rng: np.random.Generator) -> dict:
    lo, hi = bootstrap_ci(correct, rng)
    return {
        "accuracy": round(float(correct.mean()), 4) if len(correct) else None,
        "ci95": [round(lo, 4), round(hi, 4)],
        "n_correct": int(correct.sum()),
        "n": int(len(correct)),
    }


def style_axes(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    ax.grid(axis="y", color="#e6e5e0", linewidth=0.8)
    ax.set_axisbelow(True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["main", "hard"], default="main", help="benchmark set")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    suffix = "" if args.set == "main" else f"_{args.set}"
    res_dir = REPO / "benchmarks" / ("smoke" if args.smoke else f"results{suffix}")
    fig_dir = res_dir if args.smoke else REPO / "benchmarks" / f"figures{suffix}"
    fig_dir.mkdir(parents=True, exist_ok=True)

    pw = np.load(res_dir / "pairwise.npz")
    names = [str(x) for x in pw["domain"]]
    table = pd.read_csv(REPO / "benchmarks" / f"data{suffix}" / "dataset.csv").set_index("domain")
    table = table.loc[names].reset_index()
    fold = table["fold"].to_numpy()
    sf = table["superfamily"].to_numpy()
    n = len(names)
    tm, seq_score, seq_id = pw["tm"], pw["seq_score"], pw["seq_id"]
    rng = np.random.default_rng(SEED)

    # ---- predictions -------------------------------------------------------------
    pred = {
        "structure": {
            "fold": nearest_neighbour(tm, fold),
            "superfamily": nearest_neighbour(tm, sf),
        },
        "sequence": {
            "fold": nearest_neighbour(seq_score, fold),
            "superfamily": nearest_neighbour(seq_score, sf),
        },
        "majority": {"fold": majority_baseline(fold), "superfamily": majority_baseline(sf)},
    }
    truth = {"fold": fold, "superfamily": sf}
    correct = {m: {lv: (pred[m][lv] == truth[lv]).astype(float) for lv in truth} for m in pred}

    # ---- subsets -----------------------------------------------------------------
    same_sf = (sf[:, None] == sf[None, :]) & ~np.eye(n, dtype=bool)
    same_fold = (fold[:, None] == fold[None, :]) & ~np.eye(n, dtype=bool)
    has_sf_partner = same_sf.any(axis=1)
    max_id_same_sf = np.where(same_sf, seq_id, -1).max(axis=1)
    max_id_same_fold = np.where(same_fold, seq_id, -1).max(axis=1)
    remote = has_sf_partner & (max_id_same_sf < REMOTE_ID)

    summary: dict = {
        "n_domains": n,
        "n_folds": int(len(set(fold))),
        "n_superfamilies": int(len(set(sf))),
        "bootstrap": {"resamples": N_BOOT, "seed": SEED, "method": "percentile, over queries"},
        "primary_metric": "TM-score normalized by query length, leave-one-out 1-NN",
        "all_queries": {},
        "remote_homologs": {
            "definition": (
                f"query has >=1 same-superfamily domain in the set, and < {REMOTE_ID:.0%} "
                "sequence identity to every one of them"
            ),
            "n": int(remote.sum()),
        },
        "structure_minus_sequence": {},
    }
    for m in pred:
        summary["all_queries"][m] = {lv: acc_entry(correct[m][lv], rng) for lv in truth}
        summary["remote_homologs"][m] = {lv: acc_entry(correct[m][lv][remote], rng) for lv in truth}
    for lv in truth:
        a, b = correct["structure"][lv], correct["sequence"][lv]
        lo, hi = paired_diff_ci(a, b, rng)
        ra, rb = a[remote], b[remote]
        rlo, rhi = paired_diff_ci(ra, rb, rng)
        summary["structure_minus_sequence"][lv] = {
            "all": {
                "diff": round(float(a.mean() - b.mean()), 4),
                "ci95": [round(lo, 4), round(hi, 4)],
            },
            "remote": {
                "diff": round(float(ra.mean() - rb.mean()), 4),
                "ci95": [round(rlo, 4), round(rhi, 4)],
            },
        }

    # ---- sensitivity: other TM normalizations (secondary) --------------------------
    variants = {
        "query_length (primary)": tm,
        "reference_length": tm.T,
        "shorter_length (max of the two)": np.maximum(tm, tm.T),
        "longer_length (min of the two)": np.minimum(tm, tm.T),
    }
    summary["sensitivity_tm_normalization"] = {
        k: {
            "fold": round(float((nearest_neighbour(v, fold) == fold).mean()), 4),
            "superfamily": round(float((nearest_neighbour(v, sf) == sf).mean()), 4),
        }
        for k, v in variants.items()
    }

    # ---- pair-level statistics (how well does TM > 0.5 separate folds?) -------------
    iu = np.triu_indices(n, k=1)
    tm_max = np.maximum(tm, tm.T)[iu]
    rel = np.where(
        same_sf[iu],
        "same superfamily",
        np.where(same_fold[iu], "same fold, different superfamily", "different fold"),
    )
    pairs = pd.DataFrame({"relation": rel, "tm_max": tm_max, "seq_id": seq_id[iu]})
    pair_stats = (
        pairs.groupby("relation")
        .agg(
            n_pairs=("tm_max", "size"),
            median_tm=("tm_max", "median"),
            frac_tm_above_0_5=("tm_max", lambda x: float((x > 0.5).mean())),
            median_seq_id=("seq_id", "median"),
        )
        .round(4)
    )
    pair_stats.to_csv(res_dir / "pair_stats.csv")
    summary["pair_stats"] = pair_stats.reset_index().to_dict(orient="records")

    # ---- per-query and per-fold tables --------------------------------------------
    nn_struct = np.argmax(np.where(np.eye(n, dtype=bool), -np.inf, tm), axis=1)
    nn_seq = np.argmax(np.where(np.eye(n, dtype=bool), -np.inf, seq_score), axis=1)
    best_same_fold_tm = np.where(same_fold, tm, -1).max(axis=1)
    per_query = pd.DataFrame(
        {
            "domain": names,
            "fold": fold,
            "superfamily": sf,
            "length": table["n_residues"],
            "struct_nn": [names[j] for j in nn_struct],
            "struct_nn_tm": tm[np.arange(n), nn_struct].round(4),
            "struct_pred_fold": pred["structure"]["fold"],
            "struct_pred_sf": pred["structure"]["superfamily"],
            "seq_nn": [names[j] for j in nn_seq],
            "seq_nn_score": seq_score[np.arange(n), nn_seq],
            "seq_pred_fold": pred["sequence"]["fold"],
            "seq_pred_sf": pred["sequence"]["superfamily"],
            "best_same_fold_tm": best_same_fold_tm.round(4),
            "max_seq_id_same_sf": np.where(has_sf_partner, max_id_same_sf, np.nan).round(4),
            "max_seq_id_same_fold": max_id_same_fold.round(4),
            "remote_homolog": remote,
        }
    )
    per_query.to_csv(res_dir / "per_query.csv", index=False)
    per_fold = (
        per_query.assign(struct_ok=correct["structure"]["fold"], seq_ok=correct["sequence"]["fold"])
        .groupby("fold")
        .agg(
            n=("domain", "size"),
            structure_acc=("struct_ok", "mean"),
            sequence_acc=("seq_ok", "mean"),
        )
        .round(4)
    )
    per_fold.insert(0, "fold_name", table.groupby("fold")["fold_name"].first())
    per_fold.to_csv(res_dir / "per_fold.csv")

    # Failures: wrong fold by structure, most confident (highest TM to the wrong neighbour) first.
    fails = per_query[per_query["struct_pred_fold"] != per_query["fold"]].copy()
    fold_names = table.groupby("fold")["fold_name"].first().to_dict()
    lengths = dict(zip(names, table["n_residues"], strict=True))
    fails["true_fold_name"] = fails["fold"].map(fold_names)
    fails["pred_fold_name"] = fails["struct_pred_fold"].map(fold_names)
    fails["nn_length"] = fails["struct_nn"].map(lengths)
    fails = fails.sort_values("struct_nn_tm", ascending=False)
    fails[
        [
            "domain",
            "length",
            "fold",
            "true_fold_name",
            "struct_nn",
            "nn_length",
            "struct_pred_fold",
            "pred_fold_name",
            "struct_nn_tm",
            "best_same_fold_tm",
        ]
    ].to_csv(res_dir / "failures.csv", index=False)
    summary["n_structure_fold_errors"] = int(len(fails))
    sf_fails = per_query[per_query["struct_pred_sf"] != per_query["superfamily"]].copy()
    sf_fails["nn_superfamily"] = sf_fails["struct_pred_sf"]
    sf_fails["seq_correct_sf"] = sf_fails["seq_pred_sf"] == sf_fails["superfamily"]
    sf_fails.sort_values("struct_nn_tm", ascending=False)[
        [
            "domain",
            "length",
            "fold",
            "superfamily",
            "struct_nn",
            "nn_superfamily",
            "struct_nn_tm",
            "max_seq_id_same_sf",
            "seq_correct_sf",
        ]
    ].to_csv(res_dir / "superfamily_failures.csv", index=False)
    summary["n_structure_superfamily_errors"] = int(len(sf_fails))
    # Leave-one-out majority is always wrong when classes are (nearly) balanced: removing the
    # query makes its own class smaller than the others. The share of the largest class is
    # the more intuitive "chance" reference, so report it too.
    summary["largest_class_fraction"] = {
        "fold": round(float(pd.Series(fold).value_counts().iloc[0] / n), 4),
        "superfamily": round(float(pd.Series(sf).value_counts().iloc[0] / n), 4),
    }
    # Superfamily accuracy counted only over queries that have a same-superfamily partner in
    # the set (a query whose superfamily is unique here can never be classified correctly).
    summary["queries_with_superfamily_partner"] = {
        "n": int(has_sf_partner.sum()),
        **{
            m: acc_entry(correct[m]["superfamily"][has_sf_partner], rng)
            for m in ("structure", "sequence")
        },
    }
    # Abstention at the literature threshold (pre-specified, not tuned): answer only when the
    # nearest neighbour's TM-score is above 0.5.
    top_tm = tm[np.arange(n), nn_struct]
    confident = top_tm > 0.5
    summary["abstain_below_tm_0_5"] = {
        "coverage": round(float(confident.mean()), 4),
        "n_answered": int(confident.sum()),
        "fold_accuracy_answered": acc_entry(correct["structure"]["fold"][confident], rng),
        "fold_accuracy_abstained": acc_entry(correct["structure"]["fold"][~confident], rng),
        "superfamily_accuracy_answered": acc_entry(
            correct["structure"]["superfamily"][confident], rng
        ),
    }
    high_id = np.triu(seq_id, 1) >= 0.4
    summary["pairs_with_seq_id_ge_40pct"] = int(high_id.sum())

    # Accuracy vs sequence identity (queries that have a same-superfamily partner).
    bins = [0, 0.15, 0.20, 0.25, 0.30, 1.01]
    labels = ["<15%", "15-20%", "20-25%", "25-30%", "≥30%"]
    sub = per_query[has_sf_partner].copy()
    sub["bin"] = pd.cut(sub["max_seq_id_same_sf"], bins, labels=labels, right=False)
    sub["s_ok"] = correct["structure"]["superfamily"][has_sf_partner]
    sub["q_ok"] = correct["sequence"]["superfamily"][has_sf_partner]
    by_bin = (
        sub.groupby("bin", observed=False)
        .agg(n=("domain", "size"), structure=("s_ok", "mean"), sequence=("q_ok", "mean"))
        .round(4)
    )
    by_bin.to_csv(res_dir / "accuracy_vs_identity.csv")
    summary["accuracy_vs_identity_superfamily"] = (
        by_bin.reset_index().astype({"bin": str}).to_dict(orient="records")
    )

    (res_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    # ---- figures -----------------------------------------------------------------
    plt.rcParams.update({"font.size": 10, "axes.titleweight": "bold", "text.color": INK})

    # 1. headline accuracy bars
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
    for ax, (title, block) in zip(
        axes,
        [("All queries", summary["all_queries"]), ("Remote homologs", summary["remote_homologs"])],
        strict=True,
    ):
        x = np.arange(2)
        if block["structure"]["fold"]["n"] == 0:  # e.g. the smoke subset has no remote homologs
            ax.set_title(f"{title} (n=0)")
            continue
        for k, (m, col) in enumerate(
            [("structure", C_STRUCT), ("sequence", C_SEQ), ("majority", C_MAJ)]
        ):
            accs = [block[m][lv]["accuracy"] for lv in ("fold", "superfamily")]
            cis = [block[m][lv]["ci95"] for lv in ("fold", "superfamily")]
            err = np.array(
                [
                    [a - c[0] for a, c in zip(accs, cis, strict=True)],
                    [c[1] - a for a, c in zip(accs, cis, strict=True)],
                ]
            )
            ax.bar(
                x + (k - 1) * 0.26,
                accs,
                0.24,
                color=col,
                label="majority (leave-one-out)" if m == "majority" else m,
                yerr=err,
                capsize=3,
                ecolor=MUTED,
            )
            for xi, a, c in zip(x + (k - 1) * 0.26, accs, cis, strict=True):
                ax.text(xi, c[1] + 0.02, f"{100 * a:.1f}%", ha="center", fontsize=7, color=INK)
        n_q = block["structure"]["fold"]["n"]
        ax.set_title(f"{title} (n={n_q})")
        ax.set_xticks(x, ["Fold (CATH topology)", "Superfamily"])
        ax.set_ylim(0, 1.12)
        style_axes(ax)
    axes[0].set_ylabel("Leave-one-out 1-NN accuracy")
    handles, labs = axes[0].get_legend_handles_labels()
    fig.legend(handles, labs, frameon=False, loc="lower center", ncol=3, fontsize=8)
    fig.suptitle("Fold classification accuracy with 95% bootstrap CIs", fontsize=11)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(fig_dir / "accuracy.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # 2. accuracy vs sequence identity
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    x = np.arange(len(labels))
    ax.plot(x, by_bin["structure"], "-o", color=C_STRUCT, lw=2, ms=7, label="structure (TM-score)")
    ax.plot(x, by_bin["sequence"], "-o", color=C_SEQ, lw=2, ms=7, label="sequence (Smith-Waterman)")
    ax.set_xticks(x, [f"{lab}\n(n={k})" for lab, k in zip(labels, by_bin["n"], strict=True)])
    ax.set_xlabel("Highest sequence identity to any same-superfamily domain in the set")
    ax.set_ylabel("Superfamily accuracy (1-NN)")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("Structure keeps working where sequence similarity fades")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(fig_dir / "accuracy_vs_identity.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # 3. confusion matrix (fold level, structure); unreadable for many folds, so skipped then
    if len(set(fold)) <= 30:
        folds_sorted = sorted(set(fold))
        k = len(folds_sorted)
        pos = {f: i for i, f in enumerate(folds_sorted)}
        cm = np.zeros((k, k), dtype=int)
        for t, p in zip(fold, pred["structure"]["fold"], strict=True):
            cm[pos[t], pos[p]] += 1
        fig, ax = plt.subplots(figsize=(8.5, 7.5))
        ax.imshow(cm / cm.sum(axis=1, keepdims=True), cmap="Blues", vmin=0, vmax=1)
        short = [f"{f} {fold_names[f][:28]}" for f in folds_sorted]
        ax.set_xticks(range(k), folds_sorted, rotation=90, fontsize=7)
        ax.set_yticks(range(k), short, fontsize=7)
        for i in range(k):
            for j in range(k):
                if cm[i, j]:
                    frac = cm[i, j] / cm[i].sum()
                    ax.text(
                        j,
                        i,
                        cm[i, j],
                        ha="center",
                        va="center",
                        fontsize=7,
                        color="white" if frac > 0.6 else INK,
                    )
        ax.set_xlabel("Predicted fold (structure 1-NN)")
        ax.set_ylabel("True fold (CATH)")
        ax.set_title("Fold confusion matrix (counts; color = row fraction)")
        fig.tight_layout()
        fig.savefig(fig_dir / "confusion_matrix.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # 4. TM-score distributions by relationship
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    edges = np.linspace(0, 1, 41)
    for rel_name, col in [
        ("different fold", C_MAJ),
        ("same fold, different superfamily", C_THIRD),
        ("same superfamily", C_STRUCT),
    ]:
        vals = pairs.loc[pairs["relation"] == rel_name, "tm_max"]
        ax.hist(
            vals,
            edges,
            density=True,
            histtype="step",
            lw=2,
            color=col,
            label=f"{rel_name} (n={len(vals)})",
        )
    ax.axvline(0.5, color=INK, lw=1, ls="--")
    ax.text(0.51, ax.get_ylim()[1] * 0.92, "TM = 0.5", fontsize=8, color=INK)
    ax.set_xlabel("TM-score (normalized by the shorter protein)")
    ax.set_ylabel("Density of pairs")
    ax.set_title("TM-score separates folds, but not perfectly")
    ax.legend(frameon=False, fontsize=8)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(fig_dir / "tm_distributions.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    for m in ("structure", "sequence", "majority"):
        b = summary["all_queries"][m]
        print(f"{m:9s} fold {b['fold']['accuracy']} {b['fold']['ci95']}", end="  ")
        print(f"superfamily {b['superfamily']['accuracy']} {b['superfamily']['ci95']}")
    r = summary["remote_homologs"]
    rs, rq = r["structure"]["superfamily"]["accuracy"], r["sequence"]["superfamily"]["accuracy"]
    print(f"remote homologs (n={r['n']}): superfamily acc structure {rs} vs sequence {rq}")
    print(f"wrote {res_dir}/summary.json and figures in {fig_dir}")


if __name__ == "__main__":
    main()
