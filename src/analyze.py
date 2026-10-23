"""Scoring analysis: does the protocol rank known binders above decoys?

The headline question is not "what is the best score" but "can this protocol
tell a published alpha-synuclein interactor from an unrelated approved drug".
That is a ranking question, so it is answered with ranking statistics.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from scipy import stats  # noqa: E402

REFERENCE = "reference"
DECOY = "decoy"

# Vina scores are more negative for better predicted binding, so ranking is
# ascending and the sign is flipped wherever a "higher is better" score is
# needed (ROC, enrichment).
SCORE = "best_score"
EFFICIENCY = "ligand_efficiency"


def add_efficiency(df: pd.DataFrame) -> pd.DataFrame:
    """Ligand efficiency: score per heavy atom.

    Vina scores scale with molecular size, so a large molecule can out-score a
    small one without binding any better. Efficiency is the cheapest control
    for that confound.
    """
    out = df.copy()
    out[EFFICIENCY] = (out[SCORE] / out["heavy_atoms"]).round(4)
    return out


def ranking_stats(df: pd.DataFrame, column: str = SCORE) -> dict:
    """Mann-Whitney U, rank-biserial effect size and ROC AUC for one column."""
    refs = df.loc[df["label"] == REFERENCE, column].to_numpy()
    decoys = df.loc[df["label"] == DECOY, column].to_numpy()
    if len(refs) == 0 or len(decoys) == 0:
        raise ValueError("need both reference and decoy rows")

    # "less" because a more negative score is a better predicted binder.
    result = stats.mannwhitneyu(refs, decoys, alternative="less")
    n_pairs = len(refs) * len(decoys)
    auc = 1.0 - result.statistic / n_pairs

    return {
        "column": column,
        "n_reference": int(len(refs)),
        "n_decoy": int(len(decoys)),
        "reference_median": round(float(np.median(refs)), 3),
        "decoy_median": round(float(np.median(decoys)), 3),
        "u_statistic": float(result.statistic),
        "p_value": float(result.pvalue),
        "p_value_is_one_sided": True,
        "roc_auc": round(float(auc), 3),
        "rank_biserial": round(float(2 * auc - 1), 3),
    }


def enrichment(df: pd.DataFrame, column: str = SCORE, top_fractions=(0.05, 0.10, 0.20)) -> pd.DataFrame:
    """How many references land in the top X% of the ranked library."""
    ranked = df.sort_values(column, ascending=True).reset_index(drop=True)
    n_total = len(ranked)
    n_refs = int((ranked["label"] == REFERENCE).sum())
    baseline = n_refs / n_total

    rows = []
    for fraction in top_fractions:
        k = max(1, int(round(fraction * n_total)))
        hits = int((ranked.head(k)["label"] == REFERENCE).sum())
        rows.append(
            {
                "top_fraction": fraction,
                "n_screened": k,
                "references_found": hits,
                "references_total": n_refs,
                "hit_rate": round(hits / k, 3),
                "baseline_rate": round(baseline, 3),
                "enrichment_factor": round((hits / k) / baseline, 2) if baseline else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_distributions(df: pd.DataFrame, out_path: Path) -> Path:
    """Score and ligand-efficiency distributions, reference vs decoy."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    rng = np.random.default_rng(0)

    for ax, column, title in (
        (axes[0], SCORE, "Vina score (kcal/mol)"),
        (axes[1], EFFICIENCY, "Ligand efficiency (kcal/mol per heavy atom)"),
    ):
        groups = [df.loc[df["label"] == g, column].to_numpy() for g in (REFERENCE, DECOY)]
        # tick_labels replaced the old `labels` argument in matplotlib 3.9.
        ax.boxplot(groups, tick_labels=["reference", "decoy"], widths=0.5, showfliers=False)
        for i, values in enumerate(groups, start=1):
            jitter = rng.normal(0, 0.055, size=len(values))
            ax.scatter(i + jitter, values, s=26, alpha=0.75, zorder=3,
                       color="#B23A48" if i == 1 else "#4C6A92",
                       edgecolor="white", linewidth=0.5)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(column.replace("_", " "))
        ax.grid(axis="y", alpha=0.25)

    stat = ranking_stats(df, SCORE)
    fig.suptitle(
        f"Known alpha-synuclein interactors vs. unrelated drugs  "
        f"(ROC AUC {stat['roc_auc']}, one-sided p = {stat['p_value']:.3f})",
        fontsize=12,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def plot_ranked_library(df: pd.DataFrame, out_path: Path) -> Path:
    """Every compound in score order, references highlighted.

    Drawn as points rather than bars. Vina scores span a narrow range well away
    from zero, so bars rising from a zero baseline would spend most of the plot
    on empty space and make a 4 kcal/mol spread look like a rounding error.
    """
    ranked = df.sort_values(SCORE, ascending=True).reset_index(drop=True)
    is_ref = ranked["label"] == REFERENCE
    scores = ranked[SCORE]
    low, high = scores.min(), scores.max()
    headroom = (high - low) * 0.22

    fig, ax = plt.subplots(figsize=(12, 5.4))
    ax.vlines(range(len(ranked)), high + headroom * 0.15, scores,
              color="#C9D1DC", linewidth=1.0, zorder=1)
    ax.scatter(np.flatnonzero(~is_ref), scores[~is_ref], s=42, color="#4C6A92",
               edgecolor="white", linewidth=0.6, zorder=3, label="decoy")
    ax.scatter(np.flatnonzero(is_ref), scores[is_ref], s=72, color="#B23A48",
               edgecolor="white", linewidth=0.8, zorder=4,
               label="published alpha-synuclein interactor")

    ax.set_xlim(-1.5, len(ranked) + 0.5)
    ax.set_ylim(high + headroom * 0.15, low - headroom)
    ax.set_xlabel("compound (ranked best to worst)")
    ax.set_ylabel("Vina score (kcal/mol)")
    ax.set_title("Docked library, ranked by Vina score")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

    for i in np.flatnonzero(is_ref):
        ax.annotate(ranked.loc[i, "name"], (i, scores[i]),
                    textcoords="offset points", xytext=(0, 9), rotation=90,
                    ha="center", va="bottom", fontsize=8, color="#7A1F2B")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def plot_size_confound(df: pd.DataFrame, out_path: Path) -> Path:
    """Score against heavy-atom count, to show the size dependence directly."""
    fig, ax = plt.subplots(figsize=(6.4, 5))
    for label, color in ((DECOY, "#4C6A92"), (REFERENCE, "#B23A48")):
        sub = df[df["label"] == label]
        ax.scatter(sub["heavy_atoms"], sub[SCORE], s=38, alpha=0.8,
                   color=color, edgecolor="white", linewidth=0.5, label=label)

    slope, intercept, r, p, _ = stats.linregress(df["heavy_atoms"], df[SCORE])
    xs = np.linspace(df["heavy_atoms"].min(), df["heavy_atoms"].max(), 50)
    ax.plot(xs, slope * xs + intercept, "--", color="#444", linewidth=1.2,
            label=f"fit: r² = {r**2:.2f}, p = {p:.1e}")

    ax.set_xlabel("heavy atoms")
    ax.set_ylabel("Vina score (kcal/mol)")
    ax.set_title("Docking score vs. molecular size")
    ax.invert_yaxis()
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path
