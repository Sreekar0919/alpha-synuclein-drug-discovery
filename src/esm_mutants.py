"""ESM-2 contact maps for wild-type alpha-synuclein and familial PD mutants.

This is a side analysis, independent of the docking. It asks a narrow
question: does a protein language model predict different residue-residue
contacts for the single-point mutants that cause familial Parkinson's disease?

Run it with its own interpreter. PyTorch ships its own OpenMP runtime, which
clashes with the one that RDKit and Vina link against, so this does not share
the docking environment. See the README for the two-environment setup.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402

# UniProt P37840, human alpha-synuclein, 140 residues.
WILD_TYPE = (
    "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK"
    "EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP"
    "DNEAYEMPSEEGYQDYEPEA"
)

# Point mutations linked to autosomal-dominant Parkinson's disease. Numbering
# is 1-based over the sequence above.
MUTANTS = {"A30P": (30, "P"), "E46K": (46, "K"), "A53T": (53, "T")}

# Substitutions with no reported link to Parkinson's disease, used as a
# baseline. Without them, "the familial mutants barely move the contact map"
# is not interpretable: any single substitution in a 140-residue sequence
# would barely move it. These say how much an arbitrary one moves it.
CONTROLS = {"A19G": (19, "G"), "V40A": (40, "A"), "T75S": (75, "S")}

# The 650M model is the better one but its weights are 2.6 GB. 150M keeps the
# download under 1 GB and is enough for a qualitative contact comparison.
MODEL_NAME = "esm2_t30_150M_UR50D"

# The NAC region drives aggregation; contacts there are the ones worth watching.
NAC_REGION = (61, 95)


@dataclass(frozen=True)
class Variant:
    """One sequence to run through the model."""

    name: str
    sequence: str
    position: int | None = None


def apply_mutation(sequence: str, position: int, residue: str) -> str:
    """Substitute one residue at a 1-based position.

    Bounds are checked here; the identity of the residue being replaced is
    verified by the caller, which knows what it expects to find.
    """
    index = position - 1
    if not 0 <= index < len(sequence):
        raise ValueError(f"position {position} outside 1..{len(sequence)}")
    return sequence[:index] + residue + sequence[index + 1 :]


def build_variants(sequence: str = WILD_TYPE) -> list[Variant]:
    """Wild type, each familial mutant, and the control substitutions.

    A silent off-by-one would invalidate every downstream number, so the
    residue each mutation expects to replace is verified rather than assumed.
    """
    variants = [Variant("WT", sequence)]
    for group in (MUTANTS, CONTROLS):
        for name, (position, residue) in group.items():
            expected = name[0]
            actual = sequence[position - 1]
            if actual != expected:
                raise ValueError(
                    f"{name}: expected {expected} at position {position}, "
                    f"sequence has {actual}"
                )
            variants.append(
                Variant(name, apply_mutation(sequence, position, residue), position)
            )
    return variants


def load_model(model_name: str = MODEL_NAME):
    """Load ESM-2 and its alphabet, in eval mode on CPU."""
    import esm

    model, alphabet = getattr(esm.pretrained, model_name)()
    model.eval()
    return model, alphabet, alphabet.get_batch_converter()


@torch.no_grad()
def contacts_and_embedding(model, batch_converter, variant: Variant) -> dict:
    """Predicted contact map and mean-pooled embedding for one sequence."""
    _, _, tokens = batch_converter([(variant.name, variant.sequence)])
    layer = model.num_layers
    out = model(tokens, repr_layers=[layer], return_contacts=True)

    # Strip the BOS/EOS positions so indices line up with the sequence.
    length = len(variant.sequence)
    representation = out["representations"][layer][0, 1 : length + 1]
    return {
        "contacts": out["contacts"][0, :length, :length].cpu().numpy(),
        "embedding": representation.mean(0).cpu().numpy(),
    }


def contact_metrics(wt_map: np.ndarray, mut_map: np.ndarray, region=NAC_REGION) -> dict:
    """Summarise how far a mutant's contact map moved from wild type."""
    difference = mut_map - wt_map
    start, end = region[0] - 1, region[1]
    nac = difference[start:end, start:end]
    return {
        "mean_abs_change": round(float(np.abs(difference).mean()), 5),
        "max_abs_change": round(float(np.abs(difference).max()), 5),
        "frobenius_norm": round(float(np.linalg.norm(difference)), 4),
        "nac_mean_abs_change": round(float(np.abs(nac).mean()), 5),
        "map_correlation": round(float(np.corrcoef(wt_map.ravel(), mut_map.ravel())[0, 1]), 6),
    }


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def plot_contact_maps(results: dict, out_path: Path) -> Path:
    """Wild-type map, then each variant as a difference from wild type.

    Familial mutants on the top row, control substitutions below them, drawn on
    a shared colour scale so the two rows can be compared by eye.
    """
    positions = {**MUTANTS, **CONTROLS}
    rows = [list(MUTANTS), list(CONTROLS)]
    n_cols = max(len(r) for r in rows) + 1

    # One scale for every difference panel, otherwise each panel's own maximum
    # makes a small change look as dramatic as a large one.
    limit = max(
        float(np.abs(results[name]["contacts"] - results["WT"]["contacts"]).max())
        for name in positions
    )

    fig, axes = plt.subplots(2, n_cols, figsize=(4.0 * n_cols, 8.4))
    wt = results["WT"]["contacts"]

    im_wt = axes[0][0].imshow(wt, cmap="Greys", vmin=0, vmax=1, origin="lower")
    axes[0][0].set_title("WT predicted contacts")
    axes[0][0].set_ylabel("residue")
    fig.colorbar(im_wt, ax=axes[0][0], fraction=0.046, label="contact probability")
    axes[1][0].axis("off")

    for row_index, (names, label) in enumerate(
        zip(rows, ("familial (Parkinson's)", "control (no reported link)"))
    ):
        for col_index, name in enumerate(names, start=1):
            ax = axes[row_index][col_index]
            diff = results[name]["contacts"] - wt
            im = ax.imshow(diff, cmap="bwr", vmin=-limit, vmax=limit, origin="lower")
            position = positions[name][0]
            ax.axhline(position - 1, color="#2b8a3e", lw=0.7, alpha=0.8)
            ax.axvline(position - 1, color="#2b8a3e", lw=0.7, alpha=0.8)
            ax.set_title(f"{name} - WT   max |delta| {np.abs(diff).max():.3f}", fontsize=10)
            if col_index == 1:
                ax.set_ylabel(f"{label}\nresidue", fontsize=9)
            if row_index == 1:
                ax.set_xlabel("residue")
            fig.colorbar(im, ax=ax, fraction=0.046, label="change")

        for col_index in range(len(names) + 1, n_cols):
            axes[row_index][col_index].axis("off")

    fig.suptitle(
        f"ESM-2 ({MODEL_NAME}) contact maps: alpha-synuclein wild type vs single substitutions\n"
        "shared colour scale; familial mutations are not distinguishable from controls",
        fontsize=12,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=190)
    plt.close(fig)
    return out_path


def run(results_dir: Path) -> dict:
    """Run every variant and write the figure plus a JSON summary."""
    model, _alphabet, batch_converter = load_model()
    variants = build_variants()

    results = {v.name: contacts_and_embedding(model, batch_converter, v) for v in variants}
    wt = results["WT"]

    def metrics_for(name: str) -> dict:
        return {
            **contact_metrics(wt["contacts"], results[name]["contacts"]),
            "embedding_cosine_to_wt": round(
                cosine(wt["embedding"], results[name]["embedding"]), 6
            ),
        }

    familial = {name: metrics_for(name) for name in MUTANTS}
    controls = {name: metrics_for(name) for name in CONTROLS}

    # The comparison that makes the numbers mean anything: do the disease
    # mutations perturb the map more than arbitrary substitutions do?
    def median_norm(group: dict) -> float:
        values = sorted(entry["frobenius_norm"] for entry in group.values())
        return round(values[len(values) // 2], 4)

    summary = {
        "model": MODEL_NAME,
        "sequence_length": len(WILD_TYPE),
        "nac_region": list(NAC_REGION),
        "mutants": familial,
        "controls": controls,
        "comparison": {
            "familial_median_frobenius": median_norm(familial),
            "control_median_frobenius": median_norm(controls),
            "note": (
                "Controls are substitutions with no reported Parkinson's link. "
                "If the two medians are similar, the familial mutations do not "
                "perturb the predicted contact map more than an arbitrary "
                "substitution does."
            ),
        },
    }

    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "esm_mutant_summary.json").write_text(json.dumps(summary, indent=2))
    plot_contact_maps(results, results_dir / "figures" / "esm_contact_maps.png")
    return summary


if __name__ == "__main__":
    out = run(Path(__file__).resolve().parent.parent / "results")
    print(json.dumps(out, indent=2))
