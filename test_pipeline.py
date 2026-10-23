#!/usr/bin/env python
"""Self-checks for the parts of the pipeline that could fail silently.

    /opt/anaconda3/envs/asyn/bin/python test_pipeline.py

These cover the logic where a bug would produce plausible-looking numbers
rather than an error: the symmetry-aware RMSD, the ranking statistics, salt
stripping, and the Lipinski filter. Anything needing network or a long
docking run is left out on purpose.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from src import analyze, dock, library, receptor  # noqa: E402


def test_box_covers_ligand() -> None:
    coords = np.array([[0.0, 0.0, 0.0], [10.0, 2.0, 2.0]])
    box = receptor.box_from_ligand(coords, padding=5.0)
    assert box.center == (5.0, 1.0, 1.0), box.center
    assert box.size == (20.0, 12.0, 12.0), box.size

    low = np.array(box.center) - np.array(box.size) / 2
    high = np.array(box.center) + np.array(box.size) / 2
    assert np.all(coords >= low) and np.all(coords <= high)


def test_largest_fragment_drops_counterion() -> None:
    from rdkit import Chem

    # Thioflavin T as deposited: cation plus a chloride counterion.
    salt = Chem.MolFromSmiles("CC1=CC2=C(C=C1)[N+](=C(S2)C3=CC=C(C=C3)N(C)C)C.[Cl-]")
    assert len(Chem.GetMolFrags(salt)) == 2
    stripped = dock.largest_fragment(salt)
    assert stripped.GetNumHeavyAtoms() == 20, stripped.GetNumHeavyAtoms()
    assert not any(a.GetSymbol() == "Cl" for a in stripped.GetAtoms())


def test_pose_rmsd_is_zero_against_itself_and_symmetry_aware() -> None:
    """A pose compared with itself must score 0, and the minimum over copies
    must pick the nearest copy rather than the first one."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1O"))
    params = AllChem.ETKDGv3()
    params.randomSeed = 7
    AllChem.EmbedMolecule(mol, params)
    mol = Chem.RemoveHs(mol)

    # A copy translated far away, plus the molecule itself.
    far = Chem.Mol(mol)
    conf = far.GetConformer()
    for i in range(far.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + 50.0, pos.y, pos.z))

    rmsds = dock.pose_rmsds(mol, [far, mol])
    assert len(rmsds) == 1
    assert rmsds[0] < 1e-6, f"self-RMSD should be 0, got {rmsds[0]}"

    only_far = dock.pose_rmsds(mol, [far])
    assert only_far[0] > 40.0, only_far


def test_ranking_stats_direction() -> None:
    """More negative scores for references must give AUC above 0.5."""
    good = pd.DataFrame(
        {
            "label": ["reference"] * 4 + ["decoy"] * 6,
            "best_score": [-9.0, -8.5, -8.2, -8.0, -6.0, -5.8, -5.5, -5.2, -5.0, -4.8],
            "heavy_atoms": [20] * 10,
        }
    )
    stats = analyze.ranking_stats(good, analyze.SCORE)
    assert stats["roc_auc"] == 1.0, stats
    assert stats["p_value"] < 0.01, stats

    flipped = good.copy()
    flipped["best_score"] = -flipped["best_score"]
    assert analyze.ranking_stats(flipped, analyze.SCORE)["roc_auc"] == 0.0


def test_ligand_efficiency_normalises_size() -> None:
    df = pd.DataFrame(
        {
            "label": ["reference", "decoy"],
            "best_score": [-10.0, -10.0],
            "heavy_atoms": [20, 40],
        }
    )
    out = analyze.add_efficiency(df)
    assert out[analyze.EFFICIENCY].tolist() == [-0.5, -0.25]


def test_enrichment_counts_and_baseline() -> None:
    df = pd.DataFrame(
        {
            "label": ["reference"] * 2 + ["decoy"] * 8,
            "best_score": [-9.0, -8.0] + [-5.0] * 8,
            "heavy_atoms": [20] * 10,
        }
    )
    top = analyze.enrichment(df, top_fractions=(0.2,)).iloc[0]
    assert top["n_screened"] == 2
    assert top["references_found"] == 2
    assert top["baseline_rate"] == 0.2
    assert top["enrichment_factor"] == 5.0


def test_descriptors_ignore_counterions() -> None:
    """Descriptors must be computed on the fragment that actually gets docked.

    Counting a chloride or two sodiums inflates the heavy-atom count, which is
    the denominator of ligand efficiency, so a salt would be silently
    under-credited in the very column used to argue about the size confound.
    """
    from rdkit import Chem

    tht_salt = "CC1=CC2=C(C=C1)[N+](=C(S2)C3=CC=C(C=C3)N(C)C)C.[Cl-]"
    described = library.describe(tht_salt)
    docked = dock.largest_fragment(Chem.MolFromSmiles(tht_salt))

    assert described["heavy_atoms"] == docked.GetNumHeavyAtoms() == 20, described
    # The salt has 21 heavy atoms; the cation that is docked has 20.
    assert Chem.MolFromSmiles(tht_salt).GetNumHeavyAtoms() == 21

    congo_red = (
        "C1=CC=C2C(=C1)C(=CC(=C2N)N=NC3=CC=C(C=C3)C4=CC=C(C=C4)N=NC5=C(C6=CC=CC=C6"
        "C(=C5)S(=O)(=O)[O-])N)S(=O)(=O)[O-].[Na+].[Na+]"
    )
    assert (
        library.describe(congo_red)["heavy_atoms"]
        == dock.largest_fragment(Chem.MolFromSmiles(congo_red)).GetNumHeavyAtoms()
    )


def test_lipinski_filter() -> None:
    aspirin = library.describe("CC(=O)OC1=CC=CC=C1C(=O)O")
    assert library.lipinski_violations(aspirin) == 0, aspirin

    # A deliberately oversized, greasy molecule should violate at least twice.
    huge = library.describe("C" * 60)
    assert library.lipinski_violations(huge) >= 2, huge


def test_mutation_guard_rejects_wrong_residue() -> None:
    """build_variants must refuse a sequence whose residues do not match.

    Imported lazily: esm_mutants pulls in torch, which lives in a separate
    environment, so this check is skipped where torch is absent.
    """
    sys.path.insert(0, str(Path(__file__).parent / "src"))
    try:
        import esm_mutants
    except ModuleNotFoundError as exc:
        if exc.name not in {"torch", "esm"}:
            raise  # a real broken import, not the expected missing dependency
        print("  (skipped mutation guard: torch not in this environment)")
        return

    assert len(esm_mutants.WILD_TYPE) == 140

    # Must scramble a position the guard actually inspects. Changing residue 1
    # leaves 30, 46 and 53 intact, so build_variants would happily accept it
    # and this test would pass while checking nothing.
    position = esm_mutants.MUTANTS["A30P"][0]
    scrambled = (
        esm_mutants.WILD_TYPE[: position - 1] + "G" + esm_mutants.WILD_TYPE[position:]
    )
    assert scrambled != esm_mutants.WILD_TYPE
    assert len(scrambled) == len(esm_mutants.WILD_TYPE)

    try:
        esm_mutants.build_variants(scrambled)
    except ValueError:
        pass  # expected: position 30 no longer holds alanine
    else:
        raise AssertionError("build_variants accepted a mismatched sequence")


def main() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test in tests:
        test()
        print(f"  ok  {test.__name__}")
    print(f"\n{len(tests)} checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
