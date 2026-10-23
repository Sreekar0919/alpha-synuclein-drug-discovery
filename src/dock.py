"""Ligand preparation and AutoDock Vina docking."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

RDLogger.DisableLog("rdApp.*")

EMBED_SEED = 0xF00D
VINA_SEED = 42
EXHAUSTIVENESS = 16
N_POSES = 9


def largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    """Drop counterions and keep the heaviest covalent fragment.

    Several library members are deposited as salts (Thioflavin T chloride,
    Congo red disodium); docking the counterion would be meaningless.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(frags, key=lambda f: f.GetNumHeavyAtoms()) if len(frags) > 1 else mol


def embed_3d(smiles: str, seed: int = EMBED_SEED) -> Chem.Mol:
    """Parse SMILES and generate one minimised 3-D conformer."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES: {smiles}")

    mol = Chem.AddHs(largest_fragment(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError(f"3-D embedding failed for: {smiles}")

    # -1 means MMFF has no parameters for this molecule, 1 means it ran out of
    # iterations. Either way the ligand would dock from a worse geometry, so
    # say so rather than letting it pass as a clean minimisation.
    status = AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    if status != 0:
        reason = "no MMFF parameters" if status == -1 else "did not converge"
        print(f"      note: MMFF minimisation {reason} for {smiles}")
    return mol


def to_pdbqt(mol: Chem.Mol) -> str:
    """Convert an embedded RDKit molecule to a PDBQT string via Meeko."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    setups = MoleculePreparation()(mol)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(f"Meeko could not write PDBQT: {err}")
    return pdbqt


def dock_one(
    vina,
    pdbqt: str,
    exhaustiveness: int = EXHAUSTIVENESS,
    n_poses: int = N_POSES,
) -> dict:
    """Dock a single prepared ligand into an already-mapped receptor."""
    vina.set_ligand_from_string(pdbqt)
    vina.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    energies = np.atleast_2d(vina.energies(n_poses=n_poses))
    return {
        "best_score": float(energies[0][0]),
        "scores": [float(row[0]) for row in energies],
        "poses_pdbqt": vina.poses(n_poses=n_poses),
    }


def make_vina(receptor_pdbqt: Path, box, cpu: int = 0, seed: int = VINA_SEED):
    """Create a Vina instance with the receptor loaded and maps computed."""
    from vina import Vina

    vina = Vina(sf_name="vina", cpu=cpu, seed=seed, verbosity=0)
    vina.set_receptor(str(receptor_pdbqt))
    vina.compute_vina_maps(center=list(box.center), box_size=list(box.size))
    return vina


def pose_to_mol(poses_pdbqt: str) -> Chem.Mol:
    """Turn Vina's ranked poses into one RDKit molecule with a conformer each.

    Parse errors are left to propagate. Swallowing them here only moved the
    failure somewhere less informative, where a Meeko API change surfaced as a
    missing dictionary key instead of the actual error.
    """
    from meeko import PDBQTMolecule, RDKitMolCreate

    pdbqt_mol = PDBQTMolecule(poses_pdbqt, is_dlg=False, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)

    mol = next((m for m in mols if m is not None), None)
    if mol is None:
        raise ValueError("Meeko returned no usable molecule for the docked poses")
    return Chem.RemoveHs(mol)


def load_crystal_copies(copy_paths: dict[str, Path], smiles: str) -> list[Chem.Mol]:
    """Load every deposited ligand copy with bond orders assigned.

    Every copy has to load. The redocking RMSD is a minimum over these, so a
    copy quietly dropped can only inflate the reported distance and could turn
    a correct pose into a spurious failure.
    """
    template = largest_fragment(Chem.MolFromSmiles(smiles))
    mols, failed = [], []
    for key, path in copy_paths.items():
        raw = Chem.MolFromPDBFile(str(path), removeHs=True, sanitize=False)
        if raw is None:
            failed.append(f"{key}: RDKit could not read {path.name}")
            continue
        try:
            mols.append(AllChem.AssignBondOrdersFromTemplate(template, raw))
        except ValueError as exc:
            failed.append(f"{key}: {exc}")

    if failed:
        raise ValueError(
            f"{len(failed)} of {len(copy_paths)} crystal ligand copies failed to load:\n  "
            + "\n  ".join(failed)
        )
    return mols


def pose_rmsds(docked: Chem.Mol, crystal_mols: list[Chem.Mol]) -> list[float]:
    """Per-pose heavy-atom RMSD to the closest crystallographic copy.

    Two things make a single RMSD number the wrong answer here. An amyloid
    fibril repeats the same binding site every ~4.8 A along its axis, so
    several equivalent copies sit inside any practical box and a pose one
    layer over is the same binding mode; the minimum over copies handles that.
    And Vina returns a ranked list, so the top pose and the best pose it found
    are different claims -- one value per pose keeps them apart.

    CalcRMS maps atoms through the molecular graph and does not superimpose,
    so internal symmetry is handled and the comparison stays in place.
    """
    out = []
    for i in range(docked.GetNumConformers()):
        single = Chem.Mol(docked)
        single.RemoveAllConformers()
        single.AddConformer(docked.GetConformer(i), assignId=True)
        out.append(min(float(rdMolAlign.CalcRMS(single, c)) for c in crystal_mols))
    return out


def redock_control(
    receptor_pdbqt: Path,
    box,
    smiles: str,
    crystal_mols: list[Chem.Mol],
    exhaustiveness: int = EXHAUSTIVENESS,
    n_poses: int = N_POSES,
) -> dict:
    """Redock the co-structure ligand and measure how close it gets.

    A protocol that cannot rebuild a pose someone has already measured has no
    business being trusted on compounds nobody has measured, so this runs
    before the screen and its result frames everything after it.
    """
    mol = embed_3d(smiles)
    vina = make_vina(receptor_pdbqt, box)
    result = dock_one(vina, to_pdbqt(mol), exhaustiveness=exhaustiveness, n_poses=n_poses)

    rmsds = pose_rmsds(pose_to_mol(result["poses_pdbqt"]), crystal_mols)
    best_index = int(min(range(len(rmsds)), key=lambda i: rmsds[i]))
    return {
        "best_score": result["best_score"],
        "scores": result["scores"],
        "top_pose_rmsd": round(rmsds[0], 2),
        "best_pose_rmsd": round(rmsds[best_index], 2),
        "best_pose_rank": best_index + 1,
        "best_pose_score": round(result["scores"][best_index], 2),
        "all_pose_rmsds": [round(r, 2) for r in rmsds],
        "poses_pdbqt": result["poses_pdbqt"],
        # The conventional pass mark for redocking is 2 A on the top pose.
        "passes_2A_top_pose": bool(rmsds[0] <= 2.0),
    }
