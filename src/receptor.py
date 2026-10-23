"""Receptor preparation for alpha-synuclein fibril docking.

Uses PDB 7YNM: a cryo-EM alpha-synuclein fibril with Thioflavin-T (ligand code
TFX) resolved in the binding channel. Picking a structure that already has a
small molecule in it means the docking box is placed on an experimentally
observed site rather than a guess, and it gives us a redocking control.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import requests

PDB_ID = "7YNM"
LIGAND_CODE = "TFX"  # Thioflavin T, as deposited in 7YNM
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Vina needs a box big enough to let a drug-like ligand rotate, but small
# enough that the search stays on the site we chose. A fibril cannot be made
# to hold only one site: equivalent sites repeat every ~4.8 A along the axis
# (the beta-strand rise), so any usable box contains several of them. That is
# handled when scoring poses, not here -- see dock.redock_rmsd.
BOX_PADDING = 8.0


@dataclass(frozen=True)
class Box:
    """Docking box in Angstroms."""

    center: tuple[float, float, float]
    size: tuple[float, float, float]


def fetch_pdb(pdb_id: str, dest_dir: Path) -> Path:
    """Download a PDB entry, skipping the request if it is already on disk."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / f"{pdb_id.lower()}.pdb"
    if path.exists():
        return path

    response = requests.get(RCSB_URL.format(pdb_id=pdb_id), timeout=60)
    response.raise_for_status()
    path.write_text(response.text)
    return path


def _parse_coords(lines: list[str]) -> np.ndarray:
    return np.array(
        [[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])] for ln in lines]
    )


def read_reference_ligand(pdb_path: Path, ligand_code: str = LIGAND_CODE) -> dict:
    """Return the copy of the co-structure ligand nearest the middle of the stack.

    A fibril structure contains one ligand copy per beta-strand layer. The
    copies at the ends of the stack are missing neighbours on one side, so the
    central copy is the one with a complete binding environment.
    """
    het = [
        ln
        for ln in pdb_path.read_text().splitlines()
        if ln.startswith("HETATM") and ln[17:20].strip() == ligand_code
    ]
    if not het:
        raise ValueError(f"no {ligand_code} atoms found in {pdb_path}")

    copies: dict[str, list[str]] = {}
    for ln in het:
        copies.setdefault(f"{ln[21]}{ln[22:26].strip()}", []).append(ln)

    centroids = {key: _parse_coords(lines).mean(axis=0) for key, lines in copies.items()}
    stack_center = np.mean(list(centroids.values()), axis=0)
    best = min(centroids, key=lambda k: np.linalg.norm(centroids[k] - stack_center))

    return {
        "copy_id": best,
        "n_copies": len(copies),
        "lines": copies[best],
        "coords": _parse_coords(copies[best]),
        "centroid": centroids[best],
        # Every copy is related to every other by fibril symmetry, so a docked
        # pose that lands on any of them has reproduced the observed mode.
        "all_copies": {key: _parse_coords(lines) for key, lines in copies.items()},
        "all_copy_lines": copies,
    }


def box_from_ligand(coords: np.ndarray, padding: float = BOX_PADDING) -> Box:
    """Build a docking box around a set of ligand coordinates."""
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)
    size = span + 2 * padding
    return Box(
        center=tuple(round(float(v), 3) for v in center),
        size=tuple(round(float(v), 3) for v in size),
    )


def site_residues(pdb_path: Path, coords: np.ndarray, cutoff: float = 5.0) -> list[str]:
    """List protein residues with any atom within `cutoff` of the ligand."""
    hits: dict[str, None] = {}
    for ln in pdb_path.read_text().splitlines():
        if not ln.startswith("ATOM"):
            continue
        atom = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        if np.min(np.linalg.norm(coords - atom, axis=1)) <= cutoff:
            hits[f"{ln[17:20].strip()}{ln[22:26].strip()}"] = None
    return sorted(hits, key=lambda r: int(r[3:]))


def write_apo_receptor(pdb_path: Path, out_pdb: Path) -> Path:
    """Strip ligands, waters and ions, keeping only protein atoms."""
    kept = [ln for ln in pdb_path.read_text().splitlines() if ln.startswith("ATOM")]
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(kept) + "\nEND\n")
    return out_pdb


def write_holo_receptor(
    pdb_path: Path, out_pdb: Path, copy_lines: dict, target_copy: str
) -> Path:
    """Protein plus every co-structure ligand copy except one.

    Thioflavin-T binds amyloid as a stack of molecules spaced one beta-strand
    apart, so part of the wall of any one site is made of the neighbouring dye
    molecules rather than protein. Removing all of them, as the apo receptor
    does, deletes real binding surface. This variant leaves the neighbours in
    place and empties a single slot, which is the fairer test of whether the
    docking setup can rebuild the observed pose.
    """
    kept = [ln for ln in pdb_path.read_text().splitlines() if ln.startswith("ATOM")]
    for key, lines in copy_lines.items():
        if key != target_copy:
            kept.extend(lines)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(kept) + "\nEND\n")
    return out_pdb


def pdb_to_pdbqt(in_pdb: Path, out_pdbqt: Path, obabel: str = "obabel") -> Path:
    """Convert a receptor PDB to PDBQT (adds hydrogens at pH 7.4 and charges)."""
    result = subprocess.run(
        [obabel, str(in_pdb), "-O", str(out_pdbqt), "-xr", "-p", "7.4",
         "--partialcharge", "gasteiger"],
        capture_output=True,
        text=True,
    )
    # check=True alone would raise without showing why Open Babel refused.
    if result.returncode != 0:
        raise RuntimeError(
            f"obabel failed on {in_pdb.name} (exit {result.returncode}):\n{result.stderr}"
        )
    return out_pdbqt


def prepare(data_dir: Path, pdb_id: str = PDB_ID) -> dict:
    """Download the fibril, derive the box from the co-structure ligand, and
    write a receptor PDBQT with that ligand removed."""
    pdb_path = fetch_pdb(pdb_id, data_dir)
    ref = read_reference_ligand(pdb_path)
    box = box_from_ligand(ref["coords"])

    apo_pdb = write_apo_receptor(pdb_path, data_dir / "receptor.pdb")
    receptor_pdbqt = pdb_to_pdbqt(apo_pdb, data_dir / "receptor.pdbqt")

    holo_pdb = write_holo_receptor(
        pdb_path, data_dir / "receptor_holo.pdb", ref["all_copy_lines"], ref["copy_id"]
    )
    holo_pdbqt = pdb_to_pdbqt(holo_pdb, data_dir / "receptor_holo.pdbqt")

    ref_pdb = data_dir / "reference_ligand.pdb"
    ref_pdb.write_text("\n".join(ref["lines"]) + "\nEND\n")

    # One file per symmetry copy, so the redocking control can be scored
    # against all of them.
    copies_dir = data_dir / "reference_copies"
    copies_dir.mkdir(exist_ok=True)
    copy_paths = {}
    for key, lines in ref["all_copy_lines"].items():
        path = copies_dir / f"{key}.pdb"
        path.write_text("\n".join(lines) + "\nEND\n")
        copy_paths[key] = path

    return {
        "pdb_id": pdb_id,
        "pdb_path": pdb_path,
        "receptor_pdbqt": receptor_pdbqt,
        "receptor_holo_pdbqt": holo_pdbqt,
        "reference_ligand_pdb": ref_pdb,
        "reference_copy_paths": copy_paths,
        "reference_coords": ref["coords"],
        "reference_copy": ref["copy_id"],
        "n_ligand_copies": ref["n_copies"],
        "box": box,
        "site_residues": site_residues(pdb_path, ref["coords"]),
        "n_copies_in_box": _copies_in_box(ref["all_copies"], box),
    }


def _copies_in_box(all_copies: dict, box: Box) -> int:
    """Count symmetry copies whose centroid falls inside the docking box."""
    low = np.array(box.center) - np.array(box.size) / 2
    high = np.array(box.center) + np.array(box.size) / 2
    return sum(
        bool(np.all(coords.mean(axis=0) >= low) and np.all(coords.mean(axis=0) <= high))
        for coords in all_copies.values()
    )
