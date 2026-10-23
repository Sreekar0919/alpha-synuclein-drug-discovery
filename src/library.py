"""Compound library: reference alpha-synuclein interactors plus decoys.

SMILES are not written by hand. Each compound is resolved by name against
PubChem and the returned molecular formula is checked against the formula
recorded here, so a bad name lookup fails loudly instead of silently docking
the wrong molecule.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

PUBCHEM = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}"
    "/property/SMILES,ConnectivitySMILES,MolecularFormula/JSON"
)

LIPINSKI = {"mw_max": 500.0, "logp_max": 5.0, "hbd_max": 5, "hba_max": 10}


@dataclass(frozen=True)
class Compound:
    """One library member before any RDKit work."""

    name: str
    formula: str
    label: str  # "reference" or "decoy"
    note: str = ""


# Compounds with published evidence of binding alpha-synuclein or modulating
# its aggregation. These are the comparison set: a docking protocol with any
# discriminative power should rank them above arbitrary drugs.
REFERENCES = [
    Compound("Thioflavin T", "C17H19ClN2S", "reference",
             "resolved in the 7YNM binding channel; the redocking control"),
    Compound("anle138b", "C16H11BrN2O2", "reference",
             "diphenylpyrazole aggregation inhibitor, clinical trials in PD/MSA"),
    Compound("epigallocatechin gallate", "C22H18O11", "reference",
             "redirects alpha-synuclein into off-pathway oligomers"),
    Compound("SynuClean-D", "C13H5F3N4O5", "reference",
             "small-molecule aggregation inhibitor (Pujols et al., PNAS 2018)"),
    Compound("curcumin", "C21H20O6", "reference",
             "binds alpha-synuclein and slows fibrillation in vitro"),
    Compound("baicalein", "C15H10O5", "reference",
             "flavonoid inhibitor of alpha-synuclein fibrillation"),
    Compound("fasudil", "C14H17N3O2S", "reference",
             "binds the C-terminal Tyr residues by NMR"),
    Compound("Congo red", "C32H22N6Na2O6S2", "reference",
             "classical amyloid-binding dye"),
]

# Approved drugs from unrelated therapeutic areas. These are *presumed*
# non-binders: absence of a published alpha-synuclein interaction is not the
# same as a measured negative, which is stated as a limitation in the README.
DECOY_NAMES = [
    ("atenolol", "C14H22N2O3"), ("metformin", "C4H11N5"),
    ("ibuprofen", "C13H18O2"), ("naproxen", "C14H14O3"),
    ("aspirin", "C9H8O4"), ("acetaminophen", "C8H9NO2"),
    ("amoxicillin", "C16H19N3O5S"), ("ciprofloxacin", "C17H18FN3O3"),
    ("loratadine", "C22H23ClN2O2"), ("cetirizine", "C21H25ClN2O3"),
    ("omeprazole", "C17H19N3O3S"), ("ranitidine", "C13H22N4O3S"),
    ("simvastatin", "C25H38O5"), ("atorvastatin", "C33H35FN2O5"),
    ("lisinopril", "C21H31N3O5"), ("losartan", "C22H23ClN6O"),
    ("amlodipine", "C20H25ClN2O5"), ("metoprolol", "C15H25NO3"),
    ("furosemide", "C12H11ClN2O5S"), ("hydrochlorothiazide", "C7H8ClN3O4S2"),
    ("warfarin", "C19H16O4"), ("clopidogrel", "C16H16ClNO2S"),
    ("sertraline", "C17H17Cl2N"), ("fluoxetine", "C17H18F3NO"),
    ("citalopram", "C20H21FN2O"), ("diazepam", "C16H13ClN2O"),
    ("lorazepam", "C15H10Cl2N2O2"), ("haloperidol", "C21H23ClFNO2"),
    ("risperidone", "C23H27FN4O2"), ("olanzapine", "C17H20N4S"),
    ("gabapentin", "C9H17NO2"), ("pregabalin", "C8H17NO2"),
    ("tramadol", "C16H25NO2"), ("prednisone", "C21H26O5"),
    ("dexamethasone", "C22H29FO5"), ("albuterol", "C13H21NO3"),
    ("montelukast", "C35H36ClNO3S"), ("sildenafil", "C22H30N6O4S"),
    ("tadalafil", "C22H19N3O4"), ("caffeine", "C8H10N4O2"),
    ("nicotine", "C10H14N2"), ("theophylline", "C7H8N4O2"),
    ("allopurinol", "C5H4N4O"), ("colchicine", "C22H25NO6"),
    ("methotrexate", "C20H22N8O5"), ("doxycycline", "C22H24N2O8"),
    ("propranolol", "C16H21NO2"), ("verapamil", "C27H38N2O4"),
    ("ketoconazole", "C26H28Cl2N4O4"), ("indomethacin", "C19H16ClNO4"),
]

DECOYS = [Compound(name, formula, "decoy") for name, formula in DECOY_NAMES]


def fetch_smiles(name: str, timeout: int = 30) -> dict:
    """Resolve a compound name to SMILES and formula via PubChem."""
    url = PUBCHEM.format(name=urllib.parse.quote(name))
    with urllib.request.urlopen(url, timeout=timeout) as handle:  # noqa: S310
        payload = json.load(handle)
    props = payload["PropertyTable"]["Properties"][0]
    smiles = props.get("SMILES") or props.get("ConnectivitySMILES")
    if not smiles:
        raise ValueError(f"PubChem returned no SMILES for {name!r}")
    return {"smiles": smiles, "formula": props["MolecularFormula"]}


def resolve(compounds: list[Compound], cache_path: Path, pause: float = 0.3) -> list[dict]:
    """Resolve every compound, verifying the formula PubChem returns.

    Results are cached to disk so reruns do not hammer PubChem and so the
    docking run is reproducible even if a name lookup later changes.
    """
    cache: dict = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    resolved, mismatches = [], []

    for compound in compounds:
        if compound.name not in cache:
            cache[compound.name] = fetch_smiles(compound.name)
            time.sleep(pause)

        entry = cache[compound.name]
        if entry["formula"] != compound.formula:
            mismatches.append(
                f"{compound.name}: expected {compound.formula}, PubChem gave {entry['formula']}"
            )
            continue

        resolved.append(
            {
                "name": compound.name,
                "label": compound.label,
                "note": compound.note,
                "smiles": entry["smiles"],
                "formula": entry["formula"],
            }
        )

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True))

    if mismatches:
        raise ValueError("PubChem formula mismatch:\n  " + "\n  ".join(mismatches))
    return resolved


def desalt(mol):
    """Keep the heaviest covalent fragment, dropping counterions.

    Several library members are deposited as salts (Thioflavin T chloride,
    Congo red disodium). Only the heaviest fragment is docked, so descriptors
    have to be computed on that same fragment or the heavy-atom count used as
    the ligand-efficiency denominator counts atoms that were never docked.
    """
    from rdkit import Chem

    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(frags, key=lambda f: f.GetNumHeavyAtoms()) if len(frags) > 1 else mol


def describe(smiles: str) -> dict:
    """Compute the descriptors used for Lipinski filtering.

    Computed on the desalted molecule, matching what actually gets docked.
    """
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    mol = desalt(mol)

    return {
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Crippen.MolLogP(mol), 2),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
    }


def lipinski_violations(desc: dict) -> int:
    """Count Lipinski rule-of-five violations."""
    return sum(
        [
            desc["mw"] > LIPINSKI["mw_max"],
            desc["logp"] > LIPINSKI["logp_max"],
            desc["hbd"] > LIPINSKI["hbd_max"],
            desc["hba"] > LIPINSKI["hba_max"],
        ]
    )


def build(cache_path: Path, max_violations: int = 1) -> list[dict]:
    """Resolve references plus decoys and annotate with descriptors.

    Nothing is dropped here: the Lipinski verdict is recorded per compound so
    the notebook can show what the filter would remove and why.
    """
    rows = []
    for entry in resolve(REFERENCES + DECOYS, cache_path):
        desc = describe(entry["smiles"])
        violations = lipinski_violations(desc)
        rows.append(
            {
                **entry,
                **desc,
                "lipinski_violations": violations,
                "passes_lipinski": violations <= max_violations,
            }
        )
    return rows
