# Alpha-Synuclein Drug Discovery (exploratory)

Notebooks (July–September 2024) exploring computational routes to compounds and peptides that could
bind **alpha-synuclein**, the protein that aggregates in Parkinson's disease. They combine protein
language models, ligand screening, molecular docking and ML scoring.

| Notebook | Approach |
|---|---|
| `esm_contact_maps.ipynb` | ESM-2 (3B) residue contact maps for wild-type vs. mutated sequences, physicochemical properties, and a small CNN over one-hot sequences |
| `esm_structure_viz.ipynb` | ESM3 structure prediction and 3-D visualization (py3Dmol) |
| `ligand_screening_docking.ipynb` | Lipinski rule-of-five filtering of a SMILES library (RDKit), descriptors, AutoDock Vina docking, MD trajectory analysis (mdtraj), random-forest and Keras models to rank candidates |
| `esm_docking_setup.ipynb` | ESM protein embeddings, generated SMILES, PDBQT preparation (Open Babel) and Vina docking |
| `peptide_docking.ipynb` | Peptide inhibitors: build from sequence (RDKit), dock against the receptor, random-forest regression on peptide features, generate new candidates |

**Status:** exploratory learning notebooks, with no validated hits. The peptide model is trained on docking scores, and none of the
models use experimental binding data.

## Setup

```bash
conda install -c conda-forge rdkit openbabel pymol-open-source
conda install -c bioconda autodock-vina
pip install fair-esm esm torch mdtraj scikit-learn tensorflow py3Dmol
```

The receptor structure (PDB) and the ligand library CSV are not included.
