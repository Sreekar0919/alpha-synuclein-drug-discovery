# Archived notebooks — these never ran successfully

These five notebooks are the original 2024 exploration that this project grew
out of. They are kept for history. **None of them ever ran end to end**, and
their outputs were cleared, so nothing in them should be read as a result.

What was wrong with them:

| Notebook | Why it did not run |
|---|---|
| `esm_contact_maps.ipynb` | Imported `esm.pretrained` / `FastaBatchedDataset` from the wrong package. PyPI has two projects called `esm`: EvolutionaryScale's `esm` and Meta's `fair-esm`. The notebook installed the first and called the second's API. |
| `esm_structure_viz.ipynb` | Used ESM3, which is gated behind an access-approved Hugging Face token; the call returned HTTP 401. |
| `ligand_screening_docking.ipynb` | `pip install vina` fails to build on macOS arm64; the conda install that would have worked was cancelled. The docking cells never executed. |
| `esm_docking_setup.ipynb` | Same Vina install failure, plus the same wrong-`esm` import. |
| `peptide_docking.ipynb` | Same Vina failure. The random-forest model in it is also trained on docking scores as labels, so even had it run it would only have learned to imitate Vina, not to predict binding. |

The working pipeline that replaced them lives in [`src/`](../../src) and
[`run_pipeline.py`](../../run_pipeline.py), with results in
[`results/`](../../results) and a walkthrough in
[`notebooks/results.ipynb`](../results.ipynb).
