#!/usr/bin/env python
"""End-to-end alpha-synuclein docking run.

    python run_pipeline.py

Downloads the fibril structure, prepares the receptor, validates the setup by
redocking the co-structure ligand, then screens a small library against the
same site. Tables and figures are written to results/.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from src import analyze, dock, library, receptor  # noqa: E402

ROOT = Path(__file__).parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"

THT = "Thioflavin T"


def _vina_version() -> str:
    """Report the installed Vina version rather than hard-coding one."""
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version("vina")
    except PackageNotFoundError:
        return "unknown"


def run_controls(prep: dict, tht_smiles: str, exhaustiveness: int) -> dict:
    """Redock Thioflavin T into the site it was resolved in.

    Two receptors are tested. The apo receptor has every dye molecule stripped
    out and is what the screen uses. The holo receptor keeps the neighbouring
    dye copies and empties one slot, because Thioflavin T stacks along the
    fibril and part of each site is built from its own neighbours.
    """
    crystal = dock.load_crystal_copies(prep["reference_copy_paths"], tht_smiles)
    print(f"      {len(crystal)} crystallographic ThT copies loaded for comparison")

    controls = {}
    for name, key in (("apo", "receptor_pdbqt"), ("holo", "receptor_holo_pdbqt")):
        result = dock.redock_control(
            prep[key], prep["box"], tht_smiles, crystal, exhaustiveness=exhaustiveness
        )
        (RESULTS / f"tht_redock_{name}.pdbqt").write_text(result.pop("poses_pdbqt", ""))
        controls[name] = result
        print(
            f"      {name:<5} score {result['best_score']:>6.2f}  "
            f"top-pose RMSD {result['top_pose_rmsd']:>5.2f} A  "
            f"best-of-{len(result['all_pose_rmsds'])} RMSD {result['best_pose_rmsd']:>5.2f} A "
            f"(rank {result['best_pose_rank']})"
        )
    return controls


def write_analysis(df: pd.DataFrame, summary: dict) -> dict:
    """Recompute statistics and figures for a scored library.

    Shared by the full run and by --reanalyze so the two cannot drift apart
    and report different numbers from the same scores.
    """
    df = analyze.add_efficiency(df).sort_values(analyze.SCORE).reset_index(drop=True)
    df.to_csv(RESULTS / "docking_results.csv", index=False)

    score_stats = analyze.ranking_stats(df, analyze.SCORE)
    eff_stats = analyze.ranking_stats(df, analyze.EFFICIENCY)
    enrich = analyze.enrichment(df)
    enrich.to_csv(RESULTS / "enrichment.csv", index=False)

    summary = {
        **summary,
        "ranking_by_score": score_stats,
        "ranking_by_ligand_efficiency": eff_stats,
        "enrichment": enrich.to_dict(orient="records"),
        "best_overall": df.iloc[0][["name", "label", "best_score"]].to_dict(),
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2))

    analyze.plot_distributions(df, FIGURES / "reference_vs_decoy.png")
    analyze.plot_ranked_library(df, FIGURES / "ranked_library.png")
    analyze.plot_size_confound(df, FIGURES / "size_confound.png")
    return summary


def reanalyze() -> int:
    """Redo descriptors, statistics and figures from the saved docking scores.

    Descriptors are recomputed from SMILES rather than reused, so a fix to the
    descriptor code reaches the saved results without a 50-minute redock. The
    docking scores themselves are untouched.
    """
    df = pd.read_csv(RESULTS / "docking_results.csv")
    summary = json.loads((RESULTS / "summary.json").read_text())

    keep = ["name", "label", "note", "smiles", "formula", "best_score", "pose_score_spread"]
    rebuilt = []
    for _, row in df.iterrows():
        desc = library.describe(row["smiles"])
        violations = library.lipinski_violations(desc)
        rebuilt.append({
            **{k: row[k] for k in keep if k in df.columns},
            **desc,
            "lipinski_violations": violations,
            "passes_lipinski": violations <= 1,
        })

    summary = write_analysis(pd.DataFrame(rebuilt), summary)
    print(f"reanalysed {len(rebuilt)} compounds from saved docking scores")
    print(f"  ROC AUC (score):             {summary['ranking_by_score']['roc_auc']}")
    print(f"  ROC AUC (ligand efficiency): {summary['ranking_by_ligand_efficiency']['roc_auc']}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exhaustiveness", type=int, default=dock.EXHAUSTIVENESS)
    parser.add_argument("--limit", type=int, default=None,
                        help="dock only the first N compounds (smoke test)")
    parser.add_argument("--skip-controls", action="store_true")
    parser.add_argument("--reanalyze", action="store_true",
                        help="recompute descriptors, statistics and figures from the "
                             "docking scores already in results/, without redocking")
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)

    if args.reanalyze:
        return reanalyze()

    print("[1/5] preparing receptor ...")
    prep = receptor.prepare(DATA)
    box = prep["box"]
    print(f"      {prep['pdb_id']}: {prep['n_ligand_copies']} copies of "
          f"{receptor.LIGAND_CODE}; box built on copy {prep['reference_copy']}")
    print(f"      box center {box.center}, size {box.size}")
    print(f"      {prep['n_copies_in_box']} symmetry-equivalent sites fall inside the box")
    print(f"      site residues: {', '.join(prep['site_residues'])}")

    print("[2/5] building compound library ...")
    rows = library.build(DATA / "smiles_cache.json")
    tht_smiles = next(r["smiles"] for r in rows if r["name"] == THT)
    if args.limit:
        # References come first in the library, so a plain head() slice would
        # leave no decoys and the run would dock everything and only then fail
        # in the statistics. Interleave so a smoke test keeps both groups.
        refs = [r for r in rows if r["label"] == analyze.REFERENCE]
        decoys = [r for r in rows if r["label"] == analyze.DECOY]
        interleaved = [r for pair in zip(refs, decoys) for r in pair]
        interleaved += refs[len(decoys):] + decoys[len(refs):]
        rows = interleaved[: args.limit]
    n_ref = sum(r["label"] == analyze.REFERENCE for r in rows)
    print(f"      {len(rows)} compounds ({n_ref} reference, {len(rows) - n_ref} decoy)")

    print("[3/5] redocking control ...")
    controls = {} if args.skip_controls else run_controls(prep, tht_smiles, args.exhaustiveness)

    print("[4/5] screening ...")
    vina = dock.make_vina(prep["receptor_pdbqt"], box)
    records, failures = [], []
    started = time.time()

    for i, row in enumerate(rows, start=1):
        try:
            mol = dock.embed_3d(row["smiles"])
            result = dock.dock_one(vina, dock.to_pdbqt(mol), exhaustiveness=args.exhaustiveness)
        except Exception as exc:  # one bad ligand must not kill the screen
            failures.append({"name": row["name"], "error": f"{type(exc).__name__}: {exc}"})
            print(f"      [{i:>2}/{len(rows)}] {row['name']:<28} FAILED ({type(exc).__name__})")
            continue

        records.append({
            **row,
            "best_score": result["best_score"],
            "pose_score_spread": round(max(result["scores"]) - min(result["scores"]), 3),
        })
        print(f"      [{i:>2}/{len(rows)}] {row['name']:<28} {result['best_score']:>7.2f}")

    elapsed = time.time() - started
    print(f"      docked {len(records)}/{len(rows)} in {elapsed/60:.1f} min")

    # Without this, a systematic breakage (a Meeko or Vina API change, say)
    # fails every ligand and then surfaces as a missing-column error in pandas.
    if not records:
        raise SystemExit(
            f"every one of the {len(rows)} ligands failed to dock; first error: "
            f"{failures[0]['error'] if failures else 'unknown'}"
        )

    print("[5/5] analysing ...")
    summary = {
        "receptor": {
            "pdb_id": prep["pdb_id"],
            "ligand_code": receptor.LIGAND_CODE,
            "box_built_on_copy": prep["reference_copy"],
            "n_ligand_copies": prep["n_ligand_copies"],
            "n_symmetry_copies_in_box": prep["n_copies_in_box"],
            "box_center": list(box.center),
            "box_size": list(box.size),
            "site_residues": prep["site_residues"],
        },
        "docking": {
            "engine": f"AutoDock Vina {_vina_version()} (python API)",
            "exhaustiveness": args.exhaustiveness,
            "n_poses": dock.N_POSES,
            "seed": dock.VINA_SEED,
            "n_docked": len(records),
            "n_failed": len(failures),
            "failures": failures,
            "runtime_minutes": round(elapsed / 60, 2),
        },
        "redock_control": controls,
    }
    summary = write_analysis(pd.DataFrame(records), summary)
    score_stats = summary["ranking_by_score"]
    eff_stats = summary["ranking_by_ligand_efficiency"]

    print("\n--- result ---")
    if controls:
        apo = controls["apo"]
        print(f"redock (apo, screen setup):  top pose {apo['top_pose_rmsd']} A, "
              f"best of {len(apo['all_pose_rmsds'])} {apo['best_pose_rmsd']} A "
              f"-> {'PASS' if apo['passes_2A_top_pose'] else 'FAIL'} at the 2 A mark")
    if failures:
        print(f"ligands that failed to dock: {len(failures)} "
              f"({', '.join(f['name'] for f in failures)})")
    print(f"reference median score:      {score_stats['reference_median']} kcal/mol")
    print(f"decoy median score:          {score_stats['decoy_median']} kcal/mol")
    print(f"ROC AUC (score):             {score_stats['roc_auc']}  p = {score_stats['p_value']:.4f}")
    print(f"ROC AUC (ligand efficiency): {eff_stats['roc_auc']}  p = {eff_stats['p_value']:.4f}")
    print(f"wrote {RESULTS}/docking_results.csv, summary.json and 3 figures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
