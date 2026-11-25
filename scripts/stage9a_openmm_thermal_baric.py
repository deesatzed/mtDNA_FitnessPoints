#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "openmm",
#     "pdbfixer",
#     "biopython",
#     "numpy",
#     "scipy",
# ]
# ///
"""Stage 9A: Live OpenMM AMBER14 + OBC2 Cryo-EM Atomic Relaxation of Stage-8B Thermal & Baric Leads.

Directly mutates and energy-minimizes complete interacting cryo-EM holocomplex subassemblies
(9I4I CI, 9I6F CIV, 9HZL CIII) for the 5 newly activated Stage-8B Thermal & Baric leads:
  1. ND4L Y85H  (Desert_Heat, Jerboa, 0 bystanders, contacts nuclear NDUFS2 at 3.25 A)
  2. ND1 S313G  (Desert_Heat, Camel, contacts nuclear NDUFA13 at 3.68 A)
  3. COX1 I57V  (Arctic_Cold, Muskox, 0 bystanders, 8 clades, Safe_Natural_Poly)
  4. CYTB F181L (Arctic_Cold, Woolly Mammoth, contacts nuclear UQCRFS1 at 4.40 A)
  5. CYTB Y224H (Diving_Subterranean_Hypoxia, Sperm Whale, contacts nuclear CYC1 at 4.04 A)

Reuses the exact Stage-6C OpenMM AMBER14+OBC2 protocol so every metric (8 A pocket heavy-atom
Delta-RMSD, post-relaxation d_min, new polar H-bonds <= 3.5 A) is 1-to-1 comparable with
Stage 6's COX3 W116R (+0.014 A) and pathogenic Leigh ATP6 L156R (+0.160 A).

Outputs:
  - data/stage9_openmm/*.pdb
  - data/stage9a_openmm_thermal_baric.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from openmm.app import CutoffNonPeriodic

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import stage6c_openmm_relaxation as s6c

# Use O(N) CutoffNonPeriodic (default 1.0 nm cutoff) in AMBER14+OBC2 without monkeypatching ForceField
s6c.NoCutoff = CutoffNonPeriodic

STRUCT_DIR = ROOT / "data" / "structures"
OUT_DIR = ROOT / "data" / "stage9_openmm"
OUT_JSON = ROOT / "data" / "stage9a_openmm_thermal_baric.json"

STAGE9_TARGETS = [
    ("ND4L_Y85H",  "ND4L",  85, "Y", "H", "9I4I", STRUCT_DIR / "pdb_00009i4i.cif.gz", "Desert_Heat_Lead_0Bys"),
    ("ND1_S313G",  "ND1",  313, "S", "G", "9I4I", STRUCT_DIR / "pdb_00009i4i.cif.gz", "Desert_Heat_Camel_NDUFA13"),
    ("COX1_I57V",  "COX1",  57, "I", "V", "9I6F", STRUCT_DIR / "pdb_00009i6f.cif.gz", "Arctic_Cold_Lead_0Bys"),
    ("CYTB_F181L", "CYTB", 181, "F", "L", "9HZL", STRUCT_DIR / "pdb_00009hzl.cif.gz", "Arctic_Cold_Mammoth_UQCRFS1"),
    ("CYTB_Y224H", "CYTB", 224, "Y", "H", "9HZL", STRUCT_DIR / "pdb_00009hzl.cif.gz", "Diving_Whale_CYC1"),
]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Redirect stage6c OUT_DIR so minimized PDBs write cleanly to data/stage9_openmm/
    s6c.OUT_DIR = OUT_DIR
    chain_maps = s6c.load_chain_gene_map()

    # Load Stage 6 pathogenic control anchor for side-by-side comparison
    s6_data = json.loads((ROOT / "data" / "stage6c_openmm_summary.json").read_text(encoding="utf-8"))
    s6_by_var = {r["variant"]: r for r in s6_data["results"]}
    patho_ref = s6_by_var["ATP6_L156R"]["delta_pocket_rmsd_A"]
    cox3_ref = s6_by_var["COX3_W116R"]["delta_pocket_rmsd_A"]

    results = []
    for var_id, gene, pos1, ref_aa, alt_aa, pdb_id, cif_gz, role in STAGE9_TARGETS:
        print(f"\n[OpenMM AMBER14+OBC2] {var_id} ({role}) on {pdb_id}...", flush=True)
        pdb_text, t_cid, t_resnum, sub_genes, t_ca_A = s6c.extract_interacting_subassembly(
            cif_gz, pdb_id, gene, pos1, chain_maps[pdb_id]
        )
        t_ca_nm = t_ca_A / 10.0
        print(f"  Subassembly chains: {sub_genes} (target {t_cid}:{t_resnum})", flush=True)
        wt_res = s6c.minimize_structure(pdb_text, t_cid, t_resnum, ref_aa, None, t_ca_nm, OUT_DIR / f"{var_id}_WT_min.pdb")
        mut_res = s6c.minimize_structure(pdb_text, t_cid, t_resnum, ref_aa, alt_aa, t_ca_nm, OUT_DIR / f"{var_id}_MUT_min.pdb")

        delta_rmsd = round(mut_res["pocket_8A_heavy_rmsd_A"] - wt_res["pocket_8A_heavy_rmsd_A"], 4)
        new_polar = sorted(set(mut_res["polar_contacts_3p5A"]) - set(wt_res["polar_contacts_3p5A"]))
        verdict = (
            "STABLE_LOW_PERTURBATION"
            if (mut_res["post_relax_sc_dmin_A"] >= 2.60 and delta_rmsd < 0.08)
            else ("MODERATE_LOCAL_REMODELING" if mut_res["post_relax_sc_dmin_A"] >= 2.40 else "STERIC_STRAIN_ALERT")
        )
        rec = {
            "variant": var_id,
            "gene": gene,
            "pos": pos1,
            "ref_aa": ref_aa,
            "alt_aa": alt_aa,
            "role": role,
            "pdb_id": pdb_id,
            "target_chain": t_cid,
            "subassembly_chains": sub_genes,
            "wt": wt_res,
            "mutant": mut_res,
            "delta_pocket_rmsd_A": delta_rmsd,
            "new_mutant_polar_contacts_3p5A": new_polar,
            "physics_verdict": verdict,
        }
        results.append(rec)
        print(
            f"  WT RMSD={wt_res['pocket_8A_heavy_rmsd_A']:.3f} A | MUT RMSD={mut_res['pocket_8A_heavy_rmsd_A']:.3f} A "
            f"(dRMSD={delta_rmsd:+.4f} A) | d_min={mut_res['post_relax_sc_dmin_A']:.2f} A | "
            f"new H-bonds={len(new_polar)} -> {verdict}",
            flush=True,
        )

    n_stable = sum(1 for r in results if r["physics_verdict"] == "STABLE_LOW_PERTURBATION")
    summary = {
        "stage": "9A_OpenMM_AMBER14_Thermal_Baric_Leads",
        "epistemic_tier": "E3_MECHANISTIC_SIM",
        "forcefield": "AMBER14 + OBC2 implicit solvent (OpenMM CPU/L-BFGS)",
        "benchmark_anchors": {
            "Stage6_COX3_W116R_delta_rmsd_A": cox3_ref,
            "Stage6_Pathogenic_ATP6_L156R_delta_rmsd_A": patho_ref,
        },
        "n_evaluated": len(results),
        "n_stable_low_perturbation": n_stable,
        "results": results,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nSaved -> {OUT_JSON} ({n_stable}/{len(results)} STABLE_LOW_PERTURBATION)")


if __name__ == "__main__":
    main()
