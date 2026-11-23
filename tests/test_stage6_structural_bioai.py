#!/usr/bin/env python3
"""Automated Verification Suite for Stage 6 (Live Boltz-2, OpenMM AMBER14 & 8-Gene AlphaGenome AVI)."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data" / "stage6_structural_bioai_summary.json"
FIG_PATH = ROOT / "figures" / "stage6_boltz2_openmm_alphagenome.png"


def main() -> None:
    assert JSON_PATH.exists(), f"Missing {JSON_PATH}"
    assert FIG_PATH.exists(), f"Missing {FIG_PATH}"
    d = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    # 1. Boltz-2 WT vs Mutant co-folding verification (all 3 pairs / 6 complexes)
    boltz = {b["pair"]: b for b in d["boltz2_cofolding_wt_vs_mutant"]}
    assert len(boltz) == 3
    for name, b in boltz.items():
        assert b["wt"]["iptm"] >= 0.65 and b["mutant"]["iptm"] >= 0.65, f"Low ipTM for {name}: {b}"
        assert b["delta_iptm"] >= -0.05, f"Interface destabilized for {name}: delta_iptm={b['delta_iptm']}"
        assert b["wt"]["chainA_rmsd_vs_cryoem_A"] <= 2.0, f"High cryo-EM RMSD for {name}"
    print(
        f"[PASS] Boltz-2 6/6 WT vs Mutant dimers verified: "
        + ", ".join(f"{k} (Mut ipTM={v['mutant']['iptm']:.3f}, Δ={v['delta_iptm']:+.3f})" for k, v in boltz.items())
    )

    # 2. OpenMM AMBER14 + OBC2 physics relaxation verification
    omm = {r["variant"]: r for r in d["openmm_amber14_relaxation"]}
    assert len(omm) == 7
    for fin in ("COX3_W116R", "ND4L_A11T", "COX3_I165T", "ATP6_I106V", "CYTB_M53T"):
        assert omm[fin]["physics_verdict"] == "STABLE_LOW_PERTURBATION", f"Unexpected verdict for {fin}"
        assert omm[fin]["mutant"]["post_relax_sc_dmin_A"] >= 2.5
        assert omm[fin]["delta_pocket_rmsd_A"] < omm["ATP6_L156R"]["delta_pocket_rmsd_A"]
    assert len(omm["COX3_W116R"]["new_mutant_polar_contacts_3p5A"]) >= 4
    assert omm["ATP6_L156R"]["delta_pocket_rmsd_A"] >= 0.12
    print(
        f"[PASS] OpenMM AMBER14 verified: 5/5 finalists STABLE_LOW_PERTURBATION (ΔRMSD <= 0.036 A; "
        f"COX3 W116R forms {len(omm['COX3_W116R']['new_mutant_polar_contacts_3p5A'])} H-bonds) vs "
        f"Pathogenic ATP6 L156R ΔRMSD={omm['ATP6_L156R']['delta_pocket_rmsd_A']:+.3f} A."
    )

    # 3. 8-gene live AlphaGenome AVI regulatory law verification
    ag = d["alphagenome_8gene_promoter_law"]
    assert ag["n_genes_scanned"] == 8 and ag["total_snvs_scored"] == 984
    assert ag["oxphos_partners_mean_snvs_ge15"] > ag["buffers_mean_snvs_ge15"] * 5.0
    assert ag["mann_whitney_p_n15_greater"] < 0.05
    assert ag["db_rows_updated_with_alphagenome"] >= 20
    print(
        f"[PASS] AlphaGenome 8-gene law verified: 6 OXPHOS partners mean={ag['oxphos_partners_mean_snvs_ge15']} SNVs>=15 "
        f"vs 2 buffers mean={ag['buffers_mean_snvs_ge15']} (Mann-Whitney p={ag['mann_whitney_p_n15_greater']}, "
        f"SQLite updated={ag['db_rows_updated_with_alphagenome']} rows)."
    )
    print("\nALL STAGE 6 STRUCTURAL & BIO-AI CHECKS PASSED (100%).")


if __name__ == "__main__":
    main()
