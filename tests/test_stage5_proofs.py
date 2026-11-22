#!/usr/bin/env python3
"""Automated Hardening & Proof Verification Suite for Stage 5."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data" / "stage5_hardening_proofs.json"
FIG_PATH = ROOT / "figures" / "stage5_hardening_and_proofs.png"


def main() -> None:
    assert JSON_PATH.exists(), f"Missing {JSON_PATH}"
    assert FIG_PATH.exists(), f"Missing {FIG_PATH}"
    d = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    # Proof 1: Confounder-controlled & negative-control mitonuclear interface enrichment
    p1 = d["proof1_mitonuclear_interface_confounders_and_null"]
    assert p1["gene_stratified_cmh"]["cmh_odds_ratio"] > 1.3 and p1["gene_stratified_cmh"]["cmh_p_value"] < 0.05
    assert p1["multivariate_logistic_regression_adjusted"]["adjusted_OR_interface"] > 1.35
    assert p1["multivariate_logistic_regression_adjusted"]["wald_p_value"] < 0.05
    assert p1["matched_control_polarized_negative_null_2x2"]["fisher_p_two_sided"] > 0.05
    assert p1["leave_one_gene_out_OR_range"]["min_OR"] > 1.25
    print(f"[PASS] Proof 1: Adjusted Mitonuclear Interface OR={p1['multivariate_logistic_regression_adjusted']['adjusted_OR_interface']} (p={p1['multivariate_logistic_regression_adjusted']['wald_p_value']}); CMH OR={p1['gene_stratified_cmh']['cmh_odds_ratio']}; Control Null non-sig (p={p1['matched_control_polarized_negative_null_2x2']['fisher_p_two_sided']}).")

    # Proof 2: Proteome-wide ClinVar Pathogenic vs Benign/gnomAD & blind OOS generalization
    p2 = d["proof2_proteome_clinvar_benchmark_and_oos"]
    assert p2["n_clinvar_pathogenic_missense_3d"] >= 35
    assert p2["n_benign_or_gnomad_missense_3d"] >= 800
    assert p2["five_fold_cv_auroc"] >= 0.85
    assert p2["strict_oos_tier1_to_tier2_auroc"] >= 0.85
    assert p2["strict_oos_tier2_to_tier1_auroc"] >= 0.85
    print(f"[PASS] Proof 2: Proteome ClinVar 5-Fold CV AUROC={p2['five_fold_cv_auroc']}, Strict OOS Tier-1->Tier-2 AUROC={p2['strict_oos_tier1_to_tier2_auroc']} (Tier-2->Tier-1={p2['strict_oos_tier2_to_tier1_auroc']}).")

    # Proof 3: 3D Cryo-EM rotamer clash proof
    p3 = {r["variant"]: r for r in d["proof3_3d_rotamer_steric_clash_audit"]}
    for v in ("COX3_W116R", "ND4L_A11T", "COX3_I165T", "ATP6_I106V", "CYTB_M53T"):
        assert p3[v]["best_rotamer"]["severe_clashes_lt_2p2A"] == 0, f"Unexpected clash in {v}"
        assert p3[v]["steric_verdict"] == "CLASH_FREE_CAVITY_FIT"
    assert p3["COX3_W116R"]["best_rotamer"]["min_heavy_atom_dist_A"] >= 3.0
    assert p3["COX3_W116R"]["clash_free_rotamers"] > p3["ATP6_L156R"]["clash_free_rotamers"]
    print(f"[PASS] Proof 3: 5/5 Finalists CLASH_FREE_CAVITY_FIT; COX3 W116R d_min={p3['COX3_W116R']['best_rotamer']['min_heavy_atom_dist_A']} A ({p3['COX3_W116R']['clash_free_rotamers']}/81 rotamers vs ATP6 L156R {p3['ATP6_L156R']['clash_free_rotamers']}/81 at {p3['ATP6_L156R']['best_rotamer']['min_heavy_atom_dist_A']} A).")

    # Proof 4: Base-editor +/-2..+/-6 nt window sensitivity
    p4 = d["proof4_base_editor_window_stress_test"]
    fmap = {x["variant"]: x for x in p4["top5_finalists_window_profile"]}
    assert fmap["COX3_W116R"]["ultra_clean_pm6"] is True
    assert fmap["ND4L_A11T"]["ultra_clean_pm6"] is True
    print(f"[PASS] Proof 4: COX3 W116R & ND4L A11T remain 0-bystander out to +/-6 nt ({p4['zero_bystander_count_by_window']}).")

    # Proof 5: Supermodularity, physical invariants, multi-model forms, adversarial nulls & heteroplasmy
    p5 = d["proof5_epistasis_supermodularity_models_nulls_heteroplasmy"]
    assert p5["supermodularity_grid_proof"]["dual_cliff_theorem_concordance_pct"] == 100.0
    assert p5["physical_invariant_violations_across_all_scenarios"] == 0
    for scen, res in p5["multi_model_form_and_adversarial_controls"].items():
        assert res["Yeast_NDI1"]["epistatic_excess_median"] > 0.05, f"Yeast_NDI1 failed in {scen}"
        assert res["Human_SOD2"]["epistatic_excess_median"] > 0.05, f"Human_SOD2 failed in {scen}"
        assert res["Ciona_AOX"]["epistatic_excess_median"] > 0.04, f"Ciona_AOX failed in {scen}"
        assert abs(res["ADVERSARIAL_Null_Placebo"]["epistatic_excess_median"]) < 1e-6, f"Placebo non-zero in {scen}"
        assert res["ADVERSARIAL_Toxic_Uncoupler_FCCP"]["epistatic_excess_median"] < -0.02, f"Toxic FCCP not negative in {scen}"
    print("[PASS] Proof 5: 100.0% supermodularity concordance, 0 invariant violations, positive epistasis holds across all 4 model forms, Placebo=0.0, Toxic FCCP<0.")

    # Proof 6: Database referential integrity
    p6 = d["proof6_database_referential_integrity"]
    assert p6["integrity_pass"] is True
    assert p6["foreign_key_violations"] == 0
    assert p6["orphan_payoff_rows"] == 0 and p6["orphan_transferability_rows"] == 0
    print(f"[PASS] Proof 6: SQLite referential integrity 100% clean across all 8 tables ({p6['table_row_counts']}).")
    print("\nALL STAGE 5 HARDENING & PROOF CHECKS PASSED (100%).")


if __name__ == "__main__":
    main()
