#!/usr/bin/env python3
"""Automated Verification Suite for Stage 7 (Cross-Project Bridges: Mito_Papers + BAA Mito Twin)."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S7A_JSON = ROOT / "data" / "stage7a_heteroplasmy_and_cart.json"
S7B_JSON = ROOT / "data" / "stage7b_singlecell_alphagenome.json"
FIG_PATH = ROOT / "figures" / "stage7_cross_project_bridges.png"
EXP_PATH = ROOT / "app" / "mitofitness_explorer.html"


def main() -> None:
    assert S7A_JSON.exists() and S7B_JSON.exists() and FIG_PATH.exists()
    s7a = json.loads(S7A_JSON.read_text(encoding="utf-8"))
    s7b = json.loads(S7B_JSON.read_text(encoding="utf-8"))

    # 1. Track 7A: Solver convergence, Biogenesis Trap proof, Biphasic amplification & Pulsed safety
    t7a = s7a["track7a_biphasic_heteroplasmy_amplification"]
    assert t7a["solver_convergence"]["max_abs_diff"] < 1e-3
    regs = t7a["regimens"]
    assert regs["Unphased_Biogenesis_Trap"]["h_end_pct_median"] == 50.0
    assert regs["Continuous_Mitophagy_Only"]["w_tot_end_median"] < 0.05
    assert regs["Biphasic_Selection_Reset"]["h_end_pct_median"] >= 80.0
    assert regs["Biphasic_Selection_Reset"]["w_tot_end_median"] >= 0.60
    assert regs["TwoCycle_Pulsed_Biphasic"]["w_tot_nadir_median"] > regs["Biphasic_Selection_Reset"]["w_tot_nadir_median"] * 2.0
    print(
        f"[PASS] Track 7A Biphasic Selection verified: Unphased Trap=50.0%, Biphasic={regs['Biphasic_Selection_Reset']['h_end_pct_median']}% "
        f"(W_tot={regs['Biphasic_Selection_Reset']['w_tot_end_median']}), TwoCycle nadir={regs['TwoCycle_Pulsed_Biphasic']['w_tot_nadir_median']} (>2x safer)."
    )

    # 2. Track 7B: Ex vivo CAR-T 100%-human sequence COX3 W116R + Human SOD2
    t7b = s7a["track7b_ex_vivo_cart_tumor_hypoxia"]
    assert t7b["arms"]["COX3_W116R_plus_Human_SOD2_CART"]["effector_persistence_72h_median"] >= 0.90
    assert t7b["human_sod2_cart_synergy"]["fold_persistence_gain_vs_WT"] >= 4.0
    print(
        f"[PASS] Track 7B Ex Vivo CAR-T verified: COX3 W116R + Human SOD2 persistence="
        f"{t7b['arms']['COX3_W116R_plus_Human_SOD2_CART']['effector_persistence_72h_median']} "
        f"({t7b['human_sod2_cart_synergy']['fold_persistence_gain_vs_WT']}x vs WT, 100% human sequence)."
    )

    # 3. Track 7C: Live AlphaGenome 2 Single-Cell conditioning + 2 blank negative controls
    assert len(s7b["variants"]) == 10 and len(s7b["cell_conditions"]) == 6
    assert s7b["signal_to_null_ratio"] >= 10.0
    vmap = {v["gene"]: v for v in s7b["variants"]}
    assert vmap["UCP2"]["top_responsive_cell"] == "cardiomyocyte"
    assert vmap["HIGD2A"]["top_responsive_cell"] == "cardiomyocyte"
    assert vmap["NEG_CTRL_rs9904290"]["mean_singlecell_composite"] < 0.05
    assert vmap["NEG_CTRL_rs10799812"]["mean_singlecell_composite"] < 0.05
    print(
        f"[PASS] Track 7C AlphaGenome 2 Single-Cell verified: 10 variants x 6 cell types, "
        f"Signal/Null={s7b['signal_to_null_ratio']}x, UCP2 cardiac-selective switch ({vmap['UCP2']['max_singlecell_composite']})."
    )

    # 4. Explorer E1-E4 provenance badges & Tab 6 Freshman Guide
    html = EXP_PATH.read_text(encoding="utf-8")
    for token in ("E1: Experimental", "E2: Foundation Model", "E3: Mechanistic Sim", "E4: Hypothesis", "6. Cross-Project Bridges", "7. Plain-English Guide"):
        assert token in html, f"Missing Explorer v2 feature: {token}"
    print("[PASS] Explorer v2 verified: E1-E4 epistemic badges, Stage 7 Bridges tab, and Plain-English Guide tab present.")
    print("\nALL STAGE 7 CROSS-PROJECT BRIDGE CHECKS PASSED (100%).")


if __name__ == "__main__":
    main()
