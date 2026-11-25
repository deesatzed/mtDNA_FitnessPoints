#!/usr/bin/env python3
"""Verification Suite 9: Stage 8A (Full Mito_Papers 3-State Replication, Trap A/B,
Donor-Mito Reversibility & Surgery Pre-Op Clocks) + Stage 8B (Thermal & Baric Regimes)
+ Explorer v3 Single-App Integration.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S8A = ROOT / "data" / "stage8a_replication_engraftment.json"
S8B = ROOT / "data" / "stage8b_thermal_baric.json"
FIG8 = ROOT / "figures" / "stage8_replication_thermal_and_clocks.png"
REP8 = ROOT / "data" / "stage8_report.md"
HTML = ROOT / "app" / "mitofitness_explorer.html"


def test_stage8a_replication_and_clocks() -> None:
    assert S8A.exists(), f"Missing {S8A}"
    d = json.loads(S8A.read_text(encoding="utf-8"))

    # 1. Solver convergence RK4 vs DOP853/Radau
    conv = d["solver_convergence"]
    assert conv["passed"] is True, f"Stage 8A solver failed: {conv}"
    assert conv["max_abs_delta"] < 1e-4

    # 2. Trap B (HighGlucose 25mM) & Trap A (Unphased) MUST fail while Unmasked Biphasic rescues!
    traps = d["trap_a_b_wetlab_14d_lhs10k"]
    fM_trapB = traps["HighGlucose_25mM_TrapB_Failure"]["f_M_end"]["median"]
    fM_trapA = traps["Unphased_Biogenesis_TrapA"]["f_M_end"]["median"]
    fM_unmask = traps["Biphasic_Unmasked_Galactose_PDK1"]["f_M_end"]["median"]
    fD_unmask = traps["Biphasic_Unmasked_Galactose_PDK1"]["f_D_end"]["median"]
    fM_pulsed = traps["TwoCycle_Pulsed_Unmasked_Safe"]["f_M_end"]["median"]
    nadir_7d = traps["Biphasic_Unmasked_Galactose_PDK1"]["tot_nadir"]["median"]
    nadir_2c = traps["TwoCycle_Pulsed_Unmasked_Safe"]["tot_nadir"]["median"]

    assert fM_trapB < 0.45, f"Trap B failed to stall under high glucose: {fM_trapB}"
    assert fM_trapA < 0.45, f"Trap A failed to stall under unphased mTORC1: {fM_trapA}"
    assert fM_unmask > 0.70, f"Unmasked biphasic failed to amplify: {fM_unmask}"
    assert fD_unmask < 0.05, f"Unmasked biphasic failed to purge aged D: {fD_unmask}"
    assert fM_pulsed > 0.65 and nadir_2c > 2.0 * nadir_7d, "2-cycle pulsed must protect mass nadir"

    # 3. Sobol GSA completeness
    sob = d["sobol_gsa_10240"]
    assert sob["total_evaluations"] == 10240
    assert len(sob["ranked_parameters"]) == 8

    # 4. Context-sensitive Donor-Mitochondria 60-day reversibility
    don = d["donor_mitochondria_60d_reversibility"]
    assert don["Healthy_Normoxic_Organ"]["classification"] == "AUTO_REVERSIBLE_WASHOUT"
    assert don["PostTransplant_24hHypox_ThenNormox"]["classification"] == "AUTO_REVERSIBLE_WASHOUT"
    assert don["Chronically_Hypoxic_Aged_Tissue_BCR"]["classification"] == "SELF_TARGETING_ENGRAFTMENT"
    assert don["Healthy_Normoxic_Organ"]["f_M_day60"]["median"] < 0.05
    assert don["Chronically_Hypoxic_Aged_Tissue_BCR"]["f_M_day60"]["median"] > 0.45

    # 5. Scheduled surgery clocks: mRNA nullified at 0h surprise (<0.20) but >0.95 pre-op & >0.80 delayed +2h
    cs = d["scheduled_surgery_preop_clocks"]["clock_summary"]
    assert cs["Clock1_SmallMol_Malonate"]["coverage_if_dosed_at_0h_surprise"] > 0.75
    assert cs["Clock2_BankedMito_TAT_SOD2"]["coverage_if_dosed_at_0h_surprise"] > 0.70
    assert cs["Clock3_mRNA_LNP_SOD2_UCP2"]["coverage_if_dosed_at_0h_surprise"] < 0.20
    assert cs["Clock3_mRNA_LNP_SOD2_UCP2"]["peak_minute0_coverage"] > 0.95
    assert cs["Clock3_mRNA_LNP_SOD2_UCP2"]["delayed_wave_coverage_if_dosed_at_plus2h_ER"] > 0.80
    print(f"[PASS] Stage 8A: TrapB={fM_trapB*100:.1f}% vs Unmasked={fM_unmask*100:.1f}% (D={fD_unmask*100:.1f}%), Donor Healthy={don['Healthy_Normoxic_Organ']['f_M_day60']['median']*100:.1f}% vs Injured={don['Chronically_Hypoxic_Aged_Tissue_BCR']['f_M_day60']['median']*100:.1f}%")


def test_stage8b_thermal_and_baric() -> None:
    assert S8B.exists(), f"Missing {S8B}"
    d = json.loads(S8B.read_text(encoding="utf-8"))

    # 1. 3-regime solver convergence & adversarial Null_Placebo == 0.0
    assert d["solver_convergence"]["passed"] is True
    assert d["solver_convergence"]["max_abs_delta"] < 1e-4
    for rname, err in d["adversarial_null_placebo_max_err"].items():
        assert err == 0.0, f"Null_Placebo non-zero in {rname}: {err}"

    regs = d["regimes"]
    # 2. Cold Storage 4C->37C: Arctic_Cold COX1_I57V rescues & COX3_W116R shows positive buffer rescue
    cs = regs["Cold_Storage_4C_to_37C"]["matrix"]
    assert cs["COX1_I57V"]["combos"]["Human_SOD2"]["delta_vs_wt_med"] > 0.10
    assert cs["COX3_W116R"]["combos"]["Human_SOD2"]["delta_vs_wt_med"] > cs["COX3_W116R"]["combos"]["None"]["delta_vs_wt_med"] + 0.05
    assert cs["COX3_W116R"]["combos"]["Yeast_NDI1"]["delta_vs_wt_med"] > 0.15

    # 3. Hyperthermia 37C->41.5C: Desert_Heat ND4L_Y85H + SOD2 rescues strongly, while UCP2 HURTS vs Driver Alone!
    ht = regs["Hyperthermia_37C_to_41_5C"]["matrix"]
    y85_alone = ht["ND4L_Y85H"]["combos"]["None"]["delta_vs_wt_med"]
    y85_sod2 = ht["ND4L_Y85H"]["combos"]["Human_SOD2"]["delta_vs_wt_med"]
    y85_ucp2 = ht["ND4L_Y85H"]["combos"]["Human_UCP2"]["delta_vs_wt_med"]
    assert y85_sod2 > y85_alone > 0.20, f"Desert_Heat tight-coupling failed: alone={y85_alone}, sod2={y85_sod2}"
    assert y85_ucp2 < y85_alone - 0.15, f"Expected UCP2 hyperthermic negative sign reversal: ucp2={y85_ucp2} vs alone={y85_alone}"

    # 4. Hydrostatic 150 atm Dive: Sperm Whale ATP6_I106V boosts dive ATP nadir & beats WT
    dv = regs["Hydrostatic_Dive_150atm"]
    assert dv["matrix"]["ATP6_I106V"]["combos"]["None"]["delta_vs_wt_med"] > 0.05
    assert dv["matrix"]["ATP6_I106V"]["driver_only_atp_nadir_med"] > dv["wt_atp_nadir_med"] + 0.15

    # 5. Cold-storage Sobol GSA
    assert d["cold_storage_sobol_10240"]["total_evaluations"] == 10240
    print(f"[PASS] Stage 8B: Cold4C COX1_I57V+SOD2={cs['COX1_I57V']['combos']['Human_SOD2']['delta_vs_wt_med']:+.3f} | Heat41.5C ND4L_Y85H+SOD2={y85_sod2:+.3f} vs +UCP2={y85_ucp2:+.3f} | Dive150atm ATP6_I106V={dv['matrix']['ATP6_I106V']['combos']['None']['delta_vs_wt_med']:+.3f}")


def test_stage8_artifacts_and_explorer_v3() -> None:
    assert FIG8.exists() and FIG8.stat().st_size > 40_000, f"Missing/tiny {FIG8}"
    assert REP8.exists() and REP8.stat().st_size > 1_000, f"Missing/tiny {REP8}"
    assert HTML.exists(), f"Missing {HTML}"
    txt = HTML.read_text(encoding="utf-8")
    for token in [
        "fAxis",
        "regimeSelector",
        "Cold_Storage_4C_to_37C",
        "Hyperthermia_37C_to_41_5C",
        "Hydrostatic_Dive_150atm",
        "HighGlucose_25mM_TrapB_Failure",
        "Chronically_Hypoxic_Aged_Tissue_BCR",
        "Clock3_mRNA_LNP_SOD2_UCP2",
        "stage8a_replication",
        "stage8b_thermal_baric",
    ]:
        assert token in txt, f"Explorer v3 missing Stage 8 token: {token}"
    print(f"[PASS] Stage 8 Figure ({FIG8.stat().st_size:,} B), Report & Explorer v3 ({len(txt):,} B) verified.")


if __name__ == "__main__":
    test_stage8a_replication_and_clocks()
    test_stage8b_thermal_and_baric()
    test_stage8_artifacts_and_explorer_v3()
    print("\nALL STAGE 8 VERIFICATION TESTS PASSED!")
