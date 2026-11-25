#!/usr/bin/env python3
"""Verification Suite 11: Stage 10A (Live Boltz-2 Co-Folding of ND4L Y85H + NDUFS2 and
CYTB Y224H + CYC1) + Stage 10B (Tier-3 Ecological & Pollinator Conservation Transfer across
5 Live NCBI RefSeq Keystone Species) + Explorer v5 Integration.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S10A = ROOT / "data" / "stage10a_boltz2_thermal_baric.json"
S10B = ROOT / "data" / "stage10b_ecological_transfer.json"
FIG10 = ROOT / "figures" / "stage10_boltz2_thermal_and_ecological.png"
REP10 = ROOT / "data" / "stage10_report.md"
HTML = ROOT / "app" / "mitofitness_explorer.html"


def test_stage10a_boltz2_thermal_baric() -> None:
    assert S10A.exists(), f"Missing {S10A}"
    d = json.loads(S10A.read_text(encoding="utf-8"))
    assert d["n_dimers_predicted"] == 4
    assert len(d["pairs"]) == 2
    by_pair = {p["pair"]: p for p in d["pairs"]}
    # CYTB Y224H + CYC1 docks standalone at ultra-high fidelity
    cytb = by_pair["CYTB Y224H + CYC1"]
    assert cytb["mutant"]["iptm"] >= 0.85, f"Low CYTB+CYC1 ipTM: {cytb['mutant']['iptm']}"
    assert cytb["mutant"]["chainA_rmsd_vs_cryoem_A"] <= 1.2
    assert cytb["wt_vs_mut_10A_pocket_rmsd_A"] <= 0.30
    assert cytb["verdict"] == "INTERFACE_PRESERVED_HIGH_FIDELITY"
    # ND4L Y85H + NDUFS2 is a membrane-to-peripheral contact bridged by ND1/ND2 in CI cryo-EM (WT ipTM=0.154):
    # monomer folds cleanly (<=2.2 A vs 9I4I) with tiny |delta_iptm| < 0.05!
    nd4l = by_pair["ND4L Y85H + NDUFS2"]
    assert nd4l["mutant"]["chainA_rmsd_vs_cryoem_A"] <= 2.2
    assert abs(nd4l["delta_iptm"]) < 0.05
    assert nd4l["verdict"] == "SCAFFOLD_DEPENDENT_DIMER_MONOMER_CONSERVED"
    print(
        "[PASS] Stage 10A Boltz-2 (4 Thermal/Baric Dimers): "
        + ", ".join(
            f"{p['pair']} (ipTM={p['mutant']['iptm']:.3f}, d={p['delta_iptm']:+.3f}, 10A={p['wt_vs_mut_10A_pocket_rmsd_A']:.2f}A)"
            for p in d["pairs"]
        )
    )


def test_stage10b_tier3_ecological() -> None:
    assert S10B.exists(), f"Missing {S10B}"
    d = json.loads(S10B.read_text(encoding="utf-8"))
    assert d["solver_convergence"]["passed"] is True
    assert d["solver_convergence"]["max_abs_delta"] < 1e-4

    # All 5 live RefSeq species must extract all 5 OXPHOS genes cleanly
    audit = d["species_alignment_audit"]
    assert len(audit) == 5
    for sp, info in audit.items():
        assert len(info["gene_identities_vs_human"]) == 5, f"{sp} missing genes: {info['gene_identities_vs_human']}"
        assert info["gene_identities_vs_human"]["COX1"]["pct_identity_vs_human"] > 50.0, f"{sp} COX1 id too low"
        assert info["qo_strobilurin_pocket_conservation"] == "9/9", f"{sp} Qo pocket mismatch"

    # Ecological ODE synergy & thermal sign checks
    scs = d["ecological_scenarios_lhs5k"]
    # Honeybee heat+fungicide: TightCoupled alone hurts, Combined Haplotype + Apis-BCR surges > +0.50!
    bee = scs["E1_Honeybee_42C_Fungicide_Flight"]["arms"]
    assert bee["Natural_TightCoupled"]["delta_vs_wt_med"] < 0.0, "Expected tight coupling alone to worsen Qo fungicide ROS!"
    assert bee["Combined_Haplotype_BCR"]["delta_vs_wt_med"] > 0.50, "Expected massive Honeybee positive epistasis!"
    # Coral bleaching: UCP2 hurts while TightCoupled & Combined rescue strongly
    cor = scs["E2_Coral_31_5C_Marine_Heatwave"]["arms"]
    assert cor["Uncoupler_UCP2_Only"]["delta_vs_wt_med"] < 0.0
    assert cor["Combined_Haplotype_BCR"]["delta_vs_wt_med"] > 0.30
    # Atlantic Salmon warming river hypoxia: Combined beats WT by > +0.25
    sal = scs["E3_Salmon_22C_WarmRiver_Hypoxia"]["arms"]
    assert sal["Combined_Haplotype_BCR"]["delta_vs_wt_med"] > 0.25
    print(
        f"[PASS] Stage 10B Tier-3 Ecological (5/5 RefSeq 9/9 Qo pocket): "
        f"Honeybee Combined={bee['Combined_Haplotype_BCR']['delta_vs_wt_med']:+.3f} | "
        f"Coral={cor['Combined_Haplotype_BCR']['delta_vs_wt_med']:+.3f} | "
        f"Salmon={sal['Combined_Haplotype_BCR']['delta_vs_wt_med']:+.3f}."
    )


def test_stage10_artifacts_and_explorer_v5() -> None:
    assert FIG10.exists() and FIG10.stat().st_size > 40_000
    assert REP10.exists() and REP10.stat().st_size > 800
    txt = HTML.read_text(encoding="utf-8")
    for token in ["stage10a_boltz2", "stage10b_ecological", "Acropora_millepora", "Salmo_salar", "ND4L_Y85H_NDUFS2"]:
        assert token in txt, f"Explorer v5 missing Stage 10 token: {token}"
    print(f"[PASS] Stage 10 Figure ({FIG10.stat().st_size:,} B), Report & Explorer v5 ({len(txt):,} B) verified.")


if __name__ == "__main__":
    test_stage10a_boltz2_thermal_baric()
    test_stage10b_tier3_ecological()
    test_stage10_artifacts_and_explorer_v5()
    print("\nALL STAGE 10 VERIFICATION TESTS PASSED!")
