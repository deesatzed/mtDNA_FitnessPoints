#!/usr/bin/env python3
"""Verification Suite 10: Stage 9A (Live OpenMM AMBER14 Relaxation of 5 Thermal/Baric Leads)
+ Stage 9B (Live AlphaGenome 2 Clinical Organ-Compartment Conditioning across 6 Organs)
+ Explorer v4 Integration (Hardened Trap-B Wet-Lab SOP & Surgery Dosing Table).
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S9A = ROOT / "data" / "stage9a_openmm_thermal_baric.json"
S9B = ROOT / "data" / "stage9b_alphagenome2_organs.json"
FIG9 = ROOT / "figures" / "stage9_openmm_thermal_and_organ_ag2.png"
REP9 = ROOT / "data" / "stage9_report.md"
HTML = ROOT / "app" / "mitofitness_explorer.html"


def test_stage9a_openmm_thermal_baric() -> None:
    assert S9A.exists(), f"Missing {S9A}"
    d = json.loads(S9A.read_text(encoding="utf-8"))
    assert d["n_evaluated"] == 5
    # All 5 thermal/baric leads must relax cleanly with zero clashes (d_min >= 2.50 A)
    # and pocket perturbation strictly smaller than pathogenic ATP6 L156R (+0.160 A)!
    patho_ref = d["benchmark_anchors"]["Stage6_Pathogenic_ATP6_L156R_delta_rmsd_A"]
    assert d["n_stable_low_perturbation"] >= 4, f"Too few stable leads: {d['n_stable_low_perturbation']}/5"
    by_var = {r["variant"]: r for r in d["results"]}
    for v in ["ND4L_Y85H", "ND1_S313G", "COX1_I57V", "CYTB_F181L", "CYTB_Y224H"]:
        r = by_var[v]
        assert r["mutant"]["post_relax_sc_dmin_A"] >= 2.50, f"Steric clash in {v}: {r['mutant']['post_relax_sc_dmin_A']}"
        assert abs(r["delta_pocket_rmsd_A"]) < patho_ref, f"{v} exceeds pathogenic control RMSD!"
        # Verify minimized PDBs were saved
        assert (ROOT / "data" / "stage9_openmm" / f"{v}_WT_min.pdb").exists()
        assert (ROOT / "data" / "stage9_openmm" / f"{v}_MUT_min.pdb").exists()
    print(
        f"[PASS] Stage 9A OpenMM AMBER14: {d['n_stable_low_perturbation']}/5 STABLE_LOW_PERTURBATION "
        f"(ND4L_Y85H dRMSD={by_var['ND4L_Y85H']['delta_pocket_rmsd_A']:+.4f} A, "
        f"COX1_I57V dRMSD={by_var['COX1_I57V']['delta_pocket_rmsd_A']:+.4f} A vs Pathogenic +{patho_ref:.3f} A)."
    )


def test_stage9b_alphagenome2_organs() -> None:
    assert S9B.exists(), f"Missing {S9B}"
    d = json.loads(S9B.read_text(encoding="utf-8"))
    assert len(d["organ_conditions"]) == 6
    assert len(d["variants"]) == 9
    assert d["signal_to_null_ratio"] >= 5.0, f"Signal/Null too low: {d['signal_to_null_ratio']}"
    by_g = {v["gene"]: v for v in d["variants"]}
    # Verify key clinical compartment findings & blank controls
    assert by_g["UCP2"]["top_responsive_organ"] == "cardiomyocyte"
    assert by_g["SOD2"]["top_responsive_organ"] == "lung_aec2"
    assert by_g["NEG_CTRL_rs9904290"]["mean_composite"] < 0.05
    assert by_g["NEG_CTRL_rs10799812"]["mean_composite"] < 0.05
    print(
        f"[PASS] Stage 9B Live AlphaGenome 2 (6 Clinical Organs): Signal/Null = {d['signal_to_null_ratio']}x "
        f"| SOD2 top=lung_aec2 ({by_g['SOD2']['top_responsive_composite']}) "
        f"| UCP2 top=cardiomyocyte ({by_g['UCP2']['top_responsive_composite']})."
    )


def test_stage9_artifacts_and_explorer_v4() -> None:
    assert FIG9.exists() and FIG9.stat().st_size > 40_000, f"Missing/tiny {FIG9}"
    assert REP9.exists() and REP9.stat().st_size > 800, f"Missing/tiny {REP9}"
    txt = HTML.read_text(encoding="utf-8")
    for token in [
        "stage9a_openmm",
        "stage9b_organs",
        "s9aOpenmmBody",
        "s9bOrganTable",
        "AZD7545",
        "lung_aec2",
        "pancreatic_beta",
        "ND4L_Y85H",
    ]:
        assert token in txt, f"Explorer v4 missing Stage 9 token: {token}"
    print(f"[PASS] Stage 9 Figure ({FIG9.stat().st_size:,} B), Report & Explorer v4 ({len(txt):,} B) verified.")


if __name__ == "__main__":
    test_stage9a_openmm_thermal_baric()
    test_stage9b_alphagenome2_organs()
    test_stage9_artifacts_and_explorer_v4()
    print("\nALL STAGE 9 VERIFICATION TESTS PASSED!")
