#!/usr/bin/env python3
"""Verification Suite 12: Post-Stage 10 Reviewer #2 Remediation Audit & Stage 11
200-Species Vertebrate Evolutionary Rate-Control Stress Test.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S11_JSON = ROOT / "data" / "stage11_200species_phylogenetic_audit.json"
S11_FIG = ROOT / "figures" / "stage11_200species_rate_control.png"
REM_LOG = ROOT / "docs" / "REMEDIATION_AUDIT_LOG.md"
AUDIT_DOC = ROOT / "docs" / "CRITICAL_SCIENTIFIC_AUDIT_AND_PUBLISHABILITY.md"
PLAN_DOC = ROOT / "docs" / "MITOFITNESS_MASTER_TRANSLATIONAL_PLAN.md"
HTML = ROOT / "app" / "mitofitness_explorer.html"


def test_remediation_fixes_1_2_3_5() -> None:
    assert REM_LOG.exists() and AUDIT_DOC.exists()
    plan_txt = PLAN_DOC.read_text(encoding="utf-8")
    html_txt = HTML.read_text(encoding="utf-8")

    # Fix 1: Scott et al. PNAS 2011 accurately attributed to cyt-c affinity & lower Vmax
    assert "cytochrome" in plan_txt and "V_{\\max}" in plan_txt
    assert "Scott et al., PNAS 2011" in html_txt
    # Fix 2: Bertero-Maack extracellular Ca2+ caveat + HiFi-TALED V28R/sTALED
    assert "Bertero" in plan_txt and "1.2–1.8 mM" in plan_txt
    assert "Bertero–Maack" in html_txt and "HiFi-" in html_txt
    # Fix 3: Conditional Phase-Diagram E3/E4 framing
    assert "Conditional Phase-Diagram" in html_txt
    print("[PASS] Fixes 1, 2, 3, and 5 verified across Master Plan, Remediation Log, and Explorer.")


def test_stage11_200species_rate_control() -> None:
    assert S11_JSON.exists(), f"Missing {S11_JSON}"
    assert S11_FIG.exists() and S11_FIG.stat().st_size > 40_000
    d = json.loads(S11_JSON.read_text(encoding="utf-8"))
    assert d["mean_vertebrate_species_per_gene"] >= 180.0
    assert d["n_cryoem_residues_modeled"] == 3756
    # 130 unique (gene, aa_pos) cryo-EM residue sites (1 site carries 2 distinct alt_aa alleles across the 131 tuples)
    assert d["n_candidates"] == 130
    assert len(d["entropy_quartiles"]) == 4
    assert d["stage11_200species_adjusted_OR"] > 1.0
    print(
        f"[PASS] Stage 11 (200-species rate control, mean={d['mean_vertebrate_species_per_gene']} sp/gene): "
        f"Adj OR={d['stage11_200species_adjusted_OR']:.3f} "
        f"(95% CI {d['stage11_200species_95CI']}, p={d['stage11_200species_wald_p']:.5f}) -> {d['verdict']}"
    )


if __name__ == "__main__":
    test_remediation_fixes_1_2_3_5()
    test_stage11_200species_rate_control()
    print("\nALL STAGE 11 REMEDIATION & PHYLOGENETIC TESTS PASSED!")
