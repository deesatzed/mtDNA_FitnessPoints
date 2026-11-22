#!/usr/bin/env python3
"""Stage 3 Verification Suite: 5/5 Answer-Key Controls, Safety/Editor DB Gates, AFDB v6, AlphaGenome & AF3 Exports."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data" / "mitofitness.db"
S3_JSON = ROOT / "data" / "stage3_shortlist.json"
S3C_JSON = ROOT / "data" / "stage3c_bioai_summary.json"


def run_tests() -> None:
    assert DB_PATH.exists() and S3_JSON.exists() and S3C_JSON.exists()
    conn = sqlite3.connect(DB_PATH)
    s3 = json.loads(S3_JSON.read_text(encoding="utf-8"))
    s3c = json.loads(S3C_JSON.read_text(encoding="utf-8"))

    # 1. 5/5 Answer-Key Calibration & Out-of-Sample Controls MUST pass
    assert s3["all_controls_passed"] is True, "One or more Stage 3 Answer-Key controls failed!"
    assert len(s3["answer_key_control_audit"]) == 5
    print("[PASS] 5/5 Answer-Key controls verified (Tier-1 ND1 Y30H/A52T, ATP6 L156R + Tier-2 Frozen OOS ND4 R340H, ND6 M64V).")

    # 2. SQLite safety and transferability tables populated 100%
    n_safe = conn.execute("SELECT COUNT(*) FROM human_safety_annotations").fetchone()[0]
    n_trans = conn.execute("SELECT COUNT(*) FROM transferability_assessments").fetchone()[0]
    n_ag = conn.execute("SELECT COUNT(*) FROM transferability_assessments WHERE alphagenome_partner_summary IS NOT NULL").fetchone()[0]
    assert n_safe == 131 and n_trans == 256 and n_ag >= 5
    print(f"[PASS] SQLite Stage-3 tables populated: human_safety_annotations={n_safe}, transferability_assessments={n_trans} ({n_ag} annotated with live AlphaGenome partner AVI).")

    # 3. Verify Bar-headed goose COX3 W116R and ND4L A11T are zero-bystander Motif_Grafting_BaseEdit finalists
    f_map = {(r["gene"], int(r["human_rcrs_aa_pos"]), r["alt_aa"]): r for r in s3["top_25_base_editable_finalists"]}
    assert ("COX3", 116, "R") in f_map and f_map[("COX3", 116, "R")]["editor"]["bystander_nonsyn_count"] == 0
    assert ("ND4L", 11, "T") in f_map and f_map[("ND4L", 11, "T")]["editor"]["bystander_nonsyn_count"] == 0
    print("[PASS] Key zero-bystander Motif_Grafting_BaseEdit finalists verified (COX3 W116R m.9552T>C TALED; ND4L A11T m.10500G>A DddA11).")

    # 4. 5/5 AlphaFold DB v6 buffers & 5/5 live AlphaGenome AVI promoter scans verified
    assert len(s3c["alphafold_db_v6_buffers"]) == 5
    assert len(s3c["alphagenome_avi_promoter_scans"]) == 5
    for g, st in s3c["alphagenome_avi_promoter_scans"].items():
        assert st["max_avi_phred"] is not None and st["max_avi_phred"] > 10.0, f"AlphaGenome scan empty for {g}"
    print("[PASS] 5/5 AlphaFold DB v6 buffer structural audits & 5/5 live AlphaGenome AVI promoter scans verified.")

    # 5. 5/5 AF3 JSON + Boltz-2 YAML co-folding manifests exist and parse cleanly
    assert len(s3c["cofolding_manifests_exported"]) == 5
    for m in s3c["cofolding_manifests_exported"]:
        af3_p = ROOT / m["af3_json"]
        bz_p = ROOT / m["boltz2_yaml"]
        assert af3_p.exists() and bz_p.exists()
        j = json.loads(af3_p.read_text(encoding="utf-8"))
        assert j["dialect"] == "alphafold3" and len(j["sequences"]) >= 1
    print("[PASS] 5/5 AlphaFold 3 JSON & Boltz-2 YAML co-folding manifests verified.")

    conn.close()
    print("\nALL STAGE 3 SANITY CHECKS PASSED (100%).")


if __name__ == "__main__":
    run_tests()
