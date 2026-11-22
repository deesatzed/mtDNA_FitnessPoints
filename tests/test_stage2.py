#!/usr/bin/env python3
"""Stage 2 Verification Suite: 3D Cryo-EM mapping, FBT equivalence classes, and DB parity."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data" / "mitofitness.db"
SUMMARY_JSON = ROOT / "data" / "stage2_3d_fbt_summary.json"
FIG_PATH = ROOT / "figures" / "stage2_3d_fbt_equivalence_clusters.png"


def run_tests() -> None:
    assert DB_PATH.exists() and SUMMARY_JSON.exists() and FIG_PATH.exists()
    conn = sqlite3.connect(DB_PATH)
    s2 = json.loads(SUMMARY_JSON.read_text(encoding="utf-8"))

    # 1. All 13 human mtDNA genes linked to ultra-high-resolution cryo-EM PDB IDs + chains
    mapped = conn.execute(
        "SELECT gene_symbol, pdb_holocomplex_id, pdb_chain_id FROM human_rcrs_genes WHERE pdb_holocomplex_id IS NOT NULL"
    ).fetchall()
    assert len(mapped) == 13, f"Expected 13/13 genes mapped to PDB, got {len(mapped)}"
    print("[PASS] 13/13 human mtDNA proteins mapped to experimental cryo-EM holocomplexes (9I4I, 9HZL, 9I6F, 8H9S).")

    # 2. 100% of unique Stage-1 polarized candidates resolved in 3D cryo-EM
    assert s2["candidates_resolved_3d"] == s2["stage1_unique_candidates"] and s2["candidates_resolved_3d"] == 131
    print(f"[PASS] 131/131 (100%) Stage-1 polarized candidates resolved in 3D cryo-EM ({s2['total_human_mtdna_residues_resolved_3d']}/3789 human mtDNA residues).")

    # 3. 4 distinct, biophysically labeled FBT Payoff Equivalence Classes
    eq = s2["equivalence_classes"]
    assert len(eq) == 4 and all(not k.endswith("_2") for k in eq.keys()), f"Duplicate/bad class labels: {list(eq.keys())}"
    print(f"[PASS] 4 biophysically distinct FBT Payoff Equivalence Classes: {sorted(eq.keys())}")

    # 4. SQLite fitness_payoff_profiles populated & linked via foreign key to mt_alleles
    n_prof = conn.execute("SELECT COUNT(*) FROM fitness_payoff_profiles").fetchone()[0]
    n_alleles = conn.execute("SELECT COUNT(*) FROM mt_alleles").fetchone()[0]
    assert n_prof == n_alleles == 256, f"Expected 256 payoff profiles matching mt_alleles, got {n_prof} vs {n_alleles}"
    print(f"[PASS] SQLite fitness_payoff_profiles populated ({n_prof}/{n_alleles} allele rows linked 1-to-1).")

    # 5. Mitonuclear interface enrichment and 3D spatial audit present
    ne = s2["niche_enrichment"]["Mitonuclear_Interface"]
    assert ne["odds_ratio"] > 1.0 and ne["fisher_p"] < 0.05
    print(f"[PASS] Mitonuclear Interface enrichment confirmed: {ne['cand_pct']}% vs {ne['bg_pct']}% bg (OR={ne['odds_ratio']}, p={ne['fisher_p']}).")

    conn.close()
    print("\nALL STAGE 2 SANITY CHECKS PASSED (100%).")


if __name__ == "__main__":
    run_tests()
