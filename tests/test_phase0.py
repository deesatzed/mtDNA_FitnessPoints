#!/usr/bin/env python3
"""Phase 0 Verification Suite: Database schema, NCBI Table 2 translation, and rCRS coordinate parity."""

from __future__ import annotations

import sqlite3
from pathlib import Path
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data" / "mitofitness.db"


def mutate_rcrs_and_translate(conn: sqlite3.Connection, gene: str, rcrs_pos: int, ref_nt: str, alt_nt: str) -> tuple[int, str, str]:
    """Apply a 1-based rCRS H-strand nucleotide substitution and return (aa_pos, ref_aa, alt_aa)."""
    row = conn.execute(
        "SELECT rcrs_start_nt, rcrs_end_nt, strand, cds_nt_seq, protein_aa_seq FROM human_rcrs_genes WHERE gene_symbol = ?",
        (gene,),
    ).fetchone()
    assert row is not None, f"Gene {gene} missing"
    start_nt, end_nt, strand, cds_nt, prot_aa = row
    assert start_nt <= rcrs_pos <= end_nt, f"Position {rcrs_pos} outside {gene} ({start_nt}..{end_nt})"

    if strand == "+":
        offset = rcrs_pos - start_nt
        actual_ref = cds_nt[offset]
        assert actual_ref == ref_nt, f"{gene} m.{rcrs_pos}: expected {ref_nt}, found {actual_ref}"
        mut_cds = cds_nt[:offset] + alt_nt + cds_nt[offset + 1 :]
    else:
        # On '-' strand (ND6), H-strand ref_nt/alt_nt must be complemented and indexed from end_nt
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        offset = end_nt - rcrs_pos
        actual_ref = cds_nt[offset]
        assert actual_ref == comp[ref_nt], f"{gene} (-) m.{rcrs_pos}: expected {comp[ref_nt]}, found {actual_ref}"
        mut_cds = cds_nt[:offset] + comp[alt_nt] + cds_nt[offset + 1 :]

    aa_idx = offset // 3
    ref_codon = cds_nt[aa_idx * 3 : aa_idx * 3 + 3]
    alt_codon = mut_cds[aa_idx * 3 : aa_idx * 3 + 3]
    ref_aa = prot_aa[aa_idx]
    alt_aa = "M" if (aa_idx == 0 and alt_codon in {"ATG", "ATA", "ATT", "ATC", "GTG"}) else str(Seq(alt_codon).translate(table=2))
    assert str(Seq(ref_codon).translate(table=2)) == ref_aa or (aa_idx == 0 and ref_aa == "M")
    return aa_idx + 1, ref_aa, alt_aa


def run_tests() -> None:
    assert DB_PATH.exists(), f"Database not found: {DB_PATH}"
    conn = sqlite3.connect(DB_PATH)

    # 1. Verify all 8 tables exist
    tables = {
        r[0]
        for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'").fetchall()
    }
    expected_tables = {
        "human_rcrs_genes",
        "donor_species",
        "environmental_stressors",
        "mt_alleles",
        "human_safety_annotations",
        "fitness_payoff_profiles",
        "transferability_assessments",
        "epistatic_synergies",
    }
    assert expected_tables.issubset(tables), f"Missing tables: {expected_tables - tables}"
    print(f"[PASS] Schema check: all {len(expected_tables)} relational tables present.")

    # 2. Verify all 13 genes, tier counts (6 Tier-1, 7 Tier-2), and exact strand assignments
    genes = conn.execute(
        "SELECT gene_symbol, tier, strand, expected_aa_len, length(protein_aa_seq), length(cds_nt_seq) FROM human_rcrs_genes"
    ).fetchall()
    assert len(genes) == 13, f"Expected 13 genes, got {len(genes)}"
    t1 = [g[0] for g in genes if g[1] == 1]
    t2 = [g[0] for g in genes if g[1] == 2]
    assert len(t1) == 6 and len(t2) == 7, f"Tier split mismatch: T1={t1}, T2={t2}"
    for sym, tier, strand, exp_len, obs_len, cds_len in genes:
        assert obs_len == exp_len, f"{sym}: AA length {obs_len} != {exp_len}"
        assert cds_len == exp_len * 3, f"{sym}: CDS length {cds_len} != {exp_len * 3}"
        if sym == "ND6":
            assert strand == "-", "ND6 must be on '-' (L-strand)"
        else:
            assert strand == "+", f"{sym} must be on '+' (H-strand)"
    print("[PASS] 13/13 Human rCRS protein-coding genes verified (6 Tier-1 + 7 Tier-2; exact UniProt lengths & strands).")

    # 3. Landmark rCRS coordinate-to-AA spot checks (including '-' strand ND6!)
    landmarks = [
        ("ND1", 3394, "T", "C", 30, "Y", "H", "Tibetan high-altitude / adaptive ND1 Y30H"),
        ("ND1", 3460, "G", "A", 52, "A", "T", "Pathogenic LHON control ND1 A52T"),
        ("ATP6", 8993, "T", "G", 156, "L", "R", "Pathogenic Leigh/NARP control ATP6 L156R"),
        ("ND4", 11778, "G", "A", 340, "R", "H", "Out-of-sample Tier-2 LHON control ND4 R340H"),
        ("ND6", 14484, "T", "C", 64, "M", "V", "Minus-strand LHON control ND6 M64V"),
    ]
    for gene, pos, ref_nt, alt_nt, exp_aa_pos, exp_ref, exp_alt, desc in landmarks:
        aa_pos, ref_aa, alt_aa = mutate_rcrs_and_translate(conn, gene, pos, ref_nt, alt_nt)
        assert (aa_pos, ref_aa, alt_aa) == (exp_aa_pos, exp_ref, exp_alt), (
            f"Landmark failed for {desc}: got {ref_aa}{aa_pos}{alt_aa}, expected {exp_ref}{exp_aa_pos}{exp_alt}"
        )
        print(f"  [PASS] m.{pos}{ref_nt}>{alt_nt} -> {gene} {ref_aa}{aa_pos}{alt_aa} ({desc})")

    # 4. Stressors & Homo sapiens Reference seeded
    n_stressors = conn.execute("SELECT COUNT(*) FROM environmental_stressors").fetchone()[0]
    n_ref = conn.execute("SELECT COUNT(*) FROM donor_species WHERE role = 'Reference'").fetchone()[0]
    assert n_stressors == 5 and n_ref == 1
    print(f"[PASS] Seeded {n_stressors} World Space stressors and Homo sapiens NC_012920.1 Reference row.")
    conn.close()
    print("\nALL PHASE 0 SANITY CHECKS PASSED (100%).")


if __name__ == "__main__":
    run_tests()
