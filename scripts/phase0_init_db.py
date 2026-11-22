#!/usr/bin/env python3
"""Phase 0: Initialize portable MitoFitness Atlas SQLite database & verify Human rCRS (NC_012920.1).

- Executes schema/mitofitness_schema.sql (PostgreSQL-compatible DDL) on data/mitofitness.db.
- Loads the live-fetched Human rCRS 16,569-bp mitogenome and 13 NCBI RefSeq proteins.
- Locates exact 1-based rCRS CDS coordinates and strand for all 13 genes by direct 6-frame
  translation under NCBI Translation Table 2 (Vertebrate Mitochondrial), verifying 100%
  residue-by-residue identity (0 mismatches) against official RefSeq proteins.
- Seeds World Space (W) environmental_stressors and Homo sapiens Reference species row.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
SCHEMA_PATH = ROOT / "schema" / "mitofitness_schema.sql"
DB_PATH = ROOT / "data" / "mitofitness.db"
RCRS_NUCCORE_JSON = ROOT / "data" / "raw" / "human_rcrs_nuccore.json"
RCRS_PROTEINS_JSON = ROOT / "data" / "raw" / "human_rcrs_proteins.json"

# Canonical gene metadata (UniProt accessions, OXPHOS complex, algorithmic progression tier)
GENE_METADATA = {
    "YP_003024026.1": ("ND1", "Complex_I", 1, "P03886", 318),
    "YP_003024027.1": ("ND2", "Complex_I", 1, "P03891", 347),
    "YP_003024028.1": ("COX1", "Complex_IV", 1, "P00395", 513),
    "YP_003024029.1": ("COX2", "Complex_IV", 2, "P00403", 227),
    "YP_003024030.1": ("ATP8", "Complex_V", 2, "P03928", 68),
    "YP_003024031.1": ("ATP6", "Complex_V", 1, "P00846", 226),
    "YP_003024032.1": ("COX3", "Complex_IV", 1, "P00414", 261),
    "YP_003024033.1": ("ND3", "Complex_I", 2, "P03897", 115),
    "YP_003024034.1": ("ND4L", "Complex_I", 2, "P03901", 98),
    "YP_003024035.1": ("ND4", "Complex_I", 2, "P03905", 459),
    "YP_003024036.1": ("ND5", "Complex_I", 2, "P03915", 603),
    "YP_003024037.1": ("ND6", "Complex_I", 2, "P03923", 174),
    "YP_003024038.1": ("CYTB", "Complex_III", 1, "P00156", 380),
}

# Vertebrate mitochondrial start codons (NCBI Table 2: ATG, ATA, ATT, ATC, GTG -> Met at pos 1)
VERT_MITO_STARTS = {"ATG", "ATA", "ATT", "ATC", "GTG"}

INITIAL_STRESSORS = [
    (
        "High_Altitude_Hypoxia",
        "pO2_kPa",
        "kPa",
        "Chronic hypobaric hypoxia paired with cold (e.g., Tibetan plateau, Andes, trans-Himalayan flight).",
    ),
    (
        "Diving_Subterranean_Hypoxia",
        "pO2_kPa",
        "kPa",
        "Acute intermittent hypoxia and ischemia-reperfusion without altitude cold (deep-diving cetaceans/pinnipeds, subterranean mole-rats).",
    ),
    (
        "Arctic_Cold",
        "Temp_C",
        "degC",
        "Extreme cold exposure selecting for mild OXPHOS uncoupling / proton-leak thermogenesis under normoxia.",
    ),
    (
        "Desert_Heat",
        "Temp_C",
        "degC",
        "Chronic thermal / arid stress selecting for OXPHOS thermal stability and low ROS leakage.",
    ),
    (
        "Human_Ischemia_Reperfusion",
        "pO2_kPa",
        "kPa",
        "Target human clinical/physiological stress regime: acute tissue hypoxia followed by re-oxygenation RET-ROS burst.",
    ),
]


def translate_mito_cds(cds_nt: str) -> str:
    """Translate a vertebrate mitochondrial CDS (NCBI Table 2), converting valid start codon to M."""
    aa = str(Seq(cds_nt).translate(table=2))
    if cds_nt[:3] in VERT_MITO_STARTS and aa:
        aa = "M" + aa[1:]
    return aa


def locate_cds_on_rcrs(rcrs_seq: str, target_aa: str) -> tuple[int, int, str, str]:
    """Locate exact 1-based (start_nt, end_nt, strand, cds_nt) on rCRS matching target_aa[1:].

    Note: The first codon may be an alternative vertebrate mito start codon (ATT/ATA/ATC/GTG)
    translated as 'M', so we match the exact suffix target_aa[1:] across all 6 frames and then
    verify the start codon and exact full translation.
    """
    n = len(rcrs_seq)
    aa_len = len(target_aa)
    nt_len = aa_len * 3
    suffix = target_aa[1:]

    # Search '+' strand (frames 0, 1, 2)
    for frame in range(3):
        sub = rcrs_seq[frame : frame + ((n - frame) // 3) * 3]
        prot = str(Seq(sub).translate(table=2))
        idx = prot.find(suffix)
        if idx >= 1:
            start_0 = frame + (idx - 1) * 3
            cds_nt = rcrs_seq[start_0 : start_0 + nt_len]
            if translate_mito_cds(cds_nt) == target_aa:
                return start_0 + 1, start_0 + nt_len, "+", cds_nt

    # Search '-' strand (reverse complement frames 0, 1, 2)
    rc_seq = str(Seq(rcrs_seq).reverse_complement())
    for frame in range(3):
        sub = rc_seq[frame : frame + ((n - frame) // 3) * 3]
        prot = str(Seq(sub).translate(table=2))
        idx = prot.find(suffix)
        if idx >= 1:
            rc_start_0 = frame + (idx - 1) * 3
            cds_nt = rc_seq[rc_start_0 : rc_start_0 + nt_len]
            if translate_mito_cds(cds_nt) == target_aa:
                # Convert reverse-complement 0-based [rc_start_0, rc_start_0 + nt_len)
                # back to 1-based forward rCRS coordinates [fwd_start, fwd_end]
                fwd_end = n - rc_start_0
                fwd_start = n - (rc_start_0 + nt_len) + 1
                return fwd_start, fwd_end, "-", cds_nt

    raise ValueError(f"Could not locate exact CDS match on rCRS for protein length {aa_len}")


def main() -> None:
    if not RCRS_NUCCORE_JSON.exists() or not RCRS_PROTEINS_JSON.exists():
        raise FileNotFoundError("Missing live NCBI rCRS JSON files in data/raw/.")

    nuccore_data = json.loads(RCRS_NUCCORE_JSON.read_text(encoding="utf-8"))
    rcrs_seq = nuccore_data[0]["sequence"].strip().upper()
    assert len(rcrs_seq) == 16569, f"Expected 16569 bp rCRS, got {len(rcrs_seq)}"

    prot_data = json.loads(RCRS_PROTEINS_JSON.read_text(encoding="utf-8"))
    proteins = prot_data["sequences"]
    assert len(proteins) == 13, f"Expected 13 human mtDNA proteins, got {len(proteins)}"

    # Rebuild clean database from canonical schema
    if DB_PATH.exists():
        DB_PATH.unlink()
    conn = sqlite3.connect(DB_PATH)
    conn.executescript(SCHEMA_PATH.read_text(encoding="utf-8"))

    # Insert World Space stressors
    conn.executemany(
        "INSERT INTO environmental_stressors (stressor_name, physical_variable, unit, notes) VALUES (?, ?, ?, ?)",
        INITIAL_STRESSORS,
    )

    # Insert Homo sapiens rCRS reference row in donor_species
    conn.execute(
        """
        INSERT INTO donor_species (
            scientific_name, common_name, taxonomic_clade, mitogenome_accession,
            stressor_axis, role, sister_pair_id, habitat_notes
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            "Homo sapiens",
            "Human (rCRS Reference)",
            "Primates",
            "NC_012920.1",
            "Reference",
            "Reference",
            None,
            "Revised Cambridge Reference Sequence (16,569 bp); primary recipient scaffold.",
        ),
    )

    rows = []
    for entry in proteins:
        acc = entry["header"].split()[0].lstrip(">")
        if acc not in GENE_METADATA:
            raise KeyError(f"Unexpected RefSeq protein accession: {acc}")
        gene_symbol, oxphos_complex, tier, uniprot_id, expected_aa_len = GENE_METADATA[acc]
        refseq_aa = entry["sequence"].strip().upper()
        assert len(refseq_aa) == expected_aa_len, (
            f"{gene_symbol}: length mismatch ({len(refseq_aa)} vs expected {expected_aa_len})"
        )

        start_nt, end_nt, strand, cds_nt = locate_cds_on_rcrs(rcrs_seq, refseq_aa)
        retranslated = translate_mito_cds(cds_nt)
        assert retranslated == refseq_aa, f"{gene_symbol}: Table 2 retranslation mismatch!"

        rows.append(
            (
                gene_symbol,
                oxphos_complex,
                tier,
                start_nt,
                end_nt,
                strand,
                uniprot_id,
                expected_aa_len,
                cds_nt,
                refseq_aa,
                None,
                None,
            )
        )

    # Sort by rCRS start coordinate for tidy reporting
    rows.sort(key=lambda r: r[3])
    conn.executemany(
        """
        INSERT INTO human_rcrs_genes (
            gene_symbol, oxphos_complex, tier, rcrs_start_nt, rcrs_end_nt, strand,
            uniprot_id, expected_aa_len, cds_nt_seq, protein_aa_seq,
            pdb_holocomplex_id, pdb_chain_id
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )
    conn.commit()
    conn.close()

    print(f"Initialized {DB_PATH} successfully (Human rCRS NC_012920.1 = {len(rcrs_seq)} bp).")
    print(f"{'Gene':<6} {'Complex':<12} {'Tier':<5} {'rCRS Coords':<15} {'Strand':<6} {'AA Len':<7} {'Start Codon':<11} {'UniProt'}")
    print("-" * 78)
    for r in rows:
        gene_symbol, oxphos_complex, tier, start_nt, end_nt, strand, uniprot_id, aa_len, cds_nt = r[:9]
        print(
            f"{gene_symbol:<6} {oxphos_complex:<12} T{tier:<4} {start_nt:>5}..{end_nt:<7} {strand:<6} {aa_len:<7} {cds_nt[:3]:<11} {uniprot_id}"
        )


if __name__ == "__main__":
    main()
