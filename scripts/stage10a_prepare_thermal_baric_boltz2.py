#!/usr/bin/env python3
"""Stage 10A: Prepare 4 paired WT vs. Mutant Boltz-2 YAML manifests for our top 2
Stage-9A Thermal & Baric Mitonuclear Dimers:
  1. ND4L_WT_NDUFS2   vs. ND4L_Y85H_NDUFS2  (Desert_Heat Jerboa 0-bystander lead, 9I4I chains k + Q)
  2. CYTB_WT_CYC1     vs. CYTB_Y224H_CYC1   (Deep-Sea Dive Sperm Whale lead, 9HZL chains J + H)

Uses exact Human rCRS mtDNA protein sequences (data/mitofitness.db) paired with exact
cryo-EM-resolved mature human nuclear partner sequences so Boltz-2 predictions align
1-to-1 with experimental human cryo-EM coordinates (9I4I & 9HZL).
"""

from __future__ import annotations

import gzip
import json
import sqlite3
from pathlib import Path

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data" / "mitofitness.db"
YAML_DIR = ROOT / "data" / "stage10_boltz2" / "yamls"


def extract_chain_seq(cif_gz: Path, chain_id: str) -> str:
    parser = MMCIFParser(QUIET=True)
    with gzip.open(cif_gz, "rt", encoding="utf-8") as fh:
        model = next(parser.get_structure("s", fh).get_models())
    return "".join(seq1(r.get_resname()) for r in model[chain_id] if r.id[0] == " ")


def mutate(seq: str, pos1: int, expected_ref: str, alt: str) -> str:
    assert seq[pos1 - 1] == expected_ref, f"Expected {expected_ref} at {pos1}, got {seq[pos1 - 1]}"
    return seq[: pos1 - 1] + alt + seq[pos1:]


def write_boltz_yaml(path: Path, seq_a: str, seq_b: str) -> None:
    content = (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {seq_a}\n"
        "  - protein:\n"
        "      id: B\n"
        f"      sequence: {seq_b}\n"
    )
    path.write_text(content, encoding="utf-8")


def main() -> None:
    YAML_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    mt_seqs = dict(conn.execute("SELECT gene_symbol, protein_aa_seq FROM human_rcrs_genes").fetchall())
    conn.close()

    ndufs2_seq = extract_chain_seq(ROOT / "data/structures/pdb_00009i4i.cif.gz", "Q")
    cyc1_seq = extract_chain_seq(ROOT / "data/structures/pdb_00009hzl.cif.gz", "H")

    pairs = [
        ("ND4L_WT_NDUFS2",   mt_seqs["ND4L"],                         ndufs2_seq, "9I4I", "k", "Q"),
        ("ND4L_Y85H_NDUFS2", mutate(mt_seqs["ND4L"], 85, "Y", "H"),   ndufs2_seq, "9I4I", "k", "Q"),
        ("CYTB_WT_CYC1",     mt_seqs["CYTB"],                         cyc1_seq,   "9HZL", "J", "H"),
        ("CYTB_Y224H_CYC1",  mutate(mt_seqs["CYTB"], 224, "Y", "H"),  cyc1_seq,   "9HZL", "J", "H"),
    ]

    index = []
    for name, seq_a, seq_b, pdb_id, ch_a, ch_b in pairs:
        yp = YAML_DIR / f"{name}.yaml"
        write_boltz_yaml(yp, seq_a, seq_b)
        index.append({
            "name": name,
            "yaml_path": str(yp),
            "len_A": len(seq_a),
            "len_B": len(seq_b),
            "total_aa": len(seq_a) + len(seq_b),
            "ref_pdb": pdb_id,
            "ref_chain_A": ch_a,
            "ref_chain_B": ch_b,
        })
        print(f"Wrote {yp.name}: A={len(seq_a)} aa, B={len(seq_b)} aa (total {len(seq_a)+len(seq_b)} aa)")

    (YAML_DIR.parent / "manifest_index.json").write_text(json.dumps(index, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
