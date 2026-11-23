#!/usr/bin/env python3
"""Stage 6A: Prepare 6 paired WT vs. Mutant mitonuclear dimer manifests for live Boltz-2 co-folding.

Uses exact Human rCRS mtDNA protein sequences (from data/mitofitness.db) paired with exact
cryo-EM-resolved mature human nuclear partner sequences (from PDB 9I6F, 8H9S, 9HZL) so that
Boltz-2 co-folded dimers align 1-to-1 with experimental cryo-EM coordinates for RMSD validation.
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
OUT_DIR = ROOT / "data" / "stage6_boltz2" / "inputs"


def extract_chain_seq(cif_gz: Path, chain_id: str) -> str:
    parser = MMCIFParser(QUIET=True)
    with gzip.open(cif_gz, "rt", encoding="utf-8") as fh:
        model = next(parser.get_structure("s", fh).get_models())
    return "".join(seq1(r.get_resname()) for r in model[chain_id] if r.id[0] == " ")


def make_af3_json(name: str, seq_a: str, seq_b: str) -> dict:
    return {
        "name": name,
        "dialect": "alphafold3",
        "version": 1,
        "modelSeeds": [20261122],
        "sequences": [
            {"protein": {"id": ["A"], "sequence": seq_a}},
            {"protein": {"id": ["B"], "sequence": seq_b}},
        ],
    }


def mutate(seq: str, pos1: int, expected_ref: str, alt: str) -> str:
    assert seq[pos1 - 1] == expected_ref, f"Expected {expected_ref} at {pos1}, got {seq[pos1 - 1]}"
    return seq[: pos1 - 1] + alt + seq[pos1:]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    mt_seqs = dict(conn.execute("SELECT gene_symbol, protein_aa_seq FROM human_rcrs_genes").fetchall())
    conn.close()

    # Exact mature nuclear partner chains resolved in experimental human cryo-EM
    cox6a1_seq = extract_chain_seq(ROOT / "data/structures/pdb_00009i6f.cif.gz", "G")
    atp5mj_seq = extract_chain_seq(ROOT / "data/structures/pdb_00008h9s.cif.gz", "P")
    uqcrfs1_seq = extract_chain_seq(ROOT / "data/structures/pdb_00009hzl.cif.gz", "C")

    pairs = [
        ("COX3_WT_COX6A1", mt_seqs["COX3"], cox6a1_seq, "9I6F", "C", "G"),
        ("COX3_W116R_COX6A1", mutate(mt_seqs["COX3"], 116, "W", "R"), cox6a1_seq, "9I6F", "C", "G"),
        ("ATP6_WT_ATP5MJ", mt_seqs["ATP6"], atp5mj_seq, "8H9S", "a", "P"),
        ("ATP6_I106V_ATP5MJ", mutate(mt_seqs["ATP6"], 106, "I", "V"), atp5mj_seq, "8H9S", "a", "P"),
        ("CYTB_WT_UQCRFS1", mt_seqs["CYTB"], uqcrfs1_seq, "9HZL", "b", "C"),
        ("CYTB_M53T_UQCRFS1", mutate(mt_seqs["CYTB"], 53, "M", "T"), uqcrfs1_seq, "9HZL", "b", "C"),
    ]

    manifest_index = []
    for name, seq_a, seq_b, pdb_id, ch_a, ch_b in pairs:
        obj = make_af3_json(name, seq_a, seq_b)
        p = OUT_DIR / f"{name}.json"
        p.write_text(json.dumps(obj, indent=2), encoding="utf-8")
        manifest_index.append({
            "name": name,
            "json_path": str(p),
            "len_A": len(seq_a),
            "len_B": len(seq_b),
            "total_aa": len(seq_a) + len(seq_b),
            "ref_pdb": pdb_id,
            "ref_chain_A": ch_a,
            "ref_chain_B": ch_b,
        })
        print(f"Wrote {p.name}: A={len(seq_a)} aa, B={len(seq_b)} aa (total {len(seq_a)+len(seq_b)} aa)")

    (OUT_DIR / "manifest_index.json").write_text(json.dumps(manifest_index, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
