#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polite-http",
#   "python-dotenv",
# ]
# ///
"""Stage 1A: Fetch verified live NCBI RefSeq mitogenomes & all 13 OXPHOS proteins for the 44-species panel.

Uses the official `ncbi-sequence-fetch` wrapper (`ncbi_fetch.py`) for rate-limited E-utilities access.
Caches each species in `data/raw/mitogenomes/<slug>.json` so re-runs are instant and idempotent.
"""

from __future__ import annotations

import json
import re
import sqlite3
import sys
from pathlib import Path

# Import the official rate-limited ncbi_fetch wrapper
NCBI_SKILL_SCRIPTS = Path("/Users/o2satz/.gemini/config/plugins/science/skills/ncbi_sequence_fetch/scripts")
sys.path.insert(0, str(NCBI_SKILL_SCRIPTS))
import ncbi_fetch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from species_panel import SPECIES_PANEL  # noqa: E402

CACHE_DIR = ROOT / "data" / "raw" / "mitogenomes"
DB_PATH = ROOT / "data" / "mitofitness.db"


def slugify(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def classify_mito_protein(header: str, length: int) -> str | None:
    """Map an NCBI mitochondrial RefSeq protein FASTA header to one of the 13 canonical symbols."""
    h = header.lower()
    # Order matters: check ND4L before ND4, COX3/III before II/I, ATP8 before ATP6
    if "cytochrome b" in h or "cytb" in h or "cob" in h:
        return "CYTB"
    if re.search(r"subunit\s*4l|\bnd4l\b|\bnad4l\b", h):
        return "ND4L"
    if re.search(r"subunit\s*6|\bnd6\b|\bnad6\b", h) and ("nadh" in h or "nd6" in h or "nad6" in h):
        return "ND6"
    if re.search(r"subunit\s*5|\bnd5\b|\bnad5\b", h) and ("nadh" in h or "nd5" in h or "nad5" in h):
        return "ND5"
    if re.search(r"subunit\s*4|\bnd4\b|\bnad4\b", h) and ("nadh" in h or "nd4" in h or "nad4" in h):
        return "ND4"
    if re.search(r"subunit\s*3|\bnd3\b|\bnad3\b", h) and ("nadh" in h or "nd3" in h or "nad3" in h):
        return "ND3"
    if re.search(r"subunit\s*2|\bnd2\b|\bnad2\b", h) and ("nadh" in h or "nd2" in h or "nad2" in h):
        return "ND2"
    if re.search(r"subunit\s*1|\bnd1\b|\bnad1\b", h) and ("nadh" in h or "nd1" in h or "nad1" in h):
        return "ND1"
    if ("oxidase" in h or "cox" in h or "coiii" in h or "co3" in h) and re.search(r"subunit\s*(iii|3)\b|\bcox3\b|\bcoiii\b|\bco3\b", h):
        return "COX3"
    if ("oxidase" in h or "cox" in h or "coii" in h or "co2" in h) and re.search(r"subunit\s*(ii|2)\b|\bcox2\b|\bcoii\b|\bco2\b", h):
        return "COX2"
    if ("oxidase" in h or "cox" in h or "coi" in h or "co1" in h) and re.search(r"subunit\s*(i|1)\b|\bcox1\b|\bcoi\b|\bco1\b", h):
        return "COX1"
    if ("atp" in h) and re.search(r"subunit\s*8|\batp8\b|\batpase\s*8", h):
        return "ATP8"
    if ("atp" in h) and re.search(r"subunit\s*6|\batp6\b|\batpase\s*6", h):
        return "ATP6"
    # Fallback by length for ambiguous ATP/ND headers
    if "atp" in h:
        return "ATP8" if length < 100 else "ATP6"
    return None


def find_refseq_mitogenome(names: list[str]) -> tuple[str, str, str] | None:
    """Return (accession, uid, matched_organism_name) for the species' RefSeq complete mitochondrion."""
    for org in names:
        # 1. Prefer official NCBI RefSeq complete mitochondrion (NC_*)
        term = f'"{org}"[Organism] AND mitochondrion[Filter] AND srcdb_refseq[PROP] AND 15000:22000[SLEN]'
        uids, _ = ncbi_fetch.esearch("nuccore", term, retmax=5)
        if not uids:
            # 2. Broaden to complete mitochondrial genome in RefSeq or GenBank
            term = f'"{org}"[Organism] AND mitochondrion[Filter] AND "complete genome"[Title] AND 15000:22000[SLEN]'
            uids, _ = ncbi_fetch.esearch("nuccore", term, retmax=5)
        if uids:
            summaries = ncbi_fetch.esummary("nuccore", uids)
            # Pick NC_* first if available
            best_uid = None
            best_acc = None
            for uid in uids:
                info = summaries.get(uid, {})
                acc = info.get("accessionversion") or info.get("caption") or ""
                if acc.startswith("NC_"):
                    best_uid, best_acc = uid, acc
                    break
            if best_uid is None:
                best_uid = uids[0]
                info = summaries.get(best_uid, {})
                best_acc = info.get("accessionversion") or info.get("caption") or best_uid
            return best_acc, best_uid, org
    return None


def fetch_species_proteins(acc: str, uid: str) -> dict[str, dict]:
    """Fetch all 13 linked mitochondrial proteins for a nuccore mitogenome UID."""
    prot_ids = ncbi_fetch.elink("nuccore", "protein", uid)
    proteins: dict[str, dict] = {}
    if prot_ids:
        fasta_txt = ncbi_fetch.efetch("protein", prot_ids, rettype="fasta")
        for header, seq in ncbi_fetch.parse_fasta(fasta_txt):
            gene = classify_mito_protein(header, len(seq))
            if gene and gene not in proteins:
                proteins[gene] = {
                    "accession": header.split()[0].lstrip(">"),
                    "header": header,
                    "length": len(seq),
                    "sequence": seq.upper().strip(),
                }
    # Fallback: fasta_cds_aa directly from nuccore if elink missed any protein
    if len(proteins) < 13:
        cds_txt = ncbi_fetch.efetch("nuccore", acc, rettype="fasta_cds_aa")
        for header, seq in ncbi_fetch.parse_fasta(cds_txt):
            gene = classify_mito_protein(header, len(seq))
            if gene and gene not in proteins:
                proteins[gene] = {
                    "accession": header.split()[0].lstrip(">"),
                    "header": header,
                    "length": len(seq),
                    "sequence": seq.upper().strip(),
                }
    return proteins


def main() -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON;")

    fetched_count = 0
    cached_count = 0
    incomplete = []

    for sci_name, fallbacks, common_name, clade, axis, role, pair_ids, notes in SPECIES_PANEL:
        cache_file = CACHE_DIR / f"{slugify(sci_name)}.json"
        if cache_file.exists():
            rec = json.loads(cache_file.read_text(encoding="utf-8"))
            cached_count += 1
        else:
            hit = find_refseq_mitogenome([sci_name, *fallbacks])
            if hit is None:
                print(f"[ERROR] No mitogenome found for {sci_name} (fallbacks={fallbacks})")
                incomplete.append((sci_name, 0))
                continue
            acc, uid, matched_org = hit
            proteins = fetch_species_proteins(acc, uid)
            rec = {
                "scientific_name": sci_name,
                "matched_ncbi_organism": matched_org,
                "common_name": common_name,
                "taxonomic_clade": clade,
                "stressor_axis": axis,
                "role": role,
                "sister_pair_ids": pair_ids,
                "mitogenome_accession": acc,
                "ncbi_nuccore_uid": uid,
                "habitat_notes": notes,
                "protein_count": len(proteins),
                "proteins": proteins,
            }
            cache_file.write_text(json.dumps(rec, indent=2), encoding="utf-8")
            fetched_count += 1
            print(f"[FETCH] {sci_name:<26} -> {acc:<14} ({len(proteins)}/13 proteins) [{axis}]")

        if rec["protein_count"] < 13:
            missing = sorted(
                {"ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB", "COX1", "COX2", "COX3", "ATP6", "ATP8"}
                - set(rec["proteins"].keys())
            )
            incomplete.append((sci_name, rec["protein_count"], missing))

        # Upsert into donor_species table
        conn.execute(
            """
            INSERT INTO donor_species (
                scientific_name, common_name, taxonomic_clade, mitogenome_accession,
                stressor_axis, role, sister_pair_id, habitat_notes
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(scientific_name) DO UPDATE SET
                mitogenome_accession = excluded.mitogenome_accession,
                stressor_axis = excluded.stressor_axis,
                role = excluded.role,
                sister_pair_id = excluded.sister_pair_id,
                habitat_notes = excluded.habitat_notes
            """,
            (
                sci_name,
                common_name,
                clade,
                rec["mitogenome_accession"],
                axis,
                role,
                ",".join(pair_ids),
                notes,
            ),
        )

    conn.commit()
    total_db = conn.execute("SELECT COUNT(*) FROM donor_species").fetchone()[0]
    conn.close()

    print("\n=== Stage 1A Mitogenome Fetch Summary ===")
    print(f"Newly fetched: {fetched_count} | Loaded from cache: {cached_count} | Total donor_species in DB (incl. Human ref): {total_db}")
    if incomplete:
        print(f"WARNING — Incomplete species ({len(incomplete)}):")
        for item in incomplete:
            print(" ", item)
    else:
        print("SUCCESS — 100% of wild species fetched with complete 13/13 mitochondrial proteins!")


if __name__ == "__main__":
    main()
