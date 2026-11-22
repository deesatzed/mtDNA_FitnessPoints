#!/usr/bin/env python3
"""Stage 3C: Integrate Live AlphaFold DB v6 Audits, Live AlphaGenome AVI Scans & Export AF3/Boltz-2 Manifests.

- Summarizes official AlphaFold DB v6 (pLDDT + PAE) audits for all 5 Stage-4 compensatory buffers:
  Ciona AOX (H2XKF9), Yeast NDI1 (P32340), Human SOD2 (P04179), Human UCP2 (P55851), Human PINK1 (Q9BXM7).
- Summarizes live AlphaGenome Variant Impact (AVI) promoter saturation mutagenesis scans for contacting
  nuclear OXPHOS partners (NDUFA13, UQCRFS1, ATP5MJ) and endogenous buffers (SOD2, UCP2), updating
  `transferability_assessments.alphagenome_partner_summary` in data/mitofitness.db.
- Auto-exports ready-to-run AlphaFold 3 JSON and Boltz-2 YAML co-folding manifests in
  `data/cofolding_manifests/` for the Top 5 Safe, Compatible, Base-Editable mtDNA Finalists.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data" / "mitofitness.db"
AFDB_DIR = ROOT / "data" / "afdb"
AG_AVI_DIR = ROOT / "data" / "raw" / "alphagenome_avi"
AG_GTF_JSON = ROOT / "data" / "raw" / "alphagenome_gtf" / "all_9_nuclear_genes.json"
COFOLD_DIR = ROOT / "data" / "cofolding_manifests"
OUT_JSON = ROOT / "data" / "stage3c_bioai_summary.json"
OUT_REPORT = ROOT / "data" / "stage3c_bioai_report.md"

BUFFER_DEFS = [
    ("Ciona_AOX", "H2XKF9", "Ciona intestinalis Alternative Oxidase (Complex III/IV bypass)", "132-276 (145 AA rigid catalytic core; 1-131 flexible MTS/N-tail)"),
    ("Yeast_NDI1", "P32340", "S. cerevisiae Internal NADH Dehydrogenase (Complex I bypass)", "46-513 (468 AA rigid FAD core; 1-45 MTS)"),
    ("Human_SOD2", "P04179", "Human Mitochondrial Mn-Superoxide Dismutase (ROS sink)", "24-222 (199 AA rigid core; 1-23 cleaved MTS)"),
    ("Human_UCP2", "P55851", "Human Mitochondrial Uncoupling Protein 2 (proton-leak valve)", "57-304 (248 AA dynamic transmembrane carrier)"),
    ("Human_PINK1", "Q9BXM7", "Human PTEN-induced Kinase 1 (mitophagy sensor/effector)", "268-544 (277 AA rigid kinase domain; 1-267 MTS/TM sensor)"),
]

TOP_5_COFOLD_FINALISTS = [
    {
        "slug": "finalist1_ND4L_A11T",
        "gene": "ND4L", "aa_pos": 11, "ref_aa": "A", "alt_aa": "T",
        "rcrs_edit": "m.10500G>A", "editor": "DdCBE_DddA11_HC (0 bystanders)",
        "rationale": "Tier-2 Frozen Out-of-Sample Bird (Tibetan ground tit) + Mammal (Blind mole-rat) cross-hypoxia convergence at 95.6%-conserved site.",
        "ccd_ligands": [],
    },
    {
        "slug": "finalist2_COX3_W116R",
        "gene": "COX3", "aa_pos": 116, "ref_aa": "W", "alt_aa": "R",
        "rcrs_edit": "m.9552T>C", "editor": "TALED_AT_GC (0 bystanders)",
        "rationale": "Canonical Bar-headed goose high-altitude hypoxia benchmark (97.8% invariant) — unexpectedly single-step TALED-editable on human rCRS with 0 bystanders.",
        "ccd_ligands": ["HEA"],
    },
    {
        "slug": "finalist3_COX3_I165T",
        "gene": "COX3", "aa_pos": 165, "ref_aa": "I", "alt_aa": "T",
        "rcrs_edit": "m.9700T>C", "editor": "TALED_AT_GC (1 bystander)",
        "rationale": "Snow leopard (High Altitude) + Naked mole-rat (Subterranean Anoxia) cross-hypoxia convergence; Safe_Natural_Poly observed homoplasmically in gnomAD humans.",
        "ccd_ligands": ["HEA"],
    },
    {
        "slug": "finalist4_ATP6_I106V",
        "gene": "ATP6", "aa_pos": 106, "ref_aa": "I", "alt_aa": "V",
        "rcrs_edit": "m.8842A>G", "editor": "TALED_AT_GC (0 bystanders)",
        "rationale": "3-way cross-hypoxia convergence (Plateau pika + Sperm whale + Cuvier's beaked whale; 3 vs 0 controls); Safe_Natural_Poly (gnomAD hom AF=0.0023); contacts nuclear ATP5MJ.",
        "ccd_ligands": ["ATP", "MG"],
    },
    {
        "slug": "finalist5_CYTB_M53T",
        "gene": "CYTB", "aa_pos": 53, "ref_aa": "M", "alt_aa": "T",
        "rcrs_edit": "m.14904T>C", "editor": "TALED_AT_GC (0 bystanders)",
        "rationale": "High-Altitude rule-breaker contacting nuclear Rieske Fe-S protein UQCRFS1 (<=4.5 A); TALED-editable with 0 bystanders.",
        "ccd_ligands": ["HEM", "U10"],
    },
]


def summarize_afdb_buffers() -> list[dict]:
    out = []
    for name, uid, role, domain_note in BUFFER_DEFS:
        meta_path = AFDB_DIR / uid / f"AF-{uid}-F1-metadata.json"
        meta_raw = json.loads(meta_path.read_text(encoding="utf-8"))
        meta = meta_raw[0] if isinstance(meta_raw, list) else meta_raw
        out.append({
            "buffer_name": name,
            "uniprot_id": uid,
            "role": role,
            "sequence_length": len(meta.get("uniprotSequence") or ""),
            "global_plddt": round(float(meta.get("globalMetricValue") or 0.0), 2),
            "frac_very_low_lt50": round(float(meta.get("fractionPlddtVeryLow") or 0.0), 3),
            "frac_very_high_gt90": round(float(meta.get("fractionPlddtVeryHigh") or 0.0), 3),
            "pae_rigid_domain_summary": domain_note,
            "cif_path": f"data/afdb/{uid}/AF-{uid}-F1-model_v6.cif",
        })
    return out


def summarize_alphagenome_scans() -> dict[str, dict]:
    gtf = json.loads(AG_GTF_JSON.read_text(encoding="utf-8"))
    out = {}
    for g in ["NDUFA13", "UQCRFS1", "ATP5MJ", "SOD2", "UCP2"]:
        scan_path = AG_AVI_DIR / f"{g}_promoter.json"
        records = json.loads(scan_path.read_text(encoding="utf-8"))
        # Normalize records (list of dicts)
        if isinstance(records, dict):
            records = records.get("variants") or records.get("results") or []
        phreds = []
        top_var = None
        top_phred = -1.0
        top_mod = None
        for r in records:
            ph = float(r.get("phred") or r.get("avi_phred") or r.get("score_phred") or 0.0)
            phreds.append(ph)
            if ph > top_phred:
                top_phred = ph
                top_var = r.get("variant") or r.get("variant_id")
                top_mod = r.get("top_modality") or r.get("modality") or r.get("top_feature")
        n_high = sum(1 for p in phreds if p >= 15.0)
        g_info = gtf[g]
        tss = g_info["start"] if g_info["strand"] == "+" else g_info["end"]
        out[g] = {
            "gene": g,
            "role": "Nuclear_OXPHOS_Partner" if g in {"NDUFA13", "UQCRFS1", "ATP5MJ"} else "Endogenous_Human_Buffer",
            "chromosome": g_info["chromosome"],
            "tss_hg38": tss,
            "strand": g_info["strand"],
            "scanned_window": f"{g_info['chromosome']}:{tss-20}-{tss+20}",
            "snvs_above_phred10": len(records),
            "snvs_above_phred15": n_high,
            "max_avi_phred": round(top_phred, 2) if top_phred >= 0 else None,
            "top_variant": top_var,
            "top_modality": top_mod,
            "atlas_entity_url": f"https://deepmind.google.com/science/alphagenome/atlas?q={g}&m=entity",
        }
    return out


def export_cofolding_manifests(conn: sqlite3.Connection) -> list[dict]:
    COFOLD_DIR.mkdir(parents=True, exist_ok=True)
    seqs = dict(conn.execute("SELECT gene_symbol, protein_aa_seq FROM human_rcrs_genes").fetchall())
    exported = []
    for f in TOP_5_COFOLD_FINALISTS:
        wt_seq = seqs[f["gene"]]
        idx0 = f["aa_pos"] - 1
        assert wt_seq[idx0] == f["ref_aa"], f"Mismatch for {f['slug']}: {wt_seq[idx0]} != {f['ref_aa']}"
        mut_seq = wt_seq[:idx0] + f["alt_aa"] + wt_seq[idx0 + 1 :]

        # AlphaFold 3 JSON input specification (mutant human mtDNA chain + ligands)
        af3_seq_entries = [
            {"protein": {"id": ["A"], "sequence": mut_seq}},
        ]
        for l_idx, ccd in enumerate(f["ccd_ligands"], 1):
            af3_seq_entries.append({"ligand": {"id": [f"L{l_idx}"], "ccdCodes": [ccd]}})

        af3_payload = {
            "name": f"{f['slug']}_human_rCRS",
            "dialect": "alphafold3",
            "version": 1,
            "modelSeeds": [20261122],
            "sequences": af3_seq_entries,
        }
        af3_path = COFOLD_DIR / f"{f['slug']}_af3.json"
        af3_path.write_text(json.dumps(af3_payload, indent=2), encoding="utf-8")

        # Boltz-2 YAML input specification
        yaml_lines = [
            "version: 1",
            "sequences:",
            "  - protein:",
            "      id: A",
            f"      sequence: {mut_seq}",
        ]
        for l_idx, ccd in enumerate(f["ccd_ligands"], 1):
            yaml_lines += [
                "  - ligand:",
                f"      id: L{l_idx}",
                f"      ccd: {ccd}",
            ]
        boltz_path = COFOLD_DIR / f"{f['slug']}_boltz2.yaml"
        boltz_path.write_text("\n".join(yaml_lines) + "\n", encoding="utf-8")

        exported.append({
            "slug": f["slug"],
            "mutation": f"{f['gene']} {f['ref_aa']}{f['aa_pos']}{f['alt_aa']} ({f['rcrs_edit']})",
            "editor": f["editor"],
            "af3_json": str(af3_path.relative_to(ROOT)),
            "boltz2_yaml": str(boltz_path.relative_to(ROOT)),
            "ligands": f["ccd_ligands"],
        })
    return exported


def main() -> None:
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON;")

    afdb_summary = summarize_afdb_buffers()
    ag_summary = summarize_alphagenome_scans()

    # Update transferability_assessments.alphagenome_partner_summary for alleles touching NDUFA13, UQCRFS1, ATP5MJ
    updated_rows = 0
    for row_id, partners_json in conn.execute(
        "SELECT assessment_id, interacting_nuclear_subunits_json FROM transferability_assessments"
    ).fetchall():
        partners = json.loads(partners_json or "[]")
        notes = []
        for p in partners:
            if p in ag_summary:
                info = ag_summary[p]
                notes.append(
                    f"{p} ({info['scanned_window']}): max AVI Phred={info['max_avi_phred']} ({info['snvs_above_phred15']} SNVs>=15)"
                )
        if notes:
            conn.execute(
                "UPDATE transferability_assessments SET alphagenome_partner_summary = ? WHERE assessment_id = ?",
                ("; ".join(notes), row_id),
            )
            updated_rows += 1
    conn.commit()

    cofold_exports = export_cofolding_manifests(conn)
    conn.close()

    payload = {
        "alphafold_db_v6_buffers": afdb_summary,
        "alphagenome_avi_promoter_scans": ag_summary,
        "db_transferability_rows_updated_with_alphagenome": updated_rows,
        "cofolding_manifests_exported": cofold_exports,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # Write Markdown summary report
    lines = [
        "# Stage 3C Report — AlphaFold DB v6 Structural Audit, Live AlphaGenome AVI Promoter Scans & Co-Folding Manifests",
        "",
        "> *Research-Use Framing:* All AlphaGenome AVI scores reflect predicted molecular/regulatory functional impact (never clinical diagnosis or medical advice).",
        "",
        "## 1. AlphaFold DB v6 Structural Audit of Stage-4 Compensatory Buffers (`pLDDT` + `PAE`)",
        "",
        "| Buffer | UniProt | Length (AA) | Global pLDDT | Very Low (<50) | Very High (>90) | PAE Rigid Domain & Engineering Guardrail |",
        "| :--- | :---: | :---: | :---: | :---: | :---: | :--- |",
    ]
    for b in afdb_summary:
        lines.append(
            f"| `{b['buffer_name']}` | `{b['uniprot_id']}` | {b['sequence_length']} | **{b['global_plddt']}** | {b['frac_very_low_lt50']*100:.1f}% | {b['frac_very_high_gt90']*100:.1f}% | {b['pae_rigid_domain_summary']} |"
        )

    lines += [
        "",
        "## 2. Live AlphaGenome Variant Impact (`AVI`) Core Promoter Saturation Scans (`TSS ± 20 bp`, `123` SNVs/gene)",
        "",
        "| Nuclear Gene | Role | Scanned Locus (hg38) | SNVs Phred>=10 | SNVs Phred>=15 | Max AVI Phred | Top Variant | Atlas Link |",
        "| :--- | :--- | :--- | :---: | :---: | :---: | :--- | :--- |",
    ]
    for g, st in ag_summary.items():
        lines.append(
            f"| `{g}` | {st['role']} | `{st['scanned_window']}` | {st['snvs_above_phred10']} | **{st['snvs_above_phred15']}** | **{st['max_avi_phred']}** | `{st['top_variant']}` | [Atlas]({st['atlas_entity_url']}) |"
        )

    lines += [
        "",
        f"*(Updated `{updated_rows}` matching rows in SQLite `transferability_assessments.alphagenome_partner_summary`.)*",
        "",
        "## 3. Exported AlphaFold 3 JSON & Boltz-2 YAML Co-Folding Manifests (Top 5 Actionable Finalists)",
        "",
        "| Finalist | Mutation & rCRS Edit | Base Editor | Cofactors | AF3 JSON | Boltz-2 YAML |",
        "| :--- | :--- | :--- | :---: | :--- | :--- |",
    ]
    for c in cofold_exports:
        ligs = ", ".join(c["ligands"]) if c["ligands"] else "None"
        lines.append(
            f"| `{c['slug']}` | `{c['mutation']}` | `{c['editor']}` | `{ligs}` | `{c['af3_json']}` | `{c['boltz2_yaml']}` |"
        )

    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Stage 3C complete! Updated {updated_rows} DB rows; exported {len(cofold_exports)} AF3/Boltz-2 manifest pairs.")
    for g, st in ag_summary.items():
        print(f"  AlphaGenome {g:<8}: Phred>=10: {st['snvs_above_phred10']:>3} | Phred>=15: {st['snvs_above_phred15']:>3} | max={st['max_avi_phred']} ({st['top_variant']})")


if __name__ == "__main__":
    main()
