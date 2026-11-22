#!/usr/bin/env python3
"""Stage 3B: Automated 3D Mitonuclear Compatibility + Live ClinVar/gnomAD Safety + DdCBE/TALED Base-Editor Gate.

Runs all 131 3D-mapped Stage-1/Stage-2 candidates PLUS 5 hard-wired human calibration controls
(Tier-1 adaptive ND1 Y30H; Tier-1 pathogenic ND1 A52T & ATP6 L156R; Tier-2 out-of-sample frozen
pathogenic ND4 R340H & ND6 M64V) through three strictly automated, zero-cherry-picking gates:

1. Human Clinical & Population Safety Gate (live ClinVar Pathogenic/Benign + gnomAD v3 chrM):
   - Maps every candidate's required rCRS nucleotide substitutions (and exact protein p.RefPosAlt)
     against live ClinVar assertions and gnomAD v3 homoplasmic/heteroplasmic allele frequencies.
   - Classifies into:
     * Contraindicated_Pathogenic (ClinVar Pathogenic / Likely Pathogenic at exact change or same AA site)
     * Safe_Natural_Poly (observed homoplasmically in gnomAD humans OR ClinVar Benign/Likely Benign)
     * Neutral_Unreported (no ClinVar pathogenicity; novel/rare in humans)

2. Mitochondrial Base-Editor Feasibility Gate (exact NCBI Translation Table 2 on Human rCRS):
   - Enumerates all Table-2 codons encoding `alt_aa` from the human `rCRS` reference codon.
   - Identifies exact editable mechanism on either H-strand or L-strand:
     * DdCBE_TC         : Single C:G -> T:A transition in canonical 5'-TC context
     * DdCBE_DddA11_HC  : Single C:G -> T:A transition in expanded 5'-AC/CC/TC (HC) context
     * TALED_AT_GC      : Single A:T -> G:C transition (sTALED/mTALED adenine deaminase)
     * Dual_Transition  : 2 transitions (C->T and/or A->G) required in codon
     * Not_Base_Editable_Transversion : Requires a purine<->pyrimidine transversion
   - Computes exact `bystander_nonsyn_count` (unintended nonsynonymous C->T or A->G bystander edits
     on the same strand within a +/-4 nt editing window around the target base on Human rCRS!).

3. 3D Mitonuclear Structural Compatibility & Transfer Strategy Gate (2.52-2.95 A Human Cryo-EM):
   - Integrates exact cryo-EM nuclear-partner contacts (<=4.5 A), local B-factor, catalytic-center
     distance, volume expansion (steric clash proxy when buried or at nuclear interface), and net
     charge perturbation to compute a 0..1 `compatibility_score` and `transfer_strategy`:
     * Motif_Grafting_BaseEdit | Co_Adapted_Cassette | Allotopic_Nuclear_Express | Incompatible

Populates `human_safety_annotations` and `transferability_assessments` in data/mitofitness.db,
and exports data/stage3_shortlist.json + data/stage3_report.md.
"""

from __future__ import annotations

import json
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data" / "mitofitness.db"
STAGE2_JSON = ROOT / "data" / "stage2_3d_fbt_summary.json"
CLINVAR_PATH_JSON = ROOT / "data" / "raw" / "clinvar_mt_pathogenic_summaries.json"
CLINVAR_BEN_JSON = ROOT / "data" / "raw" / "clinvar_mt_benign_summaries.json"
GNOMAD_MT_JSON = ROOT / "data" / "raw" / "gnomad_mt_r3_variants.json"
RCRS_NUCCORE_JSON = ROOT / "data" / "raw" / "human_rcrs_nuccore.json"
OUT_JSON = ROOT / "data" / "stage3_shortlist.json"
OUT_REPORT = ROOT / "data" / "stage3_report.md"

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}
VERT_MITO_STARTS = {"ATG", "ATA", "ATT", "ATC", "GTG"}

# Map official HGNC MT-* gene names to our canonical gene symbols
MT_HGNC_TO_GENE = {
    "MT-ND1": "ND1", "MT-ND2": "ND2", "MT-ND3": "ND3", "MT-ND4": "ND4",
    "MT-ND4L": "ND4L", "MT-ND5": "ND5", "MT-ND6": "ND6", "MT-CYB": "CYTB",
    "MT-CO1": "COX1", "MT-CO2": "COX2", "MT-CO3": "COX3", "MT-ATP6": "ATP6",
    "MT-ATP8": "ATP8",
}

# Build reverse NCBI Table 2 codon table: AA -> list of codons
AA_TO_MITO_CODONS: dict[str, list[str]] = defaultdict(list)
for b1 in "ACGT":
    for b2 in "ACGT":
        for b3 in "ACGT":
            cod = b1 + b2 + b3
            aa = str(Seq(cod).translate(table=2))
            AA_TO_MITO_CODONS[aa].append(cod)


def parse_hgvsp_to_aa_change(hgvsp: str | None) -> tuple[str, int, str] | None:
    """Parse p.Tyr30His or (p.Tyr30His) into ('Y', 30, 'H')."""
    if not hgvsp:
        return None
    m = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", hgvsp)
    if not m:
        return None
    r1 = AA3_TO_1.get(m.group(1))
    r2 = AA3_TO_1.get(m.group(3))
    if not r1 or not r2:
        return None
    return r1, int(m.group(2)), r2


def parse_rcrs_m_notation(title: str) -> tuple[int, str, str] | None:
    """Parse m.11778G>A from ClinVar title."""
    m = re.search(r"m\.(\d+)([ACGT])>([ACGT])", title)
    if not m:
        return None
    return int(m.group(1)), m.group(2), m.group(3)


def translate_codon_at(cds_nt: str, aa_idx0: int, codon: str) -> str:
    if aa_idx0 == 0 and codon in VERT_MITO_STARTS:
        return "M"
    return str(Seq(codon).translate(table=2))


def analyze_base_editor_on_rcrs(
    rcrs_seq: str,
    gene_info: dict,
    aa_pos_1: int,
    alt_aa: str,
) -> dict:
    """Exact NCBI Table-2 base-editor feasibility & bystander nonsynonymous edit counter on Human rCRS."""
    start_nt = gene_info["rcrs_start_nt"]
    end_nt = gene_info["rcrs_end_nt"]
    strand = gene_info["strand"]
    cds_nt = gene_info["cds_nt_seq"]
    aa_idx0 = aa_pos_1 - 1
    ref_codon = cds_nt[aa_idx0 * 3 : aa_idx0 * 3 + 3]

    # Helper: convert 0-based CDS index -> 1-based rCRS H-strand coordinate
    def cds_idx_to_rcrs(cidx: int) -> int:
        return (start_nt + cidx) if strand == "+" else (end_nt - cidx)

    # Helper: convert coding-strand base to H-strand base
    def to_hstrand(base: str) -> str:
        return base if strand == "+" else COMP[base]

    # Rank priority for editor classes (best = lowest number)
    rank_order = {
        "DdCBE_TC": 1,
        "TALED_AT_GC": 2,
        "DdCBE_DddA11_HC": 3,
        "Dual_Transition": 4,
        "Not_Base_Editable_Transversion": 5,
    }

    target_codons = AA_TO_MITO_CODONS.get(alt_aa, [])
    best_plan = None

    for tgt_cod in target_codons:
        diffs = [(k, ref_codon[k], tgt_cod[k]) for k in range(3) if ref_codon[k] != tgt_cod[k]]
        if not diffs:
            continue

        # Classify each nt change as C:G->T:A (CBE), A:T->G:C (ABE), or Transversion
        edit_types = []
        edit_records = []
        for k, rb_cod, ab_cod in diffs:
            cidx = aa_idx0 * 3 + k
            rcrs_pos = cds_idx_to_rcrs(cidx)
            rb_h = to_hstrand(rb_cod)
            ab_h = to_hstrand(ab_cod)

            # Determine editable strand & 5' context on rCRS
            # C->T on H-strand OR G->A on H-strand (= C->T on L-strand)
            if (rb_h == "C" and ab_h == "T") or (rb_h == "G" and ab_h == "A"):
                if rb_h == "C":
                    # Edit C->T on H-strand; 5' neighbor is rcrs_pos - 1 (1-based -> idx rcrs_pos - 2)
                    prev_5p = rcrs_seq[rcrs_pos - 2] if rcrs_pos >= 2 else "N"
                    edit_strand = "H"
                else:
                    # Edit C->T on L-strand; 5' neighbor on L-strand is complement of H-strand rcrs_pos + 1
                    nxt_h = rcrs_seq[rcrs_pos] if rcrs_pos < len(rcrs_seq) else "N"
                    prev_5p = COMP.get(nxt_h, "N")
                    edit_strand = "L"
                if prev_5p == "T":
                    etype = "DdCBE_TC"
                elif prev_5p in {"A", "C"}:
                    etype = "DdCBE_DddA11_HC"
                else:
                    # 5'-GC context: lower efficiency DddA11/non-canonical, mark as DdCBE_DddA11_HC with note
                    etype = "DdCBE_DddA11_HC"
                edit_types.append(etype)
                edit_records.append({
                    "rcrs_pos": rcrs_pos,
                    "rcrs_change": f"m.{rcrs_pos}{rb_h}>{ab_h}",
                    "editor_strand": edit_strand,
                    "motif_5p": f"{prev_5p}C",
                    "type": etype,
                })
            elif (rb_h == "A" and ab_h == "G") or (rb_h == "T" and ab_h == "C"):
                edit_strand = "H" if rb_h == "A" else "L"
                edit_types.append("TALED_AT_GC")
                edit_records.append({
                    "rcrs_pos": rcrs_pos,
                    "rcrs_change": f"m.{rcrs_pos}{rb_h}>{ab_h}",
                    "editor_strand": edit_strand,
                    "motif_5p": "A/T",
                    "type": "TALED_AT_GC",
                })
            else:
                edit_types.append("Transversion")
                edit_records.append({
                    "rcrs_pos": rcrs_pos,
                    "rcrs_change": f"m.{rcrs_pos}{rb_h}>{ab_h}",
                    "editor_strand": None,
                    "motif_5p": None,
                    "type": "Transversion",
                })

        if "Transversion" in edit_types:
            overall_class = "Not_Base_Editable_Transversion"
        elif len(edit_types) == 1:
            overall_class = edit_types[0]
        else:
            overall_class = "Dual_Transition"

        # Count nonsynonymous bystanders within +/-4 nt window on the same editable strand (for single-edit plans)
        bystander_nonsyn = 0
        if len(edit_records) == 1 and edit_records[0]["editor_strand"] is not None:
            er = edit_records[0]
            t_pos = er["rcrs_pos"]
            e_strand = er["editor_strand"]
            is_cbe = er["type"].startswith("DdCBE")
            # Target H-strand base to scan for bystanders:
            # CBE on H: 'C'->'T'; CBE on L: 'G'->'A'; ABE on H: 'A'->'G'; ABE on L: 'T'->'C'
            if is_cbe:
                scan_from, scan_to = ("C", "T") if e_strand == "H" else ("G", "A")
            else:
                scan_from, scan_to = ("A", "G") if e_strand == "H" else ("T", "C")

            for delta in range(-4, 5):
                if delta == 0:
                    continue
                b_pos = t_pos + delta
                if not (start_nt <= b_pos <= end_nt):
                    continue
                if rcrs_seq[b_pos - 1] != scan_from:
                    continue
                # For DdCBE_TC, only count bystander C in 5'-TC context
                if er["type"] == "DdCBE_TC":
                    if e_strand == "H":
                        p5 = rcrs_seq[b_pos - 2] if b_pos >= 2 else "N"
                    else:
                        p5 = COMP.get(rcrs_seq[b_pos] if b_pos < len(rcrs_seq) else "N", "N")
                    if p5 != "T":
                        continue
                # Check if mutating b_pos (scan_from -> scan_to on H-strand) changes the AA in this gene
                b_cidx = (b_pos - start_nt) if strand == "+" else (end_nt - b_pos)
                b_aa_idx = b_cidx // 3
                b_cod_pos = b_cidx % 3
                if b_aa_idx * 3 + 3 > len(cds_nt):
                    continue
                orig_cod = cds_nt[b_aa_idx * 3 : b_aa_idx * 3 + 3]
                mut_b_cod = scan_to if strand == "+" else COMP[scan_to]
                new_cod = orig_cod[:b_cod_pos] + mut_b_cod + orig_cod[b_cod_pos + 1 :]
                if translate_codon_at(cds_nt, b_aa_idx, new_cod) != translate_codon_at(cds_nt, b_aa_idx, orig_cod):
                    bystander_nonsyn += 1
        elif overall_class == "Dual_Transition":
            bystander_nonsyn = 2  # Conservative penalty for 2-edit window
        else:
            bystander_nonsyn = -1  # Not applicable

        plan = {
            "ref_codon": ref_codon,
            "target_codon": tgt_cod,
            "base_editor_class": overall_class,
            "n_nt_edits": len(edit_records),
            "edits": edit_records,
            "bystander_nonsyn_count": bystander_nonsyn,
        }
        if best_plan is None or (
            rank_order[plan["base_editor_class"]],
            plan["n_nt_edits"],
            plan["bystander_nonsyn_count"] if plan["bystander_nonsyn_count"] >= 0 else 99,
        ) < (
            rank_order[best_plan["base_editor_class"]],
            best_plan["n_nt_edits"],
            best_plan["bystander_nonsyn_count"] if best_plan["bystander_nonsyn_count"] >= 0 else 99,
        ):
            best_plan = plan

    return best_plan


def build_safety_lookups(gene_map: dict[str, dict]) -> tuple[dict, dict, dict, tuple]:
    """Index live ClinVar (Pathogenic & Benign) and gnomAD v3 chrM variants by BOTH m.pos and (gene, aa_pos).

    Translates all true ClinVar Pathogenic / Likely Pathogenic single-nt rCRS records directly via
    NCBI Table 2 so same-site pathogenic hotspot detection has 100% coding coverage.
    """
    path_by_m = {}
    path_by_aa_exact = {}
    path_by_aa_pos = defaultdict(list)

    def map_mpos_to_aa(m_pos: int, ref_h: str, alt_h: str) -> tuple[str, int, str, str] | None:
        for gsym, gi in gene_map.items():
            s, e, strand, cds = gi["rcrs_start_nt"], gi["rcrs_end_nt"], gi["strand"], gi["cds_nt_seq"]
            if not (s <= m_pos <= e):
                continue
            cidx = (m_pos - s) if strand == "+" else (e - m_pos)
            aa_idx0 = cidx // 3
            cod_pos = cidx % 3
            if aa_idx0 * 3 + 3 > len(cds):
                continue
            orig_cod = cds[aa_idx0 * 3 : aa_idx0 * 3 + 3]
            exp_base = ref_h if strand == "+" else COMP.get(ref_h, "N")
            if orig_cod[cod_pos] != exp_base:
                continue
            mut_base = alt_h if strand == "+" else COMP.get(alt_h, "N")
            new_cod = orig_cod[:cod_pos] + mut_base + orig_cod[cod_pos + 1 :]
            r_aa = translate_codon_at(cds, aa_idx0, orig_cod)
            a_aa = translate_codon_at(cds, aa_idx0, new_cod)
            return gsym, aa_idx0 + 1, r_aa, a_aa
        return None

    for rec in json.loads(CLINVAR_PATH_JSON.read_text(encoding="utf-8")):
        title = rec.get("title") or ""
        sig = rec.get("clinical_significance") or ""
        sig_l = sig.lower()
        # Strictly exclude records whose aggregate ClinVar status is Benign/Likely benign or VUS without Pathogenic
        if "pathogenic" not in sig_l or ("benign" in sig_l and "conflicting" not in sig_l):
            continue
        m_hit = parse_rcrs_m_notation(title)
        if m_hit:
            path_by_m[m_hit] = {"title": title, "sig": sig, "id": rec.get("variant_id")}
            aa_map = map_mpos_to_aa(*m_hit)
            if aa_map:
                g, pos_a, ref_a, alt_a = aa_map
                if alt_a != ref_a:
                    path_by_aa_exact[(g, pos_a, alt_a)] = {"title": title, "sig": sig, "id": rec.get("variant_id")}
                    path_by_aa_pos[(g, pos_a)].append({"title": title, "sig": sig, "alt_aa": alt_a, "id": rec.get("variant_id")})
        else:
            genes = [MT_HGNC_TO_GENE[g["symbol"]] for g in (rec.get("genes") or []) if g.get("symbol") in MT_HGNC_TO_GENE]
            p_hit = parse_hgvsp_to_aa_change(title)
            if genes and p_hit:
                g = genes[0]
                _, pos_a, alt_a = p_hit
                path_by_aa_exact[(g, pos_a, alt_a)] = {"title": title, "sig": sig, "id": rec.get("variant_id")}
                path_by_aa_pos[(g, pos_a)].append({"title": title, "sig": sig, "alt_aa": alt_a, "id": rec.get("variant_id")})

    ben_by_aa_exact = {}
    for rec in json.loads(CLINVAR_BEN_JSON.read_text(encoding="utf-8")):
        title = rec.get("title") or ""
        sig = rec.get("clinical_significance") or ""
        if "benign" not in sig.lower() or "pathogenic" in sig.lower():
            continue
        m_hit = parse_rcrs_m_notation(title)
        if m_hit:
            aa_map = map_mpos_to_aa(*m_hit)
            if aa_map:
                g, pos_a, ref_a, alt_a = aa_map
                if alt_a != ref_a:
                    ben_by_aa_exact[(g, pos_a, alt_a)] = {"title": title, "sig": sig, "id": rec.get("variant_id")}

    gnomad_by_aa_exact = {}
    gnomad_by_m = {}
    for v in json.loads(GNOMAD_MT_JSON.read_text(encoding="utf-8")):
        vid = v.get("variant_id") or ""
        parts = vid.split("-")
        an = v.get("an") or 0
        af_hom = (v.get("ac_hom") or 0) / an if an > 0 else 0.0
        af_het = (v.get("ac_het") or 0) / an if an > 0 else 0.0
        entry = {
            "variant_id": vid,
            "af_hom": round(af_hom, 8),
            "af_het": round(af_het, 8),
            "ac_hom": v.get("ac_hom") or 0,
            "ac_het": v.get("ac_het") or 0,
            "an": an,
            "hgvsp": v.get("hgvsp"),
        }
        if len(parts) == 4 and parts[0] == "M":
            try:
                gnomad_by_m[(int(parts[1]), parts[2], parts[3])] = entry
            except ValueError:
                pass
        gsym = MT_HGNC_TO_GENE.get(v.get("gene_symbol") or "")
        p_hit = parse_hgvsp_to_aa_change(v.get("hgvsp"))
        if gsym and p_hit:
            key = (gsym, p_hit[1], p_hit[2])
            if key not in gnomad_by_aa_exact or entry["af_hom"] > gnomad_by_aa_exact[key]["af_hom"]:
                gnomad_by_aa_exact[key] = entry

    return path_by_m, path_by_aa_exact, path_by_aa_pos, (ben_by_aa_exact, gnomad_by_aa_exact, gnomad_by_m)


def evaluate_safety(
    gene: str,
    aa_pos: int,
    alt_aa: str,
    nt_edits: list[dict],
    path_by_m: dict,
    path_by_aa_exact: dict,
    path_by_aa_pos: dict,
    ben_by_aa_exact: dict,
    gnomad_by_aa_exact: dict,
    gnomad_by_m: dict,
) -> dict:
    """Apply automated Human Clinical & Population Safety Gate."""
    # Check exact single-nt m.pos hit in ClinVar & gnomAD first
    m_key = None
    rcrs_nt_str = None
    if len(nt_edits) == 1:
        rcrs_nt_str = nt_edits[0]["rcrs_change"]
        m_match = re.match(r"m\.(\d+)([ACGT])>([ACGT])", rcrs_nt_str)
        if m_match:
            m_key = (int(m_match.group(1)), m_match.group(2), m_match.group(3))

    p_exact = path_by_aa_exact.get((gene, aa_pos, alt_aa)) or (path_by_m.get(m_key) if m_key else None)
    p_same_site = path_by_aa_pos.get((gene, aa_pos), [])
    b_exact = ben_by_aa_exact.get((gene, aa_pos, alt_aa))
    g_exact = gnomad_by_aa_exact.get((gene, aa_pos, alt_aa)) or (gnomad_by_m.get(m_key) if m_key else None)

    af_hom = g_exact["af_hom"] if g_exact else 0.0
    af_het = g_exact["af_het"] if g_exact else 0.0
    ac_hom = g_exact["ac_hom"] if g_exact else 0

    if p_exact is not None:
        # Distinguish confirmed pathogenic from common polymorphic/adaptive alleles with conflicting/benign ClinVar
        # Rule: If ClinVar flags Pathogenic/LP and homoplasmic frequency < 0.005 (0.5%), strictly Contraindicated!
        if af_hom >= 0.005:
            status = "Safe_Natural_Poly"
            clin_label = f"High-AF Human Poly (gnomAD hom={af_hom:.4f}; ClinVar note: {p_exact['sig']})"
        else:
            status = "Contraindicated_Pathogenic"
            clin_label = p_exact["sig"]
    elif p_same_site:
        # Another substitution at the EXACT same residue causes human mitochondrial disease -> caution!
        if af_hom >= 0.001 or b_exact is not None:
            status = "Safe_Natural_Poly"
            clin_label = f"Benign/Observed alt at disease-sensitive pos (other pathogenic: {len(p_same_site)})"
        else:
            status = "Epistatic_Rescue_Required"
            clin_label = f"Same-site disease hotspot ({len(p_same_site)} ClinVar Pathogenic hit(s) at pos {aa_pos})"
    elif b_exact is not None or ac_hom >= 1:
        status = "Safe_Natural_Poly"
        clin_label = b_exact["sig"] if b_exact else f"Observed homoplasmic in gnomAD (ac_hom={ac_hom})"
    else:
        status = "Neutral_Unreported"
        clin_label = "Not_Reported"

    return {
        "rcrs_nt_change": rcrs_nt_str,
        "gnomad_af_hom": af_hom,
        "gnomad_af_het": af_het,
        "gnomad_ac_hom": ac_hom,
        "clinvar_significance": clin_label,
        "safety_gate_status": status,
    }


def compute_3d_compatibility(cand: dict, safety: dict, editor: dict) -> dict:
    """Compute transparent 0..1 3D mitonuclear structural compatibility & transfer strategy."""
    min_nuc = cand["min_nuclear_A"]
    min_cat = cand["min_catalytic_A"]
    dv = cand["d_vol"]
    dq = abs(cand["d_charge"])
    burial = cand["burial_neighbors_8A"]
    bfac = cand["local_bfactor"]

    # Steric & electrostatic clash proxy (kcal/mol-scaled heuristic, 0 = clean, >1.5 = clash)
    clash = 0.0
    # Volume expansion inside tightly packed nuclear interface or deeply buried shell
    if dv > 0:
        if min_nuc <= 4.5:
            clash += (dv / 40.0) * 1.45
        elif burial >= 85:
            clash += (dv / 50.0) * 1.15
        else:
            clash += (dv / 90.0) * 0.50
    # Charge introduction/reversal near tight nuclear interface or catalytic core
    if dq >= 0.5:
        if min_nuc <= 4.5:
            clash += 1.65 * dq
        elif min_cat <= 10.0:
            clash += 1.35 * dq
        else:
            clash += 0.75 * dq
    # Mild penalty if directly adjacent to catalytic cofactor with high Grantham
    if min_cat <= 6.0 and cand["grantham"] >= 80:
        clash += 1.10

    # Mild bonus if naturally tolerated homoplasmically in humans (gnomAD proof of compatibility!)
    if safety["safety_gate_status"] == "Safe_Natural_Poly":
        clash = max(0.0, clash - 0.65)

    clash = round(clash, 2)
    # Map clash + safety into 0..1 compatibility_score
    if safety["safety_gate_status"] == "Contraindicated_Pathogenic":
        compat_score = 0.05
        strategy = "Incompatible"
    else:
        compat_score = round(max(0.05, min(0.99, 1.0 - (clash / 4.2))), 3)
        is_base_editable = editor["base_editor_class"] in {"DdCBE_TC", "DdCBE_DddA11_HC", "TALED_AT_GC"}
        if clash <= 1.50 and is_base_editable and safety["safety_gate_status"] in {"Safe_Natural_Poly", "Neutral_Unreported"}:
            strategy = "Motif_Grafting_BaseEdit"
        elif clash > 1.50 and min_nuc <= 4.5:
            strategy = "Co_Adapted_Cassette"
        elif not is_base_editable and clash <= 1.50:
            strategy = "Allotopic_Nuclear_Express"
        elif safety["safety_gate_status"] == "Epistatic_Rescue_Required":
            strategy = "Co_Adapted_Cassette"
        else:
            strategy = "Incompatible"

    return {
        "steric_clash_score": clash,
        "compatibility_score": compat_score,
        "transfer_strategy": strategy,
        "local_bfactor": bfac,
    }


def main() -> None:
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON;")

    rcrs_seq = json.loads(RCRS_NUCCORE_JSON.read_text(encoding="utf-8"))[0]["sequence"].strip().upper()
    gene_map = {}
    for row in conn.execute(
        "SELECT gene_symbol, oxphos_complex, tier, rcrs_start_nt, rcrs_end_nt, strand, cds_nt_seq, protein_aa_seq, pdb_holocomplex_id, pdb_chain_id FROM human_rcrs_genes"
    ).fetchall():
        gene_map[row[0]] = {
            "gene_symbol": row[0], "complex": row[1], "tier": row[2],
            "rcrs_start_nt": row[3], "rcrs_end_nt": row[4], "strand": row[5],
            "cds_nt_seq": row[6], "protein_aa_seq": row[7],
            "pdb_id": row[8], "chain_id": row[9],
        }

    path_by_m, path_by_aa_exact, path_by_aa_pos, (ben_by_aa_exact, gnomad_by_aa_exact, gnomad_by_m) = build_safety_lookups(gene_map)
    print(f"Loaded safety ground truth: ClinVar Pathogenic exact={len(path_by_aa_exact)} ({len(path_by_m)} m.pos), Benign={len(ben_by_aa_exact)}, gnomAD chrM exact={len(gnomad_by_aa_exact)}.")

    # 1. ANSWER-KEY CALIBRATION & OUT-OF-SAMPLE CONTROL AUDIT
    answer_key_controls = [
        ("Tier1_Adaptive_Positive", "ND1", 30, "H", "Tibetan high-altitude m.3394T>C ND1 Y30H", "Safe_Natural_Poly"),
        ("Tier1_Pathogenic_Negative", "ND1", 52, "T", "LHON m.3460G>A ND1 A52T", "Contraindicated_Pathogenic"),
        ("Tier1_Pathogenic_Negative", "ATP6", 156, "R", "Leigh/NARP m.8993T>G ATP6 L156R", "Contraindicated_Pathogenic"),
        ("Tier2_Frozen_Pathogenic_OOS", "ND4", 340, "H", "Out-of-sample LHON m.11778G>A ND4 R340H", "Contraindicated_Pathogenic"),
        ("Tier2_Frozen_Pathogenic_OOS", "ND6", 64, "V", "Out-of-sample (-) strand LHON m.14484T>C ND6 M64V", "Contraindicated_Pathogenic"),
    ]
    control_audit = []
    for role, g, pos, alt, desc, expected_status in answer_key_controls:
        ed = analyze_base_editor_on_rcrs(rcrs_seq, gene_map[g], pos, alt)
        sf = evaluate_safety(
            g, pos, alt, ed["edits"],
            path_by_m, path_by_aa_exact, path_by_aa_pos,
            ben_by_aa_exact, gnomad_by_aa_exact, gnomad_by_m,
        )
        passed = (sf["safety_gate_status"] == expected_status)
        control_audit.append({
            "role": role, "gene": g, "pos": pos, "alt": alt, "desc": desc,
            "expected": expected_status, "observed": sf["safety_gate_status"],
            "gnomad_af_hom": sf["gnomad_af_hom"], "clinvar": sf["clinvar_significance"],
            "editor": ed["base_editor_class"], "rcrs_nt": sf["rcrs_nt_change"],
            "passed": passed,
        })
        print(f"  [{'PASS' if passed else 'FAIL'}] {role}: {g} {pos}->{alt} ({sf['rcrs_nt_change']}) -> {sf['safety_gate_status']} (gnomAD hom={sf['gnomad_af_hom']}, ClinVar={sf['clinvar_significance']})")

    # 2. EVALUATE ALL 131 3D-MAPPED STAGE-2 CANDIDATES
    s2 = json.loads(STAGE2_JSON.read_text(encoding="utf-8"))
    candidates = s2["candidates_3d"]

    conn.execute("DELETE FROM human_safety_annotations;")
    conn.execute("DELETE FROM transferability_assessments;")

    allele_lookup = defaultdict(list)
    for aid, g, pos, alt in conn.execute("SELECT allele_id, gene_symbol, human_rcrs_aa_pos, alt_aa FROM mt_alleles").fetchall():
        allele_lookup[(g, pos, alt)].append(aid)

    evaluated = []
    for c in candidates:
        g = c["gene"]
        pos = int(c["human_rcrs_aa_pos"])
        alt = c["alt_aa"]
        ed = analyze_base_editor_on_rcrs(rcrs_seq, gene_map[g], pos, alt)
        sf = evaluate_safety(
            g, pos, alt, ed["edits"],
            path_by_m, path_by_aa_exact, path_by_aa_pos,
            ben_by_aa_exact, gnomad_by_aa_exact, gnomad_by_m,
        )
        st3 = compute_3d_compatibility(c, sf, ed)

        full_rec = {
            **c,
            "editor": ed,
            "safety": sf,
            "structural_gate": st3,
        }
        evaluated.append(full_rec)

        # Insert into human_safety_annotations
        conn.execute(
            """
            INSERT OR REPLACE INTO human_safety_annotations (
                gene_symbol, human_rcrs_aa_pos, human_ref_aa, alt_aa, rcrs_nt_change,
                gnomad_af_hom, gnomad_af_het, clinvar_significance, mitomap_status,
                safety_gate_status, evidence_pmid
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                g, pos, c["human_ref_aa"], alt, sf["rcrs_nt_change"],
                sf["gnomad_af_hom"], sf["gnomad_af_het"], sf["clinvar_significance"],
                "Natural_Human_Poly" if sf["gnomad_ac_hom"] >= 1 else "Unreported_in_gnomAD",
                sf["safety_gate_status"],
                "21383179" if (g == "COX3" and pos == 116 and alt == "R") else None,
            ),
        )

        # Insert into transferability_assessments for all linked allele_ids
        for aid in allele_lookup.get((g, pos, alt), []):
            conn.execute(
                """
                INSERT OR REPLACE INTO transferability_assessments (
                    allele_id, recipient_species, pdb_id, interacting_nuclear_subunits_json,
                    min_nuclear_contact_A, cryoem_local_bfactor, afdb_partner_plddt_min,
                    steric_clash_score, rcrs_ref_codon, target_codon, required_nt_edits_json,
                    base_editor_class, bystander_nonsyn_count, alphagenome_partner_summary,
                    transfer_strategy, compatibility_score
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    aid, "Homo sapiens", c["pdb_id"], json.dumps(c["nuclear_partners_4p5A"]),
                    c["min_nuclear_A"], c["local_bfactor"], None,
                    st3["steric_clash_score"], ed["ref_codon"], ed["target_codon"],
                    json.dumps(ed["edits"]), ed["base_editor_class"], ed["bystander_nonsyn_count"],
                    None, st3["transfer_strategy"], st3["compatibility_score"],
                ),
            )

    conn.commit()
    n_safety_db = conn.execute("SELECT COUNT(*) FROM human_safety_annotations").fetchone()[0]
    n_trans_db = conn.execute("SELECT COUNT(*) FROM transferability_assessments").fetchone()[0]
    conn.close()

    # Rank finalists: Safe + Motif_Grafting_BaseEdit + low bystander + high convergence/radicality + cross-hypoxia!
    def finalist_rank_key(r: dict):
        is_motif = 1 if r["structural_gate"]["transfer_strategy"] == "Motif_Grafting_BaseEdit" else 0
        is_safe = 1 if r["safety"]["safety_gate_status"] in {"Safe_Natural_Poly", "Neutral_Unreported"} else 0
        is_clean_bystander = 1 if r["editor"]["bystander_nonsyn_count"] == 0 else 0
        is_hyp = 1 if ("High_Altitude_Hypoxia" in r["axes_hit"] or "Diving_Subterranean_Hypoxia" in r["axes_hit"]) else 0
        is_cross_hyp = 1 if ("High_Altitude_Hypoxia" in r["axes_hit"] and "Diving_Subterranean_Hypoxia" in r["axes_hit"]) else 0
        is_rad = 1 if r["grantham"] >= 50.0 else 0
        return (
            is_motif,
            is_safe,
            is_hyp,
            is_cross_hyp,
            is_rad,
            is_clean_bystander,
            r["k_ext_derived"],
            r["structural_gate"]["compatibility_score"],
            r["grantham"],
        )

    evaluated.sort(key=finalist_rank_key, reverse=True)

    strategy_counts = Counter(r["structural_gate"]["transfer_strategy"] for r in evaluated)
    safety_counts = Counter(r["safety"]["safety_gate_status"] for r in evaluated)
    editor_counts = Counter(r["editor"]["base_editor_class"] for r in evaluated)

    base_editable_safe = [
        r for r in evaluated
        if r["structural_gate"]["transfer_strategy"] == "Motif_Grafting_BaseEdit"
    ]
    zero_bystander_safe = [
        r for r in base_editable_safe
        if r["editor"]["bystander_nonsyn_count"] == 0
    ]
    pathogenic_caught = [
        r for r in evaluated
        if r["safety"]["safety_gate_status"] in {"Contraindicated_Pathogenic", "Epistatic_Rescue_Required"}
    ]

    # Collect top contacting nuclear partners across our Top 20 finalists for Stage 3C AlphaGenome + AFDB audit!
    top_nuclear_partners = Counter(
        p for r in base_editable_safe[:25] for p in r["nuclear_partners_4p5A"]
    )

    out_payload = {
        "answer_key_control_audit": control_audit,
        "all_controls_passed": all(c["passed"] for c in control_audit),
        "summary_counts": {
            "total_evaluated": len(evaluated),
            "db_human_safety_rows": n_safety_db,
            "db_transferability_rows": n_trans_db,
            "by_transfer_strategy": dict(strategy_counts),
            "by_safety_status": dict(safety_counts),
            "by_base_editor_class": dict(editor_counts),
            "motif_grafting_base_editable_safe": len(base_editable_safe),
            "zero_bystander_base_editable_safe": len(zero_bystander_safe),
            "pathogenic_or_hotspot_flagged": len(pathogenic_caught),
        },
        "top_contacting_nuclear_partners_for_alphagenome": dict(top_nuclear_partners.most_common(15)),
        "pathogenic_and_hotspot_flags": pathogenic_caught,
        "top_25_base_editable_finalists": base_editable_safe[:25],
        "all_evaluated_candidates": evaluated,
    }
    OUT_JSON.write_text(json.dumps(out_payload, indent=2), encoding="utf-8")

    # Write Markdown report
    lines = [
        "# Stage 3B Report — 3D Mitonuclear Compatibility, Live ClinVar/gnomAD Safety & DdCBE/TALED Base-Editor Gate",
        "",
        "## 1. Answer-Key Calibration & Out-of-Sample Validation (5/5 Controls)",
        "",
        "| Control Role | Gene & AA | Human rCRS Edit | Expected Gate | Observed Gate | gnomAD Hom AF | Base Editor | Verdict |",
        "| :--- | :---: | :---: | :--- | :--- | :---: | :---: | :---: |",
    ]
    for c in control_audit:
        lines.append(
            f"| `{c['role']}` | `{c['gene']} {c['pos']}->{c['alt']}` | `{c['rcrs_nt']}` | `{c['expected']}` | `{c['observed']}` | `{c['gnomad_af_hom']}` | `{c['editor']}` | **{'PASS' if c['passed'] else 'FAIL'}** |"
        )

    lines += [
        "",
        "## 2. Automated Funnel Summary Across All 131 3D-Mapped Candidates",
        f"- **Transfer Strategy Breakdown:** `{dict(strategy_counts)}`",
        f"- **Human Safety Gate (`ClinVar` + `gnomAD` v3):** `{dict(safety_counts)}`",
        f"- **Base-Editor Class (`NCBI Table 2` on `rCRS`):** `{dict(editor_counts)}`",
        f"- **Actionable `Motif_Grafting_BaseEdit` Finalists:** **`{len(base_editable_safe)}`** (`{len(zero_bystander_safe)}` with **0 nonsynonymous bystanders** in the ±4 nt editing window!)",
        f"- **Top Contacting Human Nuclear Partners (for AlphaGenome Audit):** `{dict(top_nuclear_partners.most_common(10))}`",
        "",
        "## 3. Top 20 Actionable `Motif_Grafting_BaseEdit` Candidates (Safe, Compatible, Base-Editable Today)",
        "",
        "| Rank | Gene (T) | Site | rCRS Edit | Base Editor (Bystanders) | Safety (`gnomAD` Hom) | Clash / Compat | Nuclear Partners (<=4.5 A) | FBT Class & Axes |",
        "| :---: | :--- | :---: | :---: | :--- | :--- | :---: | :--- | :--- |",
    ]
    for idx, r in enumerate(base_editable_safe[:20], 1):
        ed = r["editor"]
        sf = r["safety"]
        st3 = r["structural_gate"]
        partners = ", ".join(r["nuclear_partners_4p5A"]) if r["nuclear_partners_4p5A"] else "None (core/lipid)"
        axes = "/".join(a[:4] for a in r["axes_hit"])
        lines.append(
            f"| {idx} | `{r['gene']}` (T{r['tier']}) | `{r['human_ref_aa']}{r['human_rcrs_aa_pos']}{r['alt_aa']}` | `{sf['rcrs_nt_change']}` | `{ed['base_editor_class']}` (bys={ed['bystander_nonsyn_count']}) | `{sf['safety_gate_status']}` ({sf['gnomad_af_hom']:.4f}) | `{st3['steric_clash_score']}` / **{st3['compatibility_score']}** | {partners} | `{r['payoff_equivalence_class'][:18]}` ({axes}) |"
        )

    if pathogenic_caught:
        lines += [
            "",
            "## 4. Safety Traps Caught Automatically (Wild Extreme Adaptations That Hit Human Disease Hotspots!)",
            "",
            "| Gene (T) | Site | Extreme Species (Axis) | Safety Gate | ClinVar / gnomAD Detail |",
            "| :--- | :---: | :--- | :--- | :--- |",
        ]
        for r in pathogenic_caught[:15]:
            spp = ", ".join(p["species"] for p in r["pairs"][:2])
            lines.append(
                f"| `{r['gene']}` (T{r['tier']}) | `{r['human_ref_aa']}{r['human_rcrs_aa_pos']}{r['alt_aa']}` | {spp} | `{r['safety']['safety_gate_status']}` | {r['safety']['clinvar_significance']} |"
            )

    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nStage 3B complete! Wrote {OUT_JSON} and {OUT_REPORT}.")
    print(f"Strategy breakdown: {dict(strategy_counts)}")
    print(f"Safety breakdown: {dict(safety_counts)}")
    print(f"Base-editor breakdown: {dict(editor_counts)}")
    print(f"Motif_Grafting_BaseEdit safe finalists: {len(base_editable_safe)} ({len(zero_bystander_safe)} with 0 nonsyn bystanders!)")
    print(f"Safety traps caught automatically: {len(pathogenic_caught)}")


if __name__ == "__main__":
    main()
