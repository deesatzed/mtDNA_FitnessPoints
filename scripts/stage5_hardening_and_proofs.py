#!/usr/bin/env python3
"""Stage 5: Comprehensive Scientific Hardening, Adversarial Stress-Testing & Formal Proofs.

Executes 6 non-circular mathematical, statistical, 3D structural, and dynamical proofs across
Stages 1-4 using 100% real cached datasets (44 wild vertebrate mitogenomes, 4 human cryo-EM
holocomplexes 9I4I/9HZL/9I6F/8H9S, 18,164 gnomAD v3 chrM variants, 513 ClinVar Pathogenic +
1,891 ClinVar Benign summaries, and the Stage-4 biophysical ODE engine):

  Proof 1: Confounder-Controlled & Negative-Control Mitonuclear Interface Enrichment Proof
           (Matched Control-Branch Null + Gene-Stratified CMH + Burial/Entropy-Stratified
           Permutation N=5,000 + Multivariate Logistic Regression + Leave-One-Axis/Gene-Out)
  Proof 2: Full-Proteome Human ClinVar Pathogenic vs. Benign/gnomAD Benchmark & Strict
           Out-of-Sample Tier-1 (6 genes) -> Frozen Tier-2 (7 genes) Transfer AUROC
  Proof 3: 3D Atomic Neighborhood & Backbone-Independent Side-Chain Rotamer Steric-Clash Proof
           directly on experimental cryo-EM coordinates (5 Finalists vs 2 Pathogenic Controls)
  Proof 4: Base-Editor Window Sensitivity Stress Test (+/-2, +/-3, +/-4, +/-5, +/-6 nt on rCRS)
  Proof 5: Stage-4 Formal Supermodularity Theorem Grid, Physical Invariants, 3 Alternative
           Model-Forms (Hill ROS, 30-min severe I/R, 10-min mild I/R), Adversarial Placebo &
           Toxic-Uncoupler Controls, and Wet-Lab Heteroplasmy Dose-Response Sweep (0%-100%)
  Proof 6: End-to-End SQLite Referential Integrity & Cross-Stage Parity Audit

Outputs:
  - data/stage5_hardening_proofs.json
  - data/stage5_hardening_report.md
  - figures/stage5_hardening_and_proofs.png
"""

from __future__ import annotations

import gzip
import json
import math
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.Seq import Seq
from scipy.spatial import cKDTree
from scipy.stats import chi2, fisher_exact, mannwhitneyu, qmc
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import stage1_convergent_mining as s1
import stage2_3d_fbt_clustering as s2
import stage3_structural_safety_editor_gate as s3
import stage4_epistasis_synergy_uq as s4

DB_PATH = ROOT / "data" / "mitofitness.db"
OUT_JSON = ROOT / "data" / "stage5_hardening_proofs.json"
OUT_REPORT = ROOT / "data" / "stage5_hardening_report.md"
OUT_FIG = ROOT / "figures" / "stage5_hardening_and_proofs.png"

SEED = 20261122
RNG = np.random.default_rng(SEED)


# =============================================================================
# SHARED HELPERS: Load 3D microenvironments & MSA stats for all 3,756 residues
# =============================================================================
def load_all_3d_and_msa() -> tuple[dict[tuple[str, int], dict], dict[tuple[str, int], dict], dict[str, dict], dict[str, dict]]:
    """Return:
      - env3d: (gene, aa_pos_1) -> full Stage-2 3D cryo-EM microenvironment dict
      - msa_stats: (gene, aa_pos_1) -> {consensus_aa, consensus_frac, entropy, k_ext, k_ctl, ...}
      - control_candidates: (gene, aa_pos_1) -> matched Control-lineage polarized negative controls
      - gene_info: gene -> {complex, tier, rcrs_start_nt, rcrs_end_nt, strand, cds_nt_seq, protein_seq, uniprot_id}
    """
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    gene_info = {}
    uniprot_to_mt = {}
    rcrs_seqs = {}
    for row in conn.execute("SELECT * FROM human_rcrs_genes"):
        g = row["gene_symbol"]
        gene_info[g] = dict(row)
        uniprot_to_mt[row["uniprot_id"]] = g
        rcrs_seqs[g] = row["protein_aa_seq"]
    conn.close()

    meta = json.loads(s2.META_JSON.read_text(encoding="utf-8"))
    chain_meta = s2.build_chain_metadata(meta, uniprot_to_mt)

    env3d: dict[tuple[str, int], dict] = {}
    for _cname, (pid, cpath) in s2.COMPLEX_PDBS.items():
        g_env, _ = s2.extract_3d_microenvironments(pid, cpath, chain_meta[pid], rcrs_seqs)
        env3d.update(g_env)

    # Load all 44 wild species + Human rCRS proteins and rebuild exact Stage-1 MSA stats + Control-branch negative controls
    species_records: dict[str, dict] = {}
    for sci_name, *_ in s1.SPECIES_PANEL:
        p = s1.CACHE_DIR / f"{s1.slugify(sci_name)}.json"
        if p.exists():
            species_records[sci_name] = json.loads(p.read_text(encoding="utf-8"))

    active_pairs = [
        pair for pair in s1.SISTER_PAIRS
        if pair[2] in species_records and pair[3] in species_records
    ]
    aligner = s1.make_aligner()

    msa_stats: dict[tuple[str, int], dict] = {}
    control_candidates: dict[tuple[str, int], dict] = {}

    for gene, hseq in rcrs_seqs.items():
        aligned_by_sp: dict[str, list[str]] = {}
        for sci_name, rec in species_records.items():
            prot_entry = rec["proteins"].get(gene)
            if prot_entry and prot_entry.get("sequence"):
                aligned_by_sp[sci_name] = s1.align_to_human_rcrs(aligner, hseq, prot_entry["sequence"])
        n_aligned = len(aligned_by_sp)

        for pos0 in range(len(hseq)):
            h_ref = hseq[pos0]
            col_wild = [seq[pos0] for seq in aligned_by_sp.values()]
            gap_frac = sum(1 for r in col_wild if r == "-") / max(1, n_aligned)
            if gap_frac > s1.MAX_GAP_FRACTION:
                continue
            valid_col = [r for r in ([h_ref] + col_wild) if r in s1.VALID_AA]
            if not valid_col:
                continue
            counts = Counter(valid_col)
            cons_aa, cons_cnt = counts.most_common(1)[0]
            cons_frac = cons_cnt / len(valid_col)
            ent = s1.shannon_entropy(valid_col)

            # Count derived shifts in Extreme vs Control branches across active sister pairs
            k_ext_der = 0
            k_ctl_der = 0
            ctl_radical_hits = []
            for pair_id, _axis, ext_sp, ctl_sp, *_ in active_pairs:
                if ext_sp not in aligned_by_sp or ctl_sp not in aligned_by_sp:
                    continue
                ea = aligned_by_sp[ext_sp][pos0]
                ca = aligned_by_sp[ctl_sp][pos0]
                if ea not in s1.VALID_AA or ca not in s1.VALID_AA or ea == ca:
                    continue
                if ca == cons_aa and ea != cons_aa:
                    k_ext_der += 1
                elif ea == cons_aa and ca != cons_aa:
                    k_ctl_der += 1
                    gdist = s1.grantham_distance(cons_aa, ca)
                    if gdist >= 50.0 and cons_frac >= 0.85:
                        ctl_radical_hits.append((pair_id, ca, gdist))

            key = (gene, pos0 + 1)
            msa_stats[key] = {
                "human_aa": h_ref,
                "consensus_aa": cons_aa,
                "consensus_frac": round(cons_frac, 4),
                "entropy": round(ent, 4),
                "gap_frac": round(gap_frac, 4),
                "k_ext_derived": k_ext_der,
                "k_ctl_derived": k_ctl_der,
            }

            # Mirror Stage-1 Class I & II rules symmetrically on CONTROL branches as a strict negative control!
            net_ctl_bias = k_ctl_der - k_ext_der
            is_ctl_class1 = bool(ctl_radical_hits and k_ext_der == 0)
            is_ctl_class2 = bool(k_ctl_der >= 2 and net_ctl_bias >= 2)
            if is_ctl_class1 or is_ctl_class2:
                control_candidates[key] = {
                    "is_ctl_class1": is_ctl_class1,
                    "is_ctl_class2": is_ctl_class2,
                    "k_ctl_derived": k_ctl_der,
                    "k_ext_derived": k_ext_der,
                }

    return env3d, msa_stats, control_candidates, gene_info


# =============================================================================
# PROOF 1: Confounder-Controlled & Negative-Control Mitonuclear Interface Proof
# =============================================================================
def run_proof1_interface_enrichment(
    env3d: dict[tuple[str, int], dict],
    msa_stats: dict[tuple[str, int], dict],
    control_candidates: dict[tuple[str, int], dict],
    s2_data: dict,
) -> dict:
    ext_keys = {(c["gene"], c["human_rcrs_aa_pos"]) for c in s2_data["candidates_3d"]}
    ext_axes_map = {(c["gene"], c["human_rcrs_aa_pos"]): set(c["axes_hit"]) for c in s2_data["candidates_3d"]}
    ctl_keys = set(control_candidates.keys()) & set(env3d.keys())

    # Universe = all 3D-resolved human mtDNA residues with MSA QC coverage
    universe = sorted(set(env3d.keys()) & set(msa_stats.keys()))

    def calc_2x2(pos_set: set[tuple[str, int]], univ: list[tuple[str, int]]) -> dict:
        a = sum(1 for k in univ if k in pos_set and env3d[k]["min_nuclear_A"] <= 4.5)
        b = sum(1 for k in univ if k in pos_set and env3d[k]["min_nuclear_A"] > 4.5)
        c = sum(1 for k in univ if k not in pos_set and env3d[k]["min_nuclear_A"] <= 4.5)
        d = sum(1 for k in univ if k not in pos_set and env3d[k]["min_nuclear_A"] > 4.5)
        odds, p_two = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        _, p_gt = fisher_exact([[a, b], [c, d]], alternative="greater")
        return {
            "n_candidates": a + b,
            "n_interface": a,
            "pct_interface": round(100.0 * a / max(1, a + b), 2),
            "bg_pct_interface": round(100.0 * c / max(1, c + d), 2),
            "odds_ratio": round(float(odds), 3),
            "fisher_p_two_sided": round(float(p_two), 5),
            "fisher_p_greater": round(float(p_gt), 5),
        }

    ext_2x2 = calc_2x2(ext_keys, universe)
    ctl_2x2 = calc_2x2(ctl_keys, universe)

    # Direct head-to-head Extreme vs Matched Control-Branch candidates
    a_ec = sum(1 for k in ext_keys if k in env3d and env3d[k]["min_nuclear_A"] <= 4.5)
    b_ec = sum(1 for k in ext_keys if k in env3d and env3d[k]["min_nuclear_A"] > 4.5)
    c_ec = sum(1 for k in ctl_keys if env3d[k]["min_nuclear_A"] <= 4.5)
    d_ec = sum(1 for k in ctl_keys if env3d[k]["min_nuclear_A"] > 4.5)
    or_ec, p_ec = fisher_exact([[a_ec, b_ec], [c_ec, d_ec]], alternative="two-sided")

    # Gene-Stratified Cochran-Mantel-Haenszel (CMH) Test
    genes = sorted({k[0] for k in universe})
    cmh_num = 0.0
    cmh_den = 0.0
    cmh_stat_num = 0.0
    cmh_stat_var = 0.0
    for g in genes:
        gu = [k for k in universe if k[0] == g]
        a = sum(1 for k in gu if k in ext_keys and env3d[k]["min_nuclear_A"] <= 4.5)
        b = sum(1 for k in gu if k in ext_keys and env3d[k]["min_nuclear_A"] > 4.5)
        c = sum(1 for k in gu if k not in ext_keys and env3d[k]["min_nuclear_A"] <= 4.5)
        d = sum(1 for k in gu if k not in ext_keys and env3d[k]["min_nuclear_A"] > 4.5)
        n = a + b + c + d
        if n <= 1:
            continue
        cmh_num += (a * d) / n
        cmh_den += (b * c) / n
        r1, r2, c1, c2 = a + b, c + d, a + c, b + d
        e_a = (r1 * c1) / n
        v_a = (r1 * r2 * c1 * c2) / (n * n * (n - 1))
        cmh_stat_num += (a - e_a)
        cmh_stat_var += v_a

    cmh_or = cmh_num / max(1e-12, cmh_den)
    cmh_chi2 = (abs(cmh_stat_num) - 0.5) ** 2 / max(1e-12, cmh_stat_var) if abs(cmh_stat_num) > 0.5 else 0.0
    cmh_p = float(chi2.sf(cmh_chi2, df=1))

    # Burial (quartiles) x Entropy (quartiles) Stratified Permutation Test (N=5,000)
    burials = np.array([env3d[k]["burial_neighbors_8A"] for k in universe], dtype=float)
    entropies = np.array([msa_stats[k]["entropy"] for k in universe], dtype=float)
    b_bins = np.digitize(burials, np.quantile(burials, [0.25, 0.5, 0.75]))
    e_bins = np.digitize(entropies, np.quantile(entropies, [0.25, 0.5, 0.75]))
    strata = defaultdict(list)
    for idx, k in enumerate(universe):
        strata[(int(b_bins[idx]), int(e_bins[idx]))].append(idx)

    is_iface = np.array([1 if env3d[k]["min_nuclear_A"] <= 4.5 else 0 for k in universe], dtype=int)
    is_ext = np.array([1 if k in ext_keys else 0 for k in universe], dtype=int)
    obs_iface_count = int(np.sum(is_iface[is_ext == 1]))

    n_perm = 5000
    perm_counts = np.zeros(n_perm, dtype=int)
    for idxs in strata.values():
        arr = np.array(idxs, dtype=int)
        k_in = int(np.sum(is_ext[arr]))
        if k_in == 0:
            continue
        iface_vals = is_iface[arr]
        m = len(arr)
        # Hypergeometric exact draw per stratum across n_perm permutations!
        n_good = int(np.sum(iface_vals))
        n_bad = m - n_good
        perm_counts += RNG.hypergeometric(n_good, n_bad, k_in, size=n_perm)

    strat_perm_mean = float(np.mean(perm_counts))
    strat_perm_p = float((np.sum(perm_counts >= obs_iface_count) + 1) / (n_perm + 1))

    # Multivariate Logistic Regression (exact Newton-IRLS Wald test with Gene FE + Burial + Entropy)
    # Using unpenalized IRLS so beta, SE, Wald z, and p-value are exact!
    gene_list = sorted({k[0] for k in universe})
    X_cols = [
        np.ones(len(universe)),
        is_iface.astype(float),
        (burials - burials.mean()) / (burials.std() + 1e-9),
        (entropies - entropies.mean()) / (entropies.std() + 1e-9),
    ]
    for g in gene_list[1:]:
        X_cols.append(np.array([1.0 if k[0] == g else 0.0 for k in universe]))
    X = np.column_stack(X_cols)
    y = is_ext.astype(float)

    beta = np.zeros(X.shape[1])
    for _ in range(30):
        eta = np.clip(X @ beta, -20.0, 20.0)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = np.maximum(1e-9, p * (1.0 - p))
        H = X.T @ (X * w[:, None]) + 1e-6 * np.eye(X.shape[1])
        grad = X.T @ (y - p)
        step = np.linalg.solve(H, grad)
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            break
    eta = np.clip(X @ beta, -20.0, 20.0)
    p = 1.0 / (1.0 + np.exp(-eta))
    w = np.maximum(1e-9, p * (1.0 - p))
    cov = np.linalg.inv(X.T @ (X * w[:, None]) + 1e-6 * np.eye(X.shape[1]))
    beta_iface = float(beta[1])
    se_iface = float(math.sqrt(cov[1, 1]))
    wald_z = beta_iface / se_iface
    wald_p = float(2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(wald_z) / math.sqrt(2.0)))))
    adj_or = float(math.exp(beta_iface))
    adj_or_ci95 = [round(math.exp(beta_iface - 1.96 * se_iface), 3), round(math.exp(beta_iface + 1.96 * se_iface), 3)]

    # Leave-One-Axis-Out Jackknife
    all_axes = sorted({ax for axes in ext_axes_map.values() for ax in axes})
    loao = {}
    for drop_ax in all_axes:
        sub_ext = {k for k in ext_keys if len(ext_axes_map.get(k, set()) - {drop_ax}) > 0}
        loao[f"drop_{drop_ax}"] = calc_2x2(sub_ext, universe)

    # Leave-One-Gene-Out Jackknife
    logo_ors = {}
    for drop_g in gene_list:
        sub_u = [k for k in universe if k[0] != drop_g]
        sub_e = {k for k in ext_keys if k[0] != drop_g}
        logo_ors[drop_g] = calc_2x2(sub_e, sub_u)["odds_ratio"]

    return {
        "universe_qc_3d_residues": len(universe),
        "extreme_polarized_2x2": ext_2x2,
        "matched_control_polarized_negative_null_2x2": ctl_2x2,
        "extreme_vs_control_direct_fisher": {
            "odds_ratio": round(float(or_ec), 3),
            "p_two_sided": round(float(p_ec), 5),
        },
        "gene_stratified_cmh": {
            "cmh_odds_ratio": round(float(cmh_or), 3),
            "cmh_chi2": round(float(cmh_chi2), 3),
            "cmh_p_value": round(float(cmh_p), 5),
        },
        "burial_entropy_stratified_permutation_5000": {
            "observed_interface_hits": obs_iface_count,
            "expected_stratified_hits": round(strat_perm_mean, 2),
            "enrichment_fold": round(obs_iface_count / max(1e-9, strat_perm_mean), 3),
            "empirical_p_value": round(strat_perm_p, 5),
        },
        "multivariate_logistic_regression_adjusted": {
            "covariates": ["is_mitonuclear_interface", "burial_neighbors_8A", "shannon_entropy", "12_gene_fixed_effects"],
            "adjusted_OR_interface": round(adj_or, 3),
            "adjusted_OR_95CI": adj_or_ci95,
            "wald_z": round(wald_z, 3),
            "wald_p_value": round(wald_p, 5),
            "beta_burial_z": round(float(beta[2]), 3),
            "beta_entropy_z": round(float(beta[3]), 3),
        },
        "leave_one_axis_out_jackknife": loao,
        "leave_one_gene_out_OR_range": {
            "min_OR": min(logo_ors.values()),
            "max_OR": max(logo_ors.values()),
            "per_dropped_gene": logo_ors,
        },
    }


# =============================================================================
# PROOF 2: Full-Proteome Human ClinVar Pathogenic vs. Benign/gnomAD Calibration
# =============================================================================
def run_proof2_proteome_pathogenicity_benchmark(
    env3d: dict[tuple[str, int], dict],
    msa_stats: dict[tuple[str, int], dict],
    gene_info: dict[str, dict],
    s3_data: dict,
) -> tuple[dict, dict]:
    """Train zero-leakage comparative+3D model on human ClinVar Pathogenic vs Benign/gnomAD,
    verify 5-fold CV AUROC and strict out-of-sample Tier-1 -> Frozen Tier-2 transfer AUROC."""
    rcrs_seq = json.loads(s3.RCRS_NUCCORE_JSON.read_text(encoding="utf-8"))[0]["sequence"]

    # Map every rCRS 1-based coordinate -> (gene, aa_pos_1, cds_idx, strand)
    rcrs_to_cds = {}
    for g, gi in gene_info.items():
        st, en, strand, cds = gi["rcrs_start_nt"], gi["rcrs_end_nt"], gi["strand"], gi["cds_nt_seq"]
        for cidx in range(len(cds)):
            rpos = (st + cidx) if strand == "+" else (en - cidx)
            rcrs_to_cds[rpos] = (g, cidx // 3 + 1, cidx, strand, cds)

    def m_to_missense(rpos: int, ref_h: str, alt_h: str) -> tuple[str, int, str, str] | None:
        if rpos not in rcrs_to_cds:
            return None
        g, aa_pos, cidx, strand, cds = rcrs_to_cds[rpos]
        if rcrs_seq[rpos - 1] != ref_h:
            return None
        aa_idx0 = aa_pos - 1
        cod_pos = cidx % 3
        if aa_idx0 * 3 + 3 > len(cds):
            return None
        orig_cod = cds[aa_idx0 * 3 : aa_idx0 * 3 + 3]
        alt_cod_base = alt_h if strand == "+" else s3.COMP[alt_h]
        new_cod = orig_cod[:cod_pos] + alt_cod_base + orig_cod[cod_pos + 1 :]
        ref_aa = s3.translate_codon_at(cds, aa_idx0, orig_cod)
        alt_aa = s3.translate_codon_at(cds, aa_idx0, new_cod)
        if ref_aa == alt_aa or alt_aa == "*" or ref_aa == "*":
            return None
        return (g, aa_pos, ref_aa, alt_aa)

    # 1. Collect ClinVar Pathogenic / Likely Pathogenic missense variants
    path_raw = json.loads(s3.CLINVAR_PATH_JSON.read_text(encoding="utf-8"))
    path_set: set[tuple[str, int, str, str]] = set()
    for rec in path_raw:
        sig = (rec.get("clinical_significance") or "").lower()
        if "pathogenic" not in sig or "benign" in sig or "conflicting" in sig:
            continue
        title = rec.get("title") or ""
        m = s3.parse_rcrs_m_notation(title)
        if m:
            hit = m_to_missense(*m)
            if hit:
                path_set.add(hit)
                continue
        # Fallback: protein HGVS + MT-* gene
        p_hit = s3.parse_hgvsp_to_aa_change(title)
        genes = [s3.MT_HGNC_TO_GENE.get(g.get("symbol", "")) for g in (rec.get("genes") or [])]
        genes = [g for g in genes if g]
        if p_hit and len(genes) == 1:
            ref_a, pos1, alt_a = p_hit
            if ref_a != alt_a and alt_a != "*":
                path_set.add((genes[0], pos1, ref_a, alt_a))

    # 2. Collect ClinVar Benign / LB + Common Homoplasmic gnomAD v3 (af_hom >= 0.001) missense variants
    ben_raw = json.loads(s3.CLINVAR_BEN_JSON.read_text(encoding="utf-8"))
    ben_set: set[tuple[str, int, str, str]] = set()
    for rec in ben_raw:
        sig = (rec.get("clinical_significance") or "").lower()
        if "benign" not in sig or "pathogenic" in sig or "conflicting" in sig:
            continue
        title = rec.get("title") or ""
        m = s3.parse_rcrs_m_notation(title)
        if m:
            hit = m_to_missense(*m)
            if hit:
                ben_set.add(hit)

    gnomad_raw = json.loads(s3.GNOMAD_MT_JSON.read_text(encoding="utf-8"))
    for v in gnomad_raw:
        an = v.get("an") or 0
        ac_hom = v.get("ac_hom") or 0
        if an <= 0 or (ac_hom / an) < 0.001:
            continue
        vid = v.get("variant_id", "")
        parts = vid.split("-")
        if len(parts) == 4 and len(parts[2]) == 1 and len(parts[3]) == 1:
            hit = m_to_missense(int(parts[1]), parts[2], parts[3])
            if hit:
                ben_set.add(hit)

    # Remove any ambiguous overlap from benign set
    ben_set -= path_set

    def featurize(item: tuple[str, int, str, str]) -> list[float] | None:
        g, pos1, ref_a, alt_a = item
        if (g, pos1) not in env3d or (g, pos1) not in msa_stats:
            return None
        e3 = env3d[(g, pos1)]
        ms = msa_stats[(g, pos1)]
        gdist = s1.grantham_distance(ref_a, alt_a)
        d_chg = abs(s1.CHARGE_PH74.get(alt_a, 0.0) - s1.CHARGE_PH74.get(ref_a, 0.0))
        v_ref = s1.GRANTHAM_PROPS[ref_a][2] if ref_a in s1.GRANTHAM_PROPS else 120.0
        v_alt = s1.GRANTHAM_PROPS[alt_a][2] if alt_a in s1.GRANTHAM_PROPS else 120.0
        d_vol = abs(v_alt - v_ref)
        return [
            ms["consensus_frac"],
            ms["entropy"],
            min(50.0, e3["min_catalytic_A"]),
            min(30.0, e3["min_nuclear_A"]),
            float(e3["burial_neighbors_8A"]),
            gdist,
            d_chg,
            d_vol,
        ]

    feat_names = [
        "consensus_frac", "shannon_entropy", "min_catalytic_A",
        "min_nuclear_A", "burial_neighbors_8A", "grantham", "abs_d_charge", "abs_d_vol"
    ]

    X_rows, y_rows, tiers, items_kept = [], [], [], []
    for item in sorted(path_set):
        fv = featurize(item)
        if fv is not None:
            X_rows.append(fv)
            y_rows.append(1)
            tiers.append(gene_info[item[0]]["tier"])
            items_kept.append(item)
    for item in sorted(ben_set):
        fv = featurize(item)
        if fv is not None:
            X_rows.append(fv)
            y_rows.append(0)
            tiers.append(gene_info[item[0]]["tier"])
            items_kept.append(item)

    X_all = np.array(X_rows, dtype=float)
    y_all = np.array(y_rows, dtype=int)
    tiers_arr = np.array(tiers)

    # Univariate Mann-Whitney separation
    univar = {}
    for j, fn in enumerate(feat_names):
        p_vals = X_all[y_all == 1, j]
        b_vals = X_all[y_all == 0, j]
        _, p_u = mannwhitneyu(p_vals, b_vals, alternative="two-sided")
        univar[fn] = {
            "pathogenic_median": round(float(np.median(p_vals)), 3),
            "benign_median": round(float(np.median(b_vals)), 3),
            "mann_whitney_p": float(f"{p_u:.3e}"),
        }

    # 5-Fold Stratified Cross-Validation AUROC
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    oof_probs = np.zeros(len(y_all))
    fold_aucs = []
    for tr, va in skf.split(X_all, y_all):
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, C=1.0))
        clf.fit(X_all[tr], y_all[tr])
        pr = clf.predict_proba(X_all[va])[:, 1]
        oof_probs[va] = pr
        fold_aucs.append(float(roc_auc_score(y_all[va], pr)))
    cv_auc = float(roc_auc_score(y_all, oof_probs))
    fpr_cv, tpr_cv, _ = roc_curve(y_all, oof_probs)

    # Strict Out-of-Sample Generalization: Train ONLY on Tier-1 (6 genes) -> Test blindly on Frozen Tier-2 (7 genes)!
    m_t1 = tiers_arr == 1
    m_t2 = tiers_arr == 2
    clf_t1 = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, C=1.0))
    clf_t1.fit(X_all[m_t1], y_all[m_t1])
    t2_probs = clf_t1.predict_proba(X_all[m_t2])[:, 1]
    t2_oos_auc = float(roc_auc_score(y_all[m_t2], t2_probs))
    fpr_t2, tpr_t2, _ = roc_curve(y_all[m_t2], t2_probs)

    # Reverse Out-of-Sample: Train ONLY on Tier-2 -> Test on Tier-1
    clf_t2 = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, C=1.0))
    clf_t2.fit(X_all[m_t2], y_all[m_t2])
    t1_oos_auc = float(roc_auc_score(y_all[m_t1], clf_t2.predict_proba(X_all[m_t1])[:, 1]))

    # Fit full model to score our 131 wild extremophile candidates & 5 finalists
    clf_full = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, C=1.0))
    clf_full.fit(X_all, y_all)
    wild_scores = {}
    for c in s3_data["all_evaluated_candidates"]:
        fv = featurize((c["gene"], c["human_rcrs_aa_pos"], c["human_ref_aa"], c["alt_aa"]))
        if fv is not None:
            prob = float(clf_full.predict_proba([fv])[0, 1])
            wild_scores[f"{c['gene']}_{c['human_ref_aa']}{c['human_rcrs_aa_pos']}{c['alt_aa']}"] = round(prob, 4)

    finalist_keys = ["COX3_W116R", "ND4L_A11T", "COX3_I165T", "ATP6_I106V", "CYTB_M53T"]
    finalist_probs = {k: wild_scores.get(k) for k in finalist_keys}

    # Score pathogenic controls too for contrast
    ctrl_probs = {}
    for g, pos1, ra, aa in [("ATP6", 156, "L", "R"), ("ND1", 52, "A", "T"), ("ND4", 340, "R", "H"), ("ND6", 64, "M", "V"), ("ND1", 30, "Y", "H")]:
        fv = featurize((g, pos1, ra, aa))
        if fv is not None:
            ctrl_probs[f"{g}_{ra}{pos1}{aa}"] = round(float(clf_full.predict_proba([fv])[0, 1]), 4)

    summary = {
        "n_clinvar_pathogenic_missense_3d": int(np.sum(y_all == 1)),
        "n_benign_or_gnomad_missense_3d": int(np.sum(y_all == 0)),
        "n_tier1_variants": int(np.sum(m_t1)),
        "n_tier2_frozen_variants": int(np.sum(m_t2)),
        "five_fold_cv_auroc": round(cv_auc, 4),
        "five_fold_aucs": [round(x, 4) for x in fold_aucs],
        "strict_oos_tier1_to_tier2_auroc": round(t2_oos_auc, 4),
        "strict_oos_tier2_to_tier1_auroc": round(t1_oos_auc, 4),
        "univariate_feature_separation": univar,
        "wild_131_median_pathogenicity_prob": round(float(np.median(list(wild_scores.values()))), 4),
        "pathogenic_median_prob": round(float(np.median(oof_probs[y_all == 1])), 4),
        "benign_median_prob": round(float(np.median(oof_probs[y_all == 0])), 4),
        "finalist_pathogenicity_probs": finalist_probs,
        "answer_key_control_probs": ctrl_probs,
    }
    plot_data = {
        "fpr_cv": fpr_cv, "tpr_cv": tpr_cv, "cv_auc": cv_auc,
        "fpr_t2": fpr_t2, "tpr_t2": tpr_t2, "t2_auc": t2_oos_auc,
    }
    return summary, plot_data


# =============================================================================
# PROOF 3: 3D Atomic Contact Network & Full Side-Chain Rotamer Clash Scan
# =============================================================================
def _rot_matrix(axis: np.ndarray, theta_rad: float) -> np.ndarray:
    axis = axis / (np.linalg.norm(axis) + 1e-12)
    a = math.cos(theta_rad / 2.0)
    b, c, d = -axis * math.sin(theta_rad / 2.0)
    return np.array([
        [a*a + b*b - c*c - d*d, 2*(b*c + a*d), 2*(b*d - a*c)],
        [2*(b*c - a*d), a*a + c*c - b*b - d*d, 2*(c*d + a*b)],
        [2*(b*d + a*c), 2*(c*d - a*b), a*a + d*d - b*b - c*c],
    ])


def _place_next_atom(a1: np.ndarray, a2: np.ndarray, a3: np.ndarray, bond_len: float, bond_angle_deg: float, dihedral_deg: float) -> np.ndarray:
    """NeRF (Natural Extension Reference Frame) exact internal-to-Cartesian atom placement."""
    bc = a3 - a2
    bc = bc / (np.linalg.norm(bc) + 1e-12)
    ab = a2 - a1
    n = np.cross(ab, bc)
    n_norm = np.linalg.norm(n)
    if n_norm < 1e-8:
        n = np.cross(np.array([1.0, 0.0, 0.0]), bc)
        n_norm = np.linalg.norm(n)
    n = n / n_norm
    m = np.cross(n, bc)
    theta = math.radians(180.0 - bond_angle_deg)
    phi = math.radians(dihedral_deg)
    d_local = np.array([
        bond_len * math.cos(theta),
        bond_len * math.sin(theta) * math.cos(phi),
        bond_len * math.sin(theta) * math.sin(phi),
    ])
    return a3 + d_local[0] * bc + d_local[1] * m + d_local[2] * n


def generate_sidechain_rotamers(N: np.ndarray, CA: np.ndarray, C: np.ndarray, CB: np.ndarray, alt_aa: str) -> list[tuple[str, np.ndarray]]:
    """Build realistic 3D heavy-atom side-chain rotamer library using exact NeRF bond geometry."""
    rotamers: list[tuple[str, np.ndarray]] = []
    if alt_aa == "V":
        # Val: CG1, CG2 on chi1 in {-60, 60, 180} (30-deg fine grid around canonical)
        for chi1 in (-60, 60, 180, -75, -45, 45, 75, 165, -165):
            cg1 = _place_next_atom(N, CA, CB, 1.52, 110.5, chi1)
            cg2 = _place_next_atom(N, CA, CB, 1.52, 110.5, chi1 + 120.0)
            rotamers.append((f"chi1={chi1}", np.array([CB, cg1, cg2])))
    elif alt_aa == "T":
        # Thr: OG1, CG2 on chi1
        for chi1 in (-60, 60, 180, -75, -45, 45, 75, 165, -165):
            og1 = _place_next_atom(N, CA, CB, 1.43, 109.5, chi1)
            cg2 = _place_next_atom(N, CA, CB, 1.52, 110.5, chi1 + 120.0)
            rotamers.append((f"chi1={chi1}", np.array([CB, og1, cg2])))
    elif alt_aa == "R":
        # Arg: CG, CD, NE, CZ, NH1, NH2 across 81 canonical Dunbrack rotamers (chi1 x chi2 x chi3 x chi4)
        for chi1 in (-65, -175, 65):
            cg = _place_next_atom(N, CA, CB, 1.52, 114.1, chi1)
            for chi2 in (180, -65, 65):
                cd = _place_next_atom(CA, CB, cg, 1.52, 111.5, chi2)
                for chi3 in (180, -65, 65):
                    ne = _place_next_atom(CB, cg, cd, 1.46, 112.0, chi3)
                    for chi4 in (180, -85, 85):
                        cz = _place_next_atom(cg, cd, ne, 1.33, 124.0, chi4)
                        nh1 = _place_next_atom(cd, ne, cz, 1.33, 120.0, 0.0)
                        nh2 = _place_next_atom(cd, ne, cz, 1.33, 120.0, 180.0)
                        rotamers.append((f"chi=({chi1},{chi2},{chi3},{chi4})", np.array([CB, cg, cd, ne, cz, nh1, nh2])))
    else:
        rotamers.append(("CB_only", np.array([CB])))
    return rotamers


def run_proof3_3d_rotamer_clash_audit(gene_info: dict[str, dict]) -> list[dict]:
    """Run exact experimental cryo-EM 5.0 A contact audit + 3D rotamer steric clash scan
    for all 5 Stage-3 finalists + 2 human pathogenic controls."""
    targets = [
        ("COX3", 116, "W", "R", "Complex_IV", "Stage3_Finalist_#1_BarHeadedGoose"),
        ("ND4L", 11, "A", "T", "Complex_I", "Stage3_Finalist_#2_BirdMammalHypoxia"),
        ("COX3", 165, "I", "T", "Complex_IV", "Stage3_Finalist_#3_SnowLeopardNMR"),
        ("ATP6", 106, "I", "V", "Complex_V", "Stage3_Finalist_#4_3WayHypoxia"),
        ("CYTB", 53, "M", "T", "Complex_III", "Stage3_Finalist_#5_UQCRFS1_Interface"),
        ("ATP6", 156, "L", "R", "Complex_V", "Pathogenic_Control_Leigh_NARP"),
        ("ND1", 52, "A", "T", "Complex_I", "Pathogenic_Control_LHON"),
    ]

    meta = json.loads(s2.META_JSON.read_text(encoding="utf-8"))
    uniprot_to_mt = {gi["uniprot_id"]: g for g, gi in gene_info.items()}
    chain_meta_all = s2.build_chain_metadata(meta, uniprot_to_mt)

    # Cache parsed structures per complex
    parser = MMCIFParser(QUIET=True)
    structures = {}
    for cname, (pid, cpath) in s2.COMPLEX_PDBS.items():
        with gzip.open(cpath, "rt", encoding="utf-8") as fh:
            structures[cname] = (pid, next(parser.get_structure(pid, fh).get_models()))

    results = []
    for gene, pos1, ref_aa, alt_aa, cname, role in targets:
        pid, model = structures[cname]
        cmeta = chain_meta_all[pid]
        # Find the chain for this mtDNA gene
        target_chain_id = None
        for cid, info in cmeta.items():
            if info["gene_symbol"] == gene and cid in model:
                target_chain_id = cid
                break
        assert target_chain_id is not None, f"Chain not found for {gene} in {pid}"

        # Find exact residue matching rCRS pos1 (note: Stage 2 maps exact sequence numbering; verify resseq)
        target_res = None
        chain = model[target_chain_id]
        for res in chain:
            if res.id[0] == " " and res.id[1] == pos1:
                target_res = res
                break
        if target_res is None:
            # fallback by 1-based residue index if offset
            std_res = [r for r in chain if r.id[0] == " "]
            target_res = std_res[pos1 - 1]

        N = target_res["N"].coord.astype(float)
        CA = target_res["CA"].coord.astype(float)
        C = target_res["C"].coord.astype(float)
        if "CB" in target_res:
            CB = target_res["CB"].coord.astype(float)
        else:
            # Reconstruct virtual CB for Gly/Ala if missing
            CB = _place_next_atom(C, N, CA, 1.53, 110.5, -122.0)

        # Collect ALL other heavy atoms in the holocomplex (excluding target residue itself & immediate peptide-bonded N/C)
        env_coords = []
        env_labels = []
        target_sc_atoms = [a.coord.astype(float) for a in target_res.get_atoms() if a.element != "H"]
        target_center = np.mean(target_sc_atoms, axis=0)

        for ch in model:
            ch_id = ch.id
            ch_gene = cmeta.get(ch_id, {}).get("gene_symbol", ch_id)
            is_mt = cmeta.get(ch_id, {}).get("is_mtdna", False)
            for res in ch:
                if ch_id == target_chain_id and res.id == target_res.id:
                    continue
                resname = res.get_resname().strip()
                if resname == "HOH":
                    continue
                is_het = res.id[0] != " "
                for atom in res.get_atoms():
                    if atom.element == "H":
                        continue
                    # Skip immediate covalent backbone atoms of i-1 (C, O) and i+1 (N)
                    if ch_id == target_chain_id and abs(res.id[1] - target_res.id[1]) == 1 and atom.get_name() in {"N", "C", "O", "CA"}:
                        continue
                    coord = atom.coord.astype(float)
                    if np.linalg.norm(coord - target_center) <= 14.0:
                        env_coords.append(coord)
                        env_labels.append({
                            "chain": ch_id,
                            "subunit": ch_gene,
                            "is_mtdna": is_mt,
                            "is_het": is_het,
                            "resname": resname,
                            "resseq": res.id[1],
                            "atom": atom.get_name(),
                        })

        env_arr = np.array(env_coords, dtype=float)
        tree = cKDTree(env_arr)

        # Exact 5.0 A WT atomic neighborhood (grouped by residue/ligand)
        wt_contacts = {}
        for a_coord in target_sc_atoms:
            idxs = tree.query_ball_point(a_coord, r=5.0)
            for idx in idxs:
                d = float(np.linalg.norm(a_coord - env_arr[idx]))
                lb = env_labels[idx]
                ckey = f"{lb['subunit']}({lb['chain']}):{lb['resname']}{lb['resseq']}"
                if ckey not in wt_contacts or d < wt_contacts[ckey]["min_dist_A"]:
                    wt_contacts[ckey] = {
                        "contact": ckey,
                        "type": "Lipid/Ligand" if lb["is_het"] else ("mtDNA" if lb["is_mtdna"] else "Nuclear_Partner"),
                        "min_dist_A": round(d, 2),
                    }
        wt_contact_list = sorted(wt_contacts.values(), key=lambda x: x["min_dist_A"])

        # WT side-chain minimum distance to environment
        wt_sc_only = np.array([a.coord.astype(float) for a in target_res.get_atoms() if a.element != "H" and a.get_name() not in {"N", "CA", "C", "O"}], dtype=float)
        wt_dmin = float(np.min(tree.query(wt_sc_only)[0])) if len(wt_sc_only) > 0 else 4.0

        # Run 3D rotamer steric clash scan for mutant alt_aa
        rotamers = generate_sidechain_rotamers(N, CA, C, CB, alt_aa)
        best_rot = None
        clash_free_count = 0
        for rname, rcoords in rotamers:
            # Evaluate distances of new sidechain atoms (beyond CB) against all environment heavy atoms
            probe = rcoords[1:] if len(rcoords) > 1 else rcoords
            dists, _ = tree.query(probe, k=min(8, len(env_arr)))
            if dists.ndim == 1:
                dists = dists[:, None]
            d_min = float(np.min(dists))
            severe_clashes = int(np.sum(dists < 2.20))  # heavy-atom steric overlap
            mild_contacts = int(np.sum((dists >= 2.20) & (dists < 2.60)))
            # Soft repulsive Lennard-Jones proxy (sum of ((3.0/d)^6 - 1) for d < 3.0)
            close = dists[dists < 3.0]
            lj_rep = float(np.sum((3.0 / np.maximum(1.2, close)) ** 6 - 1.0)) if len(close) > 0 else 0.0
            if severe_clashes == 0:
                clash_free_count += 1
            cand = {
                "rotamer": rname,
                "min_heavy_atom_dist_A": round(d_min, 2),
                "severe_clashes_lt_2p2A": severe_clashes,
                "mild_contacts_2p2_2p6A": mild_contacts,
                "soft_lj_repulsion": round(lj_rep, 2),
            }
            if best_rot is None or (cand["severe_clashes_lt_2p2A"], cand["soft_lj_repulsion"], -cand["min_heavy_atom_dist_A"]) < (
                best_rot["severe_clashes_lt_2p2A"], best_rot["soft_lj_repulsion"], -best_rot["min_heavy_atom_dist_A"]
            ):
                best_rot = cand

        # Steric volume change & cavity verdict
        d_vol = round(s1.GRANTHAM_PROPS[alt_aa][2] - s1.GRANTHAM_PROPS[ref_aa][2], 1)
        lipid_contacts = [c["contact"] for c in wt_contact_list if c["type"] == "Lipid/Ligand"]
        nuc_contacts = [c["contact"] for c in wt_contact_list if c["type"] == "Nuclear_Partner"]

        results.append({
            "variant": f"{gene}_{ref_aa}{pos1}{alt_aa}",
            "role": role,
            "pdb_id": pid,
            "chain_id": target_chain_id,
            "delta_volume_A3": d_vol,
            "wt_min_env_dist_A": round(wt_dmin, 2),
            "rotamers_scanned": len(rotamers),
            "clash_free_rotamers": clash_free_count,
            "best_rotamer": best_rot,
            "steric_verdict": "CLASH_FREE_CAVITY_FIT" if best_rot["severe_clashes_lt_2p2A"] == 0 else "STERIC_CLASH_CONSTRAINED",
            "n_contacts_5A": len(wt_contact_list),
            "lipid_ligand_contacts_5A": lipid_contacts[:6],
            "nuclear_partner_contacts_5A": nuc_contacts[:6],
            "top_contacts_5A": wt_contact_list[:8],
        })

    return results


# =============================================================================
# PROOF 4: Base-Editor Window Sensitivity Stress Test (+/-2 to +/-6 nt)
# =============================================================================
def run_proof4_editor_window_stress_test(gene_info: dict[str, dict], s3_data: dict) -> dict:
    """Stress-test all 42 Motif_Grafting_BaseEdit candidates across +/-2, +/-3, +/-4, +/-5, +/-6 nt."""
    rcrs_seq = json.loads(s3.RCRS_NUCCORE_JSON.read_text(encoding="utf-8"))[0]["sequence"]
    windows = [2, 3, 4, 5, 6]

    editable = [
        c for c in s3_data["all_evaluated_candidates"]
        if c["structural_gate"]["transfer_strategy"] == "Motif_Grafting_BaseEdit"
    ]

    def count_bystanders_for_window(gene: str, aa_pos_1: int, edit_rec: dict, half_w: int) -> int:
        gi = gene_info[gene]
        start_nt, end_nt, strand, cds_nt = gi["rcrs_start_nt"], gi["rcrs_end_nt"], gi["strand"], gi["cds_nt_seq"]
        t_pos = edit_rec["rcrs_pos"]
        e_strand = edit_rec["editor_strand"]
        etype = edit_rec["type"]
        if e_strand is None:
            return 99
        is_cbe = etype.startswith("DdCBE")
        if is_cbe:
            scan_from, scan_to = ("C", "T") if e_strand == "H" else ("G", "A")
        else:
            scan_from, scan_to = ("A", "G") if e_strand == "H" else ("T", "C")

        cnt = 0
        for delta in range(-half_w, half_w + 1):
            if delta == 0:
                continue
            b_pos = t_pos + delta
            if not (start_nt <= b_pos <= end_nt):
                continue
            if rcrs_seq[b_pos - 1] != scan_from:
                continue
            if etype == "DdCBE_TC":
                if e_strand == "H":
                    p5 = rcrs_seq[b_pos - 2] if b_pos >= 2 else "N"
                else:
                    p5 = s3.COMP.get(rcrs_seq[b_pos] if b_pos < len(rcrs_seq) else "N", "N")
                if p5 != "T":
                    continue
            b_cidx = (b_pos - start_nt) if strand == "+" else (end_nt - b_pos)
            b_aa_idx = b_cidx // 3
            b_cod_pos = b_cidx % 3
            if b_aa_idx * 3 + 3 > len(cds_nt):
                continue
            orig_cod = cds_nt[b_aa_idx * 3 : b_aa_idx * 3 + 3]
            mut_b_cod = scan_to if strand == "+" else s3.COMP[scan_to]
            new_cod = orig_cod[:b_cod_pos] + mut_b_cod + orig_cod[b_cod_pos + 1 :]
            if s3.translate_codon_at(cds_nt, b_aa_idx, new_cod) != s3.translate_codon_at(cds_nt, b_aa_idx, orig_cod):
                cnt += 1
        return cnt

    zero_by_window = {f"pm{w}nt": 0 for w in windows}
    per_candidate = []
    for c in editable:
        ed = c["editor"]
        recs = ed.get("edits") or []
        if len(recs) != 1:
            continue
        w_map = {}
        for w in windows:
            b_cnt = count_bystanders_for_window(c["gene"], c["human_rcrs_aa_pos"], recs[0], w)
            w_map[f"pm{w}nt"] = b_cnt
            if b_cnt == 0:
                zero_by_window[f"pm{w}nt"] += 1
        per_candidate.append({
            "variant": f"{c['gene']}_{c['human_ref_aa']}{c['human_rcrs_aa_pos']}{c['alt_aa']}",
            "tier": c["tier"],
            "editor_class": ed["base_editor_class"],
            "rcrs_edit": recs[0]["rcrs_change"],
            "bystanders_by_window": w_map,
            "ultra_clean_pm6": w_map["pm6nt"] == 0,
        })

    ultra_clean = [x for x in per_candidate if x["ultra_clean_pm6"]]
    finalists_win = [
        x for x in per_candidate
        if x["variant"] in {"COX3_W116R", "ND4L_A11T", "COX3_I165T", "ATP6_I106V", "CYTB_M53T"}
    ]

    return {
        "n_single_step_base_editable_tested": len(per_candidate),
        "zero_bystander_count_by_window": zero_by_window,
        "ultra_clean_zero_bystander_at_pm6nt": ultra_clean,
        "top5_finalists_window_profile": finalists_win,
    }


# =============================================================================
# PROOF 5: Stage-4 Supermodularity, Invariants, Multi-Model Forms, Nulls & Heteroplasmy
# =============================================================================
def simulate_ir_general(
    params_vec: dict[str, np.ndarray],
    driver: dict | None,
    buffer: dict | None,
    t_isch: float = 20.0,
    t_rep: float = 40.0,
    dt: float = 0.25,
    ros_law: str = "linear_threshold",
    heteroplasmy: float = 1.0,
) -> dict[str, np.ndarray]:
    """Extended Stage-4 ODE simulator matching stage4_epistasis_synergy_uq.py identically at baseline,
    while supporting alternative ROS laws, ischemia durations, heteroplasmy dosing (0..1), and
    physical invariant tracking."""
    n = len(params_vec["f_load"])
    scale = params_vec["drv_scale"] * heteroplasmy

    if driver is not None and heteroplasmy > 0.0:
        km_o2 = params_vec["Km_O2_wt"] * (1.0 - scale * (1.0 - driver["km_o2_mult"]))
        atp_o_gain = scale * driver["atp_o_gain"]
        ret_mult = 1.0 + scale * (driver["ret_sensitivity"] - 1.0)
        leak_d = scale * driver["proton_leak_delta"]
    else:
        km_o2 = params_vec["Km_O2_wt"]
        atp_o_gain = np.zeros(n)
        ret_mult = np.ones(n)
        leak_d = np.zeros(n)

    if buffer is not None:
        aox_drain = buffer["aox_drain"]
        ndi1_bypass = buffer["ndi1_bypass"]
        theta_boost = buffer["theta_gsh_boost"]
        ucp_leak = buffer["ucp_leak"]
        mito_boost = buffer["mitophagy_boost"]
        buf_cost = buffer["atp_o_cost"]
    else:
        aox_drain = ndi1_bypass = theta_boost = ucp_leak = mito_boost = buf_cost = 0.0

    theta_eff = params_vec["theta_GSH"] + theta_boost

    def po2_fn(t_now: float) -> float:
        if t_now < t_isch:
            return 0.24
        if t_now < t_isch + 3.0:
            return 0.24 + (14.0 - 0.24) * ((t_now - t_isch) / 3.0)
        return 14.0

    def rhs_and_ep(t_now: float, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        W, Q, Succ, ATP = y[0], y[1], y[2], y[3]
        pO2 = po2_fn(t_now)
        g_iv = pO2 / (np.maximum(0.10, km_o2) + pO2)
        leak = np.clip(0.12 + leak_d + ucp_leak, 0.04, 0.95)
        dpsi = np.clip(0.60 + 0.55 * Q * g_iv - 0.65 * leak - 0.15 * aox_drain * np.maximum(0.0, Q - 0.55), 0.25, 1.35)
        isch_flag = 1.0 if pO2 < 1.0 else 0.0

        dSucc = 0.11 * params_vec["succ_surge"] * isch_flag * (1.45 - Succ) - 0.16 * (1.0 - isch_flag) * g_iv * Succ
        q_in = (0.34 * params_vec["f_load"] + 0.85 * (1.0 - isch_flag) * g_iv * Succ) * (1.0 - Q)
        aox_flux = aox_drain * 0.95 * np.maximum(0.0, Q - 0.48) * (pO2 / (0.60 + pO2))
        q_out = (0.45 * g_iv * np.maximum(0.15, W) + aox_flux) * Q
        dQ = q_in - q_out

        ret_factor = (
            params_vec["k_RET"]
            * ret_mult
            * (1.0 - 0.58 * ndi1_bypass)
            * (1.0 - 0.62 * np.minimum(1.0, ucp_leak / 0.25))
            * np.maximum(0.0, Q - 0.40) ** 1.4
            * g_iv
            * (0.25 + 1.45 * Succ)
            * (dpsi ** 2)
        )
        eR = (0.55 / np.maximum(0.35, W)) * (1.0 + 1.80 * ret_factor)
        ep = eR * (params_vec["f_load"] ** 2)

        if ros_law == "linear_threshold":
            ep_excess = np.maximum(0.0, ep - theta_eff)
        elif ros_law == "hill_cooperative":
            ratio = np.maximum(0.0, ep / np.maximum(1e-6, theta_eff))
            ep_excess = 1.15 * (ratio ** 4) / (1.0 + ratio ** 4) * np.maximum(0.0, ratio - 0.65)
        else:
            raise ValueError(ros_law)

        dW = -params_vec["k_dmg_ir"] * ep_excess * W + (0.008 + 0.042 * mito_boost) * (0.95 - W) * g_iv * (1.0 - np.minimum(1.0, ep_excess))
        coupling = np.maximum(0.05, (1.0 + atp_o_gain) * (1.0 - 0.55 * (leak - 0.12)) * (1.0 - buf_cost * (0.35 + 0.65 * np.maximum(0.0, Q - 0.45))))
        atp_prod = 1.18 * np.maximum(0.05, W) * g_iv * coupling
        atp_use = 0.52 * params_vec["f_load"] * (ATP / (0.18 + ATP))
        dATP = 0.35 * (atp_prod - atp_use)
        return np.stack([dW, dQ, dSucc, dATP], axis=0), ep

    y = np.zeros((4, n), dtype=np.float64)
    y[0, :] = 0.92
    y[1, :] = 0.35
    y[2, :] = 0.10
    y[3, :] = 0.95

    min_atp = np.full(n, 99.0)
    max_ep = np.zeros(n)
    inv_violations = 0

    total_t = t_isch + t_rep
    n_steps = int(round(total_t / dt))
    atp_win_start = min(8.0, 0.4 * t_isch)

    for step in range(n_steps):
        t = step * dt
        k1, ep1 = rhs_and_ep(t, y)
        if atp_win_start <= t <= t_isch:
            min_atp = np.minimum(min_atp, y[3])
        if t >= t_isch:
            max_ep = np.maximum(max_ep, ep1)

        k2, _ = rhs_and_ep(t + 0.5 * dt, np.clip(y + 0.5 * dt * k1, 0.0, 2.0))
        k3, _ = rhs_and_ep(t + 0.5 * dt, np.clip(y + 0.5 * dt * k2, 0.0, 2.0))
        k4, _ = rhs_and_ep(t + dt, np.clip(y + dt * k3, 0.0, 2.0))
        y_next = y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        if np.any(~np.isfinite(y_next)) or np.any(y_next[0] < -1e-6) or np.any(y_next[0] > 1.0 + 1e-6) or np.any(y_next[1] < -1e-6) or np.any(y_next[1] > 1.0 + 1e-6) or np.any(y_next[2] < -1e-6) or np.any(y_next[3] < -1e-6):
            inv_violations += 1
        y = np.clip(y_next, 0.0, 2.0)

    surv_isch = 1.0 / (1.0 + np.exp(-14.0 * (min_atp - params_vec["theta_ATP"])))
    if ros_law == "linear_threshold":
        surv_rep = 1.0 / (1.0 + np.exp(4.0 * (max_ep - theta_eff)))
    else:
        ratio_pk = np.maximum(0.0, max_ep / np.maximum(1e-6, theta_eff))
        surv_rep = 1.0 / (1.0 + ratio_pk ** 5)
    viab = surv_isch * surv_rep * (y[0] / 0.92)
    return {
        "viability": viab,
        "W_end": y[0],
        "min_atp": min_atp,
        "max_ep": max_ep,
        "invariant_violations": inv_violations,
    }


def run_proof5_epistasis_hardening() -> tuple[dict, dict]:
    """Run formal supermodularity proof, invariant check, 3 alt model forms, adversarial controls,
    and heteroplasmy dose-response sweep."""
    # 1. Exact Supermodularity Grid Proof on nominal parameters for COX3_W116R + Yeast_NDI1
    p_nom = {k: np.array([v]) for k, v in s4.NOMINAL_P.items()}
    drv = s4.DRIVERS["COX3_W116R"]
    buf_ndi = s4.BUFFERS["Yeast_NDI1"]
    grid_vals = np.linspace(0.0, 1.0, 11)
    supermod_cells = 0
    dual_cliff_cells = 0
    dual_cliff_supermod = 0
    max_delta2 = -1.0

    v_00 = float(simulate_ir_general(p_nom, None, None)["viability"][0])
    for d_frac in grid_vals[1:]:
        r_d0 = simulate_ir_general(p_nom, drv, None, heteroplasmy=float(d_frac))
        v_d0 = float(r_d0["viability"][0])
        atp_d0 = float(r_d0["min_atp"][0])
        ep_d0 = float(r_d0["max_ep"][0])
        for b_frac in grid_vals[1:]:
            buf_scaled = {
                "aox_drain": buf_ndi["aox_drain"] * b_frac,
                "ndi1_bypass": buf_ndi["ndi1_bypass"] * b_frac,
                "theta_gsh_boost": buf_ndi["theta_gsh_boost"] * b_frac,
                "ucp_leak": buf_ndi["ucp_leak"] * b_frac,
                "mitophagy_boost": buf_ndi["mitophagy_boost"] * b_frac,
                "atp_o_cost": buf_ndi["atp_o_cost"] * b_frac,
            }
            v_0b = float(simulate_ir_general(p_nom, None, buf_scaled)["viability"][0])
            r_db = simulate_ir_general(p_nom, drv, buf_scaled, heteroplasmy=float(d_frac))
            v_db = float(r_db["viability"][0])
            ep_db = float(r_db["max_ep"][0])
            d2 = v_db - v_d0 - v_0b + v_00
            if d2 > 0:
                supermod_cells += 1
            if d2 > max_delta2:
                max_delta2 = d2
            if atp_d0 > s4.NOMINAL_P["theta_ATP"] and ep_d0 > s4.NOMINAL_P["theta_GSH"] and ep_db < ep_d0:
                dual_cliff_cells += 1
                if d2 > 0:
                    dual_cliff_supermod += 1

    # 2. Multi-Model-Form & Adversarial Buffer Stress Test (N=3,000 LHS per scenario)
    n_lhs = 3000
    sampler = qmc.LatinHypercube(d=len(s4.PARAM_NAMES), seed=SEED + 5)
    U = sampler.random(n_lhs)
    p_lhs = {name: s4.NOMINAL_P[name] * (0.70 + 0.60 * U[:, j]) for j, name in enumerate(s4.PARAM_NAMES)}

    test_buffers = {
        "Yeast_NDI1": s4.BUFFERS["Yeast_NDI1"],
        "Human_SOD2": s4.BUFFERS["Human_SOD2"],
        "Ciona_AOX": s4.BUFFERS["Ciona_AOX"],
        "Human_UCP2": s4.BUFFERS["Human_UCP2"],
        "ADVERSARIAL_Null_Placebo": {
            "aox_drain": 0.0, "ndi1_bypass": 0.0, "theta_gsh_boost": 0.0,
            "ucp_leak": 0.0, "mitophagy_boost": 0.0, "atp_o_cost": 0.0,
        },
        "ADVERSARIAL_Toxic_Uncoupler_FCCP": {
            "aox_drain": 0.0, "ndi1_bypass": 0.0, "theta_gsh_boost": 0.0,
            "ucp_leak": 0.85, "mitophagy_boost": 0.0, "atp_o_cost": 0.65,
        },
    }

    scenarios = {
        "Baseline_20mIsch_40mRep_LinearROS": dict(t_isch=20.0, t_rep=40.0, ros_law="linear_threshold"),
        "AltModel_M1_Hill_Cooperative_ROS": dict(t_isch=20.0, t_rep=40.0, ros_law="hill_cooperative"),
        "AltProtocol_M2_Prolonged_30mIsch_60mRep": dict(t_isch=30.0, t_rep=60.0, ros_law="linear_threshold"),
        "AltProtocol_M3_Short_10mIsch_40mRep": dict(t_isch=10.0, t_rep=40.0, ros_law="linear_threshold"),
    }

    total_inv_violations = 0
    scenario_results = {}
    for sname, skw in scenarios.items():
        r_wt = simulate_ir_general(p_lhs, None, None, **skw)
        r_drv = simulate_ir_general(p_lhs, drv, None, **skw)
        total_inv_violations += r_wt["invariant_violations"] + r_drv["invariant_violations"]
        dv = r_drv["viability"] - r_wt["viability"]
        s_out = {}
        for bname, bcfg in test_buffers.items():
            r_buf = simulate_ir_general(p_lhs, None, bcfg, **skw)
            r_pair = simulate_ir_general(p_lhs, drv, bcfg, **skw)
            total_inv_violations += r_buf["invariant_violations"] + r_pair["invariant_violations"]
            bv = r_buf["viability"] - r_wt["viability"]
            pv = r_pair["viability"] - r_wt["viability"]
            e_add = pv - (dv + bv)
            denom = np.maximum(0.025, np.maximum(0.0, dv) + np.maximum(0.0, bv))
            syn = pv / denom
            s_out[bname] = {
                "viability_median": round(float(np.median(r_pair["viability"])), 4),
                "synergy_median": round(float(np.median(syn)), 3),
                "epistatic_excess_median": round(float(np.median(e_add)), 4),
                "pct_lhs_positive": round(float(100.0 * np.mean(e_add > 0.0)), 1),
            }
        scenario_results[sname] = s_out

    # 3. Wet-Lab TALED Heteroplasmy Dose-Response Sweep (0% -> 100%, N=2,000 LHS)
    n_het = 2000
    U_h = qmc.LatinHypercube(d=len(s4.PARAM_NAMES), seed=SEED + 99).random(n_het)
    p_het = {name: s4.NOMINAL_P[name] * (0.70 + 0.60 * U_h[:, j]) for j, name in enumerate(s4.PARAM_NAMES)}
    het_grid = [round(x, 2) for x in np.linspace(0.0, 1.0, 21)]
    r_wt_h = simulate_ir_general(p_het, None, None)

    top4_bufs = ["Yeast_NDI1", "Human_SOD2", "Ciona_AOX", "Human_UCP2"]
    r_buf_h = {bk: simulate_ir_general(p_het, None, s4.BUFFERS[bk]) for bk in top4_bufs}

    het_curves = {bk: {"heteroplasmy_pct": [], "viability_med": [], "epistatic_excess_med": [], "pct_pos": []} for bk in top4_bufs}
    het_driver_alone = []
    for h in het_grid:
        r_d = simulate_ir_general(p_het, drv, None, heteroplasmy=h)
        dv = r_d["viability"] - r_wt_h["viability"]
        het_driver_alone.append(round(float(np.median(r_d["viability"])), 4))
        for bk in top4_bufs:
            r_p = simulate_ir_general(p_het, drv, s4.BUFFERS[bk], heteroplasmy=h)
            bv = r_buf_h[bk]["viability"] - r_wt_h["viability"]
            pv = r_p["viability"] - r_wt_h["viability"]
            e_add = pv - (dv + bv)
            het_curves[bk]["heteroplasmy_pct"].append(int(round(h * 100)))
            het_curves[bk]["viability_med"].append(round(float(np.median(r_p["viability"])), 4))
            het_curves[bk]["epistatic_excess_med"].append(round(float(np.median(e_add)), 4))
            het_curves[bk]["pct_pos"].append(round(float(100.0 * np.mean(e_add > 0.0)), 1))

    # Extract exact minimum editing heteroplasmy thresholds (E_add >= 0.05 & Viability >= 0.25)
    min_het_thresholds = {}
    for bk in top4_bufs:
        hc = het_curves[bk]
        h_syn05 = next((p for p, e in zip(hc["heteroplasmy_pct"], hc["epistatic_excess_med"]) if e >= 0.05), None)
        h_viab25 = next((p for p, v in zip(hc["heteroplasmy_pct"], hc["viability_med"]) if v >= 0.25), None)
        min_het_thresholds[bk] = {
            "min_heteroplasmy_pct_for_Eadd_0p05": h_syn05,
            "min_heteroplasmy_pct_for_Viability_0p25": h_viab25,
        }

    summary = {
        "supermodularity_grid_proof": {
            "grid_cells_evaluated": 100,
            "positive_cross_partial_cells": supermod_cells,
            "dual_cliff_theorem_cells": dual_cliff_cells,
            "dual_cliff_theorem_concordance_pct": round(100.0 * dual_cliff_supermod / max(1, dual_cliff_cells), 2),
            "max_cross_partial_delta2V": round(max_delta2, 4),
        },
        "physical_invariant_violations_across_all_scenarios": total_inv_violations,
        "multi_model_form_and_adversarial_controls": scenario_results,
        "heteroplasmy_min_thresholds": min_het_thresholds,
        "heteroplasmy_driver_alone_curve": het_driver_alone,
        "heteroplasmy_curves": het_curves,
    }
    return summary, {"het_grid_pct": [int(round(x * 100)) for x in het_grid], "het_driver_alone": het_driver_alone, "het_curves": het_curves, "scenario_results": scenario_results}


# =============================================================================
# PROOF 6: End-to-End SQLite Referential Integrity & Cross-Stage Parity Audit
# =============================================================================
def run_proof6_db_integrity() -> dict:
    conn = sqlite3.connect(DB_PATH)
    fk_issues = conn.execute("PRAGMA foreign_key_check;").fetchall()
    tables = [
        "environmental_stressors", "donor_species", "human_rcrs_genes",
        "mt_alleles", "fitness_payoff_profiles", "human_safety_annotations",
        "transferability_assessments", "epistatic_synergies",
    ]
    counts = {t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0] for t in tables}

    # Check 1-to-1 parity between mt_alleles, fitness_payoff_profiles, transferability_assessments
    orphan_payoff = conn.execute(
        "SELECT COUNT(*) FROM mt_alleles a LEFT JOIN fitness_payoff_profiles p ON a.allele_id = p.allele_id WHERE p.allele_id IS NULL"
    ).fetchone()[0]
    orphan_trans = conn.execute(
        "SELECT COUNT(*) FROM mt_alleles a LEFT JOIN transferability_assessments t ON a.allele_id = t.allele_id WHERE t.allele_id IS NULL"
    ).fetchone()[0]
    conn.close()
    return {
        "foreign_key_violations": len(fk_issues),
        "table_row_counts": counts,
        "orphan_payoff_rows": orphan_payoff,
        "orphan_transferability_rows": orphan_trans,
        "integrity_pass": (len(fk_issues) == 0 and orphan_payoff == 0 and orphan_trans == 0),
    }


# =============================================================================
# RENDER 4-PANEL HARDENING & PROOFS FIGURE
# =============================================================================
def render_stage5_figure(p2_plot: dict, p3_res: list[dict], p5_plot: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.6))

    # Panel A: Full-Proteome Human Pathogenicity ROC (5-Fold CV + Strict Tier-1 -> Frozen Tier-2 OOS)
    ax = axes[0, 0]
    ax.plot(p2_plot["fpr_cv"], p2_plot["tpr_cv"], color="#0284c7", lw=2.4,
            label=f"5-Fold Proteome CV (AUROC = {p2_plot['cv_auc']:.3f})")
    ax.plot(p2_plot["fpr_t2"], p2_plot["tpr_t2"], color="#059669", lw=2.3, ls="--",
            label=f"Strict OOS: Tier-1 → Frozen Tier-2 (AUROC = {p2_plot['t2_auc']:.3f})")
    ax.plot([0, 1], [0, 1], color="#94a3b8", lw=1.1, ls=":")
    ax.set_xlabel("False Positive Rate (ClinVar Benign / Common gnomAD)")
    ax.set_ylabel("True Positive Rate (ClinVar Pathogenic / LP)")
    ax.set_title("(A) Proof 2: Proteome-Wide ClinVar Calibration & OOS Transfer", fontsize=10.5, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8.5, frameon=True)
    ax.grid(True, alpha=0.25)

    # Panel B: 3D Cryo-EM Rotamer Steric Clash Audit (5 Finalists vs 2 Pathogenic Controls)
    ax = axes[0, 1]
    labels = [r["variant"].replace("_", " ") for r in p3_res]
    dmins = [r["best_rotamer"]["min_heavy_atom_dist_A"] for r in p3_res]
    clashes = [r["best_rotamer"]["severe_clashes_lt_2p2A"] for r in p3_res]
    colors = ["#059669" if (r["role"].startswith("Stage3") and c == 0) else ("#d97706" if r["role"].startswith("Stage3") else "#dc2626")
              for r, c in zip(p3_res, clashes)]
    y = np.arange(len(labels))
    bars = ax.barh(y, dmins, color=colors, edgecolor="#1e293b", height=0.62)
    ax.axvline(2.20, color="#dc2626", lw=1.6, ls="--", label="Severe Steric Clash Threshold (2.20 Å)")
    for idx, (dm, cl) in enumerate(zip(dmins, clashes)):
        tag = f"{dm:.2f} Å (0 clashes)" if cl == 0 else f"{dm:.2f} Å ({cl} clashes!)"
        ax.text(dm + 0.06, idx, tag, va="center", fontsize=8.0, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, max(dmins) + 1.45)
    ax.set_xlabel("Best-Rotamer Min Heavy-Atom Distance in Experimental Cryo-EM (Å)")
    ax.set_title("(B) Proof 3: 3D Rotamer Steric-Clash Audit (9I4I/9HZL/9I6F/8H9S)", fontsize=10.5, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8.0, frameon=True)
    ax.grid(True, axis="x", alpha=0.25)

    # Panel C: Multi-Model-Form & Adversarial Null Epistasis Proof (COX3 W116R)
    ax = axes[1, 0]
    scen_keys = list(p5_plot["scenario_results"].keys())
    scen_short = ["Baseline\n(20m/40m)", "Hill ROS\n(n=4)", "Severe I/R\n(30m/60m)", "Mild I/R\n(10m/40m)"]
    buf_plot = ["Yeast_NDI1", "Human_SOD2", "Ciona_AOX", "ADVERSARIAL_Null_Placebo", "ADVERSARIAL_Toxic_Uncoupler_FCCP"]
    buf_labels = ["Yeast NDI1", "Human SOD2", "Ciona AOX", "Null Placebo (Ctrl)", "Toxic FCCP (Ctrl)"]
    buf_cols = ["#0284c7", "#059669", "#7c3aed", "#64748b", "#dc2626"]
    x = np.arange(len(scen_keys))
    w = 0.15
    for bi, (bk, bl, bc) in enumerate(zip(buf_plot, buf_labels, buf_cols)):
        vals = [p5_plot["scenario_results"][sk][bk]["epistatic_excess_median"] for sk in scen_keys]
        ax.bar(x + (bi - 2) * w, vals, width=w, label=bl, color=bc, edgecolor="#0f172a", linewidth=0.5)
    ax.axhline(0.0, color="#0f172a", lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(scen_short, fontsize=8.5)
    ax.set_ylabel("Median Epistatic Excess $E_{add}$ (N=3,000 LHS)")
    ax.set_ylim(-0.18, 0.21)
    ax.set_title("(C) Proof 5A: Robustness Across Model Forms & Adversarial Nulls", fontsize=10.5, fontweight="bold")
    ax.legend(fontsize=7.6, loc="upper left", ncol=2, frameon=True)
    ax.grid(True, axis="y", alpha=0.25)

    # Panel D: Wet-Lab TALED Heteroplasmy Dose-Response Curve (0% -> 100%)
    ax = axes[1, 1]
    hx = p5_plot["het_grid_pct"]
    ax.plot(hx, p5_plot["het_driver_alone"], color="#94a3b8", lw=2.0, ls="--", label="COX3 W116R alone (unbuffered)")
    curve_styles = [
        ("Yeast_NDI1", "#0284c7", "COX3 W116R + Yeast NDI1"),
        ("Human_SOD2", "#059669", "COX3 W116R + Human SOD2"),
        ("Ciona_AOX", "#7c3aed", "COX3 W116R + Ciona AOX"),
        ("Human_UCP2", "#d97706", "COX3 W116R + Human UCP2"),
    ]
    for bk, col, lbl in curve_styles:
        ax.plot(hx, p5_plot["het_curves"][bk]["viability_med"], color=col, lw=2.3, marker="o", ms=3.5, label=lbl)
    ax.axvline(60, color="#0f172a", lw=1.2, ls=":", label="Typical Wet-Lab TALED Target (60%)")
    ax.set_xlabel("COX3 W116R Base-Editing Heteroplasmy (%)")
    ax.set_ylabel("Median Acute I/R Viability (N=2,000 LHS)")
    ax.set_title("(D) Proof 5B: Wet-Lab Heteroplasmy Dose-Response Thresholds", fontsize=10.5, fontweight="bold")
    ax.legend(fontsize=7.8, loc="upper left", frameon=True)
    ax.grid(True, alpha=0.25)

    plt.tight_layout(pad=1.5)
    plt.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    plt.close()


# =============================================================================
# MAIN
# =============================================================================
def main() -> None:
    print("[0/6] Loading experimental 3D cryo-EM microenvironments & 45-species MSA...")
    env3d, msa_stats, control_candidates, gene_info = load_all_3d_and_msa()
    s2_data = json.loads(s2.OUT_JSON.read_text(encoding="utf-8"))
    s3_data = json.loads(s3.OUT_JSON.read_text(encoding="utf-8"))

    print("[1/6] Running Proof 1: Confounder-Controlled & Negative-Control Mitonuclear Interface Proof...")
    p1 = run_proof1_interface_enrichment(env3d, msa_stats, control_candidates, s2_data)
    print(f"      Extreme={p1['extreme_polarized_2x2']['pct_interface']}% (OR={p1['extreme_polarized_2x2']['odds_ratio']}) vs "
          f"Matched Control Null={p1['matched_control_polarized_negative_null_2x2']['pct_interface']}% (OR={p1['matched_control_polarized_negative_null_2x2']['odds_ratio']}) | "
          f"CMH OR={p1['gene_stratified_cmh']['cmh_odds_ratio']} (p={p1['gene_stratified_cmh']['cmh_p_value']}) | "
          f"Multivar Adj OR={p1['multivariate_logistic_regression_adjusted']['adjusted_OR_interface']} (p={p1['multivariate_logistic_regression_adjusted']['wald_p_value']})")

    print("[2/6] Running Proof 2: Full-Proteome ClinVar Pathogenic vs Benign/gnomAD Calibration & OOS Transfer...")
    p2, p2_plot = run_proof2_proteome_pathogenicity_benchmark(env3d, msa_stats, gene_info, s3_data)
    print(f"      n_Path={p2['n_clinvar_pathogenic_missense_3d']}, n_Benign/gnomAD={p2['n_benign_or_gnomad_missense_3d']} | "
          f"5-Fold CV AUROC={p2['five_fold_cv_auroc']} | Strict OOS Tier-1->Tier-2 AUROC={p2['strict_oos_tier1_to_tier2_auroc']}")

    print("[3/6] Running Proof 3: 3D Atomic Contact Network & Rotamer Steric-Clash Audit...")
    p3 = run_proof3_3d_rotamer_clash_audit(gene_info)
    for r in p3:
        br = r["best_rotamer"]
        print(f"      {r['variant']:12s} ({r['role'][:20]:20s}): best d_min={br['min_heavy_atom_dist_A']:.2f} A, "
              f"clashes(<2.2A)={br['severe_clashes_lt_2p2A']}, clash_free_rots={r['clash_free_rotamers']}/{r['rotamers_scanned']} -> {r['steric_verdict']}")

    print("[4/6] Running Proof 4: Base-Editor Window Sensitivity Stress Test (+/-2 to +/-6 nt)...")
    p4 = run_proof4_editor_window_stress_test(gene_info, s3_data)
    print("      Zero-bystander counts by window:", p4["zero_bystander_count_by_window"])

    print("[5/6] Running Proof 5: Supermodularity Grid, Invariants, Multi-Model Forms, Adversarial Nulls & Heteroplasmy Sweep...")
    p5, p5_plot = run_proof5_epistasis_hardening()
    print(f"      Dual-cliff supermodularity concordance: {p5['supermodularity_grid_proof']['dual_cliff_theorem_concordance_pct']}% | "
          f"Physical invariant violations: {p5['physical_invariant_violations_across_all_scenarios']} | "
          f"Min heteroplasmy thresholds: {p5['heteroplasmy_min_thresholds']}")

    print("[6/6] Running Proof 6: End-to-End SQLite Referential Integrity Audit & Rendering Figure...")
    p6 = run_proof6_db_integrity()
    render_stage5_figure(p2_plot, p3, p5_plot)

    full = {
        "proof1_mitonuclear_interface_confounders_and_null": p1,
        "proof2_proteome_clinvar_benchmark_and_oos": p2,
        "proof3_3d_rotamer_steric_clash_audit": p3,
        "proof4_base_editor_window_stress_test": p4,
        "proof5_epistasis_supermodularity_models_nulls_heteroplasmy": p5,
        "proof6_database_referential_integrity": p6,
    }
    OUT_JSON.write_text(json.dumps(full, indent=2), encoding="utf-8")

    # Write detailed Stage 5 Markdown Report
    lines = [
        "# Stage 5 Report — Scientific Hardening, Adversarial Stress-Tests & Formal Proofs",
        "",
        "This audit directly stress-tests every statistical, 3D structural, genetic, base-editing, and biophysical claim across Stages 1–4 against negative controls, confounders, alternative model formulations, and adversarial null inputs.",
        "",
        "## Proof 1: Confounder-Controlled & Negative-Control Mitonuclear Interface Enrichment",
        f"- **Extreme-Polarized Candidates (`n = {p1['extreme_polarized_2x2']['n_candidates']}`):** `{p1['extreme_polarized_2x2']['pct_interface']}%` at mitonuclear interface (`<= 4.5 Å`) vs `{p1['extreme_polarized_2x2']['bg_pct_interface']}%` background (**`OR = {p1['extreme_polarized_2x2']['odds_ratio']}`**, Fisher one-sided `p = {p1['extreme_polarized_2x2']['fisher_p_greater']}`).",
        f"- **Matched Control-Branch Negative Null (`n = {p1['matched_control_polarized_negative_null_2x2']['n_candidates']}`):** `{p1['matched_control_polarized_negative_null_2x2']['pct_interface']}%` at interface (**`OR = {p1['matched_control_polarized_negative_null_2x2']['odds_ratio']}`**, `p = {p1['matched_control_polarized_negative_null_2x2']['fisher_p_two_sided']}`).",
        f"- **Gene-Stratified Cochran–Mantel–Haenszel (CMH) Test:** `OR_CMH = {p1['gene_stratified_cmh']['cmh_odds_ratio']}` (`p = {p1['gene_stratified_cmh']['cmh_p_value']}`).",
        f"- **3D Burial × Vertebrate Entropy Stratified Permutation (`N = 5,000`):** Observed = `{p1['burial_entropy_stratified_permutation_5000']['observed_interface_hits']}` vs Expected = `{p1['burial_entropy_stratified_permutation_5000']['expected_stratified_hits']}` (`{p1['burial_entropy_stratified_permutation_5000']['enrichment_fold']}x`, empirical `p = {p1['burial_entropy_stratified_permutation_5000']['empirical_p_value']}`).",
        f"- **Multivariate Logistic Regression (controlling for 3D burial, Shannon entropy & 12 gene fixed effects):** Adjusted **`OR = {p1['multivariate_logistic_regression_adjusted']['adjusted_OR_interface']}`** (`95% CI: {p1['multivariate_logistic_regression_adjusted']['adjusted_OR_95CI']}`, Wald `p = {p1['multivariate_logistic_regression_adjusted']['wald_p_value']}`).",
        f"- **Leave-One-Gene-Out Jackknife Stability:** `OR` remains in `[{p1['leave_one_gene_out_OR_range']['min_OR']}, {p1['leave_one_gene_out_OR_range']['max_OR']}]` across all 13 genes.",
        "",
        "## Proof 2: Full-Proteome Human ClinVar Pathogenic vs. Benign/gnomAD Benchmark & OOS Transfer",
        f"- **Dataset (100% live ClinVar + gnomAD v3 chrM mapped onto 3D cryo-EM):** `{p2['n_clinvar_pathogenic_missense_3d']}` Pathogenic/LP vs `{p2['n_benign_or_gnomad_missense_3d']}` Benign/Common-gnomAD missense variants across all 13 human mtDNA genes.",
        f"- **Zero-Leakage 5-Fold CV AUROC (comparative + 3D features only):** **`{p2['five_fold_cv_auroc']}`** (folds: `{p2['five_fold_aucs']}`).",
        f"- **Strict Out-of-Sample Tier-1 (6 genes) → Frozen Tier-2 (7 genes) Transfer AUROC:** **`{p2['strict_oos_tier1_to_tier2_auroc']}`** (reverse Tier-2 → Tier-1: **`{p2['strict_oos_tier2_to_tier1_auroc']}`**).",
        f"- **Predicted Pathogenicity Calibration:** Pathogenic median = `{p2['pathogenic_median_prob']}`, Benign median = `{p2['benign_median_prob']}`, Wild 131 candidates median = `{p2['wild_131_median_pathogenicity_prob']}`.",
        f"- **Top-5 Finalist Risk Scores:** `{p2['finalist_pathogenicity_probs']}` vs Pathogenic Controls `{p2['answer_key_control_probs']}`.",
        "",
        "## Proof 3: 3D Cryo-EM Atomic Neighborhood & Rotamer Steric-Clash Audit",
        "",
        "| Variant | Role | PDB | `dV (Å³)` | Rotamers Clash-Free | Best `d_min (Å)` | Severe Clashes (`<2.2 Å`) | Verdict | Key 5 Å Lipid / Nuclear Contacts |",
        "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- |",
    ]
    for r in p3:
        br = r["best_rotamer"]
        ct = ", ".join((r["lipid_ligand_contacts_5A"][:2] + r["nuclear_partner_contacts_5A"][:2])) or "mtDNA core"
        lines.append(
            f"| `{r['variant']}` | {r['role']} | `{r['pdb_id']}` | `{r['delta_volume_A3']:+.1f}` | `{r['clash_free_rotamers']}/{r['rotamers_scanned']}` | `{br['min_heavy_atom_dist_A']}` | `{br['severe_clashes_lt_2p2A']}` | **{r['steric_verdict']}** | `{ct}` |"
        )

    lines += [
        "",
        "## Proof 4: Base-Editor Window Stress Test (`±2 nt` to `±6 nt` on Human `rCRS`)",
        f"- **Zero-Bystander Survivors across Windows (`n = {p4['n_single_step_base_editable_tested']}` single-step editable):** `{p4['zero_bystander_count_by_window']}`.",
        "",
        "| Finalist | Editor | rCRS Edit | `±2 nt` | `±3 nt` | `±4 nt` | `±5 nt` | `±6 nt` |",
        "| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |",
    ]
    for f in p4["top5_finalists_window_profile"]:
        bw = f["bystanders_by_window"]
        lines.append(f"| `{f['variant']}` | `{f['editor_class']}` | `{f['rcrs_edit']}` | `{bw['pm2nt']}` | `{bw['pm3nt']}` | `{bw['pm4nt']}` | `{bw['pm5nt']}` | `{bw['pm6nt']}` |")

    lines += [
        "",
        "## Proof 5: Stage-4 Formal Supermodularity, Alternative Model Forms, Adversarial Nulls & Heteroplasmy Thresholds",
        f"- **Supermodularity Theorem Concordance:** `{p5['supermodularity_grid_proof']['dual_cliff_theorem_concordance_pct']}%` of dual-cliff operating points satisfy $\\Delta^2 V > 0$ (`max Δ²V = +{p5['supermodularity_grid_proof']['max_cross_partial_delta2V']}`).",
        f"- **Physical Invariant Audit:** **`{p5['physical_invariant_violations_across_all_scenarios']}`** violations across all trajectories (bounds `0 <= W, Q_red <= 1`, `Succ, ATP, EP >= 0`, zero `NaN/Inf`).",
        "",
        "| Model Form / Stress Protocol | `Yeast NDI1` (`E_add`, `%+`) | `Human SOD2` (`E_add`, `%+`) | `Ciona AOX` (`E_add`, `%+`) | `Null Placebo` (Adversarial) | `Toxic FCCP` (Adversarial) |",
        "| :--- | :---: | :---: | :---: | :---: | :---: |",
    ]
    for sk, sv in p5["multi_model_form_and_adversarial_controls"].items():
        lines.append(
            f"| `{sk}` | `+{sv['Yeast_NDI1']['epistatic_excess_median']}` (`{sv['Yeast_NDI1']['pct_lhs_positive']}%`) | "
            f"`+{sv['Human_SOD2']['epistatic_excess_median']}` (`{sv['Human_SOD2']['pct_lhs_positive']}%`) | "
            f"`+{sv['Ciona_AOX']['epistatic_excess_median']}` (`{sv['Ciona_AOX']['pct_lhs_positive']}%`) | "
            f"`{sv['ADVERSARIAL_Null_Placebo']['epistatic_excess_median']}` | "
            f"`{sv['ADVERSARIAL_Toxic_Uncoupler_FCCP']['epistatic_excess_median']}` |"
        )

    lines += [
        "",
        "### Wet-Lab Minimum Editing Heteroplasmy Thresholds (`COX3 W116R` `TALED`)",
    ]
    for bk, th in p5["heteroplasmy_min_thresholds"].items():
        h_e = f"{th['min_heteroplasmy_pct_for_Eadd_0p05']}%" if th['min_heteroplasmy_pct_for_Eadd_0p05'] is not None else ">100% (not reached at 0.05)"
        h_v = f"{th['min_heteroplasmy_pct_for_Viability_0p25']}%" if th['min_heteroplasmy_pct_for_Viability_0p25'] is not None else ">100%"
        lines.append(f"- **`{bk}`:** `E_add >= +0.05` achieved at **`{h_e}`** heteroplasmy; `Viability >= 0.25` achieved at **`{h_v}`** heteroplasmy.")

    lines += [
        "",
        "## Proof 6: End-to-End SQLite Referential Integrity",
        f"- **Foreign-key violations:** `{p6['foreign_key_violations']}` | **Orphan payoff rows:** `{p6['orphan_payoff_rows']}` | **Orphan transferability rows:** `{p6['orphan_transferability_rows']}` | **PASS:** `{p6['integrity_pass']}`.",
    ]
    OUT_REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nStage 5 Hardening & Proofs complete! Wrote {OUT_JSON}, {OUT_REPORT}, and {OUT_FIG}.")


if __name__ == "__main__":
    main()
