#!/usr/bin/env python3
"""Stage 10B: Tier-3 Ecological & Pollinator Thermal-Hypoxia Conservation Transfer.

Bridges MitoFitness Atlas's extremophile OXPHOS discovery engine and Mito_Papers' Apis-BCR
to 5 ecologically keystone non-human species facing climate warming, marine heatwaves,
riverine hypoxia, and agricultural Complex III fungicide stress — using 100% live NCBI
RefSeq mitogenome GenBank records (exact phylum-native NCBI /translation= CDS annotations):
  1. Apis mellifera      (Western Honeybee, NC_001566.1 — Table 5; Mito_Papers Paper III bridge)
  2. Bombus terrestris   (Buff-Tailed Bumblebee, KT368150.1 / NC_045071.1 — Table 5; alpine pollinator)
  3. Acropora millepora  (Great Barrier Reef Stony Coral, NC_003522.1 — Table 4; coral bleaching)
  4. Salmo salar         (Atlantic Salmon, NC_001960.1 — Table 2; warming-river cardiac hypoxia)
  5. Oncorhynchus mykiss (Rainbow Trout, NC_001717.1 — Table 2; cold-water salmonid CT_max)

Biosafety & Governance Principle:
  Strictly NON-GENE-DRIVE conservation biology — identifies naturally occurring wild thermal/hypoxia
  haplotypes for assisted natural gene flow / selective nursery breeding + nutritional/probiotic
  Apis-BCR buffering (zero synthetic gene drives released into wild populations).

Outputs:
  - data/stage10b_ecological_mitogenomes.json
  - data/stage10b_ecological_transfer.json
"""

from __future__ import annotations

import json
import math
import sqlite3
import time
import urllib.request
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from scipy.integrate import solve_ivp
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data" / "mitofitness.db"
CACHE_JSON = ROOT / "data" / "stage10b_ecological_mitogenomes.json"
OUT_JSON = ROOT / "data" / "stage10b_ecological_transfer.json"
SEED = 20261126

R_GAS = 8.314462618e-3  # kJ/(mol*K)

ECO_SPECIES = [
    ("Apis_mellifera",      "Western Honeybee",          "NC_001566.1", "Insecta (Hymenoptera)",   "Heatwave 42C flight + Complex III strobilurin fungicide (Mito_Papers bridge)"),
    ("Bombus_terrestris",   "Buff-Tailed Bumblebee",     "NC_045179.1", "Insecta (Hymenoptera)",   "Alpine/temperate pollinator heatwaves + high-altitude hypoxia"),
    ("Acropora_millepora",  "Great Barrier Reef Coral",  "NC_003522.1", "Cnidaria (Anthozoa)",     "Marine heatwave bleaching (27C -> 31.5C coral mitochondrial ROS)"),
    ("Salmo_salar",         "Atlantic Salmon",           "NC_001960.1", "Teleostei (Salmonidae)",  "Warming-river cardiac hypoxia (12C -> 22C CT_max collapse)"),
    ("Oncorhynchus_mykiss", "Rainbow Trout",             "NC_001717.1", "Teleostei (Salmonidae)",  "Cold-water salmonid thermal-oxygen ceiling"),
]

# Order COX3 before COX1 and enforce exact-name or length bounds so "subunit i" never matches "subunit iii"!
GENE_SPECS = [
    ("COX3", {"cox3", "coiii", "co3", "cytochrome c oxidase subunit iii", "cytochrome c oxidase subunit 3"}, 220, 310),
    ("COX1", {"cox1", "coi", "co1", "cytochrome c oxidase subunit i", "cytochrome c oxidase subunit 1"}, 450, 620),
    ("CYTB", {"cytb", "cob", "cyb", "cytochrome b"}, 340, 440),
    ("ATP6", {"atp6", "atpase6", "atp synthase f0 subunit 6", "atp synthase subunit 6"}, 190, 270),
    ("ND4L", {"nd4l", "nad4l", "nadh dehydrogenase subunit 4l"}, 75, 130),
]

KEY_SITES = [
    ("COX3", 116, "W", "R", "High_Altitude_Hypoxia (Bar-headed goose O2 affinity)"),
    ("COX1",  57, "I", "V", "Arctic_Cold (Muskox 0-bystander cold fluidity)"),
    ("ND4L",  85, "Y", "H", "Desert_Heat (Jerboa 0-bystander tight coupling)"),
    ("ND4L",  11, "A", "T", "Cross_Hypoxia (Bird+Mammal Complex I E-channel)"),
    ("CYTB", 181, "F", "L", "Arctic_Cold (Woolly Mammoth UQCRFS1 contact)"),
    ("CYTB", 224, "Y", "H", "Diving_Baric (Sperm Whale CYC1 contact)"),
    ("ATP6", 106, "I", "V", "Diving_Hypoxia (Sperm Whale + Pika ATP5MJ)"),
    ("ATP6", 195, "I", "V", "Arctic_Cold (Polar Bear thermogenic dial)"),
]

# Bos taurus / Human rCRS CYTB Complex III Qo strobilurin-binding pocket (Mito_Papers Paper III 1SQB)
QO_POCKET_HUMAN_POS = [128, 131, 142, 270, 271, 272, 273, 274, 278]  # F129/Y132/G143/PEWY/F275/Y279 (1-based human CYTB)


def fetch_genbank_proteins(acc: str) -> dict[str, str]:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gbwithparts&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "MitoFitnessAtlas/1.0"})
    text = urllib.request.urlopen(req, timeout=30).read().decode("utf-8", errors="replace")
    import io
    rec = SeqIO.read(io.StringIO(text), "genbank")
    prots = {}
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        tr = feat.qualifiers.get("translation", [])
        if not tr:
            continue
        names = set()
        for qk in ("gene", "product", "gene_synonym"):
            for v in feat.qualifiers.get(qk, []):
                names.add(v.strip().lower())
        seq = tr[0].strip("*")
        for std_g, aliases, min_len, max_len in GENE_SPECS:
            if std_g in prots or not (min_len <= len(seq) <= max_len):
                continue
            if names & aliases:
                prots[std_g] = seq
                break
    return prots


def align_and_map(human_seq: str, eco_seq: str) -> tuple[float, dict[int, str]]:
    """Return (% identity, {human_1based_pos: eco_aa_or_gap})."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    # Sanitize rare non-standard residues (U/X/B/Z)
    valid = set("ARNDCQEGHILKMFPSTWYV")
    h_clean = "".join(c if c in valid else "X" for c in human_seq)
    e_clean = "".join(c if c in valid else "X" for c in eco_seq)
    aln = aligner.align(h_clean, e_clean)[0]
    h_blocks, e_blocks = aln.aligned
    pos_map: dict[int, str] = {}
    matches = 0
    aligned_len = 0
    for (hs, he), (es, ee) in zip(h_blocks, e_blocks):
        for offset in range(he - hs):
            h_idx = hs + offset + 1  # 1-based
            e_aa = eco_seq[es + offset]
            pos_map[h_idx] = e_aa
            aligned_len += 1
            if human_seq[h_idx - 1] == e_aa:
                matches += 1
    pct_id = round(100.0 * matches / max(1, len(human_seq)), 2)
    return pct_id, pos_map


# =============================================================================
# ECOLOGICAL THERMAL-HYPOXIA-FUNGICIDE ARRHENIUS ODE (3 SCENARIOS)
# =============================================================================
ECO_PARAM_NAMES = [
    "km_o2_wt", "ea_oxphos", "k_leak_heat", "k_ret_ros", "theta_gsh", "f_load", "qo_sens",
]
ECO_NOMINAL = {
    "km_o2_wt": 1.20, "ea_oxphos": 42.0, "k_leak_heat": 0.26, "k_ret_ros": 1.90,
    "theta_gsh": 0.54, "f_load": 0.72, "qo_sens": 1.00,
}

ECO_INTERVENTIONS = {
    "WT_Vulnerable":          {"km_mult": 1.00, "heat_res": 0.00, "atp_gain": 0.00, "qo_clear": 0.00, "theta_boost": 0.00, "ucp_leak": 0.00, "atp_cost": 0.00},
    "Null_Placebo":           {"km_mult": 1.00, "heat_res": 0.00, "atp_gain": 0.00, "qo_clear": 0.00, "theta_boost": 0.00, "ucp_leak": 0.00, "atp_cost": 0.00},
    "Uncoupler_UCP2_Only":    {"km_mult": 1.00, "heat_res": 0.00, "atp_gain": 0.00, "qo_clear": 0.00, "theta_boost": 0.06, "ucp_leak": 0.24, "atp_cost": 0.10},
    "Natural_TightCoupled":   {"km_mult": 0.90, "heat_res": 0.44, "atp_gain": 0.14, "qo_clear": 0.00, "theta_boost": 0.00, "ucp_leak": -0.07, "atp_cost": 0.00},
    "Natural_HighO2_Affinity":{"km_mult": 0.46, "heat_res": 0.16, "atp_gain": 0.12, "qo_clear": 0.00, "theta_boost": 0.00, "ucp_leak": -0.03, "atp_cost": 0.00},
    "Nutritional_BCR_Buffer": {"km_mult": 1.00, "heat_res": 0.10, "atp_gain": 0.02, "qo_clear": 0.72, "theta_boost": 0.44, "ucp_leak": 0.00, "atp_cost": 0.01},
    "Combined_Haplotype_BCR": {"km_mult": 0.60, "heat_res": 0.44, "atp_gain": 0.15, "qo_clear": 0.72, "theta_boost": 0.44, "ucp_leak": -0.06, "atp_cost": 0.01},
}


def eco_env(t_min: float, scenario: str) -> tuple[float, float, float, float, float]:
    """Return (T_c, T_ref_c, pO2_kPa, f_mult, qo_block)."""
    if scenario == "E1_Honeybee_42C_Fungicide_Flight":
        # 35 C normal flight ref -> 42.5 C heatwave foraging + 65% strobilurin Complex III Qo block (60 min)
        r = min(1.0, t_min / 10.0)
        return 35.0 + 7.5 * r, 35.0, 14.0, 1.0 + 0.55 * r, 0.65
    elif scenario == "E2_Coral_31_5C_Marine_Heatwave":
        # 27 C coral reef ref -> 31.5 C sustained marine heatwave (+4.5 C) for 120 min scaled stress
        r = min(1.0, t_min / 20.0)
        return 27.0 + 4.5 * r, 27.0, 6.5 - 2.0 * r, 0.90 + 0.25 * r, 0.0
    elif scenario == "E3_Salmon_22C_WarmRiver_Hypoxia":
        # 12 C cold-water salmon ref -> 22 C warm river (+10 C!) + dissolved O2 crash 8.0 -> 0.85 kPa
        r = min(1.0, t_min / 15.0)
        return 12.0 + 10.0 * r, 12.0, 8.0 - 7.15 * r, 0.85 + 0.45 * r, 0.0
    else:
        raise ValueError(scenario)


def eco_ode_rhs(t_min: float, y: np.ndarray, p: dict[str, np.ndarray], iv: dict, scenario: str) -> np.ndarray:
    W = np.clip(y[0], 0.02, 1.2)
    Q = np.clip(y[1], 0.02, 0.99)
    A = np.clip(y[2], 0.01, 1.3)

    T_c, T_ref_c, pO2, f_mult, qo_blk_raw = eco_env(t_min, scenario)
    T_k = T_c + 273.15
    T_ref_k = T_ref_c + 273.15
    dT = max(0.0, T_c - T_ref_c)

    # Arrhenius enzyme rate & exponential thermal proton leak above species' native T_ref
    k_temp = np.exp(np.clip(-(p["ea_oxphos"] / R_GAS) * (1.0 / T_k - 1.0 / T_ref_k), -2.5, 1.2))
    base_leak = np.clip(0.14 + iv["ucp_leak"], 0.05, 0.65)
    g_leak = base_leak * np.exp(p["k_leak_heat"] * dT * (1.0 - 0.55 * iv["heat_res"]))
    if iv["ucp_leak"] > 0.0:
        g_leak = g_leak + iv["ucp_leak"] * 0.55 * (dT / 5.0)

    # Complex III Qo fungicide blockade (mitigated by Apis-BCR p-coumaric CYP9Q detox!)
    qo_eff = qo_blk_raw * p["qo_sens"] * (1.0 - iv["qo_clear"])
    c3_pass = np.clip(1.0 - qo_eff, 0.12, 1.0)

    km_eff = p["km_o2_wt"] * iv["km_mult"]
    g_o2 = pO2 / (km_eff + pO2)
    v_etc = 1.18 * W * k_temp * c3_pass * g_o2 * Q

    v_in = 0.72 * k_temp * (1.0 - Q)
    dQ = v_in - 0.88 * v_etc

    pmf = np.clip((0.76 + 0.55 * Q) / (1.0 + 1.35 * g_leak), 0.15, 1.35)
    ros_qo = 12.0 * p["k_ret_ros"] * (qo_eff + 0.15) * (np.maximum(0.0, Q - 0.44) ** 2) * (pmf ** 2) * g_o2 * k_temp
    ros_heat = 0.95 * (dT / 5.0) * (1.0 - 0.50 * iv["heat_res"]) * W
    ros_excess = np.maximum(0.0, (ros_qo + ros_heat) - (p["theta_gsh"] + iv["theta_boost"]))

    dW = -0.028 * ros_excess * W

    adp_avail = np.clip(1.25 - A, 0.0, 1.25) / 0.35
    adp_factor = adp_avail / (1.0 + adp_avail)
    coup = np.clip((1.0 + iv["atp_gain"] - iv["atp_cost"]) / (1.0 + 1.55 * g_leak), 0.10, 1.30)
    dA = 2.60 * v_etc * coup * adp_factor * (1.0 - 0.15 * np.minimum(1.0, ros_excess)) - p["f_load"] * f_mult * (A / (0.22 + A))
    return np.stack([dW, dQ, dA], axis=0)


def sim_eco_batch(p: dict[str, np.ndarray], iv_name: str, scenario: str) -> dict[str, np.ndarray]:
    iv = ECO_INTERVENTIONS[iv_name]
    n = len(p["km_o2_wt"])
    y = np.zeros((3, n), dtype=np.float64)
    y[0, :], y[1, :], y[2, :] = 0.95, 0.35, 0.95
    t_end = 120.0 if scenario == "E2_Coral_31_5C_Marine_Heatwave" else 60.0
    dt = 0.15
    n_steps = int(round(t_end / dt))
    min_a = np.full(n, 99.0)
    for k in range(n_steps):
        t = k * dt
        min_a = np.minimum(min_a, y[2])
        k1 = eco_ode_rhs(t, y, p, iv, scenario)
        k2 = eco_ode_rhs(t + 0.5 * dt, y + 0.5 * dt * k1, p, iv, scenario)
        k3 = eco_ode_rhs(t + 0.5 * dt, y + 0.5 * dt * k2, p, iv, scenario)
        k4 = eco_ode_rhs(t + dt - 1e-12, y + dt * k3, p, iv, scenario)
        y = y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        y[0] = np.clip(y[0], 0.02, 1.2)
        y[1] = np.clip(y[1], 0.02, 0.99)
        y[2] = np.clip(y[2], 0.01, 1.3)
    viab = 0.50 * y[0] + 0.30 * np.minimum(1.0, y[2]) + 0.20 * np.minimum(1.0, min_a)
    return {"W_end": y[0], "ATP_end": y[2], "ATP_nadir": min_a, "viability": viab}


def verify_eco_solver() -> dict:
    p_nom = {k: np.array([v], dtype=np.float64) for k, v in ECO_NOMINAL.items()}
    max_d = 0.0
    for sc, knots in [
        ("E1_Honeybee_42C_Fungicide_Flight", [0.0, 10.0, 60.0]),
        ("E2_Coral_31_5C_Marine_Heatwave", [0.0, 20.0, 120.0]),
        ("E3_Salmon_22C_WarmRiver_Hypoxia", [0.0, 15.0, 60.0]),
    ]:
        rk = sim_eco_batch(p_nom, "Combined_Haplotype_BCR", sc)
        iv = ECO_INTERVENTIONS["Combined_Haplotype_BCR"]

        def rhs(t, y_sc):
            return eco_ode_rhs(t, y_sc[:, None], p_nom, iv, sc)[:, 0]

        for m in ("DOP853", "Radau"):
            yc = np.array([0.95, 0.35, 0.95], dtype=np.float64)
            for a, b in zip(knots[:-1], knots[1:]):
                yc = solve_ivp(rhs, [a, b], yc, method=m, rtol=1e-10, atol=1e-12).y[:, -1]
            max_d = max(max_d, abs(float(rk["W_end"][0]) - float(yc[0])), abs(float(rk["ATP_end"][0]) - float(yc[2])))
    return {"max_abs_delta": float(max_d), "passed": bool(max_d < 1e-4)}


def main() -> None:
    print("=== Stage 10B: Tier-3 Ecological & Pollinator Thermal-Hypoxia Transfer ===")
    conn = sqlite3.connect(DB_PATH)
    human_seqs = dict(conn.execute("SELECT gene_symbol, protein_aa_seq FROM human_rcrs_genes").fetchall())
    conn.close()

    # 1. Fetch or load live NCBI RefSeq mitogenomes
    if CACHE_JSON.exists():
        eco_cache = json.loads(CACHE_JSON.read_text(encoding="utf-8"))
    else:
        eco_cache = {}
        for sp_key, common, acc, clade, stress in ECO_SPECIES:
            print(f"  Fetching live NCBI RefSeq {acc} ({common})...", flush=True)
            prots = fetch_genbank_proteins(acc)
            eco_cache[sp_key] = {"common": common, "accession": acc, "clade": clade, "stress": stress, "proteins": prots}
            print(f"    Extracted {len(prots)}/5 OXPHOS proteins: {sorted(prots.keys())}", flush=True)
            time.sleep(0.35)
        CACHE_JSON.write_text(json.dumps(eco_cache, indent=2), encoding="utf-8")

    # 2. Cross-Phylum Alignment onto Human Cryo-EM Coordinates & Key Fitness-Point Audit
    alignment_audit = {}
    for sp_key, common, acc, clade, stress in ECO_SPECIES:
        prots = eco_cache[sp_key]["proteins"]
        gene_ids = {}
        site_calls = {}
        qo_conserved = 0
        for g in ["COX1", "COX3", "CYTB", "ATP6", "ND4L"]:
            if g not in prots:
                continue
            pid, pmap = align_and_map(human_seqs[g], prots[g])
            gene_ids[g] = {"pct_identity_vs_human": pid, "aa_len": len(prots[g])}
            for sg, pos, href, ealt, desc in KEY_SITES:
                if sg == g:
                    eco_aa = pmap.get(pos, "-")
                    status = "PRE_ADAPTED_EXTREMOPHILE_ALLELE" if eco_aa == ealt else ("CONSERVED_HUMAN_WT" if eco_aa == href else "LINEAGE_SPECIFIC")
                    site_calls[f"{sg}_{href}{pos}{ealt}"] = {"eco_aa": eco_aa, "status": status}
            if g == "CYTB":
                for qpos in QO_POCKET_HUMAN_POS:
                    if pmap.get(qpos, "-") == human_seqs["CYTB"][qpos - 1]:
                        qo_conserved += 1
        alignment_audit[sp_key] = {
            "common_name": common,
            "refseq_accession": acc,
            "clade": clade,
            "primary_ecological_stress": stress,
            "gene_identities_vs_human": gene_ids,
            "qo_strobilurin_pocket_conservation": f"{qo_conserved}/{len(QO_POCKET_HUMAN_POS)}",
            "fitness_point_sites": site_calls,
        }
        print(
            f"  [Align] {common:25s} ({acc}) | COX1={gene_ids.get('COX1',{}).get('pct_identity_vs_human',0):.1f}% "
            f"| CYTB={gene_ids.get('CYTB',{}).get('pct_identity_vs_human',0):.1f}% | Qo pocket={qo_conserved}/9"
        )

    # 3. Run 3 Ecological Stress ODE Scenarios (N=5,000 LHS)
    conv = verify_eco_solver()
    print(f"\n[ODE] Ecological 3-scenario solver convergence: max |Delta| = {conv['max_abs_delta']:.3e} (passed={conv['passed']})")
    assert conv["passed"]

    sampler = qmc.LatinHypercube(d=len(ECO_PARAM_NAMES), seed=SEED)
    u = sampler.random(5000)
    p_lhs = {name: ECO_NOMINAL[name] * (0.70 + 0.60 * u[:, i]) for i, name in enumerate(ECO_PARAM_NAMES)}

    scenarios_out = {}
    for sc, label in [
        ("E1_Honeybee_42C_Fungicide_Flight", "Apis mellifera: 42.5 C Heatwave + 65% Complex III Strobilurin Flight"),
        ("E2_Coral_31_5C_Marine_Heatwave",   "Acropora millepora: 27 C -> 31.5 C Marine Heatwave Bleaching"),
        ("E3_Salmon_22C_WarmRiver_Hypoxia",  "Salmo salar: 12 C -> 22 C Warming River + Dissolved O2 Crash"),
    ]:
        wt = sim_eco_batch(p_lhs, "WT_Vulnerable", sc)
        pl = sim_eco_batch(p_lhs, "Null_Placebo", sc)
        assert float(np.max(np.abs(pl["viability"] - wt["viability"]))) == 0.0
        arms = {}
        print(f"\n  --- {label} (WT med = {np.median(wt['viability']):.4f}) ---")
        for iv_name in ["WT_Vulnerable", "Uncoupler_UCP2_Only", "Natural_TightCoupled", "Natural_HighO2_Affinity", "Nutritional_BCR_Buffer", "Combined_Haplotype_BCR"]:
            r = sim_eco_batch(p_lhs, iv_name, sc)
            dv = r["viability"] - wt["viability"]
            arms[iv_name] = {
                "viab_med": round(float(np.median(r["viability"])), 4),
                "delta_vs_wt_med": round(float(np.median(dv)), 4),
                "delta_vs_wt_p05": round(float(np.percentile(dv, 5)), 4),
                "delta_vs_wt_p95": round(float(np.percentile(dv, 95)), 4),
                "atp_nadir_med": round(float(np.median(r["ATP_nadir"])), 4),
                "w_end_med": round(float(np.median(r["W_end"])), 4),
            }
            print(f"    {iv_name:24s} | V={arms[iv_name]['viab_med']:.4f} (dV={arms[iv_name]['delta_vs_wt_med']:+.4f}) | ATP_nadir={arms[iv_name]['atp_nadir_med']:.3f}")
        scenarios_out[sc] = {"label": label, "wt_viability_med": arms["WT_Vulnerable"]["viab_med"], "arms": arms}

    payload = {
        "stage": "10B_Tier3_Ecological_Pollinator_Conservation_Transfer",
        "epistemic_tier": "E1_GENOMICS_plus_E3_SIM",
        "biosafety_governance": "Strictly non-gene-drive conservation biology (assisted natural haplotype selection & nutritional Apis-BCR buffering only).",
        "solver_convergence": conv,
        "species_alignment_audit": alignment_audit,
        "ecological_scenarios_lhs5k": scenarios_out,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nSaved -> {OUT_JSON}")


if __name__ == "__main__":
    main()
