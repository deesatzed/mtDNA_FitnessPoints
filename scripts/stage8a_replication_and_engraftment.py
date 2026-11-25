#!/usr/bin/env python3
"""Stage 8A: Full Mito_Papers 3-State Mitochondrial Replication, Trap A/B Unmasking,
Context-Sensitive Donor-Mitochondria Reversibility & Scheduled-Surgery Pre-Op Clocks.

Upgrades Stage 7A ([W_wt, W_mut]) by integrating the 4 missing replication physics
mechanisms from Mito_Papers (mitorestorer/simulator/sim_bioenergetic_reset.py):
  1. Trap B — Reverse-ATPase Phantom Voltage Masking under high glucose (25 mM DMEM/RPMI):
     stalled WT and D_aged mitochondria burn glycolytic ATP via reverse Complex V to fake
     DeltaPsi_m and blind PINK1 unless glycolytic unmasking (galactose / PDK1 AZD7545) is on.
  2. Trap A — POLG/TFAM Replicative Speed Advantage of shorter deletion-damaged mtDNA
     (k_rep_adv = 1.32x) in aged patients (0%, 15%, 30%, 45% pre-existing D_aged).
  3. mTORC1–ULK1(Ser757) Mitophagy-vs-Biogenesis Antagonism & Induction Lag (tau_lag = 6-13.5 h).
  4. 3-State [W_host, D_aged, M_donor] Donor-Mitochondria Context-Sensitive Reversibility:
     auto-dilutes in healthy normoxic tissue (zero POLG speed advantage + mild exo turnover)
     vs. selectively expands inside chronically hypoxic/aged injured zones under cyclic BCR.

Also computes the 4-Clock Scheduled-Surgery Pre-Op Lead-Time Efficacy Curves (-48 h to +12 h)
linking Clock 1 (malonate/MitoSNO), Clock 2 (TAT-SOD2 / banked donor mitochondria),
Clock 3 (transient SOD2+UCP2 mRNA-LNP), and Clock 3 + tau_lag conditioning.

Outputs:
  - data/stage8a_replication_engraftment.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "data" / "stage8a_replication_engraftment.json"
SEED = 20261124

# =============================================================================
# 3-STATE MITOCHONDRIAL REPLICATION MODEL: y = [W_host, D_aged, M_eng]
# =============================================================================
REP_PARAM_NAMES = [
    "Km_wt",        # WT human Complex IV O2 half-saturation (kPa)
    "Km_mut",       # COX3 W116R Complex IV O2 half-saturation (kPa)
    "k_rep_adv",    # Trap A: POLG replicative speed advantage of shorter deletion mtDNA (D_aged)
    "k_mito_base",  # Basal mitophagy rate (/day)
    "k_mito_sel",   # PINK1/Urolithin A depolarization-gated selective mitophagy gain (/day)
    "k_bio_max",    # PGC-1a/NAD+ biogenesis replication rate (/day)
    "pink1_eff",    # Host PINK1 kinase efficiency (1.0 = WT, lower in risk genotypes)
    "s_bias_err",   # WT vs COX3_W116R depolarization discrimination fidelity
]

REP_NOMINAL = {
    "Km_wt": 1.25,
    "Km_mut": 0.53,
    "k_rep_adv": 1.32,
    "k_mito_base": 0.05,
    "k_mito_sel": 0.92,
    "k_bio_max": 0.62,
    "pink1_eff": 1.00,
    "s_bias_err": 1.00,
}


def rep3_ode_rhs_vec(
    t_day: float,
    y: np.ndarray,
    p: dict[str, np.ndarray],
    regimen: str,
    pO2_hypox: float = 1.10,
    k_exo_clear: float = 0.0,
) -> np.ndarray:
    """Vectorized 3-state organelle population ODE: y = [W_host, D_aged, M_eng].

    Matches Mito_Papers' exact biophysical gating equations:
      - voltage_unmasking = 1 / (1 + exp(10 * (glycolytic_flux - 0.35)))
      - ulk1_activation   = 1 - 0.72 * mTORC1_anabolic
      - Trap A: bio_rate_D carries k_rep_adv * (1 - 0.45 * voltage_unmasking)
      - Trap B: selective_mito for both D_aged and hypoxia-stalled W_host requires
                voltage_unmasking (reverse-ATPase blockade) AND ulk1_activation!
    """
    W = np.maximum(1e-6, y[0])
    D = np.maximum(1e-6, y[1])
    M = np.maximum(1e-6, y[2])
    tot = W + D + M

    # Regimen schedules returning (pO2, glycolytic_flux, mTORC1, u_mito, pgc1a_drive)
    if regimen == "Adversarial_Null_Placebo":
        pO2, glyc, mtor, u_mito, pgc = 14.0, 0.85, 0.75, 0.0, 0.25
    elif regimen == "HighGlucose_25mM_TrapB_Failure":
        # Classic wet-lab pitfall: mild hypoxia + Urolithin A in standard 25 mM glucose DMEM/RPMI!
        # High glycolytic_flux (0.90) powers reverse Complex V -> phantom DeltaPsi masks WT & D!
        if t_day < 7.0:
            pO2, glyc, mtor, u_mito, pgc = pO2_hypox, 0.90, 0.70, 1.00, 0.25
        else:
            pO2, glyc, mtor, u_mito, pgc = 14.0, 0.90, 0.85, 0.05, 1.00
    elif regimen == "Unphased_Biogenesis_TrapA":
        # Taking NAD+/PGC-1a + anabolic mTORC1 ON continuously without Stage-1 purge
        pO2, glyc, mtor, u_mito, pgc = 14.0, 0.75, 0.88, 0.50, 1.15
    elif regimen == "Biphasic_Unmasked_Galactose_PDK1":
        # Stage 1 (Days 0-7): Mild hypoxia + galactose/AZD7545 glycolytic unmasking + low mTORC1 + UroA
        # Stage 2 (Days 7-14): Normoxia + PGC-1a/NAD+ biogenesis expansion
        if t_day < 7.0:
            # Include 6-hour lag-compressed induction sigmoid (tau_lag = 0.25 d)
            cat_depth = 1.0 / (1.0 + math.exp(-(t_day - 0.25) / 0.08))
            pO2 = pO2_hypox
            glyc = 0.45 * (1.0 - cat_depth) + 0.14 * cat_depth
            mtor = 0.65 * (1.0 - cat_depth) + 0.18 * cat_depth
            u_mito = 0.30 + 0.70 * cat_depth
            pgc = 0.10
        else:
            pO2, glyc, mtor, u_mito, pgc = 14.0, 0.40, 0.88, 0.05, 1.15
    elif regimen == "TwoCycle_Pulsed_Unmasked_Safe":
        # 2 x (3.5 d unmasked mild purge + 3.5 d normoxic biogenesis) -> safe mass floor!
        cyc = t_day % 7.0
        if cyc < 3.5:
            cat_depth = 1.0 / (1.0 + math.exp(-(cyc - 0.25) / 0.08))
            pO2 = pO2_hypox
            glyc = 0.45 * (1.0 - cat_depth) + 0.14 * cat_depth
            mtor = 0.65 * (1.0 - cat_depth) + 0.20 * cat_depth
            u_mito = 0.25 + 0.65 * cat_depth
            pgc = 0.16
        else:
            pO2, glyc, mtor, u_mito, pgc = 14.0, 0.40, 0.88, 0.05, 1.15
    elif regimen == "Donor_Healthy_Normoxic_60d":
        # Healthy transplanted organ or normoxic tissue after acute event: normal fed physiology
        pO2, glyc, mtor, u_mito, pgc = 14.0, 0.65, 0.75, 0.10, 0.40
    elif regimen == "Donor_PostTransplant_24hHypox_ThenNormox":
        # 24 h peri-transplant low-O2 window followed by 59 d healthy normoxia
        if t_day < 1.0:
            pO2, glyc, mtor, u_mito, pgc = 1.20, 0.35, 0.30, 0.40, 0.20
        else:
            pO2, glyc, mtor, u_mito, pgc = 14.0, 0.65, 0.75, 0.10, 0.40
    elif regimen == "Donor_ChronicHypox_Aged_CyclicBCR":
        # Chronically hypoxic / aged damaged tissue (COPD muscle, ischemic limb, failing RV)
        # under cyclic weekly BCR (2.5 d unmasked purge : 4.5 d normoxic/mild-hypoxic build)
        cyc = t_day % 7.0
        pO2 = 1.35  # Persistent tissue hypoxia (e.g., ~10 mmHg local tissue pO2)
        if cyc < 2.5:
            cat_depth = 1.0 / (1.0 + math.exp(-(cyc - 0.25) / 0.08))
            glyc = 0.40 * (1.0 - cat_depth) + 0.15 * cat_depth
            mtor = 0.60 * (1.0 - cat_depth) + 0.20 * cat_depth
            u_mito = 0.25 + 0.70 * cat_depth
            pgc = 0.14
        else:
            glyc, mtor, u_mito, pgc = 0.40, 0.85, 0.08, 1.05
    else:
        raise ValueError(f"Unknown regimen: {regimen}")

    # 1. Mito_Papers Trap B & mTORC1-ULK1(Ser757) gates
    voltage_unmasking = 1.0 / (1.0 + math.exp(10.0 * (glyc - 0.35)))
    ulk1_activation = 1.0 - 0.72 * mtor

    # 2. Complex IV O2 conductance & depolarization probabilities
    g_wt = pO2 / (p["Km_wt"] + pO2)
    g_mut = pO2 / (p["Km_mut"] + pO2)
    dep_wt = 1.0 / (1.0 + np.exp(6.0 * (g_wt - 0.58)))
    dep_mut = 1.0 / (1.0 + np.exp(6.0 * (g_mut - 0.58)))
    # Deletion-damaged mtDNA (D_aged) is constitutively bioenergetically defective (dep_D ~ 0.94)
    dep_d = 0.94

    eff_dep_wt = dep_wt * p["s_bias_err"]
    eff_dep_mut = dep_mut / np.maximum(0.5, p["s_bias_err"])

    # Selective mitophagy requires BOTH voltage unmasking (Trap B off) AND ULK1 active (low mTORC1)
    gate_sel = p["pink1_eff"] * u_mito * voltage_unmasking * ulk1_activation
    clear_W = p["k_mito_base"] + p["k_mito_sel"] * gate_sel * eff_dep_wt
    # Basal homeostatic quality control clears defective D slightly faster than healthy W in normoxia,
    # while unmasked PINK1/UroA accelerates D clearance dramatically
    clear_D = p["k_mito_base"] * 1.18 + p["k_mito_sel"] * 1.25 * gate_sel * dep_d
    clear_M = p["k_mito_base"] + k_exo_clear + p["k_mito_sel"] * gate_sel * eff_dep_mut

    # 3. Carrying-capacity-limited POLG/PGC-1a biogenesis (Trap A: D replicates k_rep_adv faster!)
    cap_factor = np.maximum(0.0, 1.0 - (tot / 1.05) ** 2)
    anabolic_drive = p["k_bio_max"] * pgc * (0.25 + 0.95 * mtor) * cap_factor
    maint = 0.055 * cap_factor

    bio_W = anabolic_drive + maint
    bio_M = anabolic_drive + maint  # Full-length 16,569 bp edit has NO POLG speed advantage over WT!
    bio_D = anabolic_drive * p["k_rep_adv"] * (1.0 - 0.45 * voltage_unmasking) + maint

    # 4. Stress-coupled ROS conversion W -> D (tiny ~0.08%/d in healthy normoxia; rises under hypoxic stall)
    k_dmg_W = 0.0008 + 0.014 * dep_wt
    k_dmg_M = 0.0005 + 0.003 * dep_mut

    dW = (bio_W - clear_W - k_dmg_W) * W
    dD = (bio_D - clear_D) * D + k_dmg_W * W + k_dmg_M * M
    dM = (bio_M - clear_M - k_dmg_M) * M
    return np.stack([dW, dD, dM], axis=0)


def simulate_rep3_batch(
    p: dict[str, np.ndarray],
    regimen: str,
    w0: float,
    d0: float,
    m0: float,
    t_end: float = 14.0,
    dt: float = 0.02,
    pO2_hypox: float = 1.10,
    k_exo_clear: float = 0.0,
) -> dict[str, np.ndarray]:
    n = len(p["Km_wt"])
    y = np.zeros((3, n), dtype=np.float64)
    y[0, :] = w0
    y[1, :] = d0
    y[2, :] = m0

    n_steps = int(round(t_end / dt))
    min_tot = np.full(n, 99.0)
    for k in range(n_steps):
        t = k * dt
        tot = y[0] + y[1] + y[2]
        min_tot = np.minimum(min_tot, tot)
        # Keep k4 strictly inside the current half-open interval when a discrete switch sits at (k+1)*dt
        t_mid = t + 0.5 * dt
        t_end_step = t + dt - 1e-12
        k1 = rep3_ode_rhs_vec(t, y, p, regimen, pO2_hypox, k_exo_clear)
        k2 = rep3_ode_rhs_vec(t_mid, np.maximum(1e-6, y + 0.5 * dt * k1), p, regimen, pO2_hypox, k_exo_clear)
        k3 = rep3_ode_rhs_vec(t_mid, np.maximum(1e-6, y + 0.5 * dt * k2), p, regimen, pO2_hypox, k_exo_clear)
        k4 = rep3_ode_rhs_vec(t_end_step, np.maximum(1e-6, y + dt * k3), p, regimen, pO2_hypox, k_exo_clear)
        y = np.maximum(1e-6, y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4))

    tot = y[0] + y[1] + y[2]
    return {
        "W_end": y[0],
        "D_end": y[1],
        "M_end": y[2],
        "tot_end": tot,
        "tot_nadir": min_tot,
        "f_M_end": y[2] / np.maximum(1e-9, tot),
        "f_D_end": y[1] / np.maximum(1e-9, tot),
        "f_W_end": y[0] / np.maximum(1e-9, tot),
    }


def simulate_rep3_trace(
    regimen: str,
    w0: float,
    d0: float,
    m0: float,
    t_end: float,
    dt: float = 0.02,
    sample_stride: int = 25,
    k_exo_clear: float = 0.0,
) -> dict:
    p = {k: np.array([v], dtype=np.float64) for k, v in REP_NOMINAL.items()}
    y = np.array([[w0], [d0], [m0]], dtype=np.float64)
    n_steps = int(round(t_end / dt))
    t_out, fM_out, fD_out, tot_out = [], [], [], []
    for k in range(n_steps + 1):
        t = k * dt
        if k % sample_stride == 0 or k == n_steps:
            tot = float((y[0] + y[1] + y[2])[0])
            t_out.append(round(t, 2))
            fM_out.append(round(float(y[2][0] / tot), 4))
            fD_out.append(round(float(y[1][0] / tot), 4))
            tot_out.append(round(tot, 4))
        if k < n_steps:
            k1 = rep3_ode_rhs_vec(t, y, p, regimen, 1.10, k_exo_clear)
            k2 = rep3_ode_rhs_vec(t + 0.5 * dt, np.maximum(1e-6, y + 0.5 * dt * k1), p, regimen, 1.10, k_exo_clear)
            k3 = rep3_ode_rhs_vec(t + 0.5 * dt, np.maximum(1e-6, y + 0.5 * dt * k2), p, regimen, 1.10, k_exo_clear)
            k4 = rep3_ode_rhs_vec(t + dt - 1e-12, np.maximum(1e-6, y + dt * k3), p, regimen, 1.10, k_exo_clear)
            y = np.maximum(1e-6, y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4))
    return {"t_days": t_out, "f_M": fM_out, "f_D": fD_out, "tot_mass": tot_out}


def verify_solver_convergence() -> dict:
    """Verify RK4 (dt=0.02 d) against scipy DOP853 & stiff Radau (rtol=1e-10, atol=1e-12)."""
    p_nom = {k: np.array([v], dtype=np.float64) for k, v in REP_NOMINAL.items()}
    # Aged cell baseline: 40% W_host, 20% D_aged, 40% M_edit (total 0.92)
    w0, d0, m0 = 0.92 * 0.40, 0.92 * 0.20, 0.92 * 0.40
    rk4 = simulate_rep3_batch(p_nom, "Biphasic_Unmasked_Galactose_PDK1", w0, d0, m0, t_end=14.0, dt=0.02)

    def rhs_scalar(t, y):
        dy = rep3_ode_rhs_vec(
            t,
            np.array([[y[0]], [y[1]], [y[2]]], dtype=np.float64),
            p_nom,
            "Biphasic_Unmasked_Galactose_PDK1",
        )
        return dy[:, 0]

    # Integrate piece-wise across Day-7 switch for exact event handling
    def solve_piecewise(method: str):
        sol1 = solve_ivp(rhs_scalar, [0.0, 7.0], [w0, d0, m0], method=method, rtol=1e-10, atol=1e-12)
        sol2 = solve_ivp(rhs_scalar, [7.0, 14.0], sol1.y[:, -1], method=method, rtol=1e-10, atol=1e-12)
        y_end = sol2.y[:, -1]
        return float(y_end[2] / np.sum(y_end)), float(y_end[1] / np.sum(y_end))

    dop_m, dop_d = solve_piecewise("DOP853")
    rad_m, rad_d = solve_piecewise("Radau")
    rk4_m, rk4_d = float(rk4["f_M_end"][0]), float(rk4["f_D_end"][0])
    max_diff = max(abs(rk4_m - dop_m), abs(rk4_m - rad_m), abs(rk4_d - dop_d), abs(rk4_d - rad_d))
    return {
        "rk4_fM": round(rk4_m, 7),
        "dop853_fM": round(dop_m, 7),
        "radau_fM": round(rad_m, 7),
        "rk4_fD": round(rk4_d, 7),
        "dop853_fD": round(dop_d, 7),
        "radau_fD": round(rad_d, 7),
        "max_abs_delta": float(max_diff),
        "passed": bool(max_diff < 1e-4),
    }


def sample_lhs_params(n: int, seed: int, rel_halfwidth: float = 0.30) -> dict[str, np.ndarray]:
    sampler = qmc.LatinHypercube(d=len(REP_PARAM_NAMES), seed=seed)
    u = sampler.random(n)
    p = {}
    for i, name in enumerate(REP_PARAM_NAMES):
        nom = REP_NOMINAL[name]
        lo = nom * (1.0 - rel_halfwidth)
        hi = nom * (1.0 + rel_halfwidth)
        if name == "k_rep_adv":
            lo, hi = 1.12, 1.52  # Always > 1.0 (deletion POLG advantage range)
        elif name == "pink1_eff":
            lo, hi = 0.65, 1.15
        elif name == "s_bias_err":
            lo, hi = 0.72, 1.28
        p[name] = lo + u[:, i] * (hi - lo)
    return p


def run_sobol_gsa(n_base: int = 1024) -> dict:
    """Saltelli-Sobol GSA (n_base * (k + 2) = 10,240 evals) on Aged-Patient Biphasic f_M_end."""
    k = len(REP_PARAM_NAMES)
    rng = np.random.default_rng(SEED + 801)
    A = rng.random((n_base, k))
    B = rng.random((n_base, k))

    def u_to_p(U: np.ndarray) -> dict[str, np.ndarray]:
        p = {}
        for i, name in enumerate(REP_PARAM_NAMES):
            nom = REP_NOMINAL[name]
            lo, hi = nom * 0.70, nom * 1.30
            if name == "k_rep_adv":
                lo, hi = 1.12, 1.52
            elif name == "pink1_eff":
                lo, hi = 0.65, 1.15
            elif name == "s_bias_err":
                lo, hi = 0.72, 1.28
            p[name] = lo + U[:, i] * (hi - lo)
        return p

    w0, d0, m0 = 0.92 * 0.40, 0.92 * 0.20, 0.92 * 0.40
    yA = simulate_rep3_batch(u_to_p(A), "Biphasic_Unmasked_Galactose_PDK1", w0, d0, m0)["f_M_end"]
    yB = simulate_rep3_batch(u_to_p(B), "Biphasic_Unmasked_Galactose_PDK1", w0, d0, m0)["f_M_end"]
    y_all = np.concatenate([yA, yB])
    var_y = float(np.var(y_all))

    s1, st = {}, {}
    for i, name in enumerate(REP_PARAM_NAMES):
        AB_i = A.copy()
        AB_i[:, i] = B[:, i]
        yAB_i = simulate_rep3_batch(u_to_p(AB_i), "Biphasic_Unmasked_Galactose_PDK1", w0, d0, m0)["f_M_end"]
        s1_val = float(np.mean(yB * (yAB_i - yA)) / max(1e-12, var_y))
        st_val = float(0.5 * np.mean((yA - yAB_i) ** 2) / max(1e-12, var_y))
        s1[name] = round(max(0.0, s1_val), 4)
        st[name] = round(max(0.0, st_val), 4)

    ranked = sorted(REP_PARAM_NAMES, key=lambda x: st[x], reverse=True)
    return {
        "n_base": n_base,
        "total_evaluations": n_base * (k + 2),
        "output_variance": round(var_y, 6),
        "S1": s1,
        "ST": st,
        "ranked_parameters": ranked,
    }


def compute_scheduled_surgery_preop_clocks() -> dict:
    """Compute exact protective fraction vs. dosing lead time (t_lead in hours relative to clamp release t=0).

    Negative t_lead = dosed BEFORE scheduled surgical unclamping (e.g., -18 h pre-op).
    Positive t_lead = dosed AFTER unscheduled emergency reperfusion onset (e.g., +2 h in ER).

    Uses pharmacokinetic Bateman onset/decay curves anchored to biological clocks:
      - Clock 1 (Malonate/MitoSNO): t_on = 0.05 h (3 min), t_half = 1.5 h
      - Clock 2 (TAT-SOD2 / Banked Donor Mito): t_on = 0.25 h (15 min), t_half = 48 h
      - Clock 3 (Transient SOD2+UCP2 mRNA-LNP): t_on = 4.5 h (translation+TOM/TIM import), t_half = 72 h
      - Clock 3 + tau_lag=6 h (Compressed BCR priming) vs tau_lag=13.5 h (Uncompressed fasting)
    Evaluated against TWO injury waves:
      - Wave A (Minute-0 Clamp-Release RET-ROS Burst, 0-1 h post-reperfusion)
      - Wave B (Delayed 6-48 h Wave: Blast ear, Burn stasis, Birth HIE)
    """
    lead_hours = [-48, -36, -24, -18, -12, -6, -2, -0.5, 0.0, 1.0, 2.0, 6.0, 12.0]

    def bateman_activity(dt_since_dose_h: float, tau_on: float, tau_off: float) -> float:
        if dt_since_dose_h <= 0.0:
            return 0.0
        # Normalized so peak = 1.0
        t_peak = (tau_on * tau_off / (tau_off - tau_on)) * math.log(tau_off / tau_on)
        peak = math.exp(-t_peak / tau_off) - math.exp(-t_peak / tau_on)
        val = math.exp(-dt_since_dose_h / tau_off) - math.exp(-dt_since_dose_h / tau_on)
        return max(0.0, min(1.0, val / max(1e-9, peak)))

    def wave_coverage(t_lead_h: float, tau_on: float, tau_off: float, wave_start: float, wave_end: float, tau_lag: float = 0.0) -> float:
        # Mean active fraction across the injury wave [wave_start, wave_end]
        grid = np.linspace(wave_start, wave_end, 40)
        vals = []
        for t_ev in grid:
            dt_dose = t_ev - t_lead_h
            if dt_dose <= 0.0:
                vals.append(0.0)
            else:
                act = bateman_activity(dt_dose, tau_on, tau_off)
                if tau_lag > 0.0:
                    lag_gate = 1.0 / (1.0 + math.exp(-(dt_dose - tau_lag) / 1.8))
                    act *= lag_gate
                vals.append(act)
        return round(float(np.mean(vals)), 4)

    clocks = {
        "Clock1_SmallMol_Malonate": {"tau_on": 0.06, "tau_off": 2.2, "tau_lag": 0.0, "label": "Clock 1: Malonate / MitoSNO (30s-5m)"},
        "Clock2_BankedMito_TAT_SOD2": {"tau_on": 0.25, "tau_off": 68.0, "tau_lag": 0.0, "label": "Clock 2: Banked Donor Mito / TAT-SOD2 (5-30m)"},
        "Clock3_mRNA_LNP_SOD2_UCP2": {"tau_on": 4.8, "tau_off": 84.0, "tau_lag": 0.0, "label": "Clock 3: Transient mRNA-LNP SOD2+UCP2 (2-6h)"},
        "Clock3_BCR_Compressed_6h": {"tau_on": 3.5, "tau_off": 72.0, "tau_lag": 6.0, "label": "Clock 3 + Mito_Papers Compressed Lag (tau=6h)"},
        "Clock3_BCR_Standard_13_5h": {"tau_on": 3.5, "tau_off": 72.0, "tau_lag": 13.5, "label": "Clock 3 + Standard Fasting Lag (tau=13.5h)"},
    }

    minute0_curves = {}
    delayed_curves = {}
    summary = {}
    for cname, cfg in clocks.items():
        m0 = [wave_coverage(tl, cfg["tau_on"], cfg["tau_off"], 0.0, 1.0, cfg["tau_lag"]) for tl in lead_hours]
        dw = [wave_coverage(tl, cfg["tau_on"], cfg["tau_off"], 12.0, 48.0, cfg["tau_lag"]) for tl in lead_hours]
        minute0_curves[cname] = m0
        delayed_curves[cname] = dw
        best_idx = int(np.argmax(m0))
        idx_0h = lead_hours.index(0.0)
        idx_p2h = lead_hours.index(2.0)
        summary[cname] = {
            "label": cfg["label"],
            "optimal_preop_lead_h": lead_hours[best_idx],
            "peak_minute0_coverage": m0[best_idx],
            "coverage_if_dosed_at_0h_surprise": m0[idx_0h],
            "delayed_wave_coverage_if_dosed_at_plus2h_ER": dw[idx_p2h],
        }

    return {
        "lead_hours": lead_hours,
        "minute0_clamp_release_wave_0_to_1h": minute0_curves,
        "delayed_injury_wave_12_to_48h": delayed_curves,
        "clock_summary": summary,
    }


def summarize_array(arr: np.ndarray) -> dict:
    return {
        "median": round(float(np.median(arr)), 4),
        "p05": round(float(np.percentile(arr, 5)), 4),
        "p95": round(float(np.percentile(arr, 95)), 4),
        "mean": round(float(np.mean(arr)), 4),
    }


def main() -> None:
    print("=== Stage 8A: Full Mito_Papers 3-State Replication, Donor Reversibility & Surgery Clocks ===")
    conv = verify_solver_convergence()
    print(f"[1/4] Solver convergence (RK4 vs DOP853/Radau): max |Delta| = {conv['max_abs_delta']:.3e} (passed={conv['passed']})")
    assert conv["passed"], f"Solver convergence failed: {conv}"

    # 1. N = 10,000 LHS across 5 wet-lab protocols (Aged cell: 40% W_host, 20% D_aged, 40% M_edit)
    N_LHS = 10_000
    p_lhs = sample_lhs_params(N_LHS, SEED)
    w0_aged, d0_aged, m0_aged = 0.92 * 0.40, 0.92 * 0.20, 0.92 * 0.40

    regimens = [
        "HighGlucose_25mM_TrapB_Failure",
        "Unphased_Biogenesis_TrapA",
        "Biphasic_Unmasked_Galactose_PDK1",
        "TwoCycle_Pulsed_Unmasked_Safe",
        "Adversarial_Null_Placebo",
    ]
    trap_results = {}
    for reg in regimens:
        res = simulate_rep3_batch(p_lhs, reg, w0_aged, d0_aged, m0_aged, t_end=14.0, dt=0.02)
        trap_results[reg] = {
            "f_M_end": summarize_array(res["f_M_end"]),
            "f_D_end": summarize_array(res["f_D_end"]),
            "f_W_end": summarize_array(res["f_W_end"]),
            "tot_end": summarize_array(res["tot_end"]),
            "tot_nadir": summarize_array(res["tot_nadir"]),
            "prob_fM_gt_80pct": round(float(np.mean(res["f_M_end"] >= 0.80)), 4),
            "prob_fD_lt_5pct": round(float(np.mean(res["f_D_end"] <= 0.05)), 4),
        }
        print(
            f"  {reg:34s} | f_M={trap_results[reg]['f_M_end']['median']*100:5.1f}% "
            f"| f_D={trap_results[reg]['f_D_end']['median']*100:5.1f}% "
            f"| nadir={trap_results[reg]['tot_nadir']['median']:.3f}"
        )

    # 2D Aged-Patient Deletion Grid: D_aged(0) in {0%, 15%, 30%, 45%}
    aged_grid = []
    p_grid = sample_lhs_params(2000, SEED + 77)
    for d_frac in [0.00, 0.15, 0.30, 0.45]:
        rem = 1.0 - d_frac
        w0_g = 0.92 * (0.50 * rem)
        m0_g = 0.92 * (0.50 * rem)
        d0_g = max(1e-6, 0.92 * d_frac)
        row = {"d0_frac": d_frac, "m0_frac": round(0.50 * rem, 3)}
        for reg in ["HighGlucose_25mM_TrapB_Failure", "Unphased_Biogenesis_TrapA", "Biphasic_Unmasked_Galactose_PDK1", "TwoCycle_Pulsed_Unmasked_Safe"]:
            r = simulate_rep3_batch(p_grid, reg, w0_g, d0_g, m0_g, t_end=14.0, dt=0.02)
            row[reg] = {
                "f_M_med": round(float(np.median(r["f_M_end"])), 4),
                "f_D_med": round(float(np.median(r["f_D_end"])), 4),
                "nadir_med": round(float(np.median(r["tot_nadir"])), 4),
            }
        aged_grid.append(row)

    # 2. Saltelli-Sobol GSA (10,240 evals)
    print("[2/4] Running 10,240-eval Saltelli-Sobol GSA on 3-state replication...")
    sobol = run_sobol_gsa(n_base=1024)
    print("  Top 4 Sobol S_T drivers:", [(k, sobol["ST"][k]) for k in sobol["ranked_parameters"][:4]])

    # 3. Track 8A-2: Context-Sensitive Donor-Mitochondria Reversibility vs Self-Expansion (60 days, M0 = 10%)
    print("[3/4] Running 60-day Donor-Mitochondria Context-Sensitive Reversibility (N=5,000 LHS)...")
    p_don = sample_lhs_params(5000, SEED + 202)
    donor_contexts = {
        "Healthy_Normoxic_Organ": {
            "regimen": "Donor_Healthy_Normoxic_60d",
            "w0": 0.92 * 0.88, "d0": 0.92 * 0.02, "m0": 0.92 * 0.10,
            "k_exo_clear": 0.025,  # Mild exogenous organelle turnover/dilution in healthy normoxia
            "label": "Healthy Normoxic Organ (14 kPa, 2% D)",
        },
        "PostTransplant_24hHypox_ThenNormox": {
            "regimen": "Donor_PostTransplant_24hHypox_ThenNormox",
            "w0": 0.92 * 0.85, "d0": 0.92 * 0.05, "m0": 0.92 * 0.10,
            "k_exo_clear": 0.025,
            "label": "Post-Transplant (24h Stress -> 59d Healthy)",
        },
        "Chronically_Hypoxic_Aged_Tissue_BCR": {
            "regimen": "Donor_ChronicHypox_Aged_CyclicBCR",
            "w0": 0.92 * 0.60, "d0": 0.92 * 0.30, "m0": 0.92 * 0.10,
            "k_exo_clear": 0.025,  # SAME mild exogenous turnover, yet overcomes via hypoxic PINK1 selectivity!
            "label": "Hypoxic Aged Tissue + Cyclic BCR (1.35 kPa, 30% D)",
        },
    }
    donor_results = {}
    for cname, cfg in donor_contexts.items():
        res = simulate_rep3_batch(
            p_don, cfg["regimen"], cfg["w0"], cfg["d0"], cfg["m0"],
            t_end=60.0, dt=0.05, k_exo_clear=cfg["k_exo_clear"],
        )
        trace = simulate_rep3_trace(
            cfg["regimen"], cfg["w0"], cfg["d0"], cfg["m0"],
            t_end=60.0, dt=0.05, sample_stride=20, k_exo_clear=cfg["k_exo_clear"],
        )
        fM_med = float(np.median(res["f_M_end"]))
        regime_tag = "AUTO_REVERSIBLE_WASHOUT" if fM_med < 0.06 else ("SELF_TARGETING_ENGRAFTMENT" if fM_med > 0.45 else "STABLE_PERSISTENCE")
        donor_results[cname] = {
            "label": cfg["label"],
            "m0_frac": 0.10,
            "f_M_day60": summarize_array(res["f_M_end"]),
            "f_D_day60": summarize_array(res["f_D_end"]),
            "tot_day60": summarize_array(res["tot_end"]),
            "tot_nadir": summarize_array(res["tot_nadir"]),
            "classification": regime_tag,
            "trace_60d": trace,
        }
        print(
            f"  {cname:38s} | Day-60 M_donor={fM_med*100:5.1f}% (from 10%) "
            f"| Day-60 D={donor_results[cname]['f_D_day60']['median']*100:4.1f}% -> {regime_tag}"
        )

    # 4. Track 8A-3: Scheduled Surgery Pre-Op Lead-Time Clocks
    print("[4/4] Computing 4-Clock Scheduled-Surgery Pre-Op Lead-Time Curves...")
    surgery_clocks = compute_scheduled_surgery_preop_clocks()
    for k, s in surgery_clocks["clock_summary"].items():
        print(
            f"  {k:30s} | Opt Lead={s['optimal_preop_lead_h']:+5.1f}h (cov={s['peak_minute0_coverage']:.3f}) "
            f"| Surprise 0h={s['coverage_if_dosed_at_0h_surprise']:.3f} "
            f"| Delayed +2h ER={s['delayed_wave_coverage_if_dosed_at_plus2h_ER']:.3f}"
        )

    payload = {
        "stage": "8A_Full_MitoPapers_Replication_Engraftment_and_Surgery_Clocks",
        "epistemic_tier": "E3_MECHANISTIC_SIM",
        "solver_convergence": conv,
        "trap_a_b_wetlab_14d_lhs10k": trap_results,
        "aged_patient_deletion_grid": aged_grid,
        "sobol_gsa_10240": sobol,
        "donor_mitochondria_60d_reversibility": donor_results,
        "scheduled_surgery_preop_clocks": surgery_clocks,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    print(f"Saved -> {OUT_JSON}")


if __name__ == "__main__":
    main()
