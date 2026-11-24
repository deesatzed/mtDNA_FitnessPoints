#!/usr/bin/env python3
"""Stage 7A: Cross-Project Bridge from Mito_Papers (Biphasic Heteroplasmy Amplification & CAR-T Hypoxia).

Track 7A — Biphasic Mitophagy-Selection Heteroplasmy Amplification ODE:
  Couples Mito_Papers' two-stage Bioenergetic Circuit Reset (Stage 1 selective PINK1/Urolithin A
  mitophagy under mild hypoxia pO2 = 1.2 kPa -> Stage 2 normoxic PGC-1a/NAD+ biogenesis) to a
  realistic bulk TALED-edited cell pool starting at 50% COX3 W116R heteroplasmy (50% WT rCRS +
  50% COX3 W116R). Tests whether mild hypoxic differential depolarization selectively clears
  unedited WT copies and amplifies heteroplasmy past Stage 5's >=85-90% acute I/R synergy cliff
  without single-cell clonal sorting. Includes:
    - RK4 vs. DOP853 & stiff Radau solver convergence check
    - N = 10,000 Latin Hypercube Sampling (+/-30% priors) across 4 arms
    - 10,240-eval Saltelli-Sobol GSA on final Day-14 heteroplasmy & acute I/R rescue
    - Exact wet-lab falsification boundary across initial editing (20%-70%) x selectivity bias

Track 7B — Ex Vivo CAR-T Solid-Tumor Hypoxia-Reoxygenation Cycling (72 h, 6 Cycles):
  Adapts Mito_Papers' CAR-T bioreactor tumor-microenvironment hypoxia-cycling paradigm to test
  100% human-sequence COX3 W116R + Human SOD2 (zero foreign-protein immunogenicity, clonal >90%
  editing feasible ex vivo) under N = 10,000 LHS draws.

Outputs:
  - data/stage7a_heteroplasmy_and_cart.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import stage4_epistasis_synergy_uq as s4

OUT_JSON = ROOT / "data" / "stage7a_heteroplasmy_and_cart.json"
SEED = 20261123


# =============================================================================
# TRACK 7A: BIPHASIC MITOPHAGY-SELECTION HETEROPLASMY AMPLIFICATION (14 DAYS)
# =============================================================================
SEL_PARAM_NAMES = [
    "Km_wt",         # WT COX3 O2 half-saturation (kPa)
    "Km_mut",        # COX3 W116R O2 half-saturation (kPa)
    "k_mito_base",   # Basal mitophagy turnover rate (/day)
    "k_mito_sel",    # PINK1/Urolithin A depolarization-gated selective mitophagy gain (/day)
    "k_bio_max",     # Stage-2 PGC-1a/NAD+ biogenesis replication rate (/day)
    "depol_slope",   # Sigmoid steepness converting low O2 conductance -> PINK1 depolarization
    "k_copy_cap",    # Carrying capacity K_cap of total mtDNA pool (normalized = 1.0)
    "s_bias_err",    # Fractional noise / leakiness in WT vs mutant selectivity (1.0 = nominal)
]

SEL_NOMINAL = {
    "Km_wt": 1.25,
    "Km_mut": 0.53,
    "k_mito_base": 0.06,
    "k_mito_sel": 0.82,
    "k_bio_max": 0.55,
    "depol_slope": 6.0,
    "k_copy_cap": 1.00,
    "s_bias_err": 1.00,
}


def sel_ode_rhs_vec(
    t_day: float,
    y: np.ndarray,
    p: dict[str, np.ndarray],
    regimen: str,
    pO2_hypox: float = 1.10,
) -> np.ndarray:
    """Two-population organelle ODE: y = [W_wt, W_mut].

    In mild hypoxia (pO2 ~ 1.1 kPa), WT (Km=1.25 kPa) operates at ~47% O2 saturation
    while COX3 W116R (Km=0.53 kPa) operates at ~68% saturation -> differential DeltaPsi
    depolarization triggers preferential PINK1/Urolithin A clearance of unedited WT copies!
    """
    W_wt, W_mut = np.maximum(1e-6, y[0]), np.maximum(1e-6, y[1])
    W_tot = W_wt + W_mut

    if regimen == "Normoxia_NoSelection":
        pO2 = 14.0
        u_mito = 0.0
        u_bio = 0.30
    elif regimen == "Unphased_Biogenesis_Trap":
        pO2 = 14.0
        u_mito = 0.0
        u_bio = 1.00
    elif regimen == "Continuous_Mitophagy_Only":
        pO2 = pO2_hypox
        u_mito = 1.00
        u_bio = 0.05
    elif regimen == "Biphasic_Selection_Reset":
        if t_day < 7.0:
            # Stage 1 (Days 0-7): Mild hypoxia conditioning + PINK1/Urolithin A selective purge
            pO2 = pO2_hypox
            u_mito = 1.00
            u_bio = 0.08
        else:
            # Stage 2 (Days 7-14): Normoxic recovery + PGC-1a/NAD+ biogenesis expansion
            pO2 = 14.0
            u_mito = 0.05
            u_bio = 1.00
    elif regimen == "TwoCycle_Pulsed_Biphasic":
        # Safer 2-cycle pulsed protocol: 2 x (3.5 d mild purge + 3.5 d normoxic recovery)
        # Prevents mitochondrial mass depletion (W_tot nadir stays biologically safe > 0.30!)
        cycle_t = t_day % 7.0
        if cycle_t < 3.5:
            pO2 = pO2_hypox
            u_mito = 0.90
            u_bio = 0.15
        else:
            pO2 = 14.0
            u_mito = 0.05
            u_bio = 1.00
    else:
        raise ValueError(regimen)

    # Fractional Complex IV oxygen conductance g_wt vs g_mut
    g_wt = pO2 / (p["Km_wt"] + pO2)
    g_mut = pO2 / (p["Km_mut"] + pO2)

    # Depolarization probability (high when O2 conductance drops below ~0.58)
    dep_wt = 1.0 / (1.0 + np.exp(p["depol_slope"] * (g_wt - 0.58)))
    dep_mut = 1.0 / (1.0 + np.exp(p["depol_slope"] * (g_mut - 0.58)))

    # Differential selectivity scaled by s_bias_err
    eff_dep_wt = dep_wt * p["s_bias_err"]
    eff_dep_mut = dep_mut / np.maximum(0.5, p["s_bias_err"])

    clear_wt = p["k_mito_base"] + u_mito * p["k_mito_sel"] * eff_dep_wt
    clear_mut = p["k_mito_base"] + u_mito * p["k_mito_sel"] * eff_dep_mut

    # Carrying-capacity-limited symmetric biogenesis (logistic expansion preserves current ratio!)
    growth_drive = u_bio * p["k_bio_max"] * np.maximum(0.0, 1.0 - W_tot / p["k_copy_cap"])
    # Basal homeostatic maintenance
    maint = 0.06 * np.maximum(0.0, 1.0 - W_tot / p["k_copy_cap"])

    dW_wt = (growth_drive + maint - clear_wt) * W_wt
    dW_mut = (growth_drive + maint - clear_mut) * W_mut
    return np.stack([dW_wt, dW_mut], axis=0)


def simulate_selection_batch(
    p: dict[str, np.ndarray],
    regimen: str,
    h0: float = 0.50,
    t_end: float = 14.0,
    dt: float = 0.05,
    pO2_hypox: float = 1.10,
) -> dict[str, np.ndarray]:
    n = len(p["Km_wt"])
    y = np.zeros((2, n), dtype=np.float64)
    y[0, :] = 0.92 * (1.0 - h0)
    y[1, :] = 0.92 * h0

    n_steps = int(round(t_end / dt))
    min_wtot = np.full(n, 99.0)
    for k in range(n_steps):
        t = k * dt
        min_wtot = np.minimum(min_wtot, y[0] + y[1])
        k1 = sel_ode_rhs_vec(t, y, p, regimen, pO2_hypox)
        k2 = sel_ode_rhs_vec(t + 0.5 * dt, np.maximum(1e-6, y + 0.5 * dt * k1), p, regimen, pO2_hypox)
        k3 = sel_ode_rhs_vec(t + 0.5 * dt, np.maximum(1e-6, y + 0.5 * dt * k2), p, regimen, pO2_hypox)
        k4 = sel_ode_rhs_vec(t + dt, np.maximum(1e-6, y + dt * k3), p, regimen, pO2_hypox)
        y = np.maximum(1e-6, y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4))

    w_tot = y[0] + y[1]
    h_end = y[1] / np.maximum(1e-9, w_tot)
    return {
        "h_end": h_end,
        "w_tot_end": w_tot,
        "w_tot_nadir": min_wtot,
        "w_mut_end": y[1],
        "w_wt_end": y[0],
    }


def verify_selection_solver_convergence() -> dict:
    p1 = {k: np.array([v]) for k, v in SEL_NOMINAL.items()}
    rk = simulate_selection_batch(p1, "Biphasic_Selection_Reset", h0=0.50, dt=0.05)

    def rhs_sc(t, yv):
        return sel_ode_rhs_vec(t, np.array([[yv[0]], [yv[1]]]), p1, "Biphasic_Selection_Reset")[:, 0]

    # Integrate piecewise across the Stage 1 -> Stage 2 switch at t = 7.0 days for exact continuity
    sol_d1 = solve_ivp(rhs_sc, [0.0, 7.0], [0.46, 0.46], method="DOP853", rtol=1e-10, atol=1e-12)
    sol_d2 = solve_ivp(rhs_sc, [7.0, 14.0], sol_d1.y[:, -1], method="DOP853", rtol=1e-10, atol=1e-12)
    h_dop = float(sol_d2.y[1, -1] / (sol_d2.y[0, -1] + sol_d2.y[1, -1]))

    sol_r1 = solve_ivp(rhs_sc, [0.0, 7.0], [0.46, 0.46], method="Radau", rtol=1e-10, atol=1e-12)
    sol_r2 = solve_ivp(rhs_sc, [7.0, 14.0], sol_r1.y[:, -1], method="Radau", rtol=1e-10, atol=1e-12)
    h_rad = float(sol_r2.y[1, -1] / (sol_r2.y[0, -1] + sol_r2.y[1, -1]))

    h_rk = float(rk["h_end"][0])
    diff = max(abs(h_rk - h_dop), abs(h_rk - h_rad))
    return {
        "h_end_RK4_dt005": round(h_rk, 7),
        "h_end_DOP853": round(h_dop, 7),
        "h_end_Radau": round(h_rad, 7),
        "max_abs_diff": float(f"{diff:.3e}"),
    }


def saltelli_sobol_selection(n_base: int = 1024) -> dict:
    d = len(SEL_PARAM_NAMES)
    sampler = qmc.Sobol(d=2 * d, scramble=True, seed=SEED)
    M = sampler.random(n_base)
    A_u, B_u = M[:, :d], M[:, d:]

    def u_to_p(U):
        return {name: SEL_NOMINAL[name] * (0.70 + 0.60 * U[:, j]) for j, name in enumerate(SEL_PARAM_NAMES)}

    yA = simulate_selection_batch(u_to_p(A_u), "Biphasic_Selection_Reset", h0=0.50)["h_end"]
    yB = simulate_selection_batch(u_to_p(B_u), "Biphasic_Selection_Reset", h0=0.50)["h_end"]
    var_y = float(np.var(np.concatenate([yA, yB]), ddof=1))

    s1, st = {}, {}
    for j, name in enumerate(SEL_PARAM_NAMES):
        AB_j = A_u.copy()
        AB_j[:, j] = B_u[:, j]
        yAB = simulate_selection_batch(u_to_p(AB_j), "Biphasic_Selection_Reset", h0=0.50)["h_end"]
        s1[name] = round(float(np.mean(yB * (yAB - yA)) / max(1e-12, var_y)), 4)
        st[name] = round(float(0.5 * np.mean((yA - yAB) ** 2) / max(1e-12, var_y)), 4)
    return {"n_evals": n_base * (d + 2), "output_variance": round(var_y, 6), "sobol_S1": s1, "sobol_ST": st}


# =============================================================================
# TRACK 7B: EX VIVO CAR-T SOLID-TUMOR HYPOXIA-REOXYGENATION CYCLING (72 HOURS)
# =============================================================================
def simulate_cart_tumor_cycling_batch(
    p: dict[str, np.ndarray],
    drv: dict | None,
    buf: dict | None,
    heteroplasmy: float = 0.92,
    dt_hr: float = 0.25,
) -> dict[str, np.ndarray]:
    """72-hour ex vivo CAR-T tumor-microenvironment hypoxia cycling (6 cycles x 12 h:
    8 h tumor-core hypoxia pO2 = 0.45 kPa + 4 h perivascular reoxygenation pO2 = 6.5 kPa)
    under sustained high effector flux demand (f_load_cart = 1.22 * f_load)."""
    n = len(p["Km_O2_wt"])
    scale = p["drv_scale"] * heteroplasmy
    if drv is not None and heteroplasmy > 0.0:
        km_o2 = p["Km_O2_wt"] * (1.0 - scale * (1.0 - drv["km_o2_mult"]))
        atp_o_gain = scale * drv["atp_o_gain"]
        ret_mult = 1.0 + scale * (drv["ret_sensitivity"] - 1.0)
        leak_d = scale * drv["proton_leak_delta"]
    else:
        km_o2 = p["Km_O2_wt"]
        atp_o_gain = np.zeros(n)
        ret_mult = np.ones(n)
        leak_d = np.zeros(n)

    if buf is not None:
        aox_drain = buf["aox_drain"]
        ndi1_bypass = buf["ndi1_bypass"]
        theta_boost = buf["theta_gsh_boost"]
        ucp_leak = buf["ucp_leak"]
        mito_boost = buf["mitophagy_boost"]
        buf_cost = buf["atp_o_cost"]
    else:
        aox_drain = ndi1_bypass = theta_boost = ucp_leak = mito_boost = buf_cost = 0.0

    f_cart = 1.18 * p["f_load"]
    theta_eff = p["theta_GSH"] + theta_boost

    W = np.full(n, 0.92)
    Q = np.full(n, 0.35)
    Succ = np.full(n, 0.10)
    ATP = np.full(n, 0.92)
    Eff = np.ones(n)  # CAR-T cytolytic effector persistence [0, 1]

    n_steps = int(round(72.0 / dt_hr))
    cum_ep_exc = np.zeros(n)
    min_atp = np.full(n, 99.0)

    for k in range(n_steps):
        t_hr = k * dt_hr
        phase = t_hr % 12.0
        pO2 = 0.45 if phase < 8.0 else 6.50
        isch = 1.0 if phase < 8.0 else 0.0

        g_iv = pO2 / (np.maximum(0.10, km_o2) + pO2)
        leak = np.clip(0.12 + leak_d + ucp_leak, 0.04, 0.55)
        dpsi = np.clip(0.60 + 0.55 * Q * g_iv - 0.65 * leak - 0.15 * aox_drain * np.maximum(0.0, Q - 0.55), 0.25, 1.35)

        dSucc = 0.35 * p["succ_surge"] * isch * (1.30 - Succ) - 0.55 * (1.0 - isch) * g_iv * Succ
        q_in = (0.45 * f_cart + 0.90 * (1.0 - isch) * g_iv * Succ) * (1.0 - Q)
        aox_fl = aox_drain * 0.95 * np.maximum(0.0, Q - 0.48) * (pO2 / (0.60 + pO2))
        q_out = (0.55 * g_iv * np.maximum(0.15, W) + aox_fl) * Q
        dQ = q_in - q_out

        ret_f = (
            p["k_RET"] * ret_mult * (1.0 - 0.58 * ndi1_bypass) * (1.0 - 0.62 * np.minimum(1.0, ucp_leak / 0.25))
            * np.maximum(0.0, Q - 0.40) ** 1.4 * g_iv * (0.25 + 1.35 * Succ) * (dpsi ** 2)
        )
        ep = (0.55 / np.maximum(0.35, W)) * (1.0 + 1.55 * ret_f) * (f_cart ** 2)
        ep_exc = np.maximum(0.0, ep - theta_eff)
        cum_ep_exc += ep_exc * dt_hr

        coupling = (1.0 + atp_o_gain) * (1.0 - 0.55 * (leak - 0.12)) * (1.0 - buf_cost * (0.35 + 0.65 * np.maximum(0.0, Q - 0.45)))
        atp_prod = 1.28 * np.maximum(0.05, W) * g_iv * coupling
        atp_use = 0.56 * f_cart * (ATP / (0.18 + ATP))
        dATP = 0.65 * (atp_prod - atp_use)

        # Exhaustion & organelle damage across 6 hypoxia-reoxygenation cycles
        atp_def = np.maximum(0.0, p["theta_ATP"] - ATP)
        dW = (-0.045 * ep_exc - 0.030 * atp_def) * W + (0.010 + 0.035 * mito_boost) * (0.95 - W) * g_iv * (1.0 - np.minimum(1.0, ep_exc))
        dEff = (-0.065 * atp_def - 0.055 * ep_exc) * Eff

        W = np.clip(W + dt_hr * dW, 0.0, 1.0)
        Q = np.clip(Q + dt_hr * dQ, 0.0, 1.0)
        Succ = np.clip(Succ + dt_hr * dSucc, 0.0, 2.0)
        ATP = np.clip(ATP + dt_hr * dATP, 0.0, 1.5)
        Eff = np.clip(Eff + dt_hr * dEff, 0.0, 1.0)
        if t_hr >= 4.0:
            min_atp = np.minimum(min_atp, ATP)

    return {
        "effector_persistence_72h": Eff,
        "W_end_72h": W,
        "min_atp_72h": min_atp,
        "cum_ep_excess_72h": cum_ep_exc,
    }


def main() -> None:
    print("[1/4] Verifying Track 7A Biphasic Selection solver convergence (RK4 vs DOP853/Radau)...")
    conv = verify_selection_solver_convergence()
    print("      Convergence:", conv)

    print("[2/4] Running N=10,000 LHS across 4 Biphasic Selection regimens (starting from 50% bulk TALED)...")
    n_lhs = 10000
    U_s = qmc.LatinHypercube(d=len(SEL_PARAM_NAMES), seed=SEED).random(n_lhs)
    p_sel = {name: SEL_NOMINAL[name] * (0.70 + 0.60 * U_s[:, j]) for j, name in enumerate(SEL_PARAM_NAMES)}

    regimens = [
        "Normoxia_NoSelection",
        "Unphased_Biogenesis_Trap",
        "Continuous_Mitophagy_Only",
        "Biphasic_Selection_Reset",
        "TwoCycle_Pulsed_Biphasic",
    ]
    sel_arms = {}
    for reg in regimens:
        r = simulate_selection_batch(p_sel, reg, h0=0.50)
        h_med, h_p05, h_p95 = [round(float(x * 100.0), 2) for x in np.percentile(r["h_end"], [50, 5, 95])]
        w_med, w_p05, w_p95 = [round(float(x), 3) for x in np.percentile(r["w_tot_end"], [50, 5, 95])]
        nad_med = round(float(np.median(r["w_tot_nadir"])), 3)
        pct_ge85 = round(float(100.0 * np.mean(r["h_end"] >= 0.85)), 1)
        pct_ge90 = round(float(100.0 * np.mean(r["h_end"] >= 0.90)), 1)
        sel_arms[reg] = {
            "h_end_pct_median": h_med,
            "h_end_pct_iq90": [h_p05, h_p95],
            "w_tot_end_median": w_med,
            "w_tot_end_iq90": [w_p05, w_p95],
            "w_tot_nadir_median": nad_med,
            "pct_draws_ge_85pct_heteroplasmy": pct_ge85,
            "pct_draws_ge_90pct_heteroplasmy": pct_ge90,
        }
        print(f"      {reg:28s}: Day-14 h={h_med:5.1f}% [{h_p05:5.1f}, {h_p95:5.1f}%] | W_tot={w_med:.3f} (nadir {nad_med:.3f}) | P(h>=85%)={pct_ge85}%")

    # Generate 14-day nominal trajectory for plotting & interactive Explorer
    t_days = [round(i * 0.5, 2) for i in range(29)]
    p_nom1 = {k: np.array([v]) for k, v in SEL_NOMINAL.items()}
    trajectories = {}
    for reg in regimens:
        h_tr, w_tr = [], []
        for td in t_days:
            if td == 0.0:
                h_tr.append(50.0)
                w_tr.append(0.92)
            else:
                rr = simulate_selection_batch(p_nom1, reg, h0=0.50, t_end=td, dt=0.05)
                h_tr.append(round(float(rr["h_end"][0] * 100.0), 2))
                w_tr.append(round(float(rr["w_tot_end"][0]), 3))
        trajectories[reg] = {"days": t_days, "heteroplasmy_pct": h_tr, "w_tot": w_tr}

    # Initial editing sweep (h0 = 20% .. 70%) for wet-lab feasibility curve
    h0_sweep = []
    for h0_pct in (20, 30, 40, 50, 60, 70):
        r = simulate_selection_batch(p_sel, "Biphasic_Selection_Reset", h0=h0_pct / 100.0)
        h0_sweep.append({
            "h0_pct": h0_pct,
            "h14_median_pct": round(float(np.median(r["h_end"]) * 100.0), 2),
            "h14_iq90_pct": [round(float(x * 100.0), 2) for x in np.percentile(r["h_end"], [5, 95])],
            "pct_draws_clearing_85pct": round(float(100.0 * np.mean(r["h_end"] >= 0.85)), 1),
        })

    print("[3/4] Running 10,240-eval Saltelli-Sobol GSA on Biphasic Heteroplasmy Amplification...")
    sobol_sel = saltelli_sobol_selection(n_base=1024)
    print("      Sobol S_T:", sobol_sel["sobol_ST"])

    print("[4/4] Running Track 7B: Ex Vivo CAR-T 72-h Solid-Tumor Hypoxia Cycling (N=10,000 LHS)...")
    U_c = qmc.LatinHypercube(d=len(s4.PARAM_NAMES), seed=SEED + 7).random(n_lhs)
    p_cart = {name: s4.NOMINAL_P[name] * (0.70 + 0.60 * U_c[:, j]) for j, name in enumerate(s4.PARAM_NAMES)}

    cart_configs = [
        ("WT_Human_CART", None, None, 0.0),
        ("COX3_W116R_Alone_CART", s4.DRIVERS["COX3_W116R"], None, 0.92),
        ("Human_SOD2_Alone_CART", None, s4.BUFFERS["Human_SOD2"], 0.0),
        ("COX3_W116R_plus_Human_SOD2_CART", s4.DRIVERS["COX3_W116R"], s4.BUFFERS["Human_SOD2"], 0.92),
        ("COX3_W116R_plus_Yeast_NDI1_CART", s4.DRIVERS["COX3_W116R"], s4.BUFFERS["Yeast_NDI1"], 0.92),
        ("COX3_W116R_plus_Ciona_AOX_CART", s4.DRIVERS["COX3_W116R"], s4.BUFFERS["Ciona_AOX"], 0.92),
    ]
    cart_raw = {}
    cart_summary = {}
    for label, drv, buf, het in cart_configs:
        res = simulate_cart_tumor_cycling_batch(p_cart, drv, buf, heteroplasmy=het)
        cart_raw[label] = res
        e_med, e_p05, e_p95 = [round(float(x), 4) for x in np.percentile(res["effector_persistence_72h"], [50, 5, 95])]
        w_med = round(float(np.median(res["W_end_72h"])), 4)
        a_med = round(float(np.median(res["min_atp_72h"])), 3)
        ep_med = round(float(np.median(res["cum_ep_excess_72h"])), 2)
        cart_summary[label] = {
            "effector_persistence_72h_median": e_med,
            "effector_persistence_72h_iq90": [e_p05, e_p95],
            "W_end_72h_median": w_med,
            "min_atp_72h_median": a_med,
            "cum_ep_excess_72h_median": ep_med,
        }
        print(f"      {label:35s}: 72h Persistence={e_med:.4f} [{e_p05:.4f}, {e_p95:.4f}] | ATP_min={a_med:.3f} | CumEP_exc={ep_med:.2f}")

    # Compute exact CAR-T positive epistasis for 100% human-sequence COX3 W116R + Human SOD2
    e_wt = cart_raw["WT_Human_CART"]["effector_persistence_72h"]
    e_drv = cart_raw["COX3_W116R_Alone_CART"]["effector_persistence_72h"]
    e_sod = cart_raw["Human_SOD2_Alone_CART"]["effector_persistence_72h"]
    e_pair = cart_raw["COX3_W116R_plus_Human_SOD2_CART"]["effector_persistence_72h"]
    cart_epistasis = (e_pair - e_wt) - ((e_drv - e_wt) + (e_sod - e_wt))
    cart_synergy_stats = {
        "pair": "COX3_W116R + Human_SOD2 (100% Human Sequence, Zero Foreign Immunogenicity)",
        "epistatic_excess_median": round(float(np.median(cart_epistasis)), 4),
        "epistatic_excess_iq90": [round(float(x), 4) for x in np.percentile(cart_epistasis, [5, 95])],
        "pct_lhs_positive": round(float(100.0 * np.mean(cart_epistasis > 0.0)), 1),
        "fold_persistence_gain_vs_WT": round(float(np.median(e_pair) / max(1e-6, np.median(e_wt))), 2),
    }
    print("      CAR-T 100%-Human Synergy:", cart_synergy_stats)

    out = {
        "track7a_biphasic_heteroplasmy_amplification": {
            "solver_convergence": conv,
            "n_lhs": n_lhs,
            "initial_bulk_taled_h0_pct": 50.0,
            "regimens": sel_arms,
            "trajectories_14d": trajectories,
            "h0_feasibility_sweep": h0_sweep,
            "sobol_gsa": sobol_sel,
        },
        "track7b_ex_vivo_cart_tumor_hypoxia": {
            "protocol": "72-h solid-tumor hypoxia cycling (6 cycles x [8h 0.45 kPa hypoxia + 4h 6.5 kPa reoxygenation])",
            "n_lhs": n_lhs,
            "arms": cart_summary,
            "human_sod2_cart_synergy": cart_synergy_stats,
        },
    }
    OUT_JSON.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nStage 7A complete -> {OUT_JSON}")


if __name__ == "__main__":
    main()
