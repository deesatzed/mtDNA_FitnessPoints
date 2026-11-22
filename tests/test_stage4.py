#!/usr/bin/env python3
"""Stage 4 Verification Suite: Solver convergence, trade-off check, positive epistasis UQ, Sobol GSA & DB."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data" / "mitofitness.db"
S4_JSON = ROOT / "data" / "stage4_epistasis_summary.json"
FIG_PATH = ROOT / "figures" / "stage4_epistasis_synergy_and_cliffs.png"


def run_tests() -> None:
    assert DB_PATH.exists() and S4_JSON.exists() and FIG_PATH.exists()
    conn = sqlite3.connect(DB_PATH)
    s4 = json.loads(S4_JSON.read_text(encoding="utf-8"))

    # 1. Solver convergence RK4 vs DOP853/Radau (< 1e-3)
    sc = s4["solver_convergence"]
    assert sc["max_abs_diff"] < 1e-3, f"Solver convergence error too large: {sc}"
    print(f"[PASS] Solver convergence verified (RK4 vs DOP853/Radau max |Delta| = {sc['max_abs_diff']}).")

    # 2. Realistic thermodynamic trade-off check (unbuffered COX3 W116R improves hypoxic ATP nadir BUT raises reperfusion EP!)
    wt = s4["standalone_baselines"]["WT_Human_rCRS"]
    goose = s4["standalone_baselines"]["drivers_alone"]["COX3_W116R"]
    assert goose["atp_nadir_median"] > wt["atp_nadir_median"], "Expected COX3 W116R to improve hypoxic ATP nadir"
    assert goose["ep_peak_rep_median"] > wt["ep_peak_rep_median"], "Expected COX3 W116R alone to exhibit reperfusion EP trade-off"
    print(f"[PASS] Non-monotonic trade-off verified: COX3 W116R ATP nadir {wt['atp_nadir_median']}->{goose['atp_nadir_median']} (gain) while EP peak {wt['ep_peak_rep_median']}->{goose['ep_peak_rep_median']} (cost).")

    # 3. Positive epistasis (S > 1.0 and E_add > 0) across top pairs under N=10,000 LHS
    top = s4["top_epistatic_pair"]
    assert s4["n_lhs"] == 10000
    assert top["synergy_index_median"] > 1.25 and top["epistatic_excess_median"] > 0.05 and top["prob_positive_epistasis_pct"] >= 80.0
    print(f"[PASS] Positive epistasis verified on N=10,000 LHS: {top['driver']} + {top['buffer']} -> S={top['synergy_index_median']}x, E_add=+{top['epistatic_excess_median']} ({top['prob_positive_epistasis_pct']}% LHS draws > 0).")

    # 4. Saltelli-Sobol GSA valid and non-degenerate
    sob = s4["sobol_gsa"]
    assert sob["total_evals"] == 10240 and sum(sob["sobol_ST"].values()) > 0.5
    print(f"[PASS] Saltelli-Sobol GSA verified ({sob['total_evals']} evals, sum S_T={sum(sob['sobol_ST'].values()):.3f}).")

    # 5. SQLite epistatic_synergies populated with all 25 Driver x Buffer pairs
    n_syn = conn.execute("SELECT COUNT(*) FROM epistatic_synergies").fetchone()[0]
    n_cliffs = conn.execute("SELECT COUNT(*) FROM fitness_payoff_profiles WHERE lower_survival_cliff IS NOT NULL").fetchone()[0]
    assert n_syn == 25 and n_cliffs >= 5
    print(f"[PASS] SQLite epistatic_synergies populated ({n_syn}/25 pairs) & {n_cliffs} payoff profile cliff rows updated.")

    conn.close()
    print("\nALL STAGE 4 SANITY CHECKS PASSED (100%).")


if __name__ == "__main__":
    run_tests()
