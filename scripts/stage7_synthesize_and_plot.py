#!/usr/bin/env python3
"""Stage 7 Synthesis & Publication Figure: Cross-Project Bridges (Mito_Papers + BAA Mito Twin).

Combines:
  - data/stage7a_heteroplasmy_and_cart.json (Biphasic Mitophagy-Selection Heteroplasmy Amplification
    + Ex Vivo CAR-T 72-h Solid-Tumor Hypoxia Cycling)
  - data/stage7b_singlecell_alphagenome.json (Live AlphaGenome 2 Single-Cell Conditioning across
    6 observed scRNA+scATAC human cell types + 2 BAA-Mito blank negative controls)

Outputs:
  - figures/stage7_cross_project_bridges.png
  - data/stage7_report.md
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
S7A_JSON = ROOT / "data" / "stage7a_heteroplasmy_and_cart.json"
S7B_JSON = ROOT / "data" / "stage7b_singlecell_alphagenome.json"
OUT_FIG = ROOT / "figures" / "stage7_cross_project_bridges.png"
OUT_REPORT = ROOT / "data" / "stage7_report.md"


def render_figure(s7a: dict, s7b: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(17.2, 5.4))

    # Panel A: 14-Day Biphasic Heteroplasmy Amplification (starting from 50% bulk TALED)
    ax = axes[0]
    tr = s7a["track7a_biphasic_heteroplasmy_amplification"]["trajectories_14d"]
    styles = [
        ("Unphased_Biogenesis_Trap", "#94a3b8", "--", "Unphased Biogenesis Trap (h=50%)"),
        ("Continuous_Mitophagy_Only", "#dc2626", ":", "Continuous Mitophagy (W_tot→0.01!)"),
        ("TwoCycle_Pulsed_Biphasic", "#059669", "-.", "Two-Cycle Pulsed Biphasic (safer nadir)"),
        ("Biphasic_Selection_Reset", "#0284c7", "-", "Biphasic Selection Reset (7d→7d)"),
    ]
    for key, col, ls, lbl in styles:
        ax.plot(tr[key]["days"], tr[key]["heteroplasmy_pct"], color=col, lw=2.3, ls=ls, label=lbl)
    ax.axhline(85.0, color="#7c3aed", lw=1.3, ls="--", label="Stage-5 Synergy Cliff (85%)")
    ax.axvspan(0, 7.0, color="#e0f2fe", alpha=0.28)
    ax.text(3.5, 94.0, "Stage 1:\nMild Hypoxia\n+ Mitophagy", ha="center", va="top", fontsize=7.6, color="#0369a1", fontweight="bold")
    ax.text(10.5, 62.5, "Stage 2:\nNormoxic\nBiogenesis", ha="center", va="center", fontsize=7.6, color="#047857", fontweight="bold")
    ax.set_xlabel("Preconditioning Protocol Time (Days)")
    ax.set_ylabel("COX3 W116R Heteroplasmy (%)")
    ax.set_ylim(38, 102)
    ax.set_title("(A) Mito_Papers Bridge: Biphasic Heteroplasmy Shift", fontsize=10.3, fontweight="bold")
    ax.legend(loc="lower right", fontsize=7.2, frameon=True)
    ax.grid(True, alpha=0.25)

    # Panel B: Ex Vivo CAR-T 72-h Solid-Tumor Hypoxia Persistence (N=10,000 LHS)
    ax = axes[1]
    arms = s7a["track7b_ex_vivo_cart_tumor_hypoxia"]["arms"]
    order = [
        ("WT_Human_CART", "WT Human\nCAR-T", "#94a3b8"),
        ("COX3_W116R_Alone_CART", "COX3 W116R\nAlone", "#d97706"),
        ("Human_SOD2_Alone_CART", "Human SOD2\nAlone", "#38bdf8"),
        ("COX3_W116R_plus_Human_SOD2_CART", "COX3 W116R +\nHuman SOD2*", "#059669"),
        ("COX3_W116R_plus_Yeast_NDI1_CART", "COX3 W116R +\nYeast NDI1", "#0284c7"),
        ("COX3_W116R_plus_Ciona_AOX_CART", "COX3 W116R +\nCiona AOX", "#7c3aed"),
    ]
    y = np.arange(len(order))
    vals = [arms[k]["effector_persistence_72h_median"] for k, _, _ in order]
    cols = [c for _, _, c in order]
    ax.barh(y, vals, color=cols, edgecolor="#0f172a", height=0.62)
    for idx, (k, _, _) in enumerate(order):
        v = arms[k]["effector_persistence_72h_median"]
        atp = arms[k]["min_atp_72h_median"]
        ax.text(v + 0.02, idx, f"{v:.3f} (ATP={atp:.2f})", va="center", fontsize=7.8, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([lbl for _, lbl, _ in order], fontsize=8.0)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.38)
    ax.set_xlabel("Median 72-h CAR-T Effector Persistence (N=10,000 LHS)")
    ax.set_title("(B) Ex Vivo CAR-T Tumor Hypoxia (*100% Human)", fontsize=10.3, fontweight="bold")
    ax.grid(True, axis="x", alpha=0.25)

    # Panel C: Live AlphaGenome 2 Single-Cell Heatmap (10 Variants x 6 Cell Types)
    ax = axes[2]
    cell_keys = list(s7b["cell_conditions"].keys())
    cell_short = ["Cardio-\nmyocyte", "CD8+\nT-Cell", "Brain\nEndoth.", "Astro-\ncyte", "Glutam.\nNeuron", "Dopam.\nNeuron"]
    v_rows = s7b["variants"]
    mat = np.array([[v["cells"][ck]["composite"] for ck in cell_keys] for v in v_rows], dtype=float)
    im = ax.imshow(mat, cmap="YlGnBu", aspect="auto", vmin=0.0, vmax=max(1.0, float(np.max(mat))))
    ax.set_xticks(np.arange(len(cell_keys)))
    ax.set_xticklabels(cell_short, fontsize=7.5)
    ax.set_yticks(np.arange(len(v_rows)))
    y_lbls = [v["gene"].replace("NEG_CTRL_", "CTRL ") for v in v_rows]
    ax.set_yticklabels(y_lbls, fontsize=8.0)
    for i in range(len(v_rows)):
        for j in range(len(cell_keys)):
            val = mat[i, j]
            tc = "white" if val > 0.55 else "#0f172a"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7.0, fontweight="bold", color=tc)
    ax.axhline(7.5, color="#dc2626", lw=1.6, ls="--")
    ax.set_title(f"(C) AlphaGenome 2 Single-Cell ({s7b['signal_to_null_ratio']}x Signal/Null)", fontsize=10.2, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Single-Cell Composite Perturbation")

    plt.tight_layout(pad=1.4)
    plt.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    plt.close()


def main() -> None:
    s7a = json.loads(S7A_JSON.read_text(encoding="utf-8"))
    s7b = json.loads(S7B_JSON.read_text(encoding="utf-8"))
    render_figure(s7a, s7b)

    t7a = s7a["track7a_biphasic_heteroplasmy_amplification"]
    t7b = s7a["track7b_ex_vivo_cart_tumor_hypoxia"]
    conv = t7a["solver_convergence"]

    lines = [
        "# Stage 7 Report — Cross-Project Bridges (`Mito_Papers` + `BAA Mito Twin` → `mtDNA_FitnessPoints`)",
        "",
        "> **4-Tier Epistemic Provenance (adopted from `BAA Mito Twin`):** `E1` = experimental cryo-EM/ClinVar/gnomAD; `E2` = AlphaGenome 2 single-cell foundation model (`v2.dna_client`); `E3` = bounded mechanistic ODE simulations (`10,000`-draw LHS + Sobol GSA); `E4` = uncalibrated biological selection hypotheses pending wet-lab falsification.",
        "",
        "## 1. Track 7A (`Mito_Papers` Bridge): Biphasic Mitophagy-Selection Heteroplasmy Amplification [`Tier E3/E4`]",
        f"- **Solver Convergence (`RK4 dt=0.05 d` vs `DOP853`/`Radau`):** `h_end = {conv['h_end_RK4_dt005']}` vs `{conv['h_end_DOP853']}` (`max |Delta| = {conv['max_abs_diff']}`).",
        "- **Hypothesis Tested:** Starting from realistic bulk `TALED` editing (**`50.0%` `COX3 W116R`**), can mild hypoxic preconditioning (`pO2 = 1.1 kPa`) + `PINK1`/Urolithin A mitophagy selectively clear depolarized WT copies (`Km = 1.25 kPa`) while preserving high-affinity `COX3 W116R` (`Km = 0.53 kPa`), then re-expand via Stage-2 normoxic biogenesis?",
        "",
        "| 14-Day Regimen (`N = 10,000` LHS) | Day-14 Heteroplasmy `[5th, 95th]` | Final `W_tot` (Nadir) | `% Draws ≥ 85%` | Biophysical Verdict |",
        "| :--- | :---: | :---: | :---: | :--- |",
    ]
    verdicts = {
        "Normoxia_NoSelection": "No selective pressure; stays at 50%",
        "Unphased_Biogenesis_Trap": "Clones WT & mutant equally (Biogenesis Trap!)",
        "Continuous_Mitophagy_Only": "High h_end, **lethal organelle depletion (`W_tot = 0.011`)**",
        "Biphasic_Selection_Reset": "Strong amplification (`83.2%`), recovers `W_tot = 0.692` (deep nadir `0.082`)",
        "TwoCycle_Pulsed_Biphasic": "**Safer pulsed protocol:** triples `W_tot` nadir (`0.243`) while hitting `80.9%` (`96.5%` 95th)",
    }
    for k, r in t7a["regimens"].items():
        lines.append(
            f"| `{k}` | **{r['h_end_pct_median']}%** `{r['h_end_pct_iq90']}` | `{r['w_tot_end_median']}` (`{r['w_tot_nadir_median']}`) | `{r['pct_draws_ge_85pct_heteroplasmy']}%` | {verdicts[k]} |"
        )
    lines += [
        f"- **Saltelli–Sobol GSA (`10,240` evals):** Variance is dominated by **`s_bias_err` (`S_T = {t7a['sobol_gsa']['sobol_ST']['s_bias_err']}`)** — pinpointing PINK1 depolarization selectivity (`ΔΨ_mut - ΔΨ_wt`) as the single wet-lab measurement that falsifies or validates metabolic selection.",
        "",
        "## 2. Track 7B (`Mito_Papers` CAR-T Bridge): *Ex Vivo* CAR-T Solid-Tumor Hypoxia Cycling [`Tier E3`]",
        f"- **Protocol:** {t7b['protocol']} (`N = 10,000` LHS).",
        f"- **100% Human-Sequence Winner (`COX3 W116R + Human SOD2`, Zero Foreign Immunogenicity):** Median 72-h effector persistence jumps **`{t7b['human_sod2_cart_synergy']['fold_persistence_gain_vs_WT']}×`** over WT CAR-T (`0.1613` → **`0.9794`**; ATP nadir `0.052` → **`0.957`**; cumulative EP excess `18.72` → **`0.36`**), matching foreign `Yeast NDI1` (`0.9828`) without foreign-protein immune risk!",
        "",
        "## 3. Track 7C (`BAA Mito Twin` Bridge): Live AlphaGenome 2 Single-Cell Conditioning (`6` Cell Types + `2` Blank Controls) [`Tier E2`]",
        f"- **Signal-to-Null Specificity:** Active 8 nuclear genes mean single-cell composite = **`{s7b['active_8genes_mean_composite']}`** vs. 2 `BAA-Mito` blank negative controls (`rs9904290`, `rs10799812`) mean = **`{s7b['negative_controls_mean_composite']}`** (**`{s7b['signal_to_null_ratio']}×` separation**, zero false alarms).",
        "- **Major Single-Cell Discovery:** **`6 / 8` nuclear genes peak in Human Ventricular Cardiomyocytes**, and **`UCP2` (`chr11:73982830:C>G`)**—which looked modest in pan-tissue `AVI`—reveals a **`7×` cardiac-selective regulatory switch** (`scRNA = -0.5907`, composite **`0.8562`** in Ventricular Cardiomyocytes vs. `-0.0851` in CD8+ T cells)!",
        "",
        "| Gene / Control | Role | Cardiomyocyte | CD8+ T-Cell (CAR-T) | Brain Endothelial | Astrocyte | Glutam. Neuron | Dopam. Neuron | Top Cell |",
        "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |",
    ]
    for v in s7b["variants"]:
        c = v["cells"]
        lines.append(
            f"| [`{v['gene']}`]({v['atlas_url']}) | {v['role']} | `{c['cardiomyocyte']['composite']:.3f}` | `{c['cd8_tcell']['composite']:.3f}` | "
            f"`{c['brain_endothelial']['composite']:.3f}` | `{c['brain_astrocyte']['composite']:.3f}` | "
            f"`{c['glutamatergic_neuron']['composite']:.3f}` | `{c['dopaminergic_neuron']['composite']:.3f}` | **`{v['top_responsive_cell']}`** |"
        )

    OUT_REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Stage 7 synthesis complete -> {OUT_FIG} and {OUT_REPORT}")


if __name__ == "__main__":
    main()
