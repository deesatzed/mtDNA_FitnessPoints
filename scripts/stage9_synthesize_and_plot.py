#!/usr/bin/env python3
"""Stage 9 Synthesis & Publication Figure:
  - Stage 9A: Live OpenMM AMBER14 + OBC2 Atomic Relaxation of 5 Stage-8B Thermal/Baric Leads
    (ND4L Y85H, ND1 S313G, COX1 I57V, CYTB F181L, CYTB Y224H) vs. Stage-6 Benchmarks
    (COX3 W116R +0.014 A vs. Pathogenic Leigh ATP6 L156R +0.160 A).
  - Stage 9B: Live AlphaGenome 2 Single-Cell Organ-Compartment Conditioning across 6 Clinical
    Organs (Heart Cardiomyocyte, Liver Hepatocyte, Lung AEC2, Kidney Proximal Tubule,
    Pancreatic Beta Cell, Blood CD8+ T Cell) x 7 Active Promoters + 2 Blank Negative Controls.

Outputs:
  - figures/stage9_openmm_thermal_and_organ_ag2.png
  - data/stage9_report.md
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
S9A_JSON = ROOT / "data" / "stage9a_openmm_thermal_baric.json"
S9B_JSON = ROOT / "data" / "stage9b_alphagenome2_organs.json"
OUT_FIG = ROOT / "figures" / "stage9_openmm_thermal_and_organ_ag2.png"
OUT_REPORT = ROOT / "data" / "stage9_report.md"


def render_figure(s9a: dict, s9b: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 6.2), gridspec_kw={"width_ratios": [1.0, 1.25]})

    # --- Panel A: OpenMM AMBER14 8 A Pocket Delta-RMSD (Stage 9A Leads vs Stage 6 Anchors) ---
    ax = axes[0]
    rows = []
    # Add Stage 6 benchmark anchors for direct comparison
    anchors = s9a["benchmark_anchors"]
    rows.append(("COX3 W116R\n(Stage 6 Hypoxia Ref)", anchors["Stage6_COX3_W116R_delta_rmsd_A"], "#0284c7", "6 new H-bonds"))
    for r in s9a["results"]:
        col = "#d97706" if "Desert_Heat" in r["role"] else ("#0891b2" if "Arctic_Cold" in r["role"] else "#7c3aed")
        tag = r["role"].replace("_Lead_0Bys", " (0 bys)").replace("_", " ")
        hb = f"d_min={r['mutant']['post_relax_sc_dmin_A']:.2f}A"
        if r["new_mutant_polar_contacts_3p5A"]:
            hb += f", +{len(r['new_mutant_polar_contacts_3p5A'])} H-bond"
        rows.append((f"{r['variant'].replace('_',' ')}\n({tag})", r["delta_pocket_rmsd_A"], col, hb))
    rows.append(("ATP6 L156R\n(Pathogenic Leigh Ctrl)", anchors["Stage6_Pathogenic_ATP6_L156R_delta_rmsd_A"], "#dc2626", "Pathogenic distortion"))

    y = np.arange(len(rows))
    vals = [x[1] for x in rows]
    cols = [x[2] for x in rows]
    ax.barh(y, vals, color=cols, edgecolor="#0f172a", height=0.60)
    ax.axvline(0.0, color="#0f172a", lw=1.0)
    ax.axvline(0.08, color="#dc2626", ls="--", lw=1.2, label="Low-Perturbation Threshold (+0.08 A)")

    for i, (_, val, _, note) in enumerate(rows):
        xpos = max(val, 0.0) + 0.005
        ax.text(xpos, i, f"{val:+.3f} A ({note})", va="center", ha="left", fontsize=7.8, fontweight="bold")

    ax.set_yticks(y)
    ax.set_yticklabels([x[0] for x in rows], fontsize=8.0)
    ax.invert_yaxis()
    ax.set_xlim(min(-0.04, min(vals) - 0.02), max(vals) + 0.215)
    ax.set_xlabel("OpenMM AMBER14+OBC2 8 A Pocket Heavy-Atom Delta-RMSD (Mut - WT, Angstroms)")
    ax.set_title("(A) Stage 9A: Live OpenMM Cryo-EM Stability of Thermal & Baric Leads", fontsize=10.3, fontweight="bold")
    ax.legend(loc="lower right", fontsize=7.6, frameon=True)
    ax.grid(True, axis="x", alpha=0.25)

    # --- Panel B: Live AlphaGenome 2 Clinical Organ-Compartment Heatmap ---
    ax = axes[1]
    org_keys = list(s9b["organ_conditions"].keys())
    org_labels = [
        "Heart\nCardiomyocyte\n(CABG/NMP)",
        "Liver\nHepatocyte\n(Pringle/MASH)",
        "Lung AEC2\nPneumocyte\n(IPF/Preemie)",
        "Kidney\nProx. Tubule\n(CSA-AKI)",
        "Pancreatic\nBeta Cell\n(T1D Islets)",
        "Blood CD8+\nT Cell\n(CAR-T)",
    ]
    v_rows = s9b["variants"]
    mat = np.array([[v["organs"][ok]["composite"] for ok in org_keys] for v in v_rows], dtype=float)
    vmax = min(1.2, max(0.85, float(np.percentile(mat, 97))))
    im = ax.imshow(mat, cmap="YlGnBu", aspect="auto", vmin=0.0, vmax=vmax)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            fg = "#ffffff" if val > 0.55 * vmax else "#0f172a"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8.0, fontweight="bold", color=fg)

    ax.axhline(6.5, color="#dc2626", lw=1.6, ls="--")
    ax.set_xticks(np.arange(len(org_labels)))
    ax.set_xticklabels(org_labels, fontsize=7.8)
    ylabels = [
        f"{v['gene']} ({v['role'].replace('Negative_Control_Blank','Blank Ctrl').replace('_',' ')})"
        for v in v_rows
    ]
    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=7.9)
    cbar = fig.colorbar(im, ax=ax, fraction=0.040, pad=0.02)
    cbar.set_label("AlphaGenome 2 Organ Composite Score", fontsize=8.4)
    ax.set_title(
        f"(B) Stage 9B: Live AlphaGenome 2 Across 6 Clinical Organs (Sig/Null = {s9b['signal_to_null_ratio']}x)",
        fontsize=10.3, fontweight="bold",
    )

    fig.suptitle(
        "MitoFitness Atlas Stage 9: Atomic Cryo-EM Stability of Thermal/Baric Leads & Live AlphaGenome 2 Organ Conditioning",
        fontsize=12.0, fontweight="bold", y=0.985,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=220)
    plt.close(fig)
    print(f"Saved figure -> {OUT_FIG}")


def write_report(s9a: dict, s9b: dict) -> None:
    lines = [
        "# Stage 9 Verification Report — Atomic Physics (OpenMM AMBER14) for Thermal/Baric Leads & Live AlphaGenome 2 Clinical Organ Conditioning",
        "",
        "## 1. Stage 9A — Live OpenMM `AMBER14 + OBC2` Cryo-EM Relaxation (`5` Thermal & Baric Leads)",
        f"- **Pass Rate:** **`{s9a['n_stable_low_perturbation']}/{s9a['n_evaluated']}` `STABLE_LOW_PERTURBATION`** vs. Stage-6 Pathogenic Leigh `ATP6 L156R` (`ΔRMSD = +{s9a['benchmark_anchors']['Stage6_Pathogenic_ATP6_L156R_delta_rmsd_A']:.3f} Å`).",
    ]
    for r in s9a["results"]:
        hb = ", ".join(r["new_mutant_polar_contacts_3p5A"][:3]) if r["new_mutant_polar_contacts_3p5A"] else "None"
        lines.append(
            f"- **`{r['variant']}`** (`{r['role']}`, `{r['pdb_id']}` chains `{r['subassembly_chains']}`): "
            f"WT `{r['wt']['pocket_8A_heavy_rmsd_A']:.3f} Å` → MUT `{r['mutant']['pocket_8A_heavy_rmsd_A']:.3f} Å` "
            f"(**`ΔRMSD = {r['delta_pocket_rmsd_A']:+.4f} Å`**), `d_min = {r['mutant']['post_relax_sc_dmin_A']:.2f} Å`, "
            f"new H-bonds: `{hb}` → **`{r['physics_verdict']}`**."
        )
    lines.extend([
        "",
        "## 2. Stage 9B — Live `AlphaGenome 2` Across 6 Clinical Organ Compartments (`Heart, Liver, Lung AEC2, Kidney, Pancreatic Beta, CD8+ T-Cell`)",
        f"- **Signal-to-Null Specificity:** **`{s9b['signal_to_null_ratio']}×`** (`7` active promoters mean `{s9b['active_7genes_mean_composite']:.4f}` vs. `2` `BAA-Mito` blank controls mean `{s9b['negative_controls_mean_composite']:.4f}`).",
    ])
    for v in s9b["variants"]:
        lines.append(
            f"- **`{v['gene']}`** (`{v['role']}`, `{v['variant']}`): mean `{v['mean_composite']:.4f}`, "
            f"top organ = **`{v['top_responsive_organ']}` (`{v['top_responsive_composite']:.4f}`)** "
            f"(Cardio=`{v['organs']['cardiomyocyte']['composite']:.3f}`, Liver=`{v['organs']['hepatocyte']['composite']:.3f}`, "
            f"Lung_AEC2=`{v['organs']['lung_aec2']['composite']:.3f}`, Kidney=`{v['organs']['kidney_pt']['composite']:.3f}`, "
            f"Beta=`{v['organs']['pancreatic_beta']['composite']:.3f}`, TCell=`{v['organs']['cd8_tcell']['composite']:.3f}`)."
        )
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Saved report -> {OUT_REPORT}")


def main() -> None:
    s9a = json.loads(S9A_JSON.read_text(encoding="utf-8"))
    s9b = json.loads(S9B_JSON.read_text(encoding="utf-8"))
    render_figure(s9a, s9b)
    write_report(s9a, s9b)


if __name__ == "__main__":
    main()
