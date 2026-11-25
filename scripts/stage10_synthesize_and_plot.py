#!/usr/bin/env python3
"""Stage 10 Synthesis & Publication Figure:
  - Stage 10A: Live Apple Silicon MPS Boltz-2 Diffusion Co-Folding of 2 Thermal & Baric
    Mitonuclear Dimers (4 structures: ND4L WT/Y85H + NDUFS2 vs. 9I4I; CYTB WT/Y224H + CYC1 vs. 9HZL).
  - Stage 10B: Tier-3 Ecological & Pollinator Thermal-Hypoxia Conservation Transfer across
    5 keystone species (Western Honeybee, Buff-Tailed Bumblebee, Great Barrier Reef Stony Coral,
    Atlantic Salmon, Rainbow Trout) with strict non-gene-drive biosafety governance.

Outputs:
  - data/stage10a_boltz2_thermal_baric.json
  - data/stage10_report.md
  - figures/stage10_boltz2_thermal_and_ecological.png
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser, Superimposer
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
STRUCT_DIR = ROOT / "data" / "structures"
BOLTZ_PRED_DIR = ROOT / "data" / "stage10_boltz2" / "preds" / "boltz_results_yamls" / "predictions"
S10A_JSON = ROOT / "data" / "stage10a_boltz2_thermal_baric.json"
S10B_JSON = ROOT / "data" / "stage10b_ecological_transfer.json"
OUT_FIG = ROOT / "figures" / "stage10_boltz2_thermal_and_ecological.png"
OUT_REPORT = ROOT / "data" / "stage10_report.md"

DIMER_PAIRS = [
    {
        "pair_label": "ND4L_Y85H + NDUFS2 (Desert Heat)",
        "short_label": "ND4L Y85H + NDUFS2",
        "wt_name": "ND4L_WT_NDUFS2",
        "mut_name": "ND4L_Y85H_NDUFS2",
        "mt_gene": "ND4L",
        "nuc_gene": "NDUFS2",
        "mut_pos": 85,
        "pdb_id": "9I4I",
        "cif_gz": STRUCT_DIR / "pdb_00009i4i.cif.gz",
        "ref_ch_a": "k",
        "ref_ch_b": "Q",
    },
    {
        "pair_label": "CYTB_Y224H + CYC1 (Deep-Sea Dive)",
        "short_label": "CYTB Y224H + CYC1",
        "wt_name": "CYTB_WT_CYC1",
        "mut_name": "CYTB_Y224H_CYC1",
        "mt_gene": "CYTB",
        "nuc_gene": "CYC1",
        "mut_pos": 224,
        "pdb_id": "9HZL",
        "cif_gz": STRUCT_DIR / "pdb_00009hzl.cif.gz",
        "ref_ch_a": "J",
        "ref_ch_b": "H",
    },
]


def load_ca_atoms(cif_path: Path, ch_a: str, ch_b: str, gzipped: bool = False) -> tuple[list, list]:
    parser = MMCIFParser(QUIET=True)
    if gzipped:
        with gzip.open(cif_path, "rt", encoding="utf-8") as fh:
            model = next(parser.get_structure("s", fh).get_models())
    else:
        with cif_path.open("r", encoding="utf-8") as fh:
            model = next(parser.get_structure("s", fh).get_models())
    ca_a = [r["CA"] for r in model[ch_a] if r.id[0] == " " and "CA" in r]
    ca_b = [r["CA"] for r in model[ch_b] if r.id[0] == " " and "CA" in r]
    return ca_a, ca_b


def kabsch_rmsd(atoms_ref: list, atoms_mob: list) -> float:
    n = min(len(atoms_ref), len(atoms_mob))
    sup = Superimposer()
    sup.set_atoms(atoms_ref[:n], atoms_mob[:n])
    return round(float(sup.rms), 3)


def analyze_stage10_boltz2() -> dict:
    results = []
    for cfg in DIMER_PAIRS:
        ref_a, ref_b = load_ca_atoms(cfg["cif_gz"], cfg["ref_ch_a"], cfg["ref_ch_b"], gzipped=True)
        wt_dir = BOLTZ_PRED_DIR / cfg["wt_name"]
        mut_dir = BOLTZ_PRED_DIR / cfg["mut_name"]
        wt_conf = json.loads((wt_dir / f"confidence_{cfg['wt_name']}_model_0.json").read_text(encoding="utf-8"))
        mut_conf = json.loads((mut_dir / f"confidence_{cfg['mut_name']}_model_0.json").read_text(encoding="utf-8"))

        wt_a, wt_b = load_ca_atoms(wt_dir / f"{cfg['wt_name']}_model_0.cif", "A", "B", gzipped=False)
        mut_a, mut_b = load_ca_atoms(mut_dir / f"{cfg['mut_name']}_model_0.cif", "A", "B", gzipped=False)

        wt_chA_rmsd = kabsch_rmsd(ref_a, wt_a)
        mut_chA_rmsd = kabsch_rmsd(ref_a, mut_a)
        wt_dimer_rmsd = kabsch_rmsd(ref_a + ref_b, wt_a + wt_b)
        mut_dimer_rmsd = kabsch_rmsd(ref_a + ref_b, mut_a + mut_b)
        wt_vs_mut_dimer = kabsch_rmsd(wt_a + wt_b, mut_a + mut_b)

        # 10 A pocket C-alpha RMSD around mutated residue
        wt_all = wt_a + wt_b
        mut_all = mut_a + mut_b
        site_coord = wt_a[cfg["mut_pos"] - 1].get_coord()
        wt_coords = np.array([a.get_coord() for a in wt_all])
        idx_10A = [i for i in range(len(wt_all)) if float(np.linalg.norm(wt_coords[i] - site_coord)) <= 10.0]
        pocket_rmsd = kabsch_rmsd([wt_all[i] for i in idx_10A], [mut_all[i] for i in idx_10A]) if len(idx_10A) >= 4 else 0.0

        wt_iptm = round(float(wt_conf.get("iptm", 0.0)), 4)
        mut_iptm = round(float(mut_conf.get("iptm", 0.0)), 4)
        delta_iptm = round(mut_iptm - wt_iptm, 4)

        if mut_iptm >= 0.65 and delta_iptm >= -0.05 and pocket_rmsd <= 1.2:
            verdict = "INTERFACE_PRESERVED_HIGH_FIDELITY"
        elif wt_iptm < 0.30 and mut_chA_rmsd <= 2.2 and abs(delta_iptm) < 0.05:
            verdict = "SCAFFOLD_DEPENDENT_DIMER_MONOMER_CONSERVED"
        else:
            verdict = "INTERFACE_REMODELED_REVIEW"
        rec = {
            "pair": cfg["short_label"],
            "pair_label": cfg["pair_label"],
            "pdb_id": cfg["pdb_id"],
            "wt": {
                "name": cfg["wt_name"],
                "iptm": wt_iptm,
                "ptm": round(float(wt_conf.get("ptm", 0.0)), 4),
                "complex_plddt": round(float(wt_conf.get("complex_plddt", 0.0)), 4),
                "chainA_rmsd_vs_cryoem_A": wt_chA_rmsd,
                "dimer_rmsd_vs_cryoem_A": wt_dimer_rmsd,
            },
            "mutant": {
                "name": cfg["mut_name"],
                "iptm": mut_iptm,
                "ptm": round(float(mut_conf.get("ptm", 0.0)), 4),
                "complex_plddt": round(float(mut_conf.get("complex_plddt", 0.0)), 4),
                "chainA_rmsd_vs_cryoem_A": mut_chA_rmsd,
                "dimer_rmsd_vs_cryoem_A": mut_dimer_rmsd,
            },
            "delta_iptm": delta_iptm,
            "wt_vs_mut_dimer_rmsd_A": wt_vs_mut_dimer,
            "wt_vs_mut_10A_pocket_rmsd_A": pocket_rmsd,
            "verdict": verdict,
        }
        results.append(rec)
        print(
            f"  [Boltz-2] {cfg['short_label']:20s} ({cfg['pdb_id']}) | WT ipTM={wt_iptm:.4f} -> MUT ipTM={mut_iptm:.4f} "
            f"(d={delta_iptm:+.4f}) | ChainA vs CryoEM={mut_chA_rmsd:.2f} A | WT-vs-Mut Dimer={wt_vs_mut_dimer:.2f} A "
            f"| 10A Pocket={pocket_rmsd:.2f} A -> {verdict}"
        )

    payload = {
        "stage": "10A_Live_Boltz2_Thermal_Baric_Dimers",
        "epistemic_tier": "E2_FOUNDATION_MODEL",
        "engine": "Boltz-2 (MIT License, Apple Silicon Metal MPS + ColabFold MSA)",
        "n_dimers_predicted": len(results) * 2,
        "pairs": results,
    }
    S10A_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def render_figure(s10a: dict, s10b: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(16.2, 6.0), gridspec_kw={"width_ratios": [1.0, 1.35]})

    # --- Panel A: Stage 6 + Stage 10A Live Boltz-2 ipTM Portfolio (All 5 Mitonuclear Dimers) ---
    ax = axes[0]
    s6 = json.loads((ROOT / "data" / "stage6_structural_bioai_summary.json").read_text(encoding="utf-8"))
    all_pairs = []
    for p in s6["boltz2_cofolding_wt_vs_mutant"]:
        all_pairs.append((p["pair"] + "\n(Stage 6 Hypoxia)", p["wt"]["iptm"], p["mutant"]["iptm"], p["wt_vs_mut_10A_pocket_rmsd_A"]))
    for p in s10a["pairs"]:
        tag = "Stage 10A Heat" if "ND4L" in p["pair"] else "Stage 10A Dive"
        all_pairs.append((p["pair"] + f"\n({tag})", p["wt"]["iptm"], p["mutant"]["iptm"], p["wt_vs_mut_10A_pocket_rmsd_A"]))

    y = np.arange(len(all_pairs))
    h = 0.34
    wt_v = [x[1] for x in all_pairs]
    mut_v = [x[2] for x in all_pairs]
    ax.barh(y - h / 2, wt_v, height=h, color="#94a3b8", edgecolor="#0f172a", label="WT Human Dimer (Boltz-2 ipTM)")
    ax.barh(y + h / 2, mut_v, height=h, color="#0284c7", edgecolor="#0f172a", label="Extremophile Mutant Dimer (Boltz-2 ipTM)")
    for i, (_, wv, mv, prmsd) in enumerate(all_pairs):
        if wv < 0.30:
            lbl_txt = f"{mv:.3f} (WT {wv:.3f}; needs ND1/ND2 bridge)"
        else:
            lbl_txt = f"{mv:.3f} (d={mv-wv:+.3f}, 10A={prmsd:.2f}A)"
        ax.text(mv + 0.015, i + h / 2, lbl_txt, va="center", fontsize=7.5, fontweight="bold", color="#0369a1")
    ax.set_yticks(y)
    ax.set_yticklabels([x[0] for x in all_pairs], fontsize=8.0)
    ax.invert_yaxis()
    ax.set_xlim(0.0, 1.44)
    ax.set_ylim(len(all_pairs) - 0.35, -0.95)
    ax.set_xlabel("Boltz-2 Interface Confidence (ipTM)")
    ax.set_title("(A) Complete 10-Structure Live Boltz-2 Portfolio (Hypoxia, Heat & Dive)", fontsize=10.1, fontweight="bold")
    ax.legend(loc="upper center", fontsize=7.4, frameon=True, ncol=2)
    ax.grid(True, axis="x", alpha=0.25)

    # --- Panel B: Stage 10B Tier-3 Ecological Conservation Transfer (3 Keystone Scenarios, N=5k LHS) ---
    ax = axes[1]
    scs = s10b["ecological_scenarios_lhs5k"]
    sc_order = [
        ("E1_Honeybee_42C_Fungicide_Flight", "Honeybee (Apis)\n42.5C + Fungicide\n(Mito_Papers Bridge)"),
        ("E2_Coral_31_5C_Marine_Heatwave",   "Reef Coral (Acropora)\n27C -> 31.5C\nMarine Heatwave"),
        ("E3_Salmon_22C_WarmRiver_Hypoxia",  "Atlantic Salmon\n12C -> 22C +\nRiver Hypoxia"),
    ]
    xg = np.arange(len(sc_order))
    bw = 0.18
    arm_specs = [
        ("WT_Vulnerable",         "#94a3b8", "WT Vulnerable"),
        ("Natural_TightCoupled",  "#d97706", "Natural Heat Haplotype Alone"),
        ("Nutritional_BCR_Buffer","#0284c7", "Nutritional / Apis-BCR Buffer"),
        ("Combined_Haplotype_BCR","#059669", "Combined Haplotype + BCR"),
    ]
    for idx, (akey, col, lbl) in enumerate(arm_specs):
        vals = [scs[sk]["arms"][akey]["viab_med"] for sk, _ in sc_order]
        off = (idx - 1.5) * bw
        bars = ax.bar(xg + off, vals, width=bw, color=col, edgecolor="#0f172a", label=lbl)
        for b in bars:
            vh = b.get_height()
            ax.text(b.get_x() + b.get_width() / 2, vh + 0.015, f"{vh:.2f}", ha="center", va="bottom", fontsize=7.4, fontweight="bold")

    ax.set_xticks(xg)
    ax.set_xticklabels([lbl for _, lbl in sc_order], fontsize=8.4)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Median Ecological Bioenergetic Viability (N=5,000 LHS)")
    ax.set_title("(B) Stage 10B: Tier-3 Non-Gene-Drive Conservation Transfer (5 RefSeq Species)", fontsize=10.1, fontweight="bold")
    ax.legend(loc="upper right", fontsize=7.5, frameon=True, ncol=2)
    ax.grid(True, axis="y", alpha=0.25)

    fig.suptitle(
        "MitoFitness Atlas Stage 10: Live Boltz-2 Thermal/Baric Dimer Co-Folding & Tier-3 Ecological Conservation Transfer",
        fontsize=12.0, fontweight="bold", y=0.985,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=220)
    plt.close(fig)
    print(f"Saved figure -> {OUT_FIG}")


def write_report(s10a: dict, s10b: dict) -> None:
    lines = [
        "# Stage 10 Verification Report — Live Boltz-2 Thermal/Baric Co-Folding & Tier-3 Ecological Conservation Transfer",
        "",
        "## 1. Stage 10A — Live Apple Silicon MPS Boltz-2 Thermal & Baric Dimers (`4` New Co-Folded Structures)",
    ]
    for p in s10a["pairs"]:
        lines.append(
            f"- **`{p['pair_label']}` (`{p['pdb_id']}`):** WT `ipTM = {p['wt']['iptm']:.4f}` → Mutant **`ipTM = {p['mutant']['iptm']:.4f}`** "
            f"(`Δ = {p['delta_iptm']:+.4f}`), Chain A vs. cryo-EM = **`{p['mutant']['chainA_rmsd_vs_cryoem_A']:.2f} Å`**, "
            f"WT vs. Mutant Dimer RMSD = **`{p['wt_vs_mut_dimer_rmsd_A']:.2f} Å`**, `10 Å` pocket RMSD = **`{p['wt_vs_mut_10A_pocket_rmsd_A']:.2f} Å`** → **`{p['verdict']}`**."
        )
    lines.extend([
        "",
        "## 2. Stage 10B — Tier-3 Ecological & Pollinator Conservation Transfer (Strictly Non-Gene-Drive)",
        f"- **Solver Convergence (`RK4` vs `DOP853`/`Radau`):** `max |Delta| = {s10b['solver_convergence']['max_abs_delta']:.3e}` (`passed = {s10b['solver_convergence']['passed']}`).",
        "- **Live NCBI RefSeq Mitogenome Cross-Phylum Audit (`5` Keystone Species):**",
    ])
    for sp, info in s10b["species_alignment_audit"].items():
        gids = info["gene_identities_vs_human"]
        lines.append(
            f"  - **{info['common_name']}** (`{info['refseq_accession']}`, `{info['clade']}`): "
            f"`COX1 = {gids.get('COX1',{}).get('pct_identity_vs_human',0):.1f}%`, "
            f"`CYTB = {gids.get('CYTB',{}).get('pct_identity_vs_human',0):.1f}%`, "
            f"Complex III $Q_o$ strobilurin pocket = **`{info['qo_strobilurin_pocket_conservation']}`**."
        )
    lines.append("- **Ecological Conservation Stress Scenarios (`N = 5,000` LHS):**")
    for sk, sc in s10b["ecological_scenarios_lhs5k"].items():
        a = sc["arms"]
        lines.append(
            f"  - **{sc['label']}:** WT `{a['WT_Vulnerable']['viab_med']:.3f}` | "
            f"Natural Tight-Coupled `{a['Natural_TightCoupled']['viab_med']:.3f}` (`{a['Natural_TightCoupled']['delta_vs_wt_med']:+.3f}`) | "
            f"Nutritional/BCR Buffer `{a['Nutritional_BCR_Buffer']['viab_med']:.3f}` (`{a['Nutritional_BCR_Buffer']['delta_vs_wt_med']:+.3f}`) | "
            f"**Combined Haplotype + BCR `{a['Combined_Haplotype_BCR']['viab_med']:.3f}` (`{a['Combined_Haplotype_BCR']['delta_vs_wt_med']:+.3f}`)**."
        )
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Saved report -> {OUT_REPORT}")


def main() -> None:
    s10a = analyze_stage10_boltz2()
    s10b = json.loads(S10B_JSON.read_text(encoding="utf-8"))
    render_figure(s10a, s10b)
    write_report(s10a, s10b)


if __name__ == "__main__":
    main()
