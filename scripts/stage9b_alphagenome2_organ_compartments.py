#!/usr/bin/env python3
"""Stage 9B: Live AlphaGenome 2 Single-Cell Conditioning Across Clinical Organ Compartments.

Directly tests our Master Translational Plan's 6 target human organ/cell compartments via live
`alphagenome.models.v2.dna_client` across 5 key buffer & Stage-8B thermal/baric partner promoters
(SOD2, UCP2, PINK1, NDUFS2 [contacts Desert-Heat ND4L Y85H], CYC1 [contacts Diving CYTB Y224H])
PLUS 2 matched BAA-Mito closed-chromatin negative controls (rs9904290, rs10799812):
  1. Heart Ventricular Cardiomyocyte     (UBERON:0000948_CL:2000046) — CABG, STEMI, NMP Heart
  2. Liver Hepatocyte                    (UBERON:0002107_CL:0000182) — Pringle Hepatectomy, NMP Liver, MASH
  3. Lung Type II Pneumocyte (AEC2)      (UBERON:0002048_CL:0002063) — IPF/ILD, 24-wk Preemie, EVLP Lung
  4. Kidney Cortex Proximal Tubule       (UBERON:0001225_CL:0002306) — CABG CSA-AKI, Partial Nephrectomy
  5. Pancreatic Beta Cell (Islet Type B) (UBERON:0000006_CL:0000169) — T1D iPSC-Islets (VX-880)
  6. Blood Effector Memory CD8+ T Cell   (UBERON:0000178_CL:0000913) — Ex Vivo CAR-T

Outputs:
  - data/stage9b_alphagenome2_organs.json
"""

from __future__ import annotations

import json
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import dotenv
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as dna_client_v2
from alphagenome.models.v2 import variant_scorers as variant_scorers_v2

ROOT = Path(__file__).resolve().parents[1]
dotenv.load_dotenv(ROOT / ".env", override=True)
dotenv.load_dotenv(Path.home() / ".env")

OUT_JSON = ROOT / "data" / "stage9b_alphagenome2_organs.json"

ORGAN_CONDITIONS = {
    "cardiomyocyte": ("UBERON:0000948_CL:2000046", "Heart Ventricular Cardiomyocyte (CABG/NMP)"),
    "hepatocyte": ("UBERON:0002107_CL:0000182", "Liver Hepatocyte (Pringle/MASH)"),
    "lung_aec2": ("UBERON:0002048_CL:0002063", "Lung Type II Pneumocyte AEC2 (IPF/Preemie/EVLP)"),
    "kidney_pt": ("UBERON:0001225_CL:0002306", "Kidney Proximal Tubule (CSA-AKI/Nephrectomy)"),
    "pancreatic_beta": ("UBERON:0000006_CL:0000169", "Pancreatic Beta Cell (T1D iPSC-Islets)"),
    "cd8_tcell": ("UBERON:0000178_CL:0000913", "Blood Effector Memory CD8+ T Cell (CAR-T)"),
}

VARIANTS_TO_SCORE = [
    ("SOD2",    "Endogenous_Buffer",      "chr6",  159693241, "T", "G", "Universal matrix superoxide shield"),
    ("UCP2",    "Endogenous_Buffer",      "chr11",  73982830, "C", "G", "Regulated mild proton-leak buffer (caution in beta/heat)"),
    ("PINK1",   "Endogenous_Buffer",      "chr1",   20633437, "C", "A", "Mitophagy quality-control kinase (AVI Phred 20.97; IPF AEC2)"),
    ("NDUFS2",  "Thermal_Partner_CI",     "chr1",  161197085, "G", "T", "3.25 A CI contact of Desert-Heat ND4L Y85H (AVI Phred 25.10)"),
    ("NDUFA13", "Thermal_Partner_CI",     "chr19",  19516241, "G", "T", "3.68 A CI contact of Desert-Heat Camel ND1 S313G"),
    ("UQCRFS1", "Cold_Partner_CIII",      "chr19",  29213160, "C", "T", "4.40 A CIII contact of Arctic-Cold Mammoth CYTB F181L"),
    ("CYC1",    "Baric_Partner_CIII",     "chr8",  144095032, "C", "A", "4.04 A CIII contact of Sperm-Whale CYTB Y224H (AVI Phred 23.10)"),
    ("NEG_CTRL_rs9904290",  "Negative_Control_Blank", "chr17",  45985000, "C", "T", "BAA-Mito blank negative control 1"),
    ("NEG_CTRL_rs10799812", "Negative_Control_Blank", "chr1",  215190400, "G", "A", "BAA-Mito blank negative control 2"),
]


def max_signed_abs(s: pd.Series) -> float:
    clean = s.dropna()
    if clean.empty:
        return 0.0
    return round(float(clean.loc[clean.abs().idxmax()]), 5)


def main() -> None:
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not api_key:
        raise RuntimeError("ALPHAGENOME_API_KEY not found in .env")

    model = dna_client_v2.create(api_key=api_key, address="dns:///gdmscience.googleapis.com:443")
    scorer_keys = [
        k for k in ["SC_GENE_EXPRESSION", "SC_ATAC", "RNA_SEQ", "CHIP_TF", "SPLICE_SITES"]
        if k in variant_scorers_v2.RECOMMENDED_VARIANT_SCORERS
    ]
    scorers = [variant_scorers_v2.RECOMMENDED_VARIANT_SCORERS[k] for k in scorer_keys]
    seq_len = dna_client_v2.SEQUENCE_LENGTH_1MB
    print(f"Scoring {len(VARIANTS_TO_SCORE)} loci across {len(ORGAN_CONDITIONS)} clinical organ compartments...", flush=True)

    results = []
    for gene, role, chrom, pos, ref, alt, note in VARIANTS_TO_SCORE:
        t0 = time.time()
        v_obj = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)
        interval = v_obj.reference_interval.resize(seq_len)
        v_str = f"{chrom}:{pos}:{ref}>{alt}"

        def score_organ(ck: str, curie: str, label: str) -> tuple[str, dict]:
            adatas = model.score_variant(interval=interval, variant=v_obj, condition=curie, variant_scorers=scorers)
            df = variant_scorers_v2.tidy_scores([adatas], match_gene_strand=True)
            if df is None or df.empty:
                return ck, {"label": label, "curie": curie, "sc_rna_target_lfc": 0.0, "sc_atac_max_delta": 0.0, "top_tf_delta": 0.0, "composite": 0.0}

            sc_rna = df[df["output_type"] == "SC_GENE_EXPRESSION"]
            sc_atac = df[df["output_type"] == "SC_ATAC"]
            rna = df[df["output_type"] == "RNA_SEQ"]
            sp = df[df["output_type"] == "SPLICE_SITES"]

            if role != "Negative_Control_Blank" and "gene_name" in df.columns:
                sc_rna_t = sc_rna[sc_rna["gene_name"].astype(str).str.upper() == gene.upper()]
                rna_t = rna[rna["gene_name"].astype(str).str.upper() == gene.upper()]
                sp_t = sp[sp["gene_name"].astype(str).str.upper() == gene.upper()]
            else:
                sc_rna_t = sc_rna
                rna_t = rna
                sp_t = sp

            val_scrna = max_signed_abs(sc_rna_t["raw_score"]) if not sc_rna_t.empty else max_signed_abs(sc_rna["raw_score"])
            val_scatac = max_signed_abs(sc_atac["raw_score"]) if not sc_atac.empty else 0.0
            val_rna = max_signed_abs(rna_t["raw_score"]) if not rna_t.empty else max_signed_abs(rna["raw_score"])
            val_sp = max_signed_abs(sp_t["raw_score"]) if not sp_t.empty else 0.0
            comp = round(abs(val_scrna) + 0.5 * abs(val_scatac) + 0.25 * abs(val_rna) + 0.25 * abs(val_sp), 4)
            return ck, {
                "label": label,
                "curie": curie,
                "sc_rna_target_lfc": val_scrna,
                "sc_atac_max_delta": val_scatac,
                "rna_target_lfc": val_rna,
                "splice_delta": val_sp,
                "composite": comp,
            }

        organ_scores = {}
        with ThreadPoolExecutor(max_workers=6) as ex:
            futs = [ex.submit(score_organ, ck, curie, lbl) for ck, (curie, lbl) in ORGAN_CONDITIONS.items()]
            for fut in as_completed(futs):
                ck, sc = fut.result()
                organ_scores[ck] = sc

        ordered_scores = {ck: organ_scores[ck] for ck in ORGAN_CONDITIONS}
        top_organ = max(ordered_scores.items(), key=lambda kv: kv[1]["composite"])
        mean_comp = round(float(np.mean([v["composite"] for v in ordered_scores.values()])), 4)
        rec = {
            "gene": gene,
            "role": role,
            "variant": v_str,
            "note": note,
            "atlas_url": f"https://deepmind.google.com/science/alphagenome/atlas/variant/{v_str}",
            "organs": ordered_scores,
            "mean_composite": mean_comp,
            "top_responsive_organ": top_organ[0],
            "top_responsive_composite": top_organ[1]["composite"],
        }
        results.append(rec)
        print(
            f"  {gene:20s} ({role:22s}) | mean={mean_comp:.4f} | top={top_organ[0]:15s} "
            f"({top_organ[1]['composite']:.4f}, scRNA={top_organ[1]['sc_rna_target_lfc']:+.4f}) [{time.time()-t0:.1f}s]",
            flush=True,
        )

    active_means = [r["mean_composite"] for r in results if r["role"] != "Negative_Control_Blank"]
    ctrl_means = [r["mean_composite"] for r in results if r["role"] == "Negative_Control_Blank"]
    sig_null_ratio = round(float(np.mean(active_means) / max(1e-6, float(np.mean(ctrl_means)))), 2)

    payload = {
        "stage": "9B_AlphaGenome2_Clinical_Organ_Compartments",
        "epistemic_tier": "E2_FOUNDATION_MODEL",
        "model": "alphagenome.models.v2.dna_client (SEQUENCE_LENGTH_1MB)",
        "organ_conditions": {k: {"curie": c, "label": l} for k, (c, l) in ORGAN_CONDITIONS.items()},
        "active_7genes_mean_composite": round(float(np.mean(active_means)), 4),
        "negative_controls_mean_composite": round(float(np.mean(ctrl_means)), 4),
        "signal_to_null_ratio": sig_null_ratio,
        "variants": results,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nSignal/Null specificity across 6 clinical organ compartments = {sig_null_ratio}x -> {OUT_JSON}")


if __name__ == "__main__":
    main()
