#!/usr/bin/env python3
"""Stage 7B: Cross-Project Bridge from BAA Mito Twin — Live AlphaGenome 2 Single-Cell Conditioning.

Scores the top promoter hotspot variant of all 8 nuclear OXPHOS & buffer genes from Stage 6
(COX6A1, NDUFA13, HIGD2A, ATP5MJ, COX4I1, UQCRFS1, UCP2, SOD2) PLUS 2 matched closed-chromatin
negative-control SNPs (rs9904290 & rs10799812 from BAA-Mito) using live AlphaGenome 2
(`alphagenome.models.v2.dna_client`) across 6 observed scRNA-seq + scATAC-seq human cell types
bridging Heart I/R, CAR-T, and Brain/Vascular aging:
  1. Heart Ventricular Cardiomyocyte (UBERON:0000948_CL:2000046)
  2. Blood Effector Memory CD8+ T Cell (UBERON:0000178_CL:0000913)
  3. Brain Microvascular Endothelial Cell (UBERON:0000955_CL:2000044)
  4. Cortical Astrocyte (UBERON:0009834_CL:0000127)
  5. Cortical Glutamatergic Neuron (UBERON:0009834_CL:0000679)
  6. Brain Dopaminergic Neuron (UBERON:0000955_CL:0000700)

Outputs:
  - data/stage7b_singlecell_alphagenome.json
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

OUT_JSON = ROOT / "data" / "stage7b_singlecell_alphagenome.json"

CELL_CONDITIONS = {
    "cardiomyocyte": ("UBERON:0000948_CL:2000046", "Heart Ventricular Cardiomyocyte"),
    "cd8_tcell": ("UBERON:0000178_CL:0000913", "Blood Effector Memory CD8+ T Cell (CAR-T)"),
    "brain_endothelial": ("UBERON:0000955_CL:2000044", "Brain Microvascular Endothelial Cell"),
    "brain_astrocyte": ("UBERON:0009834_CL:0000127", "Cortical Astrocyte"),
    "glutamatergic_neuron": ("UBERON:0009834_CL:0000679", "Cortical Glutamatergic Neuron"),
    "dopaminergic_neuron": ("UBERON:0000955_CL:0000700", "Brain Dopaminergic Neuron"),
}

# 8 Nuclear OXPHOS/Buffer top promoter SNVs from Stage 6 + 2 BAA-Mito Blank Negative Controls
VARIANTS_TO_SCORE = [
    ("COX6A1",  "Nuclear_OXPHOS_Partner", "chr12", 120438129, "G", "C"),
    ("NDUFA13", "Nuclear_OXPHOS_Partner", "chr19",  19516241, "G", "T"),
    ("HIGD2A",  "Nuclear_OXPHOS_Partner", "chr5",  176388760, "G", "T"),
    ("ATP5MJ",  "Nuclear_OXPHOS_Partner", "chr14", 103921529, "C", "A"),
    ("COX4I1",  "Nuclear_OXPHOS_Partner", "chr16",  85799675, "C", "G"),
    ("UQCRFS1", "Nuclear_OXPHOS_Partner", "chr19",  29213160, "C", "T"),
    ("UCP2",    "Endogenous_Buffer",      "chr11",  73982830, "C", "G"),
    ("SOD2",    "Endogenous_Buffer",      "chr6",  159693241, "T", "G"),
    ("NEG_CTRL_rs9904290",  "Negative_Control_Blank", "chr17",  45985000, "C", "T"),
    ("NEG_CTRL_rs10799812", "Negative_Control_Blank", "chr1",  215190400, "G", "A"),
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
    print(f"Scoring {len(VARIANTS_TO_SCORE)} variants across {len(CELL_CONDITIONS)} single-cell types via AlphaGenome 2 ({scorer_keys})...", flush=True)

    results = []
    for gene, role, chrom, pos, ref, alt in VARIANTS_TO_SCORE:
        t0 = time.time()
        v_obj = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)
        interval = v_obj.reference_interval.resize(seq_len)
        v_str = f"{chrom}:{pos}:{ref}>{alt}"

        def score_cell(ck: str, curie: str, label: str) -> tuple[str, dict]:
            adatas = model.score_variant(interval=interval, variant=v_obj, condition=curie, variant_scorers=scorers)
            df = variant_scorers_v2.tidy_scores([adatas], match_gene_strand=True)
            if df is None or df.empty:
                return ck, {"label": label, "curie": curie, "sc_rna_target_lfc": 0.0, "sc_atac_max_delta": 0.0, "rna_target_lfc": 0.0, "top_tf_delta": 0.0, "splice_delta": 0.0, "composite": 0.0}

            sc_rna = df[df["output_type"] == "SC_GENE_EXPRESSION"]
            sc_atac = df[df["output_type"] == "SC_ATAC"]
            rna = df[df["output_type"] == "RNA_SEQ"]
            tf = df[df["output_type"] == "CHIP_TF"]
            sp = df[df["output_type"] == "SPLICE_SITES"]

            if role != "Negative_Control_Blank" and "gene_name" in df.columns:
                sc_rna_t = sc_rna[sc_rna["gene_name"].astype(str).str.upper() == gene.upper()]
                rna_t = rna[rna["gene_name"].astype(str).str.upper() == gene.upper()]
                sp_t = sp[sp["gene_name"].astype(str).str.upper() == gene.upper()]
            else:
                sc_rna_t = sc_rna
                rna_t = rna
                sp_t = sp

            sc_rna_val = max_signed_abs(sc_rna_t["raw_score"]) if not sc_rna_t.empty else max_signed_abs(sc_rna["raw_score"])
            sc_atac_val = max_signed_abs(sc_atac["raw_score"]) if not sc_atac.empty else 0.0
            rna_val = max_signed_abs(rna_t["raw_score"]) if not rna_t.empty else max_signed_abs(rna["raw_score"])
            tf_val = max_signed_abs(tf["raw_score"]) if not tf.empty else 0.0
            sp_val = max_signed_abs(sp_t["raw_score"]) if not sp_t.empty else 0.0
            comp = round(float(abs(sc_rna_val) + 0.5 * abs(sc_atac_val) + 0.25 * abs(rna_val) + 0.25 * abs(sp_val)), 4)
            return ck, {
                "label": label,
                "curie": curie,
                "sc_rna_target_lfc": sc_rna_val,
                "sc_atac_max_delta": sc_atac_val,
                "rna_target_lfc": rna_val,
                "top_tf_delta": tf_val,
                "splice_delta": sp_val,
                "composite": comp,
            }

        cell_map = {}
        with ThreadPoolExecutor(max_workers=6) as pool:
            futs = {pool.submit(score_cell, ck, curie, lbl): ck for ck, (curie, lbl) in CELL_CONDITIONS.items()}
            for fut in as_completed(futs):
                ck, cdict = fut.result()
                cell_map[ck] = cdict

        ordered = {ck: cell_map[ck] for ck in CELL_CONDITIONS}
        composites = [ordered[ck]["composite"] for ck in CELL_CONDITIONS]
        top_cell = max(CELL_CONDITIONS.keys(), key=lambda ck: ordered[ck]["composite"])
        rec = {
            "gene": gene,
            "role": role,
            "variant": v_str,
            "mean_singlecell_composite": round(float(np.mean(composites)), 4),
            "max_singlecell_composite": round(float(np.max(composites)), 4),
            "top_responsive_cell": top_cell,
            "atlas_url": f"https://deepmind.google.com/science/alphagenome/atlas?q={v_str.replace('>', '%3E')}&m=variant",
            "cells": ordered,
        }
        results.append(rec)
        print(
            f"  [{time.time()-t0:4.1f}s] {gene:20s} ({role:22s}) | mean_comp={rec['mean_singlecell_composite']:.4f} | "
            f"max={rec['max_singlecell_composite']:.4f} ({top_cell}) | "
            f"cardio_scRNA={ordered['cardiomyocyte']['sc_rna_target_lfc']:+.4f}, Tcell_scRNA={ordered['cd8_tcell']['sc_rna_target_lfc']:+.4f}",
            flush=True,
        )

    # Compute specificity summary: 8 active loci vs 2 blank negative controls
    active_means = [r["mean_singlecell_composite"] for r in results if r["role"] != "Negative_Control_Blank"]
    neg_means = [r["mean_singlecell_composite"] for r in results if r["role"] == "Negative_Control_Blank"]
    summary = {
        "model": "AlphaGenome 2 (alphagenome.models.v2.dna_client)",
        "scorers": scorer_keys,
        "cell_conditions": {k: {"curie": v[0], "label": v[1]} for k, v in CELL_CONDITIONS.items()},
        "active_8genes_mean_composite": round(float(np.mean(active_means)), 4),
        "negative_controls_mean_composite": round(float(np.mean(neg_means)), 4),
        "signal_to_null_ratio": round(float(np.mean(active_means) / max(1e-6, np.mean(neg_means))), 2),
        "variants": results,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nStage 7B complete -> {OUT_JSON} (Signal/Null ratio = {summary['signal_to_null_ratio']}x)")


if __name__ == "__main__":
    main()
