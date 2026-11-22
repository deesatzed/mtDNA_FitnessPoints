#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polite-http",
#   "python-dotenv",
# ]
# ///
"""Stage 3A: Fetch live human mtDNA safety ground truth from ClinVar (Pathogenic & Benign) and gnomAD chrM.

Uses the official `clinvar_database` wrapper (`clinvar_api.py`) for ClinVar summaries and
the `polite-http` rate-limited client (0.166 QPS) for gnomAD's `mitochondrial_variants` GraphQL query.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from polite_http import http_client

CLINVAR_SCRIPTS = Path("/Users/o2satz/.gemini/config/plugins/science/skills/clinvar_database/scripts")
sys.path.insert(0, str(CLINVAR_SCRIPTS))
import clinvar_api  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
CLINVAR_PATH_OUT = RAW_DIR / "clinvar_mt_pathogenic_summaries.json"
CLINVAR_BENIGN_OUT = RAW_DIR / "clinvar_mt_benign_summaries.json"
GNOMAD_MT_OUT = RAW_DIR / "gnomad_mt_r3_variants.json"


def fetch_clinvar_summaries_for_query(query: str, out_path: Path, max_ids: int = 2000) -> list[dict]:
    if out_path.exists():
        data = json.loads(out_path.read_text(encoding="utf-8"))
        print(f"[CACHE] {out_path.name}: {len(data)} summaries")
        return data
    cv = clinvar_api.ClinVarClient()
    sr = cv.search_variants(query=query, retmax=max_ids)
    ids = sr["variant_ids"]
    print(f"[CLINVAR] Query matched {sr['total_count']} (fetching summaries for {len(ids)} IDs)...")
    all_summaries = []
    chunk = 100
    for i in range(0, len(ids), chunk):
        batch = ids[i : i + chunk]
        res = cv.get_interpretation_summary(batch)
        all_summaries.extend(res)
        print(f"  ClinVar summary batch {i + 1}-{min(i + chunk, len(ids))}/{len(ids)} -> total {len(all_summaries)}")
    out_path.write_text(json.dumps(all_summaries, indent=2), encoding="utf-8")
    return all_summaries


def fetch_gnomad_chrm(out_path: Path) -> list[dict]:
    if out_path.exists():
        data = json.loads(out_path.read_text(encoding="utf-8"))
        print(f"[CACHE] {out_path.name}: {len(data)} chrM variants")
        return data

    client = http_client.HttpClient("https://gnomad.broadinstitute.org", qps=0.1666)
    query = """
    query($chrom: String!, $start: Int!, $stop: Int!, $dataset: DatasetId!) {
      region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
        mitochondrial_variants(dataset: $dataset) {
          variant_id
          pos
          gene_symbol
          consequence
          hgvsc
          hgvsp
          ac_hom
          ac_het
          an
          max_heteroplasmy
          flags
        }
      }
    }
    """
    resp = client.fetch_json(
        "https://gnomad.broadinstitute.org/api",
        method="POST",
        json_body={
            "query": query,
            "variables": {"chrom": "M", "start": 1, "stop": 16569, "dataset": "gnomad_r3"},
        },
    )
    mt_vars = ((resp or {}).get("data") or {}).get("region", {}).get("mitochondrial_variants") or []
    print(f"[GNOMAD] Fetched {len(mt_vars)} live human chrM variants from gnomAD v3!")
    out_path.write_text(json.dumps(mt_vars, indent=2), encoding="utf-8")
    return mt_vars


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    mt_vars = fetch_gnomad_chrm(GNOMAD_MT_OUT)
    path_s = fetch_clinvar_summaries_for_query(
        'MT[chr] AND ("pathogenic"[clinsig] OR "likely pathogenic"[clinsig])',
        CLINVAR_PATH_OUT,
        max_ids=1000,
    )
    ben_s = fetch_clinvar_summaries_for_query(
        'MT[chr] AND ("benign"[clinsig] OR "likely benign"[clinsig])',
        CLINVAR_BENIGN_OUT,
        max_ids=2000,
    )
    print(f"\nStage 3A fetch complete: gnomAD chrM={len(mt_vars)}, ClinVar Path/LP={len(path_s)}, ClinVar Ben/LB={len(ben_s)}")


if __name__ == "__main__":
    main()
