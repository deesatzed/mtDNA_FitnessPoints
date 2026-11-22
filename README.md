# MitoFitness Atlas (`mtDNA_FitnessPoints`)

A standalone, reproducible, 100% real-data computational discovery engine and relational catalog (**SQLite / PostgreSQL-compatible**) that maps wild extremophile mitochondrial DNA (mtDNA) adaptations onto experimental human OXPHOS 3D cryo-EM structures, audits human clinical safety and base-editor feasibility, integrates **AlphaFold DB v6**, **AlphaFold 3 / Boltz-2**, and **AlphaGenome Variant Impact (`AVI`)**, and simulates non-monotonic **Driver + Nuclear Buffer** positive epistasis under the Energy Resistance Principle ($\mathrm{EP} = \dot{e}_R \cdot f^2$).

> **Research Scope:** All structural, regulatory, and biophysical ODE predictions in this repository are strictly *in silico* computational analyses intended to prioritize benchtop *in vitro* human cybrid experiments—not clinical guidance or medical advice.

---

## Architecture & 4-Stage Pipeline

All 13 human mtDNA protein-coding genes (`NC_012920.1`, rCRS, `16,569 bp`) are evaluated under a **3-Tier Algorithmic Progression** (Tier 1 Calibration Core: `ND1, ND2, CYTB, COX1, COX3, ATP6`; Tier 2 Frozen Out-of-Sample Sweep: `ND3, ND4, ND4L, ND5, ND6, COX2, ATP8`; Tier 3 species-agnostic transfer to ecological recipients).

| Phase / Stage | Core Script(s) | Live Public Data & Models | Key Verified Outcome |
| :--- | :--- | :--- | :--- |
| **Phase 0: Foundation** | `scripts/phase0_init_db.py`, `schema/mitofitness_schema.sql` | NCBI RefSeq Human `rCRS` (`NC_012920.1`, 13 CDS under NCBI Table 2) | 8-table portable relational DB (`data/mitofitness.db`), 13/13 exact UniProt lengths, 5/5 landmark parity checks |
| **Stage 1: Convergent Mining** | `scripts/species_panel.py`, `scripts/stage1_fetch_mitogenomes.py`, `scripts/stage1_convergent_mining.py` | 44 wild vertebrate RefSeq mitogenomes (572 wild + 13 human = 585 proteins; 23 matched Extreme-vs-Control pairs across 4 stressor axes) | Honest Thomas & Hahn control null exposed 1D permutation confounding; consensus polarization isolated **131** unique candidates (`65` Class I conserved-site rule-breakers incl. Bar-headed goose `COX3 W116R`, `70` Class II directional dials; `256` rows in `mt_alleles`) |
| **Stage 2: 3D Cryo-EM & FBT Clustering** | `scripts/stage2_3d_fbt_clustering.py` | 4 experimental human cryo-EM OXPHOS holocomplexes from RCSB PDB: **`9I4I`** (CI, 2.63 Å), **`9HZL`** (CIII, 2.52 Å), **`9I6F`** (CIV, 2.95 Å), **`8H9S`** (CV, 2.53 Å) + 85/85 UniProt→HGNC map | **131/131 (100%)** candidates mapped to 3D atomic coordinates; statistically significant **Mitonuclear Interface Enrichment (`39.7%` vs `31.5%`, `OR = 1.45, Fisher p = 0.044`)**; 4 FBT Payoff Equivalence Classes (`256/256` rows in `fitness_payoff_profiles`) |
| **Stage 3: Safety, Base-Editor & Bio-AI Gate** | `scripts/stage3a_fetch_safety.py`, `scripts/stage3_structural_safety_editor_gate.py`, `scripts/stage3c_bioai_cofolding_export.py` | Live **`18,164` gnomAD v3 chrM** variants + **`513` ClinVar Pathogenic/LP** + **AlphaFold DB v6** (`pLDDT`/`PAE`) + live **AlphaGenome `AVI`** promoter scans (`123` SNVs/gene) | **5/5 answer-key controls passed** (Tier-1 + Frozen Tier-2 incl. minus-strand `ND6 M64V`); caught NMR `ND6 A74S` pathogenic hotspot; isolated **42** base-editable candidates (**12** with 0 bystanders, led by `TALED`-editable **`COX3 W116R` `m.9552T>C`**); exported 5 **AF3 JSON & Boltz-2 YAML** manifests |
| **Stage 4: Epistasis ODE, UQ & GSA** | `scripts/stage4_epistasis_synergy_uq.py` | 4-state I/R ODE (`[W, Q_red, Succ, ATP]`) coupled to $\mathrm{EP} = \dot{e}_R \cdot f^2$; `N = 10,000` LHS (`±30%` priors); `10,240`-eval Saltelli-Sobol GSA | `RK4` vs `DOP853`/`Radau` converged (`max |Delta| = 1.029e-07`); proved non-monotonic ATP vs. RET-ROS trade-off and strong positive epistasis (**`COX3 W116R + Yeast NDI1`**, `S = 1.322x`, `E_add = +0.1187`, `93.5%` LHS positive; **`+ Human SOD2`**, `95.1%`; **`+ Ciona AOX`**, `+4.8 kPa` cliff expansion) |

---

## Key Reports, Artifacts & Figures

- **Design Decisions & Audit Trail:** `DESIGN_LOG.md`
- **Database Schema & SQLite Catalog:** `schema/mitofitness_schema.sql`, `data/mitofitness.db`
- **Stage Reports:**
  - `data/stage1_report.md` (`data/stage1_convergent_summary.json`)
  - `data/stage2_report.md` (`data/stage2_3d_fbt_summary.json`)
  - `data/stage3_report.md` (`data/stage3_shortlist.json`)
  - `data/stage4_report.md` (`data/stage4_epistasis_summary.json`)
- **Publication Figures:**
  - `figures/stage2_3d_fbt_equivalence_clusters.png`
  - `figures/stage4_epistasis_synergy_and_cliffs.png`
- **Co-Folding Manifests (AF3 JSON & Boltz-2 YAML):** `data/cofolding_manifests/`
- **Third-Party Data Attribution Licenses:** `.licenses/`

---

## Quickstart & Reproducibility

### 1. Environment Setup
```bash
uv sync --python /opt/homebrew/opt/python3.13
```
*(Optional for live AlphaGenome AVI re-queries: create a local `.env` with `ALPHAGENOME_API_KEY=...` — gitignored by default.)*

### 2. Verify All 5 Test Suites (Phase 0 $\rightarrow$ Stage 4)
```bash
.venv/bin/python tests/run_all_tests.py
```

### 3. Re-Run Local Analysis Pipeline from Cached Public Datasets
```bash
.venv/bin/python scripts/stage1_convergent_mining.py
.venv/bin/python scripts/stage2_3d_fbt_clustering.py
.venv/bin/python scripts/stage3_structural_safety_editor_gate.py
.venv/bin/python scripts/stage3c_bioai_cofolding_export.py
.venv/bin/python scripts/stage4_epistasis_synergy_uq.py
```
