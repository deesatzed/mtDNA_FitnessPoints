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
| **Stage 5: Hardening, 6 Proofs & Explorer** | `scripts/stage5_hardening_and_proofs.py`, `scripts/build_interactive_explorer.py` | Full-proteome ClinVar/gnomAD (`42` Path vs `935` Benign), 3D Dunbrack rotamer scans, ±2..±6 nt windows, 4 alt ODE forms + adversarial controls | Adjusted mitonuclear interface **`OR = 1.519` (`p = 0.027`)** vs non-sig Control Null; blind Tier-1→Tier-2 ClinVar **`AUROC = 0.888`**; 5/5 finalists clash-free in cryo-EM (`COX3 W116R` `3.21 Å`, `0` bystanders to `±6 nt`); `100%` supermodularity & `>=85–95%` clonal heteroplasmy threshold; standalone offline explorer (`app/mitofitness_explorer.html`) |
| **Stage 6: Boltz-2, OpenMM & 8-Gene AlphaGenome** | `scripts/stage6a_prepare_cofold_inputs.py`, `scripts/stage6c_openmm_relaxation.py`, `scripts/stage6_synthesize_and_plot.py` | Live **Boltz-2 (MIT)** on Apple Silicon MPS (6 WT/Mut dimers) + **OpenMM AMBER14+OBC2** cryo-EM relaxation + 8 live **AlphaGenome AVI** promoter scans (`984` SNVs) | Boltz-2 recovers cryo-EM at `0.85–1.06 Å` (`COX3 W116R+COX6A1` `ipTM = 0.917`, WT vs Mut `0.70 Å`); OpenMM proves 5/5 finalists `ΔRMSD ≤ 0.036 Å` (`COX3 W116R` forms 6 H-bonds) vs pathogenic `ATP6 L156R` (`+0.160 Å`); AlphaGenome 6 OXPHOS partners (`41.0` SNVs ≥15) vs 2 buffers (`3.0`), `p = 0.0357` |
| **Stage 7: Cross-Project Bridges (`Mito_Papers` + `BAA-Mito`)** | `scripts/stage7a_heteroplasmy_amplification_and_cart.py`, `scripts/stage7b_singlecell_alphagenome_buffers.py`, `scripts/stage7_synthesize_and_plot.py` | Biphasic mitophagy-selection ODE (`10,000` LHS + Sobol) + 72-h CAR-T hypoxia cycling + Live **AlphaGenome 2 `v2.dna_client`** across 6 observed single-cell types + 2 blank controls | Biphasic selection shifts `50%` → **`83.2%`** (`TwoCycle` triples safety nadir to `0.243`); 100%-human **`COX3 W116R + Human SOD2`** boosts CAR-T persistence **`6.07×`** (`0.161 → 0.979`); AlphaGenome 2 single-cell **`14.47×` Signal/Null** uncovers **`UCP2`'s `7×` cardiac-selective switch**; Explorer v2 adds `E1–E4` badges & in-app Freshman Guide |

---

## Key Reports, Artifacts & Figures

- **Design Decisions & Audit Trail:** `DESIGN_LOG.md`
- **Database Schema & SQLite Catalog:** `schema/mitofitness_schema.sql`, `data/mitofitness.db`
- **Interactive Offline Explorer & Guides:** `app/mitofitness_explorer.html` (7 tabs incl. `E1–E4` provenance & in-app Freshman Guide), `docs/EXPLORER_GUIDE.md`, `docs/PHYSICIAN_EXECUTIVE_BRIEF.md`, `docs/cross_project_bridges/`
- **Stage Reports:**
  - `data/stage1_report.md` (`data/stage1_convergent_summary.json`)
  - `data/stage2_report.md` (`data/stage2_3d_fbt_summary.json`)
  - `data/stage3_report.md` (`data/stage3_shortlist.json`)
  - `data/stage4_report.md` (`data/stage4_epistasis_summary.json`)
  - `data/stage5_hardening_report.md` (`data/stage5_hardening_proofs.json`)
  - `data/stage6_report.md` (`data/stage6_structural_bioai_summary.json`)
  - `data/stage7_report.md` (`data/stage7a_heteroplasmy_and_cart.json`, `data/stage7b_singlecell_alphagenome.json`)
- **Publication Figures:**
  - `figures/stage2_3d_fbt_equivalence_clusters.png`
  - `figures/stage4_epistasis_synergy_and_cliffs.png`
  - `figures/stage5_hardening_and_proofs.png`
  - `figures/stage6_boltz2_openmm_alphagenome.png`
  - `figures/stage7_cross_project_bridges.png`
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
