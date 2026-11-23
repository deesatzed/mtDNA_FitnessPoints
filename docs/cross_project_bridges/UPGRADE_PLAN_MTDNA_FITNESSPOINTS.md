# Cross-Project Upgrade Blueprint 1: Strengthening `mtDNA_FitnessPoints` (MitoFitness Atlas)
**Contributions from `Mito_Papers` (`mitorestorer`) and `BAA Mito Twin` (`baa_mito_twin`)**
*Prepared: 2026-11-23 | Target Repository: `/Users/o2satz/antigravity/mitorestorer/mtDNA_FinessPoints`*

---

## 1. Current State & Honest Gaps in `mtDNA_FitnessPoints`

`mtDNA_FitnessPoints` has completed **Phase 0 → Stage 6** (`7/7` automated verification suites passing), delivering:
- 44 wild extremophile mitogenomes mapped onto 4 experimental human OXPHOS cryo-EM holocomplexes (`9I4I`, `9HZL`, `9I6F`, `8H9S`), proving confounder-adjusted mitonuclear interface enrichment (`OR = 1.519, p = 0.027`) vs. a non-significant matched control null (`OR = 1.194, p = 0.136`).
- Full-proteome ClinVar/gnomAD safety calibration (`5-fold CV AUROC = 0.890`; blind Tier-1 → Frozen Tier-2 `AUROC = 0.888`) and `±6 nt` zero-bystander `TALED`/`DdCBE` screening (`6/131` ultra-clean edits, led by `COX3 W116R m.9552T>C` and `ND4L A11T m.10500G>A`).
- Live Apple Silicon MPS **Boltz-2** WT vs. Mutant co-folding (`COX3 W116R + COX6A1` `ipTM = 0.917`, `0.70 Å` dimer RMSD; `0.85–1.06 Å` Chain A vs. cryo-EM) and **OpenMM `AMBER14 + OBC2`** relaxation proving `COX3 W116R` forms a **6-bond stabilizing polar H-bond trident** (`ΔRMSD = +0.014 Å` vs. pathogenic `ATP6 L156R` `+0.160 Å`).
- Stage-4/5 acute I/R ODE proving **Driver + Nuclear Buffer positive epistasis** (`COX3 W116R + Yeast NDI1 / Human SOD2 / Ciona AOX`, `100%` supermodularity, robust across 4 model forms, `0.0000` placebo, negative toxic uncoupler).

### Three Real Bottlenecks That Our Other Two Projects Can Directly Solve
1. **The `≥85–95%` Heteroplasmy Wall (Stage 5 Proof 5B):** Bulk `TALED` base editing in human cells typically reaches `40–60%` heteroplasmy, which our own heteroplasmy sweep proved is *too low* to clear the acute ischemic ATP cliff without tedious clonal isolation.
2. **Pan-Tissue vs. Cell-Type-Specific Nuclear Regulation (Stage 6 AlphaGenome):** Our 8-gene AlphaGenome `AVI` promoter constraint law (`6` OXPHOS partners `41.0` SNVs ≥15 vs. `2` buffers `3.0`, `p = 0.0357`) used pan-tissue promoter `AVI` rather than cell-type-conditioned chromatin/expression tracks.
3. **Delivery Realism & Epistemic Labeling:** *In vivo* mitochondrial gene editing faces major delivery barriers, and our interactive Explorer (`app/mitofitness_explorer.html`) does not yet carry explicit evidence-tier provenance tags separating experimental facts from AI predictions and ODE projections.

---

## 2. Specific Upgrades Imported from `Mito_Papers` (`mitorestorer`)

### Recommendation 1A — Build Stage 7A: "Biphasic Mitophagy-Selection Amplification" of `TALED`-Edited `COX3 W116R`
- **Source Asset in `Mito_Papers`:** `simulator/sim_bioenergetic_reset.py` (the Biphasic Bioenergetic Circuit Reset ODE solving the *Biogenesis Trap* via Stage-1 selective `PINK1`/Urolithin A mitophagy → Stage-2 `NAD+`/`PGC-1α` biogenesis).
- **The Scientific Breakthrough Idea:**
  - Suppose a wet lab edits human cells with `TALED` to introduce `COX3 W116R` (`m.9552T>C`) at realistic bulk efficiency (**`45–55%` heteroplasmy**).
  - Under mild, controlled sub-lethal hypoxic preconditioning (`pO2 ~ 1.0–1.5 kPa`), wild-type human `COX3` mitochondria have higher `Km_O2` (`1.25 kPa`) and depolarize slightly *earlier* than high-affinity `COX3 W116R` mitochondria (`Km_O2 ~ 0.53 kPa`).
  - Applying `Mito_Papers`'s **Stage-1 Selective Mitophagy** (`PINK1/Parkin` + Urolithin A) during mild hypoxia selectively clears the depolarized **unedited WT** organelles, followed by **Stage-2 Biogenesis** (`PGC-1α`/NAD+) to re-expand the enriched `COX3 W116R` pool—pushing heteroplasmy from **`50%` → `≥85–90%` *without* single-cell clonal sorting**!
- **Concrete Implementation Here:**
  - Add `scripts/stage7a_heteroplasmy_amplification_ode.py` coupling `Mito_Papers`'s two-population `[W_wt, W_mut]` mitophagy/biogenesis dynamics directly into our Stage-4/5 I/R ODE, with `N = 10,000` LHS priors and a quantitative falsification threshold (minimum `ΔΨ` separation between WT and `W116R` required to cross `85%` heteroplasmy in ≤14 days).
- **Payoff:** Directly solves the #1 wet-lab critique of our Stage 5 heteroplasmy result and creates a novel, bridge-defining mechanistic hypothesis uniting both repos.

### Recommendation 1B — Add an *Ex Vivo* CAR-T Solid-Tumor Hypoxia Module
- **Source Asset in `Mito_Papers`:** `simulator/cart_bioreactor_simulator.html` and `docs/cart_bioreactor_sop.md`.
- **Why Here:** Our [PHYSICIAN_EXECUTIVE_BRIEF.md](../PHYSICIAN_EXECUTIVE_BRIEF.md) identified **ex vivo CAR-T cell engineering** as the highest-feasibility 5–10 year translational path (no systemic *in vivo* delivery needed; clonal selection to >90% heteroplasmy is routine in manufacturing bioreactors).
- **Concrete Implementation Here:**
  - Add a **sustained tumor-microenvironment hypoxia + re-oxygenation cycling** stress profile (matching `Mito_Papers`'s CAR-T bioreactor parameters) into Stage 4/5 to benchmark `COX3 W116R + Human SOD2` (100% human-sequence, zero foreign-protein immunogenicity!) specifically for CAR-T persistence.

---

## 3. Specific Upgrades Imported from `BAA Mito Twin` (`baa_mito_twin`)

### Recommendation 1C — Upgrade 8-Gene AlphaGenome `AVI` to `AlphaGenome 2 v2.dna_client` Single-Cell Conditioning
- **Source Asset in `BAA Mito Twin`:** `sync_alphagenome_baa.py` (which uses `alphagenome.models.v2.dna_client` to condition `SC_GENE_EXPRESSION`, `SC_ATAC`, and `CHIP_TF` across specific `UBERON`/`CL` cell types—boosting GWAS concordance in `BAA-Mito` Suite D from `ρ = 0.494` to **`ρ = 0.899`** with 100% negative-control specificity).
- **Concrete Implementation Here:**
  - Add `scripts/stage7b_singlecell_alphagenome_buffers.py` querying cell-type-conditioned `SC_GENE_EXPRESSION` and `SC_ATAC` across **cardiomyocytes, T-cells (CAR-T relevance), dopaminergic/glutamatergic neurons, astrocytes, and brain endothelial cells** for our 8 nuclear genes (`COX6A1`, `COX4I1`, `HIGD2A`, `NDUFA13`, `UQCRFS1`, `ATP5MJ`, `UCP2`, `SOD2`) + 2 matched closed-chromatin negative controls (mirroring `BAA-Mito`'s `rs9904290` / `rs10799812` discipline).
  - Specifically test whether `UCP2`'s localized 6-SNV promoter switch (`chr11:73982823-73982863`) is selectively open and tunable in cardiomyocytes/T-cells/neurons vs. repressed elsewhere.

### Recommendation 1D — Embed 4-Tier Epistemic Provenance (`E1`–`E4`) + Organelle-Transfer Delivery Track into Explorer
- **Source Asset in `BAA Mito Twin`:** `EPISTEMIC_TIERS` in `baa_mito_catalog.py` (`E1_EXPERIMENTAL`, `E2_FOUNDATION_MODEL`, `E3_MECHANISTIC_SIM`, `E4_UNVERIFIED_HYPOTHESIS`) and the Onraet & Zuryn (2026) intercellular mitochondrial transfer pathways (`CD38`/`LRP1`).
- **Concrete Implementation Here:**
  - Tag every column/panel in `app/mitofitness_explorer.html` with a crisp `E1 / E2 / E3 / E4` pill badge (`E1` = cryo-EM `9I6F`/ClinVar/gnomAD; `E2` = Boltz-2/AlphaFold/AlphaGenome; `E3` = OpenMM `AMBER14` & Stage-4/5 ODE; `E4` = untested *in vivo* organelle transfer / biphasic selection hypotheses).
  - Add an **Engineered Donor Mitochondria Transfer** comparison option to Tab 5 (Wet-Lab Protocol) alongside direct `TALED` editing.

---

## 4. Prioritized Execution Order & Non-Regression Rules

| Priority | Upgrade | Effort | Scientific Impact | Regression Guard |
| :---: | :--- | :---: | :---: | :--- |
| **1** | **Rec 1C:** Single-cell AlphaGenome 2 conditioning on our 8 nuclear genes + 2 negative controls | Low–Med | Very High (closes pan-tissue limitation) | Keeps all existing 7/7 tests passing |
| **2** | **Rec 1A:** Biphasic Mitophagy-Selection Heteroplasmy Amplification ODE (`50% → 90%`) | Medium | Breakthrough (solves #1 wet-lab barrier) | New `tests/test_stage7_bridges.py` |
| **3** | **Rec 1D:** `E1–E4` Epistemic Provenance badges + Tab 6 Freshman Guide inside `mitofitness_explorer.html` | Low | High (reviewer & clinician trust) | JS syntax check + offline self-containment |
| **4** | **Rec 1B:** CAR-T sustained-hypoxia stress arm for `COX3 W116R + Human SOD2` | Low | High (strongest 5–10 yr clinical bridge) | Separate scenario in Stage 7 JSON |
