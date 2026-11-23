# Cross-Project Upgrade Blueprint 3: Strengthening `BAA Mito Twin` (`baa_mito_twin`)
**Contributions from `mtDNA_FitnessPoints` (MitoFitness Atlas) and `Mito_Papers` (`mitorestorer`)**
*Prepared: 2026-11-23 | Target Repository: `/Users/o2satz/antigravity/baa_mito_twin`*

---

## 0. Mandatory Preprint Hash Guardrail

> **CRITICAL:** In `BAA Mito Twin`, `preprint/manuscript.tex` has an active Pass-1 audit (`preprint/.preprint_gate/audit_pass1.json`) awaiting `preprint-blind-auditor` and Pass-2 `reverify`.
> Decide explicitly before touching `preprint/` whether to:
> - **(Path A — Fast v1 Seal):** Freeze `preprint/` as-is, run blind review + Pass-2 seal for v1, and implement all upgrades below in **v2 / Phase 2**; **or**
> - **(Path B — Supercharged v1):** Fold Priority #1 & #2 below (real `chrM` alleles + Boltz-2/OpenMM on `APOE ε4` / `TFAM S12T`) into the codebase and manuscript *before* running `preprint-blind-auditor` (remembering to cite Boltz-2, OpenMM, PDB, ClinVar, and gnomAD in `preprint/references.bib` per `G2_SKILL_LICENSES`).

---

## 1. Current Strengths & Specific Gaps in `BAA Mito Twin`

`BAA Mito Twin` has achieved remarkable multiscale regulatory and epistemic rigor (`12/12` `BAA-GuardBench`, `1ba715f`):
- **AlphaGenome 2 (`v2.dna_client`) single-cell conditioning** across 6 brain/vascular cell types (`brain_endothelial`, `pericyte`, `astrocyte`, `microglia`, `glutamatergic`, `dopaminergic`) + GTEx v10 (`ρ = 0.899, p = 0.0022` vs. GWAS, `100%` specificity on 2 blank controls vs. `ρ = 0.494` naive pan-tissue).
- Coupled 4-cell intercellular mitochondrial transfer ODE (`Ages 30 → 80`) discovering `APOE × MAOB` positive epistasis (`+1.92 yrs`), `FOXO3` buffering (`-0.42 yrs`), the `12.25%` metabolic/vascular `chrM` heteroplasmy tipping point, and the early-intervention window (`92.7%` AUC reduction at age 30 vs. `12.9%` at age 70).
- Gold-standard **4-Tier Epistemic Governance (`E1`–`E4`)** quarantining unverified links (`Tier E4`).

### Three Real Gaps That `mtDNA_FitnessPoints` & `Mito_Papers` Directly Solve
1. **Stylized `chrM` Background in Experiment 4 / Suite B:** `BAA Mito Twin` explicitly set out to fix GWAS's omission of `chrM`, yet its `Nuclear × chrM` selector currently uses only **3 abstract placeholders** (`Haplogroup J/T`, `Haplogroup H`, `28% Heteroplasmy`)—with no real `chrM` sequences, ClinVar/gnomAD allele frequencies, or 3D OXPHOS structural scores.
2. **Missing 3D Atomic Structural Layer for Coding Missense Loci:** `BAA Mito Twin`'s own Phase-2 backlog flagged 3D structural validation of **`APOE ε4 (rs429358, p.Cys130Arg)`**, **`TFAM (rs1937, p.Ser12Thr)`**, and the **`ARF1` K-lactylation interface** as unbuilt.
3. **Continuous Rescue Sliders Ignore the "Biogenesis Trap" & Lack Full Sobol GSA:** Suite C tests continuous rescues (`CD38 / NAD+ boost`, `OPTN transmitophagy`, `MAO-B inhibition`), without testing whether **chronologically phasing** clearance before `CD38/NAD+` biogenesis rescues the steep age-60/70 drop-off (`12.9%` AUC reduction), and without a `10,000`-draw LHS + Saltelli–Sobol global sensitivity decomposition.

---

## 2. Specific Upgrades Imported from `mtDNA_FitnessPoints`

### Recommendation 3A — Replace Stylized `chrM` Placeholders with Real Calibrated Human `chrM` Alleles from `mitofitness.db`
- **Source Asset in `mtDNA_FitnessPoints`:** `data/mitofitness.db` (18,164 live `gnomAD v3 chrM` variants + 513 ClinVar Pathogenic + 1,891 ClinVar Benign + 131 3D-mapped extremophile candidates + blind Tier-1→Tier-2 pathogenicity model `AUROC = 0.888`) and our Stage-4 non-monotonic RET-ROS cliff (`EP = eR * f^2`).
- **Concrete Application to `BAA Mito Twin`:**
  - Upgrade `baa_mito_catalog.py`'s `chrM` selector (`Experiment 4` & `Suite B`) from abstract labels to **6 concrete, experimentally grounded human `chrM` alleles** with exact `rCRS` coordinates, `gnomAD` homoplasmic allele frequencies, and 3D cryo-EM distances:
    1. **Mildly Uncoupled / Adaptive Human Polymorphisms:** `ND1 Y30H` (`m.3394T>C`, `gnomAD AF_hom = 0.0091`, Tibetan high-altitude) and `ATP6 I106V` (`m.8842A>G`, `AF_hom = 0.0023`, contacts `ATP5MJ`).
    2. **Coupled Reference Baseline:** Human `rCRS` (`NC_012920.1`).
    3. **Hypoxia-Gain / High-RET-Liability Test Variant:** `COX3 W116R` (`m.9552T>C`) — directly testing in `BAA-Mito` whether an unbuffered high-RET `chrM` variant lowers or raises the `12.25%` vascular/T2D tipping point unless paired with `SOD2`/`FOXO3`!
    4. **Sub-Clinical Pathogenic / High-`eR` Controls:** `ND4 R340H` (`m.11778G>A`, LHON) and `ATP6 L156R` (`m.8993T>G`, Leigh/NARP) at sub-threshold heteroplasmy (`5%–30%`).
  - **Payoff:** Turns `BAA Mito Twin`'s `Nuclear × chrM` module from a conceptual placeholder into a genuinely sequence-resolved, ClinVar/gnomAD-calibrated mitonuclear epistasis engine!

### Recommendation 3B — Run Live Apple Silicon **Boltz-2** + **OpenMM `AMBER14`** on `BAA-Mito`'s Coding Missense & `ARF1` Lactylation Targets
- **Source Asset in `mtDNA_FitnessPoints`:** Our working Apple Silicon Metal (`mps`) **Boltz-2** workflow + `scripts/stage6c_openmm_relaxation.py` (`PDBFixer` + `OpenMM AMBER14 + OBC2`).
- **Concrete Application to `BAA Mito Twin`:**
  - Execute real WT vs. Mutant **OpenMM `AMBER14`** energy relaxation and **Boltz-2** co-folding on:
    1. **`APOE ε3` vs. `APOE ε4` (`rs429358`, `p.Cys130Arg`)** — quantifying exact `Arg130–Glu109` salt-bridge remodeling and domain-interaction pocket $\Delta\text{RMSD}$ (grounding `APOE`'s `E1/E2` structural link in `BAA-Mito`).
    2. **`TFAM` WT vs. `p.Ser12Thr` (`rs1937`)** — auditing mitochondrial targeting sequence / DNA-binding stability.
    3. **Astrocytic `ARF1` + `LRP1` interface** — providing 3D structural grounding (`ipTM`, `PAE`) for the lactate/ARF1-lactylation pathway!

---

## 3. Specific Upgrades Imported from `Mito_Papers` (`mitorestorer`)

### Recommendation 3C — Test **Biphasic Phased Rescue (Stage 1 `OPTN/STX17` Clearance → Stage 2 `CD38/NAD+` Donation)** to Rescue the Late-Life (Age 60–70) Cliff in Suite C!
- **Source Asset in `Mito_Papers`:** The core **Biphasic BCR principle** (`sim_bioenergetic_reset.py`): activating `NAD+`/biogenesis continuously without first clearing damaged organelles triggers the *Biogenesis Trap*.
- **Why This Is Huge for `BAA Mito Twin`:**
  - In `BAA Mito Twin` Suite C, continuous rescue at Age 30 cuts lifetime synaptic damage by **`92.7%`**, but waiting until Age 70 cuts it by only **`12.9%`**.
  - **Testable Hypothesis in `intercellular_mito_simulator.py`:** At Age 60–70, `D_axon` (damaged stationary axonal mitochondria) and `M_sen` (senescent microglial migrasomes) are already elevated. Turning on continuous astrocytic `CD38/NAD+` boost simultaneously with high `D_axon` wastes donor organelles inside a toxic `MAO-B+` / DAMP microenvironment.
  - Add a **Suite C2 ("Phased vs. Simultaneous Late-Life Rescue")** experiment comparing:
    - *Simultaneous continuous* (`CD38` + `OPTN` + `LRP1` from Age 60), vs.
    - *Biphasic staged* (**Phase 1 Ages 60–65:** `OPTN/STX17` transmitophagy + `LRP1`/lactate clearance + `MAO-B` clamp *first* → **Phase 2 Ages 65–80:** `CD38/NAD+` astrocytic mitochondrial donation once `D_axon` drops below a biomarker gate!).
  - Check quantitatively whether biphasic phasing meaningfully rescues the late-life `12.9%` floor!

### Recommendation 3D — Add `N = 10,000` LHS Uncertainty Quantification + Saltelli–Sobol Global Sensitivity Analysis to `BAA-GuardBench`
- **Source Asset in `Mito_Papers` (Paper III) & `mtDNA_FitnessPoints` (Stage 4):** `run_apis_bcr_uq_sobol.py` & `stage4_epistasis_synergy_uq.py` (`scipy.stats.qmc.LatinHypercube` + `Saltelli–Sobol` $S_1 / S_T$ + `RK4` vs. `DOP853`/`Radau` stiff solver parity).
- **Concrete Application to `BAA Mito Twin`:**
  - Add a **`Suite E: Global Sensitivity & Solver-Convergence Audit`** to `run_extended_experiments.py` running `10,000`-draw LHS (`±30%` priors on transfer rates) and Saltelli–Sobol decomposition on `BAA-Mito`'s 4-cell ODE—explicitly ranking which of the 10 Onraet & Zuryn (2026) intercellular fluxes (`OPTN` transmitophagy vs. `CD38` donation vs. `LRP1` lactylation vs. `SH3PXD2A` migrasomes) drives the largest share of Brain-Age Acceleration variance!

---

## 4. Prioritized Execution Order for `BAA Mito Twin`

| Priority | Upgrade | Effort | Scientific Impact | Preprint Seal Timing |
| :---: | :--- | :---: | :--- | :--- |
| **1** | **Rec 3D:** `10,000` LHS + Saltelli–Sobol GSA (`Suite E`) on the 10 Onraet & Zuryn transfer fluxes | Low (1 script) | Very High (ranks exact wet-lab priorities) | Great for v1 or v2 |
| **2** | **Rec 3A:** Plug real ClinVar/gnomAD-calibrated human `chrM` alleles (`ND1 Y30H`, `ATP6 I106V`, `COX3 W116R`, `ND4 R340H`) into Experiment 4 / Suite B | Low–Med | Very High (closes `chrM` placeholder gap) | High-value v2 / supercharged v1 |
| **3** | **Rec 3C:** Suite C2 Phased (`OPTN/LRP1` first → `CD38/NAD+` second) late-life rescue test | Low | High (directly bridges `Mito_Papers` BCR to brain aging) | Great for v2 |
| **4** | **Rec 3B:** Live Apple Silicon Boltz-2 + OpenMM `AMBER14` on `APOE ε4` (`C130R`), `TFAM` (`S12T`), & `ARF1` | Med (~10 min MPS) | High (adds real atomic 3D layer to `BAA-Mito`) | Phase 2 / v2 |
