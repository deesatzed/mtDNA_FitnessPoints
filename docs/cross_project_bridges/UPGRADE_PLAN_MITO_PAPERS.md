# Cross-Project Upgrade Blueprint 2: Strengthening `Mito_Papers` (`mitorestorer`)
**Contributions from `mtDNA_FitnessPoints` (MitoFitness Atlas) and `BAA Mito Twin` (`baa_mito_twin`)**
*Prepared: 2026-11-23 | Target Repository: `/Users/o2satz/antigravity/mitorestorer`*

---

## 0. Mandatory Preprint Seal Guardrail

> **CRITICAL:** In `Mito_Papers`, **Paper III** (`docs/arxiv_submission/paper3_apis_pollinator_smart_hive.tex`, co-authored by Wayne A. Satz, MD and Matthew Wong, MD, MPH, FAAEM) is currently **`SEALED_FOR_SUBMISSION`** (`seal_id: 41cdb7649ea580e57c24ac4a` for `bioRxiv`).
> **Do NOT modify any sealed v1 manuscript file until after `bioRxiv` submission** (or unless you intentionally choose to re-run `preprint-blind-auditor` + Pass-2 `reverify`). All upgrades below are cleanly scoped for **Paper III v2 / journal revision (e.g., *PLOS Computational Biology* / *Apidologie*), Paper IV (CAR-T), Paper V (MitoFitness Atlas), and the live interactive simulators**.

---

## 1. Current Strengths & Specific Gaps in `Mito_Papers`

`Mito_Papers` excels at **organelle-level lifecycle scheduling** ("software"):
- 10-state Biphasic Bioenergetic Circuit Reset (`BCR`) ODEs solving the *Biogenesis Trap* across Parkinson's, neuropsychiatric/metabolic cohorts, CAR-T bioreactors, and *Apis mellifera* honeybees.
- 10-failed-clinical-trial audit (`STEADY-PD III`, `SURE-PD3`, `NET-PD Creatine`, etc.) + `10,000`-draw LHS & `20,480`-run Saltelli–Sobol GSA.

### Three Gaps That `mtDNA_FitnessPoints` & `BAA Mito Twin` Directly Fill
1. **Static Structural Homology Only:** Paper III relies on pairwise sequence alignment of *Apis mellifera* `CYTB` (`YGH60699.1`) against bovine `PDB: 1SQB` (`50.1%` global identity, `9/9` conserved $Q_o$ pocket) plus a static AlphaFold monomer (`Q868N5`)—without **direct 3D ligand co-folding** or **physics-based atomic energy minimization**.
2. **Intracellular-Only Modeling in Neurodegeneration:** The Parkinson's & Alzheimer's simulators (`parkinsons_simulator.html`, `alzheimers_simulator.html`) model intracellular mitophagy/biogenesis inside a single representative neuron, without **astrocyte-to-neuron mitochondrial transfer** or **extracellular `MAO-B+` mitovesicle toxicity**.
3. **Software-Only Interventions (Fixed WT Hardware):** In hypoxic solid-tumor CAR-T (`cart_bioreactor_simulator.html`) and severe Complex I/IV failure, scheduling drugs/nutrients alone hits a ceiling set by native human WT OXPHOS kinetics (`Km_O2`, RET propensity).

---

## 2. Specific Upgrades Imported from `mtDNA_FitnessPoints`

### Recommendation 2A — Live Apple Silicon **Boltz-2** Ligand Co-Folding + **OpenMM `AMBER14`** for Paper III v2 (*Apis* `CYTB` + Azoxystrobin) & Paper I/II (`PINK1`/`SIRT1`/`GBA1`)
- **Source Asset in `mtDNA_FitnessPoints`:** Proven local **Apple Silicon Metal (`mps`) Boltz-2** batch pipeline + `scripts/stage6c_openmm_relaxation.py` (`PDBFixer` + `OpenMM AMBER14 + OBC2` implicit solvent).
- **Concrete Application to `Mito_Papers`:**
  1. **Paper III v2 / Journal Revision (*Apis mellifera*):** Directly **co-fold *Apis mellifera* `CYTB` (`YGH60699.1`) + `azoxystrobin` (SMILES) + `HEM` cofactors** in Boltz-2 (MIT License), compute exact `ligand_iptm`, `pLDDT`, and 3D protein-ligand interaction fingerprints (`compute_plif.py`), and align directly onto experimental bovine `PDB: 1SQB` (`AZO` pose RMSD!). This upgrades a sequence-homology argument into direct 3D atomic ligand-binding validation.
  2. **Paper I/II (Parkinson's `GBA1`/`LRRK2`/`PINK1`/`SIRT1`):** Run our OpenMM `AMBER14 + OBC2` WT vs. Mutant pocket relaxation script on the exact human PD variants (`GBA1 N370S/L444P`, `LRRK2 G2019S`) and quantify exact $\Delta\text{RMSD}$ and H-bond remodeling.

### Recommendation 2B — "Hardware + Software" Hybrid Arm in `cart_bioreactor_simulator.html` & `parkinsons_simulator.html`
- **Source Asset in `mtDNA_FitnessPoints`:** Calibrated kinetic shift parameters for **`COX3 W116R` (`m.9552T>C`, `0` bystanders to `±6 nt`) + `Human SOD2`** and our verified Stage-4 non-monotonic RET-ROS cliff.
- **Concrete Application to `Mito_Papers`:**
  - Add a **5th Regimen Arm ("Hardware-Upgraded `COX3 W116R + SOD2` + Phased BCR")** to `cart_bioreactor_simulator.html` (and as an exploratory research arm in `bcr_simulator.html`).
  - Quantitatively test whether installing a hypoxia-resilient `COX3` hardware upgrade + endogenous `SOD2` buffer expands the CAR-T survival envelope inside deep solid-tumor hypoxia (`0.5% O2`) beyond what nutrient/metabolic BCR scheduling alone can achieve.

### Recommendation 2C — Import `mtDNA_FitnessPoints`'s 6-Proof Adversarial Hardening Suite
- **Source Asset in `mtDNA_FitnessPoints`:** `scripts/stage5_hardening_and_proofs.py` (multi-model-form robustness testing across alternative ROS laws, adversarial `Null_Placebo = 0.0000` and `Toxic_FCCP < 0` negative controls, and physical invariant violation counters).
- **Concrete Application to `Mito_Papers`:**
  - Add an automated **Adversarial Placebo + Toxic-Intervention + Alternative-ROS-Law** test module to `mitorestorer/tests/run_all_verifications.py` so every ODE simulator in `Mito_Papers` proves zero false-positive drift and model-form invariance.

---

## 3. Specific Upgrades Imported from `BAA Mito Twin` (`baa_mito_twin`)

### Recommendation 2D — Couple Glial-Neuronal Intercellular Rescue (`Onraet & Zuryn 2026`) into Parkinson's & Alzheimer's Simulators
- **Source Asset in `BAA Mito Twin`:** `intercellular_mito_simulator.py` (coupled `Endothelial/Pericyte ↔ Astrocyte ↔ Microglia/BAM ↔ Neuron` organelle transfer ODEs, including `CD38`/`LRP1` astrocytic mitochondrial donation and `MAO-B+` synaptotoxic mitovesicles).
- **Why It Matters for `Mito_Papers`:**
  - In `Mito_Papers`'s Parkinson's audit, `MAO-B` inhibitors (selegiline/rasagiline) and single-pathway drugs fail to halt progression. `BAA Mito Twin` Suite C provides the exact multicellular explanation: **blocking extracellular `MAO-B+` mitovesicles alone cuts lifetime synaptic damage by only `4.6%–12.5%` because upstream astrocytic mitochondrial donation remains collapsed!**
  - Adding an optional **"Multicellular Glial Transfer Coupling" toggle** to `parkinsons_simulator.html` (and Paper II / future neurology papers) shows how Stage-1/Stage-2 BCR interacts with astrocytic `CD38`/`LRP1` rescue and explains why early-stage intervention vastly outperforms late-stage symptomatic `MAO-B` blockade.

### Recommendation 2E — Adopt `BAA-Mito`'s 4-Tier Epistemic Governance (`E1`–`E4`) & Single-Cell AlphaGenome 2 Links
- **Source Asset in `BAA Mito Twin`:** `EPISTEMIC_TIERS` (`E1`–`E4`) + `sync_alphagenome_baa.py` (`alphagenome.models.v2.dna_client` across 6 brain cell types).
- **Concrete Application to `Mito_Papers`:**
  - Replace any remaining static/generic gene links in `parkinsons_simulator.html` (`GBA1`, `EPAS1`, `TET2`, `LRRK2`, `PINK1`) with **cell-type-conditioned AlphaGenome 2 dopaminergic-neuron / microglia / astrocyte predictions** and explicit `E1 / E2 / E3 / E4` provenance badges (matching Paper III's strict separation of calibrated vs. uncalibrated hypotheses).

---

## 4. Prioritized Execution Order for `Mito_Papers`

| Priority | Upgrade | Target Scope | Effort | Impact |
| :---: | :--- | :--- | :---: | :--- |
| **1** | **Rec 2A (Part 1):** Boltz-2 *Apis mellifera* `CYTB` + `azoxystrobin` 3D ligand co-folding | Paper III journal revision (*PLOS Comp Biol* / *Apidologie*) | Low (5 min on MPS!) | Huge structural upgrade |
| **2** | **Rec 2B:** Add `COX3 W116R + Human SOD2` hardware arm to `cart_bioreactor_simulator.html` | Paper IV (CAR-T) & Web App | Low–Med | Unifies hardware + software |
| **3** | **Rec 2C:** Add adversarial Placebo/Toxic + alt-ROS-law checks to `run_all_verifications.py` | All `Mito_Papers` ODEs | Low | Bulletproofs peer review |
| **4** | **Rec 2D + 2E:** Add glial-transfer toggle & `E1–E4` badges to `parkinsons_simulator.html` | Next Web UI / Neurology paper | Medium | Bridges intracellular ↔ tissue scale |
