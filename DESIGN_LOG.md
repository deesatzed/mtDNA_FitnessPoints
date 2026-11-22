# MitoFitness Atlas — Brainstorming & Architecture Decision Log

**Workspace:** `/Users/o2satz/antigravity/mitorestorer/mtDNA_FinessPoints`  
**Started:** 2026-11-22  
**Context Ingested:**
1. Prior series: *MitoRestorer* Papers I–IV (`Mito_Papers` conversation: Core Control & Coral NPU [Preprints.org `202611.2000/v1`], *Ex Vivo* CAR-T Bioreactor [bioRxiv], *Apis mellifera* Flight-Muscle Complex III Fungicide Stress [bioRxiv], and Nigrostriatal QSP/PDE).
2. Background: `mtDNA Alleles and Variants.txt` (human vs. nuclear vs. *Drosophila* vs. pan-eukaryote/fungal mtDNA scales, animal "frozen accident" ~16.5 kb / 37 genes, and the mitonuclear co-evolution transferability bottleneck).
3. Target Direction: `Fitness Point mtDNA Adaptations.txt` (three-way OXPHOS thermodynamic dial, motif grafting / allotopic bypass / co-adapted modules, Fitness-Beats-Truth [FBT] payoff equivalence classes & non-monotonic cliffs, 5-table relational schema, and 4-stage bioinformatics pipeline).

---

## Core Scientific Synthesis & Lessons Carried Over from *MitoRestorer*

* **Conceptual Bridge:** *MitoRestorer* Papers I–IV solved *dynamic control of fixed OXPHOS hardware* (staying below the Energy Resistance redox cliff $\mathrm{EP} = \dot{e}_R \cdot f^2 < \theta_{\mathrm{GSH}}$). This standalone project investigates *evolutionary hardware re-tuning*—mining convergent mtDNA motifs across extreme environments and pairing them with nuclear compensatory buffers (`AOX`, `NDI1`, `SOD2`, `UCP2`, `PINK1`) to shift non-monotonic survival envelopes safely.
* **Rigor Standards Retained from *Mito_Papers*:**
  1. **Zero fabricated biochemical values:** Separate directly computed molecular/structural/population metrics from literature-anchored thermodynamic archetypes, and propagate kinetic uncertainty via `LHS` + `Sobol` global sensitivity analysis.
  2. **Biophysical grounding of FBT:** Ground Hoffman's FBT *payoff equivalence classes* and *non-monotonic operating cliffs* in multi-objective thermodynamic Pareto fronts and Picard & Murugan's Energy Resistance Principle.
  3. **Real-world engineering gates:** Add both a **human safety filter** (`MITOMAP` / `gnomAD` / `ClinVar`) and a **mitochondrial base-editor editability gate** (`DdCBE`/`DddA11` `C•G→T•A` and `TALED` `A•T→G•C` codon-transition + bystander analysis on human `rCRS`).
  4. **Two-kingdom split:** Mine **metazoans/vertebrates** for structurally homologous mtDNA point-motif grafting, and **fungi/plants/tunicates** (*Ciona* `AOX`, yeast `NDI1`) for nuclear-encoded allotopic bypass buffers.

---

## Architectural Decisions

### Decision #4 — Repository Structure: **LOCKED (Standalone)**
* **Choice:** Standalone project in `/Users/o2satz/antigravity/mitorestorer/mtDNA_FinessPoints`.
* **Rationale:** Keeps this comparative mitonuclear & evolutionary fitness-point program modular and cleanly separated from the *MitoRestorer* dynamical control repository, while conceptually extending it.

### Decision #1 — Version 1 Primary Vehicle: **LOCKED (Option A — Discovery Engine + Pilot Catalog)**
* **Choice:** Science-first reproducible Python pipeline + portable database (SQLite/DuckDB with PostgreSQL-compatible DDL) built to generate verifiable, publishable computational discoveries first, with interactive UI/platform layers optional later.
* **Rationale:** Avoids building an empty software shell ahead of validated biology; mirrors the high-rigor, bioRxiv-ready workflow that succeeded in *MitoRestorer* Papers I–III.

### Decision #2 — Target Recipient Strategy: **LOCKED (Option 2A First — Human Anchor, Built for Reuse in Separate 2B Climate/Ecological Endeavor)**
* **Choice:** Keep **2A (Human)** and **2B (Climate/Ecological: corals, pollinators)** as **two separate scientific endeavors**, executing **2A first** on a species-agnostic schema and codebase.
* **Rationale:**
  * ~80–85% of the schema, wild extreme-environment donor database, FBT payoff math, and Stage 4 Driver+Buffer synergy simulator transfers directly to 2B.
  * Crucially, **2A provides ground-truth "answer keys"** (experimental human cryo-EM OXPHOS structures, `MITOMAP`, `gnomAD` mtDNA, `ClinVar`, and human `DdCBE`/`TALED` editing rules) that allow us to catch false positives and calibrate our Stage 3 structural compatibility filter *before* applying it to data-sparser ecological species in 2B.
  * For future **2B**, the calibrated catalog can support not only allotopic expression (`AOX`) but also **assisted gene flow** (identifying wild resilient coral/bee lineages already carrying convergent fitness-point alleles).

### Decision #3 — OXPHOS Gene / Complex Scope: **LOCKED (3-Tier Algorithmic Progression)**
* **Choice:** Rather than hard-coding an either/or subset of genes, execute a **gated algorithmic progression** using one unified, config-driven pipeline:
  * **Tier 1 (Calibration Core — 6 Genes across all 4 mtDNA-encoded OXPHOS complexes):** `ND1`, `ND2` (Complex I), `CYTB` (Complex III), `COX1`, `COX3` (Complex IV), `ATP6` (Complex V). Used to calibrate objective thresholds against known human adaptive positive controls (e.g., Tibetan `ND1 Y30H`, bar-headed goose `COX3`, Arctic `ATP6`) and known pathogenic negative controls (e.g., `ND1` LHON, `ATP6` Leigh/NARP).
  * **Tier 2 (Frozen Out-of-Sample Sweep — Remaining 7 Genes $\rightarrow$ All 13 Protein-Coding Genes):** Run `ND3`, `ND4`, `ND4L`, `ND5`, `ND6`, `COX2`, `ATP8` through the **exact same frozen rules** without re-tuning. Provides an unbiased out-of-sample validation test (e.g., `ND4 m.11778G>A` LHON detection) while objective per-site alignment/resolution quality gates automatically quarantine noisy columns in short/fast-evolving subunits (`ATP8`, `ND6`).
  * **Tier 3 (Future 2B Ecological Endeavor):** Re-run the calibrated engine on a non-human target recipient (e.g., reef coral or *Apis mellifera*), where the same homology/alignment gates naturally select deeply conserved loci (`COX1`, `CYTB`, `ATP6`).

---

## Concrete Build Blueprint (Standalone `mtDNA_FinessPoints`)

### 0. Project Scaffolding & Portable Database Schema
* **Environment:** Local `uv`-managed `.venv` (`biopython`, `numpy`, `scipy`, `pandas`, `scikit-learn`, `matplotlib`, `requests`) — zero global `pip` installs.
* **Database:** Portable **SQLite** database (`data/mitofitness.db`) built from **PostgreSQL-compatible DDL** (`schema/mitofitness_schema.sql`), extending the 5 tables from `Fitness Point mtDNA Adaptations.txt` with:
  * `donor_species` (NCBI RefSeq mitogenome accession, stressor axis, sister-control pair ID, taxonomic clade),
  * `human_safety_annotations` (`rCRS` nucleotide/AA position, `gnomAD` heteroplasmy/homoplasmy AF, `ClinVar` / `MITOMAP` pathogenicity classification),
  * `base_editor_feasibility` (`DdCBE_TC`, `DdCBE_DddA11`, `TALED`, bystander-edit count in the editing window, and `Not_Directly_Base_Editable` flags inside `transferability_assessments`).

### 1. Stage 1 — Paired Extreme-vs-Control Convergent Mining (Live NCBI RefSeq)
* **Phylogenetically Paired Wild Donor Panel (~40–50 complete vertebrate mitogenomes + Human `rCRS` `NC_012920.1`):**
  * **High-Altitude Hypoxia (+ Cold):** Bar-headed goose vs. Greylag goose; Yak vs. Cattle; Plateau pika vs. Rabbit; Tibetan antelope vs. lowland bovid; Alpaca/Vicuña vs. Camel; Snow leopard vs. Lion/Tiger; Yunnan snub-nosed monkey vs. Rhesus macaque.
  * **Diving & Subterranean Hypoxia / Ischemia-Reperfusion (NO Altitude Cold — key FBT deconvolution control):** Sperm whale, Cuvier's beaked whale, Weddell seal, Elephant seal vs. terrestrial/shallow controls; Naked mole-rat & Blind mole-rat vs. Guinea pig / Mouse.
  * **Arctic / Extreme Cold (Normal $\mathrm{O}_2$):** Polar bear vs. Brown/Black bear; Arctic fox vs. Red fox; Reindeer vs. White-tailed deer; Bowhead whale & Narwhal; Woolly mammoth (*Mammuthus primigenius*) vs. African/Asian elephant.
  * **Extreme Heat / Desert (Normal $\mathrm{O}_2$):** Dromedary camel, Fennec fox, Arabian oryx, Desert rodent vs. temperate sister taxa.
* **Algorithm:** Live NCBI RefSeq mitogenome fetch $\rightarrow$ CDS extraction & translation for all 13 mtDNA genes $\rightarrow$ protein MSA mapped onto Human `rCRS` coordinates $\rightarrow$ sister-pair directional substitution calling, independent clade convergence count, physicochemical radicality scoring (Grantham / volume / hydropathy / charge / polarity), and automated column quality gating (gap fraction & entropy).

### 2. Stage 2 — FBT / Energy-Resistance Payoff Equivalence Clustering
* **Tier A (Computed Molecular + Structural Feature Space):** Cluster convergent sites using `scikit-learn` (PCA/spectral embedding + `HDBSCAN`) over observable physicochemical shifts + 3D cryo-EM microenvironment features (distance to heme/Fe-S/Q-site/proton channel, solvent accessibility, mitonuclear interface contact count, evolutionary constraint).
* **Tier B (Literature-Anchored Thermodynamic Archetypes):** Anchor equivalence classes with experimentally characterized adaptations (with live-verified PMIDs/DOIs) to test whether distinct stressors (e.g., diving hypoxia, subterranean hypoxia, and high-altitude hypoxia) collapse onto shared OXPHOS structural/thermodynamic classes.

### 3. Stage 3 — 3D Human Cryo-EM Mitonuclear Compatibility + Safety + `DdCBE`/`TALED` Gate
* **Live Experimental Structures:** Query RCSB PDB live for the highest-resolution human OXPHOS cryo-EM holocomplexes (**Complex I, Complex III, Complex IV, Complex V**), map every mtDNA chain directly onto `rCRS` numbering, and inspect top candidates interactively in the 3D Biomolecular Structure Viewer.
* **Three Hard Automated Gates:**
  1. **Mitonuclear Structural Gate:** Contact map ($\le 4.5\text{ \AA}$) to all neighboring nuclear-encoded subunits, burial/active-site proximity, local steric/electrostatic clash & packing perturbation score, and cryo-EM local B-factor/resolution filter.
  2. **Human Clinical Safety Gate:** Live cross-check via `ClinVar`, `gnomAD` mtDNA, and `MITOMAP` to quarantine pathogenic/LHON/MELAS/Leigh loci (or flag `Epistatic_Rescue_Required`).
  3. **Base-Editor Feasibility Gate:** Exact human mitochondrial genetic code (NCBI Translation Table 2!) codon-transition analysis for `DdCBE` (`5'-TC` $\rightarrow$ `C•G>T•A`), `DddA11` (`HC`), and `TALED` (`A•T>G•C`), including bystander nonsynonymous edit risk within the editing window.

### 4. Stage 4 — Epistatic Driver + Nuclear Buffer Synergy & Non-Monotonic Cliff Simulation
* **Biophysical Model:** Couple the *MitoRestorer* Energy Resistance Principle ($\mathrm{EP} = \dot{e}_R \cdot f^2$) with an OXPHOS electron-flux / proton-leak / RET-ROS thermodynamic model conditioned on human physiological world states $W$ (e.g., acute ischemia-reperfusion $p\mathrm{O}_2$ drop & re-oxygenation surge, thermal/acidic stress, high-flux bioenergetic load).
* **Synergy Discovery:** Evaluate Stage-3-passing mtDNA **drivers** paired against canonical nuclear **compensatory buffers** (*Ciona* `AOX`, yeast `NDI1`, `SOD2`, `UCP2`, `PINK1`/mitophagy), quantifying **positive epistasis ($S > 1$)** and expanded survival-cliff margins with `N = 10,000` Latin Hypercube Sampling (`LHS`) and Saltelli–Sobol global sensitivity analysis (matching *MitoRestorer* Paper III rigor).

---

## AlphaFold2 (`AFDB v6`) & AlphaGenome Integration Plan (Added 2026-11-22)

* **AlphaGenome (`AVI` + Regulatory / Tissue Expression API, `ALPHAGENOME_API_KEY` configured in `.env`):**
  * *Boundary:* Strictly applies to **nuclear chromosomes (`chr1–22, X, Y`)**, not circular `chrM`. Research-use molecular/functional framing only.
  * *Stage 3 (Mitonuclear Partner Audit):* Query tissue-specific expression (heart, brain, skeletal muscle), isoform switches (`COX4I1` vs. hypoxia-inducible `COX4I2`), and `AVI` regulatory variant vulnerability for the exact nuclear OXPHOS subunits contacting candidate mtDNA sites ($\le 4.5\text{ \AA}$).
  * *Stage 4 (Endogenous Buffer Regulatory Mapping):* Run promoter/enhancer saturation mutagenesis (`region` scan) on human endogenous nuclear buffers (`SOD2`, `UCP2`, `PINK1`, `PPARGC1A`, `NFE2L2`) plus clickable AlphaGenome Atlas deep-links.
* **AlphaFold2 (`AlphaFold DB v6` — `pLDDT` + `PAE`) & AF3 / Boltz-2 Export:**
  * *Stage 3/4 (Allotopic & Nuclear Buffer Structural Audit):* Fetch AFDB v6 models + `pLDDT`/`PAE` rigid-domain boundaries for *Ciona* `AOX` (no experimental PDB structure), yeast `NDI1`, and human `UCP2`, `SOD2`, `PINK1`.
  * *Stage 3 (Disorder Guardrail):* Cross-check unresolved/low-resolution cryo-EM loops on human OXPHOS partners against AFDB v6 `pLDDT` (`< 50` = intrinsically disordered) to prevent false-positive steric clash calls.
  * *End of Stage 3 (Co-Folding Manifests):* Auto-export validated **AlphaFold 3 JSON / Boltz-2 YAML** co-folding specs for the top 3–5 Safe, Compatible, Base-Editable finalists (mutant human mtDNA chain + contacting nuclear partners + cofactors).

---

## Phase-by-Phase Checkpoint & Sanity-Check Contract

Every phase ends with a human review checkpoint before advancing:
* **Phase 0:** Local `uv` `.venv` + SQLite/Postgres-compatible DB + verified human `rCRS` (`NC_012920.1`) 13-gene translation (Pass = 13/13 exact UniProt AA length match under NCBI Table 2).
* **Stage 1:** Paired extreme-vs-control vertebrate mitogenomes (NCBI RefSeq), Tier 1 (6 calibration genes) $\rightarrow$ Tier 2 (all 13 genes frozen). Pass = recovers known benchmarks cleanly above naive frequency; auto-filters noisy columns.
* **Stage 2:** Observable physicochemical + structural clustering (`HDBSCAN`/PCA). Pass = permutation-tested cross-stressor enrichment (reported honestly whether significant or null).
* **Stage 3:** Live human cryo-EM holocomplexes + `AFDB v6` disorder/buffer audit + `ClinVar`/`gnomAD`/`MITOMAP` safety + `DdCBE`/`TALED` codon editability + `AlphaGenome` nuclear partner audit. Pass = Tier 1 catches `ND1`/`ATP6` pathogenic controls & passes adaptive controls; Tier 2 frozen catches `ND4 m.11778G>A` LHON out-of-sample.
* **Stage 4:** Coupled Energy Resistance + OXPHOS redox/RET model with `N=10,000` `LHS` & Sobol GSA + `AlphaGenome` endogenous buffer regulatory scan + quantitative *in vitro* human cybrid falsification criteria.

---

## Execution Progress

### Phase 0 — **COMPLETE & VERIFIED (2026-11-22)**
* Isolated `uv` environment built (`.venv`, Python 3.13, `biopython`, `numpy`, `scipy`, `pandas`, `scikit-learn`, `matplotlib`, `requests`).
* `.env` protected in `.gitignore`; confirmed `ALPHAGENOME_API_KEY` is present.
* Created PostgreSQL-compatible 8-table schema (`schema/mitofitness_schema.sql`) and initialized `data/mitofitness.db`.
* Fetched live Human `rCRS` (`NC_012920.1`, `16,569 bp`) and all 13 RefSeq proteins (`YP_003024026.1`–`YP_003024038.1`) via `ncbi-sequence-fetch`.
* Automated test suite (`tests/test_phase0.py`) passed 100%: 13/13 genes match canonical UniProt lengths under NCBI Table 2 (including `ATA`/`ATT` start codons and `-` strand `ND6`), and all 5 landmark rCRS coordinate tests (`m.3394T>C -> ND1 Y30H`, `m.3460G>A -> ND1 A52T`, `m.8993T>G -> ATP6 L156R`, `m.11778G>A -> ND4 R340H`, `m.14484T>C -> ND6 M64V`) verified with zero offset errors.

### Stage 1 — **COMPLETE & VERIFIED (2026-11-22)**
* Fetched **44/44 complete wild vertebrate NCBI RefSeq mitogenomes** (`13/13` proteins each = `572` wild proteins + `13` human `rCRS` = `585` total) across **23 matched Extreme-vs-Control sister pairs** (`scripts/species_panel.py`, `scripts/stage1_fetch_mitogenomes.py`).
* **Major Methodological Catch & Upgrade (`stage1_convergent_mining.py` v2):**
  * Caught two classic convergent-evolution confounders during Stage 1 audit: (a) symmetric neutral `I<->V` homoplasy at high-entropy sites, and (b) uniform/entropy-only permutation inflation (`1.34x–1.69x` for Extreme radical shifts vs. `1.39x–1.99x` for matched Control radical shifts!).
  * Proved honestly that genome-wide sequence overlap alone cannot distinguish adaptive cross-stressor convergence from 3D structural permissiveness—confirming why **3D cryo-EM structural mapping + FBT payoff clustering** (Stages 2–3) is essential.
  * Implemented vertebrate consensus polarization (`k_ext_derived` vs `k_ctl_derived`) separating **65 Class I Conserved-Site Radical Rule-Breakers** (`41` Tier-1, `24` Tier-2 Frozen; automatically recovering canonical **Bar-headed goose `COX3 W116R`**, 97.8% invariant, Grantham `101.29`) and **70 Class II Directionally Biased Convergent Dials** (`32` Tier-1, `38` Tier-2 Frozen; `k_ext >= 2`, `net_bias >= +2`).
  * Populated `256` allele-pair rows in `mt_alleles`; `tests/test_stage1.py` passed 100%.

### Stage 2 — **COMPLETE & VERIFIED (2026-11-22)**
* Queried RCSB PDB live (`pdb-database` skill) and downloaded 4 ultra-high-resolution experimental human cryo-EM OXPHOS holocomplexes:
  * Complex I: **`9I4I`** (`2.63 Å`, 44 protein subunits + `FMN/FES/SF4/CDL`)
  * Complex III: **`9HZL`** (`2.52 Å`, 10 protein subunits + `HEM/HEC/FES/U10/CDL`)
  * Complex IV: **`9I6F`** (`2.95 Å`, 14 protein subunits incl. `HIGD2A` + `HEA/CU/MG/CDL`)
  * Complex V: **`8H9S`** (`2.53 Å`, 17 protein subunits incl. c-ring + `ATP/ADP/MG/CDL`)
* Mapped **13/13 human mtDNA proteins** (`3,756/3,789` residues) and **131/131 (100%) Stage-1 polarized candidates** directly onto 3D atomic coordinates (`scripts/stage2_3d_fbt_clustering.py`, `tests/test_stage2.py` 100% PASS).
* **Two Key 3D Biological Findings:**
  1. **Statistically Significant Mitonuclear Interface Enrichment (`OR = 1.45, Fisher p = 0.044`):** `39.7%` (`52/131`) of extreme-adaptation candidates contact a human nuclear partner subunit within `4.5 Å` vs. `31.5%` background across the human OXPHOS proteome—quantitatively proving why mitonuclear compatibility (Stage 3) is the decisive engineering bottleneck.
  2. **4 Distinct FBT Structural-Physicochemical Payoff Equivalence Classes** (`CONSERVATIVE_TRANSMEMBRANE_DIAL` n=52, `HYDROPHOBIC_BULK_INTERFACE_SEAL` n=37, `POLARITY_GAIN_CHANNEL_REMODELING` n=34, `ELECTROSTATIC_CHARGE_REMODELING` n=8) populated into `fitness_payoff_profiles` (`256/256` rows) + rendered in `figures/stage2_3d_fbt_equivalence_clusters.png`. Also transparently decomposed 3D cross-hypoxia proximity (`9.71 Å` vs `12.26 Å`, `p = 0.007` in Complex I driven by 14 shared exact sites vs `15.25 Å` for strictly non-overlapping sites).





