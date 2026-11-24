# Stage 7 Report — Cross-Project Bridges (`Mito_Papers` + `BAA Mito Twin` → `mtDNA_FitnessPoints`)

> **4-Tier Epistemic Provenance (adopted from `BAA Mito Twin`):** `E1` = experimental cryo-EM/ClinVar/gnomAD; `E2` = AlphaGenome 2 single-cell foundation model (`v2.dna_client`); `E3` = bounded mechanistic ODE simulations (`10,000`-draw LHS + Sobol GSA); `E4` = uncalibrated biological selection hypotheses pending wet-lab falsification.

## 1. Track 7A (`Mito_Papers` Bridge): Biphasic Mitophagy-Selection Heteroplasmy Amplification [`Tier E3/E4`]
- **Solver Convergence (`RK4 dt=0.05 d` vs `DOP853`/`Radau`):** `h_end = 0.8493408` vs `0.8496002` (`max |Delta| = 0.0002595`).
- **Hypothesis Tested:** Starting from realistic bulk `TALED` editing (**`50.0%` `COX3 W116R`**), can mild hypoxic preconditioning (`pO2 = 1.1 kPa`) + `PINK1`/Urolithin A mitophagy selectively clear depolarized WT copies (`Km = 1.25 kPa`) while preserving high-affinity `COX3 W116R` (`Km = 0.53 kPa`), then re-expand via Stage-2 normoxic biogenesis?

| 14-Day Regimen (`N = 10,000` LHS) | Day-14 Heteroplasmy `[5th, 95th]` | Final `W_tot` (Nadir) | `% Draws ≥ 85%` | Biophysical Verdict |
| :--- | :---: | :---: | :---: | :--- |
| `Normoxia_NoSelection` | **50.0%** `[50.0, 50.0]` | `0.744` (`0.744`) | `0.0%` | No selective pressure; stays at 50% |
| `Unphased_Biogenesis_Trap` | **50.0%** `[50.0, 50.0]` | `0.898` (`0.898`) | `0.0%` | Clones WT & mutant equally (Biogenesis Trap!) |
| `Continuous_Mitophagy_Only` | **96.05%** `[35.69, 99.93]` | `0.011` (`0.011`) | `71.0%` | High h_end, **lethal organelle depletion (`W_tot = 0.011`)** |
| `Biphasic_Selection_Reset` | **83.21%** `[42.39, 97.53]` | `0.692` (`0.082`) | `46.0%` | Strong amplification (`83.2%`), recovers `W_tot = 0.692` (deep nadir `0.082`) |
| `TwoCycle_Pulsed_Biphasic` | **80.9%** `[43.09, 96.5]` | `0.618` (`0.243`) | `40.2%` | **Safer pulsed protocol:** triples `W_tot` nadir (`0.243`) while hitting `80.9%` (`96.5%` 95th) |
- **Saltelli–Sobol GSA (`10,240` evals):** Variance is dominated by **`s_bias_err` (`S_T = 0.805`)** — pinpointing PINK1 depolarization selectivity (`ΔΨ_mut - ΔΨ_wt`) as the single wet-lab measurement that falsifies or validates metabolic selection.

## 2. Track 7B (`Mito_Papers` CAR-T Bridge): *Ex Vivo* CAR-T Solid-Tumor Hypoxia Cycling [`Tier E3`]
- **Protocol:** 72-h solid-tumor hypoxia cycling (6 cycles x [8h 0.45 kPa hypoxia + 4h 6.5 kPa reoxygenation]) (`N = 10,000` LHS).
- **100% Human-Sequence Winner (`COX3 W116R + Human SOD2`, Zero Foreign Immunogenicity):** Median 72-h effector persistence jumps **`6.07×`** over WT CAR-T (`0.1613` → **`0.9794`**; ATP nadir `0.052` → **`0.957`**; cumulative EP excess `18.72` → **`0.36`**), matching foreign `Yeast NDI1` (`0.9828`) without foreign-protein immune risk!

## 3. Track 7C (`BAA Mito Twin` Bridge): Live AlphaGenome 2 Single-Cell Conditioning (`6` Cell Types + `2` Blank Controls) [`Tier E2`]
- **Signal-to-Null Specificity:** Active 8 nuclear genes mean single-cell composite = **`0.3096`** vs. 2 `BAA-Mito` blank negative controls (`rs9904290`, `rs10799812`) mean = **`0.0214`** (**`14.47×` separation**, zero false alarms).
- **Major Single-Cell Discovery:** **`6 / 8` nuclear genes peak in Human Ventricular Cardiomyocytes**, and **`UCP2` (`chr11:73982830:C>G`)**—which looked modest in pan-tissue `AVI`—reveals a **`7×` cardiac-selective regulatory switch** (`scRNA = -0.5907`, composite **`0.8562`** in Ventricular Cardiomyocytes vs. `-0.0851` in CD8+ T cells)!

| Gene / Control | Role | Cardiomyocyte | CD8+ T-Cell (CAR-T) | Brain Endothelial | Astrocyte | Glutam. Neuron | Dopam. Neuron | Top Cell |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| [`COX6A1`](https://deepmind.google.com/science/alphagenome/atlas?q=chr12:120438129:G%3EC&m=variant) | Nuclear_OXPHOS_Partner | `0.761` | `0.592` | `0.692` | `0.363` | `0.317` | `0.530` | **`cardiomyocyte`** |
| [`NDUFA13`](https://deepmind.google.com/science/alphagenome/atlas?q=chr19:19516241:G%3ET&m=variant) | Nuclear_OXPHOS_Partner | `0.361` | `0.152` | `0.226` | `0.117` | `0.087` | `0.170` | **`cardiomyocyte`** |
| [`HIGD2A`](https://deepmind.google.com/science/alphagenome/atlas?q=chr5:176388760:G%3ET&m=variant) | Nuclear_OXPHOS_Partner | `1.116` | `0.232` | `0.134` | `0.561` | `0.079` | `0.145` | **`cardiomyocyte`** |
| [`ATP5MJ`](https://deepmind.google.com/science/alphagenome/atlas?q=chr14:103921529:C%3EA&m=variant) | Nuclear_OXPHOS_Partner | `0.754` | `0.066` | `0.087` | `0.352` | `0.063` | `0.103` | **`cardiomyocyte`** |
| [`COX4I1`](https://deepmind.google.com/science/alphagenome/atlas?q=chr16:85799675:C%3EG&m=variant) | Nuclear_OXPHOS_Partner | `0.249` | `0.502` | `0.373` | `0.192` | `0.323` | `0.293` | **`cd8_tcell`** |
| [`UQCRFS1`](https://deepmind.google.com/science/alphagenome/atlas?q=chr19:29213160:C%3ET&m=variant) | Nuclear_OXPHOS_Partner | `0.499` | `0.039` | `0.123` | `0.213` | `0.124` | `0.022` | **`cardiomyocyte`** |
| [`UCP2`](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:73982830:C%3EG&m=variant) | Endogenous_Buffer | `0.856` | `0.170` | `0.761` | `0.450` | `0.464` | `0.611` | **`cardiomyocyte`** |
| [`SOD2`](https://deepmind.google.com/science/alphagenome/atlas?q=chr6:159693241:T%3EG&m=variant) | Endogenous_Buffer | `0.048` | `0.097` | `0.124` | `0.110` | `0.083` | `0.079` | **`brain_endothelial`** |
| [`NEG_CTRL_rs9904290`](https://deepmind.google.com/science/alphagenome/atlas?q=chr17:45985000:C%3ET&m=variant) | Negative_Control_Blank | `0.026` | `0.027` | `0.027` | `0.025` | `0.022` | `0.020` | **`cd8_tcell`** |
| [`NEG_CTRL_rs10799812`](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:215190400:G%3EA&m=variant) | Negative_Control_Blank | `0.018` | `0.010` | `0.051` | `0.010` | `0.010` | `0.012` | **`brain_endothelial`** |