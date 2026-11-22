# Stage 5 Report — Scientific Hardening, Adversarial Stress-Tests & Formal Proofs

This audit directly stress-tests every statistical, 3D structural, genetic, base-editing, and biophysical claim across Stages 1–4 against negative controls, confounders, alternative model formulations, and adversarial null inputs.

## Proof 1: Confounder-Controlled & Negative-Control Mitonuclear Interface Enrichment
- **Extreme-Polarized Candidates (`n = 130`):** `39.23%` at mitonuclear interface (`<= 4.5 Å`) vs `31.18%` background (**`OR = 1.425`**, Fisher one-sided `p = 0.0341`).
- **Matched Control-Branch Negative Null (`n = 360`):** `35.0%` at interface (**`OR = 1.194`**, `p = 0.13562`).
- **Gene-Stratified Cochran–Mantel–Haenszel (CMH) Test:** `OR_CMH = 1.467` (`p = 0.0482`).
- **3D Burial × Vertebrate Entropy Stratified Permutation (`N = 5,000`):** Observed = `51` vs Expected = `42.78` (`1.192x`, empirical `p = 0.06679`).
- **Multivariate Logistic Regression (controlling for 3D burial, Shannon entropy & 12 gene fixed effects):** Adjusted **`OR = 1.519`** (`95% CI: [1.049, 2.199]`, Wald `p = 0.02696`).
- **Leave-One-Gene-Out Jackknife Stability:** `OR` remains in `[1.356, 1.502]` across all 13 genes.

## Proof 2: Full-Proteome Human ClinVar Pathogenic vs. Benign/gnomAD Benchmark & OOS Transfer
- **Dataset (100% live ClinVar + gnomAD v3 chrM mapped onto 3D cryo-EM):** `42` Pathogenic/LP vs `935` Benign/Common-gnomAD missense variants across all 13 human mtDNA genes.
- **Zero-Leakage 5-Fold CV AUROC (comparative + 3D features only):** **`0.8899`** (folds: `[0.8764, 0.9335, 0.8656, 0.8543, 0.9291]`).
- **Strict Out-of-Sample Tier-1 (6 genes) → Frozen Tier-2 (7 genes) Transfer AUROC:** **`0.8882`** (reverse Tier-2 → Tier-1: **`0.9378`**).
- **Predicted Pathogenicity Calibration:** Pathogenic median = `0.1769`, Benign median = `0.0075`, Wild 131 candidates median = `0.0551`.
- **Top-5 Finalist Risk Scores:** `{'COX3_W116R': 0.4294, 'ND4L_A11T': 0.2956, 'COX3_I165T': 0.026, 'ATP6_I106V': 0.0132, 'CYTB_M53T': 0.0582}` vs Pathogenic Controls `{'ATP6_L156R': 0.4979, 'ND1_A52T': 0.0233, 'ND4_R340H': 0.276, 'ND6_M64V': 0.0171, 'ND1_Y30H': 0.1662}`.

## Proof 3: 3D Cryo-EM Atomic Neighborhood & Rotamer Steric-Clash Audit

| Variant | Role | PDB | `dV (Å³)` | Rotamers Clash-Free | Best `d_min (Å)` | Severe Clashes (`<2.2 Å`) | Verdict | Key 5 Å Lipid / Nuclear Contacts |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- |
| `COX3_W116R` | Stage3_Finalist_#1_BarHeadedGoose | `9I6F` | `-46.0` | `13/81` | `3.21` | `0` | **CLASH_FREE_CAVITY_FIT** | `mtDNA core` |
| `ND4L_A11T` | Stage3_Finalist_#2_BirdMammalHypoxia | `9I4I` | `+30.0` | `5/9` | `2.33` | `0` | **CLASH_FREE_CAVITY_FIT** | `mtDNA core` |
| `COX3_I165T` | Stage3_Finalist_#3_SnowLeopardNMR | `9I6F` | `-50.0` | `9/9` | `3.61` | `0` | **CLASH_FREE_CAVITY_FIT** | `COX6A1(G):PHE43` |
| `ATP6_I106V` | Stage3_Finalist_#4_3WayHypoxia | `8H9S` | `-27.0` | `9/9` | `3.56` | `0` | **CLASH_FREE_CAVITY_FIT** | `ATP5MJ(P):MET31, ATP5MJ(P):GLY27` |
| `CYTB_M53T` | Stage3_Finalist_#5_UQCRFS1_Interface | `9HZL` | `-44.0` | `6/9` | `3.75` | `0` | **CLASH_FREE_CAVITY_FIT** | `UQCRFS1(C):MET140` |
| `ATP6_L156R` | Pathogenic_Control_Leigh_NARP | `8H9S` | `+13.0` | `4/81` | `2.26` | `0` | **CLASH_FREE_CAVITY_FIT** | `mtDNA core` |
| `ND1_A52T` | Pathogenic_Control_LHON | `9I4I` | `+30.0` | `9/9` | `3.4` | `0` | **CLASH_FREE_CAVITY_FIT** | `NDUFS7(E):PEE301` |

## Proof 4: Base-Editor Window Stress Test (`±2 nt` to `±6 nt` on Human `rCRS`)
- **Zero-Bystander Survivors across Windows (`n = 42` single-step editable):** `{'pm2nt': 31, 'pm3nt': 13, 'pm4nt': 12, 'pm5nt': 9, 'pm6nt': 6}`.

| Finalist | Editor | rCRS Edit | `±2 nt` | `±3 nt` | `±4 nt` | `±5 nt` | `±6 nt` |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `ND4L_A11T` | `DdCBE_DddA11_HC` | `m.10500G>A` | `0` | `0` | `0` | `0` | `0` |
| `COX3_I165T` | `TALED_AT_GC` | `m.9700T>C` | `0` | `1` | `1` | `1` | `3` |
| `ATP6_I106V` | `TALED_AT_GC` | `m.8842A>G` | `0` | `0` | `0` | `0` | `1` |
| `COX3_W116R` | `TALED_AT_GC` | `m.9552T>C` | `0` | `0` | `0` | `0` | `0` |
| `CYTB_M53T` | `TALED_AT_GC` | `m.14904T>C` | `0` | `0` | `0` | `1` | `2` |

## Proof 5: Stage-4 Formal Supermodularity, Alternative Model Forms, Adversarial Nulls & Heteroplasmy Thresholds
- **Supermodularity Theorem Concordance:** `100.0%` of dual-cliff operating points satisfy $\Delta^2 V > 0$ (`max Δ²V = +0.3853`).
- **Physical Invariant Audit:** **`0`** violations across all trajectories (bounds `0 <= W, Q_red <= 1`, `Succ, ATP, EP >= 0`, zero `NaN/Inf`).

| Model Form / Stress Protocol | `Yeast NDI1` (`E_add`, `%+`) | `Human SOD2` (`E_add`, `%+`) | `Ciona AOX` (`E_add`, `%+`) | `Null Placebo` (Adversarial) | `Toxic FCCP` (Adversarial) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `Baseline_20mIsch_40mRep_LinearROS` | `+0.1117` (`92.7%`) | `+0.1053` (`95.1%`) | `+0.0814` (`92.7%`) | `0.0` | `-0.1021` |
| `AltModel_M1_Hill_Cooperative_ROS` | `+0.103` (`92.0%`) | `+0.1035` (`96.8%`) | `+0.0782` (`91.7%`) | `0.0` | `-0.0994` |
| `AltProtocol_M2_Prolonged_30mIsch_60mRep` | `+0.0931` (`90.3%`) | `+0.0826` (`91.2%`) | `+0.0616` (`87.7%`) | `0.0` | `-0.0512` |
| `AltProtocol_M3_Short_10mIsch_40mRep` | `+0.1302` (`97.4%`) | `+0.1401` (`100.0%`) | `+0.1025` (`97.8%`) | `0.0` | `-0.156` |

### Wet-Lab Minimum Editing Heteroplasmy Thresholds (`COX3 W116R` `TALED`)
- **`Yeast_NDI1`:** `E_add >= +0.05` achieved at **`90%`** heteroplasmy; `Viability >= 0.25` achieved at **`95%`** heteroplasmy.
- **`Human_SOD2`:** `E_add >= +0.05` achieved at **`85%`** heteroplasmy; `Viability >= 0.25` achieved at **`95%`** heteroplasmy.
- **`Ciona_AOX`:** `E_add >= +0.05` achieved at **`95%`** heteroplasmy; `Viability >= 0.25` achieved at **`100%`** heteroplasmy.
- **`Human_UCP2`:** `E_add >= +0.05` achieved at **`>100% (not reached at 0.05)`** heteroplasmy; `Viability >= 0.25` achieved at **`100%`** heteroplasmy.

## Proof 6: End-to-End SQLite Referential Integrity
- **Foreign-key violations:** `0` | **Orphan payoff rows:** `0` | **Orphan transferability rows:** `0` | **PASS:** `True`.