# Stage 9 Verification Report — Atomic Physics (OpenMM AMBER14) for Thermal/Baric Leads & Live AlphaGenome 2 Clinical Organ Conditioning

## 1. Stage 9A — Live OpenMM `AMBER14 + OBC2` Cryo-EM Relaxation (`5` Thermal & Baric Leads)
- **Pass Rate:** **`5/5` `STABLE_LOW_PERTURBATION`** vs. Stage-6 Pathogenic Leigh `ATP6 L156R` (`ΔRMSD = +0.160 Å`).
- **`ND4L_Y85H`** (`Desert_Heat_Lead_0Bys`, `9I4I` chains `['NDUFS2(Q)', 'ND2(i)', 'ND4L(k)']`): WT `0.513 Å` → MUT `0.496 Å` (**`ΔRMSD = -0.0170 Å`**), `d_min = 3.61 Å`, new H-bonds: `None` → **`STABLE_LOW_PERTURBATION`**.
- **`ND1_S313G`** (`Desert_Heat_Camel_NDUFA13`, `9I4I` chains `['NDUFA3(U)', 'NDUFA13(W)', 'ND3(j)', 'ND1(s)']`): WT `0.421 Å` → MUT `0.471 Å` (**`ΔRMSD = +0.0500 Å`**), `d_min = 4.00 Å`, new H-bonds: `None` → **`STABLE_LOW_PERTURBATION`**.
- **`COX1_I57V`** (`Arctic_Cold_Lead_0Bys`, `9I6F` chains `['COX1(A)', 'COX7C(L)']`): WT `0.591 Å` → MUT `0.511 Å` (**`ΔRMSD = -0.0800 Å`**), `d_min = 3.45 Å`, new H-bonds: `None` → **`STABLE_LOW_PERTURBATION`**.
- **`CYTB_F181L`** (`Arctic_Cold_Mammoth_UQCRFS1`, `9HZL` chains `['CYTB(J)', 'UQCRFS1(P)', 'CYTB(V)']`): WT `0.380 Å` → MUT `0.384 Å` (**`ΔRMSD = +0.0040 Å`**), `d_min = 3.43 Å`, new H-bonds: `None` → **`STABLE_LOW_PERTURBATION`**.
- **`CYTB_Y224H`** (`Diving_Whale_CYC1`, `9HZL` chains `['UQCRB(F)', 'CYC1(H)', 'CYTB(J)']`): WT `0.483 Å` → MUT `0.514 Å` (**`ΔRMSD = +0.0310 Å`**), `d_min = 2.80 Å`, new H-bonds: `NE2--J:THR23.OG1(2.80A)` → **`STABLE_LOW_PERTURBATION`**.

## 2. Stage 9B — Live `AlphaGenome 2` Across 6 Clinical Organ Compartments (`Heart, Liver, Lung AEC2, Kidney, Pancreatic Beta, CD8+ T-Cell`)
- **Signal-to-Null Specificity:** **`9.62×`** (`7` active promoters mean `0.2245` vs. `2` `BAA-Mito` blank controls mean `0.0233`).
- **`SOD2`** (`Endogenous_Buffer`, `chr6:159693241:T>G`): mean `0.1084`, top organ = **`lung_aec2` (`0.1508`)** (Cardio=`0.048`, Liver=`0.099`, Lung_AEC2=`0.151`, Kidney=`0.113`, Beta=`0.143`, TCell=`0.097`).
- **`UCP2`** (`Endogenous_Buffer`, `chr11:73982830:C>G`): mean `0.3998`, top organ = **`cardiomyocyte` (`0.8562`)** (Cardio=`0.856`, Liver=`0.576`, Lung_AEC2=`0.229`, Kidney=`0.258`, Beta=`0.309`, TCell=`0.170`).
- **`PINK1`** (`Endogenous_Buffer`, `chr1:20633437:C>A`): mean `0.0793`, top organ = **`cd8_tcell` (`0.1167`)** (Cardio=`0.046`, Liver=`0.083`, Lung_AEC2=`0.100`, Kidney=`0.071`, Beta=`0.060`, TCell=`0.117`).
- **`NDUFS2`** (`Thermal_Partner_CI`, `chr1:161197085:G>T`): mean `0.3591`, top organ = **`cd8_tcell` (`1.5594`)** (Cardio=`0.089`, Liver=`0.206`, Lung_AEC2=`0.147`, Kidney=`0.105`, Beta=`0.049`, TCell=`1.559`).
- **`NDUFA13`** (`Thermal_Partner_CI`, `chr19:19516241:G>T`): mean `0.2055`, top organ = **`cardiomyocyte` (`0.3608`)** (Cardio=`0.361`, Liver=`0.200`, Lung_AEC2=`0.201`, Kidney=`0.165`, Beta=`0.154`, TCell=`0.152`).
- **`UQCRFS1`** (`Cold_Partner_CIII`, `chr19:29213160:C>T`): mean `0.1843`, top organ = **`cardiomyocyte` (`0.4992`)** (Cardio=`0.499`, Liver=`0.223`, Lung_AEC2=`0.065`, Kidney=`0.186`, Beta=`0.094`, TCell=`0.039`).
- **`CYC1`** (`Baric_Partner_CIII`, `chr8:144095032:C>A`): mean `0.2354`, top organ = **`lung_aec2` (`0.3277`)** (Cardio=`0.089`, Liver=`0.283`, Lung_AEC2=`0.328`, Kidney=`0.244`, Beta=`0.219`, TCell=`0.251`).
- **`NEG_CTRL_rs9904290`** (`Negative_Control_Blank`, `chr17:45985000:C>T`): mean `0.0252`, top organ = **`pancreatic_beta` (`0.0346`)** (Cardio=`0.026`, Liver=`0.018`, Lung_AEC2=`0.024`, Kidney=`0.021`, Beta=`0.035`, TCell=`0.027`).
- **`NEG_CTRL_rs10799812`** (`Negative_Control_Blank`, `chr1:215190400:G>A`): mean `0.0215`, top organ = **`hepatocyte` (`0.0444`)** (Cardio=`0.018`, Liver=`0.044`, Lung_AEC2=`0.031`, Kidney=`0.013`, Beta=`0.013`, TCell=`0.010`).
