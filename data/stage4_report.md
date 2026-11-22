# Stage 4 Report — Driver + Nuclear Buffer Positive Epistasis, Non-Monotonic Cliffs & Prospective Falsification

> *Computational Scope:* Strictly *in silico* biophysical simulation parameterized with `N = 10,000` Latin Hypercube uncertainty priors (`+/-30%`) and `10,240`-eval Saltelli-Sobol GSA to prioritize benchtop *in vitro* human cybrid experiments.

## 1. Numerical Verification & Standalone Trade-Off Proof (No Free Lunch)
- **Solver Convergence (`RK4 dt=0.25 min` vs `DOP853`/`Radau`, `rtol=1e-10`):** `W_end = 0.9290026` vs `0.9290027` / `0.9290027` (`max |Delta| = 1.029e-07`).
- **Human WT (`rCRS`) Ischemia-Reperfusion Baseline:** Viability = `0.0012` (Ischemic ATP nadir = `0.099`, Reperfusion EP peak = `1.155`).

| Standalone Arm | Median I/R Viability | Ischemic ATP Nadir | Reperfusion EP Peak | Biophysical Trade-Off |
| :--- | :---: | :---: | :---: | :--- |
| `Human WT (rCRS)` | `0.0012` | `0.099` | `1.155` | Starves during severe ischemia (`ATP < 0.52`) |
| Driver: `COX3_W116R` | `0.1142` | **`0.802`** | `1.486` | Rescues ischemic ATP nadir, **amplifies reperfusion EP/RET-ROS cliff!** |
| Driver: `ND4L_A11T` | `0.0136` | **`0.327`** | `1.42` | Rescues ischemic ATP nadir, **amplifies reperfusion EP/RET-ROS cliff!** |
| Driver: `COX3_I165T` | `0.0227` | **`0.346`** | `1.322` | Rescues ischemic ATP nadir, **amplifies reperfusion EP/RET-ROS cliff!** |
| Driver: `ATP6_I106V` | `0.0035` | **`0.18`** | `1.235` | Rescues ischemic ATP nadir, **amplifies reperfusion EP/RET-ROS cliff!** |
| Driver: `CYTB_M53T` | `0.0085` | **`0.274`** | `1.371` | Rescues ischemic ATP nadir, **amplifies reperfusion EP/RET-ROS cliff!** |
| Buffer: `Ciona_AOX` | `0.002` | `0.094` | **`0.938`** | Clamps reperfusion EP/ROS, still fails ischemic ATP floor alone |
| Buffer: `Yeast_NDI1` | `0.0021` | `0.092` | **`0.868`** | Clamps reperfusion EP/ROS, still fails ischemic ATP floor alone |
| Buffer: `Human_SOD2` | `0.0023` | `0.098` | **`1.152`** | Clamps reperfusion EP/ROS, still fails ischemic ATP floor alone |
| Buffer: `Human_UCP2` | `0.0018` | `0.074` | **`0.755`** | Clamps reperfusion EP/ROS, still fails ischemic ATP floor alone |
| Buffer: `Human_PINK1` | `0.0015` | `0.098` | **`1.148`** | Clamps reperfusion EP/ROS, still fails ischemic ATP floor alone |

## 2. Top 10 Positive-Epistasis Driver + Nuclear Buffer Pairs (`N = 10,000` LHS)

| Rank | Driver (`rCRS` Edit) | Nuclear Buffer (`AFDB`) | Median I/R Viability `[5th, 95th]` | Synergy Index `S` `[5th, 95th]` | Epistatic Excess `E_add` | `% LHS > 0` | EP Peak (`Drv` -> `Pair`) | Cliff Expansion |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | `COX3_W116R` (`m.9552T>C`) | `Yeast_NDI1` (`P32340`) | **0.4112** `[0.0005, 0.9653]` | **1.322x** `[0.019, 8.122]` | **+0.1187** | `93.5%` | `1.486` -> **`1.044`** | `+3.5 kPa` |
| 2 | `COX3_W116R` (`m.9552T>C`) | `Human_SOD2` (`P04179`) | **0.3703** `[0.0002, 0.9849]` | **1.385x** `[0.008, 4.667]` | **+0.1149** | `95.1%` | `1.486` -> **`1.474`** | `+2.3 kPa` |
| 3 | `COX3_W116R` (`m.9552T>C`) | `Ciona_AOX` (`H2XKF9`) | **0.3095** `[0.0002, 0.9555]` | **1.237x** `[0.008, 4.623]` | **+0.0849** | `92.5%` | `1.486` -> **`1.198`** | `+4.8 kPa` |
| 4 | `COX3_W116R` (`m.9552T>C`) | `Human_UCP2` (`P55851`) | **0.2888** `[0.0008, 0.9743]` | **1.104x** `[0.028, 12.01]` | **+0.0529** | `77.8%` | `1.486` -> **`0.859`** | `+4.6 kPa` |
| 5 | `COX3_W116R` (`m.9552T>C`) | `Human_PINK1` (`Q9BXM7`) | **0.1673** `[0.0, 0.9459]` | **1.105x** `[-0.0, 1.649]` | **+0.0377** | `91.0%` | `1.486` -> **`1.469`** | `-0.8 kPa` |
| 6 | `COX3_I165T` (`m.9700T>C`) | `Human_SOD2` (`P04179`) | **0.0495** `[0.0001, 0.9756]` | **1.048x** `[0.004, 2.723]` | **+0.0147** | `92.8%` | `1.322` -> **`1.318`** | `+3.8 kPa` |
| 7 | `ND4L_A11T` (`m.10500G>A`) | `Human_SOD2` (`P04179`) | **0.0336** `[0.0001, 0.9657]` | **1.025x** `[0.002, 2.809]` | **+0.0110** | `85.3%` | `1.42` -> **`1.411`** | `+2.8 kPa` |
| 8 | `CYTB_M53T` (`m.14904T>C`) | `Human_SOD2` (`P04179`) | **0.0191** `[0.0001, 0.9506]` | **0.686x** `[0.002, 2.289]` | **+0.0054** | `85.0%` | `1.371` -> **`1.364`** | `+3.1 kPa` |
| 9 | `COX3_I165T` (`m.9700T>C`) | `Ciona_AOX` (`H2XKF9`) | **0.0342** `[0.0001, 0.9366]` | **0.792x** `[0.005, 1.957]` | **+0.0021** | `78.1%` | `1.322` -> **`1.072`** | `+6.7 kPa` |
| 10 | `ND4L_A11T` (`m.10500G>A`) | `Yeast_NDI1` (`P32340`) | **0.0257** `[0.0002, 0.9287]` | **0.69x** `[0.006, 2.428]` | **+0.0021** | `76.8%` | `1.42` -> **`1.003`** | `+4.0 kPa` |

## 3. Saltelli-Sobol Global Sensitivity Analysis (`10,240` Evals on `COX3_W116R + Yeast_NDI1`)

| Biophysical Parameter | First-Order $S_1$ | Total-Order $S_T$ | Experimental Measurement Needed |
| :--- | :---: | :---: | :--- |
| `f_load` | `0.2267` | **0.5794** | Baseline cellular respiration flux demand (Seahorse OCR) |
| `drv_scale` | `0.2345` | **0.4599** | Actual heteroplasmy editing efficiency & biochemical effect size of TALED/DdCBE edit |
| `Km_O2_wt` | `0.0678` | **0.2061** | Oroboros O2k high-resolution O2 affinity (p50) in WT vs edited cybrids |
| `theta_GSH` | `0.0189` | **0.1713** | MitoSOX / GSH:GSSG redox buffering capacity under acute oxidative burst |
| `k_RET` | `0.0155` | **0.0916** | Amplex Red H2O2/O2 flux during succinate-driven RET |
| `theta_ATP` | `0.0253` | **0.0676** | CellTiter-Glo / Seahorse ATP-production viability threshold under hypoxia |
| `succ_surge` | `0.0126` | **0.0227** | LC-MS matrix succinate accumulation after 20-min anoxia |
| `k_dmg_ir` | `-0.0001` | **0.0** | TMRM depolarization / mtDNA lesion accumulation post-reperfusion |

## 4. Prospective *In Vitro* Human Cybrid Falsification Protocol (4-Arm Benchtop Test)
To falsify or confirm our top-ranked positive-epistasis prediction (**`COX3_W116R (m.9552T>C)` + `Yeast_NDI1`**) in human `HEK293T` or `iPSC`-cardiomyocyte cybrids:
1. **4-Arm Factorial Design:** (i) Unedited WT `rCRS`, (ii) `TALED`/`DdCBE` single-edit Driver alone (`>=60%` heteroplasmy verified by amplicon NGS with 0 bystander edits), (iii) Lentiviral/AAV Buffer alone (`Ciona AOX` rigid core `132-276` with human Cox8 MTS, or endogenous `SOD2`/`UCP2`), (iv) **Combined Driver + Buffer**.
2. **Assay:** 20-min severe hypoxia (`0.5% O2`) followed by 40-min re-oxygenation (`21% O2`) in an Oroboros Oxygraph-2k / Seahorse XF with simultaneous Amplex Red `H2O2/O2` and `TMRM` membrane-potential readouts.
3. **Quantitative Falsification Criteria:**
   - **Trade-off falsification:** The FBT non-monotonic trade-off hypothesis is falsified if the unbuffered Driver (`COX3 W116R` or `ND4L A11T`) fails to retain `>=15%` higher hypoxic ATP synthesis OR fails to exhibit a `>=20%` higher reperfusion `H2O2/O2` surge than WT.
   - **Positive-epistasis falsification:** The synergy prediction (`S > 1`, `E_add > 0`) is falsified if the combined `Driver + Buffer` arm fails to suppress reperfusion `H2O2/O2` by `>=35%` relative to Driver-alone while maintaining `>=85%` of the Driver's hypoxic ATP gain (i.e., `E_add <= 0` at `alpha = 0.05`).
