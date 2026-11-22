# Stage 2 Report — 3D Cryo-EM Microenvironment Mapping & FBT Payoff Equivalence Classes

- **Experimental Human Cryo-EM Atlas:** Complex I (`9I4I`, `2.63 A`), Complex III (`9HZL`, `2.52 A`), Complex IV (`9I6F`, `2.95 A`), Complex V (`8H9S`, `2.53 A`).
- **3D Coverage:** `3756` / `3789` human mtDNA residues resolved in 3D; **`131` / `131`** Stage-1 candidates mapped directly to atomic coordinates (`0` unresolved flexible N/C-termini).
- **Unsupervised FBT Clustering:** `KMeans (k=4, silhouette=0.179)` + `HDBSCAN (0 dense cores)` over 10 observable physicochemical + 3D structural features (`PC1+PC2 = 38.3%` variance).
- **SQLite `fitness_payoff_profiles` Populated:** `256` rows linked to `mt_alleles`.

## 1. 3D Cryo-EM Micro-Niche Distribution (Candidates vs. Human mtDNA Background)

| 3D OXPHOS Micro-Niche | Candidates (%) | Human Background (%) | Class I / Class II | Odds Ratio | Fisher Exact p |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `Mitonuclear_Interface` | 52 (39.7%) | 1183 (31.5%) | 25 / 29 | **1.45** | `0.04427` |
| `Cardiolipin_Lipid_Interface` | 6 (4.6%) | 339 (9.0%) | 2 / 4 | **0.47** | `0.0856` |
| `Catalytic_Cofactor_Shell` | 14 (10.7%) | 508 (13.5%) | 9 / 6 | **0.76** | `0.43416` |
| `Intra_mtDNA_Interface` | 18 (13.7%) | 487 (13.0%) | 7 / 11 | **1.07** | `0.79092` |
| `Buried_Core` | 0 (0.0%) | 14 (0.4%) | 0 / 0 | **0.00** | `1.0` |
| `Solvent_Lipid_Periphery` | 41 (31.3%) | 1225 (32.6%) | 22 / 20 | **0.94** | `0.77677` |

## 2. 3D Spatial Proximity Test (High-Altitude vs. Diving/Subterranean Hypoxia Within Each Complex)

| Complex (PDB) | Altitude Sites | Diving/Sub Sites | Observed Median 3D NN Dist (A) | Random Expected (A) | p (closer than random) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `Complex_I_9I4I_Altitude_to_DivingSub` | 34 | 43 | **9.71 A** | 12.26 A | `0.007` |
| `Complex_III_9HZL_Altitude_to_DivingSub` | 9 | 8 | **11.04 A** | 13.08 A | `0.2449` |
| `Complex_IV_9I6F_Altitude_to_DivingSub` | 19 | 14 | **11.3 A** | 14.08 A | `0.061` |

## 3. Discovered FBT Payoff Equivalence Classes (Cross-Stressor Breakdown)

| Payoff Equivalence Class | n (T1/T2) | Mean Grantham | Mean Dist Nuclear (A) | Mean Dist Catalytic (A) | Mean Dist Lipid (A) | Stressor Axes Present | Key Exemplars |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- |
| `CONSERVATIVE_TRANSMEMBRANE_DIAL` | 52 (21/31) | 49.7 | 6.48 | 80.75 | 13.45 | Arct:9, Dese:36, Divi:42, High:28 | COX1:I57V (Arct/Dese/Divi/High); CYTB:I115V (Arct/Dese/Divi/High); ND1:S313G (Dese/Divi/High) |
| `POLARITY_GAIN_CHANNEL_REMODELING` | 34 (25/9) | 98.3 | 9.1 | 42.37 | 12.45 | Arct:2, Dese:7, Divi:14, High:14 | ND5:L552S (Dese/Divi); ND4L:A11T (Divi/High); COX1:A39T (Divi) |
| `HYDROPHOBIC_BULK_INTERFACE_SEAL` | 37 (19/18) | 93.7 | 6.67 | 71.36 | 14.02 | Arct:3, Dese:15, Divi:11, High:17 | ND5:H264Y (Arct/Dese/Divi); ATP6:A180V (Divi); ND6:G15C (Dese/Divi) |
| `ELECTROSTATIC_CHARGE_REMODELING` | 8 (4/4) | 118.4 | 6.15 | 59.4 | 13.47 | Arct:1, Divi:3, High:4 | ND5:K581Q (Divi); ND6:W158C (High); ND4:Y334D (High) |
