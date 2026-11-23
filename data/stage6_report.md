# Stage 6 Report — Deep Structural & Bio-AI Validation (Live Boltz-2, OpenMM AMBER14 & 8-Gene AlphaGenome AVI)

> *Research Scope & Attribution:* Co-folding executed locally on Apple Silicon Metal (`mps`) using **Boltz-2** (MIT License) with live ColabFold MSA; physics-based atomic energy relaxation executed via **OpenMM + PDBFixer** (`AMBER14 + OBC2` implicit solvent) directly on experimental human cryo-EM holocomplexes (`9I6F`, `9I4I`, `9HZL`, `8H9S`); regulatory promoter impact scored via live **AlphaGenome Variant Impact (`AVI`)** API (research-use only).

## 1. Live Boltz-2 Diffusion Co-Folding: WT vs. Mutant Mitonuclear Dimers (6 Structures)

| Mitonuclear Dimer | PDB | WT `ipTM` | Mutant `ipTM` (`Δ`) | WT / Mut `pLDDT` | Chain A RMSD vs Cryo-EM (WT / Mut) | WT vs Mut Dimer RMSD | 10 Å Pocket RMSD | Verdict |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **`COX3_W116R + COX6A1`** | `9I6F` | `0.9216` | **`0.9172`** (`-0.0043`) | `0.921` / `0.926` | `1.21 Å` / `1.03 Å` | **`0.70 Å`** | `0.31 Å` | `INTERFACE_PRESERVED_HIGH_CONFIDENCE` |
| **`ATP6_I106V + ATP5MJ`** | `8H9S` | `0.6536` | **`0.6770`** (`+0.0234`) | `0.916` / `0.920` | `1.00 Å` / `1.06 Å` | **`1.10 Å`** | `0.33 Å` | `INTERFACE_PRESERVED_HIGH_CONFIDENCE` |
| **`CYTB_M53T + UQCRFS1`** | `9HZL` | `0.7563` | **`0.7693`** (`+0.0129`) | `0.917` / `0.920` | `0.92 Å` / `0.85 Å` | **`1.15 Å`** | `0.08 Å` | `INTERFACE_PRESERVED_HIGH_CONFIDENCE` |

## 2. Physics-Based OpenMM AMBER14 + OBC2 Cryo-EM Energy Relaxation (5 Finalists vs. 2 Pathogenic Controls)

| Variant | Role | Cryo-EM Subassembly | WT 8 Å RMSD | Mutant 8 Å RMSD (`Δ`) | Mut `d_min` | New Stabilizing Polar H-Bonds (`≤3.5 Å`) | Verdict |
| :--- | :--- | :--- | :---: | :---: | :---: | :--- | :--- |
| **`COX3_W116R`** | Stage3_Finalist_1 | `COX1(A), COX3(C), COX6B1(H), HIGD2A(N)` | `0.654 Å` | `0.668 Å` (**`+0.014 Å`**) | `2.93 Å` | `NE--C:SER104.OG(3.38A), NE--C:SER105.N(3.46A), NH1--C:SER105.OG(3.50A) (+3 more)` | **STABLE_LOW_PERTURBATION** |
| **`ND4L_A11T`** | Stage3_Finalist_2 | `ND4L(k), ND6(m)` | `0.391 Å` | `0.368 Å` (**`-0.023 Å`**) | `2.64 Å` | `OG1--k:ASN7.O(2.64A)` | **STABLE_LOW_PERTURBATION** |
| **`COX3_I165T`** | Stage3_Finalist_3 | `COX3(C), COX6A1(G)` | `0.344 Å` | `0.368 Å` (**`+0.024 Å`**) | `3.22 Å` | `OG1--C:GLN161.O(3.22A)` | **STABLE_LOW_PERTURBATION** |
| **`ATP6_I106V`** | Stage3_Finalist_4 | `ATP6(N), ATP5MJ(P)` | `0.277 Å` | `0.298 Å` (**`+0.021 Å`**) | `3.45 Å` | `None` | **STABLE_LOW_PERTURBATION** |
| **`CYTB_M53T`** | Stage3_Finalist_5 | `UQCRFS1(C), CYTB(J), CYTB(V)` | `0.506 Å` | `0.542 Å` (**`+0.036 Å`**) | `3.22 Å` | `OG1--J:LEU49.O(3.22A)` | **STABLE_LOW_PERTURBATION** |
| **`ATP6_L156R`** | Pathogenic_Control_Leigh | `ATP5MC1(2), ATP5MC1(3), ATP6(N)` | `0.406 Å` | `0.566 Å` (**`+0.160 Å`**) | `3.23 Å` | `NE--N:PHE214.O(3.37A), NH2--N:PHE214.O(3.32A)` | **MODERATE_LOCAL_REMODELING** |
| **`ND1_A52T`** | Pathogenic_Control_LHON | `NDUFS7(E), ND1(s)` | `0.405 Å` | `0.444 Å` (**`+0.039 Å`**) | `3.20 Å` | `OG1--s:PRO48.O(3.20A)` | **STABLE_LOW_PERTURBATION** |

## 3. Expanded Live AlphaGenome `AVI` 8-Gene Promoter Constraint Law (`984` SNVs Scored)
- **Nuclear OXPHOS Partners (`n = 6`):** Mean **`41.0`** SNVs with Phred ≥ 15 (mean max Phred = **`27.49`**).
- **Endogenous Human Compensatory Buffers (`n = 2`):** Mean **`3.0`** SNVs with Phred ≥ 15 (mean max Phred = **`15.54`**).
- **One-Sided Mann–Whitney U Test (OXPHOS Partners > Buffers):** `p = 0.0357` (exact minimum attainable p for `n = 6 vs 2` is `1/28 = 0.0357` — complete non-overlapping rank separation!).
- **SQLite Updated:** `28` rows in `transferability_assessments` now carry live AlphaGenome partner regulatory annotations (including Complex IV `COX6A1` & `HIGD2A`).

| Gene | Role | Promoter Window (hg38) | SNVs Phred ≥ 15 (`/123`) | Max AVI Phred | Top Modality |
| :--- | :--- | :--- | :---: | :---: | :--- |
| [`COX6A1`](https://deepmind.google.com/science/alphagenome/atlas?q=COX6A1&m=entity) | Nuclear_OXPHOS_Partner | `chr12:120438092-120438132` | **51** | **34.87** | `Start Lost` |
| [`NDUFA13`](https://deepmind.google.com/science/alphagenome/atlas?q=NDUFA13&m=entity) | Nuclear_OXPHOS_Partner | `chr19:19516204-19516244` | **38** | **32.61** | `Start Lost` |
| [`HIGD2A`](https://deepmind.google.com/science/alphagenome/atlas?q=HIGD2A&m=entity) | Nuclear_OXPHOS_Partner | `chr5:176388730-176388770` | **26** | **26.39** | `Cactus` |
| [`ATP5MJ`](https://deepmind.google.com/science/alphagenome/atlas?q=ATP5MJ&m=entity) | Nuclear_OXPHOS_Partner | `chr14:104679803-104679843` | **49** | **24.24** | `Cactus` |
| [`COX4I1`](https://deepmind.google.com/science/alphagenome/atlas?q=COX4I1&m=entity) | Nuclear_OXPHOS_Partner | `chr16:85799674-85799714` | **43** | **24.04** | `Cactus` |
| [`UQCRFS1`](https://deepmind.google.com/science/alphagenome/atlas?q=UQCRFS1&m=entity) | Nuclear_OXPHOS_Partner | `chr19:29213878-29213918` | **39** | **22.78** | `Cactus` |
| [`UCP2`](https://deepmind.google.com/science/alphagenome/atlas?q=UCP2&m=entity) | Endogenous_Buffer | `chr11:73982823-73982863` | **6** | **19.53** | `PRO-cap` |
| [`SOD2`](https://deepmind.google.com/science/alphagenome/atlas?q=SOD2&m=entity) | Endogenous_Buffer | `chr6:159693221-159693261` | **0** | **11.54** | `PRO-cap` |