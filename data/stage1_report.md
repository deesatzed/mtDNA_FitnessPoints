# Stage 1 Report (v2 — Rate-Controlled & Directionally Polarized) — Paired Extreme-vs-Control Convergent Mining

- **Wild Vertebrate Panel:** 44 complete NCBI RefSeq mitogenomes (23 matched Extreme-vs-Control pairs across 4 orthogonal stressor axes) + Human `rCRS` (`NC_012920.1`).
- **Objective Alignment QC (`gap_fraction <= 0.15`):** `3774` / `3789` rCRS AA columns passed (`ATP8` auto-trimmed to `83.8%`, core genes `99.4%–100%`).
- **Global Branch-Polarity Sanity Check (Thomas & Hahn Control):**
  - All derived changes from vertebrate consensus: Extreme = `2795` vs Control = `2431` (ratio = **`1.15x`**).
  - Ultra-conserved radical shifts (`entropy <= 0.65`, `consensus >= 85%`, `Grantham >= 75`): Extreme = **`78`** vs Control = `84` (ratio = **`0.929x`**).
- **Clean Candidate Yields (after eliminating symmetric neutral `I<->V` wobble):**
  - **Class I (Conserved-Site Radical Rule-Breakers):** **`65`** sites (`Tier 1 = 41`, `Tier 2 Frozen = 24`).
  - **Class II (Directionally Biased Convergent Dials, `k_ext >= 2`, `net_bias >= +2`):** **`70`** sites (`Tier 1 = 32`, `Tier 2 Frozen = 38`).
  - **SQLite `mt_alleles` rows populated:** `256`.

## 1. Benchmark Sanity Check
- **[RECOVERED — Class I]** Bar-headed goose COX3 W116R (Scott et al. 2011 PNAS): Human `COX3` pos `116` (Consensus=`W` [97.8%], Shift=`W->R`, Grantham=`101.29`, Entropy=`0.154` bits).

## 2. Entropy-Stratified (Rate-Controlled) FBT Cross-Stressor Permutation Tests (10,000 Draws)
*Controls for site-level evolutionary rate across 10 entropy quantile bins so invariant/variable rate heterogeneity cannot create fake enrichment.*

| Comparison | Axis A Sites | Axis B Sites | Observed | Rate-Controlled Expected | Fold Enrichment | Exact p-value |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `Derived_All_Altitude_vs_DivingSubterranean` | 594 | 917 | **282** | 283.36 | **0.995x** | `0.58084` |
| `Derived_Radical50_Altitude_vs_DivingSubterranean` | 321 | 486 | **123** | 92.07 | **1.336x** | `0.0001` |
| `Derived_Radical50_Altitude_vs_ArcticCold` | 321 | 102 | **40** | 23.61 | **1.694x** | `0.0002` |
| `Derived_Radical50_DivingSubterranean_vs_ArcticCold` | 486 | 102 | **51** | 34.46 | **1.48x** | `0.0003` |
| `Derived_Radical50_ArcticCold_vs_DesertHeat` | 102 | 273 | **33** | 20.34 | **1.623x** | `0.0011` |
| `Control_Null_All_AltitudeCtl_vs_DivingCtl` | 564 | 807 | **244** | 236.12 | **1.033x** | `0.21378` |
| `Control_Null_Radical50_AltitudeCtl_vs_DivingCtl` | 282 | 411 | **91** | 65.26 | **1.394x** | `0.0002` |
| `Control_Null_Radical50_AltitudeCtl_vs_ArcticCtl` | 282 | 103 | **36** | 18.11 | **1.988x** | `0.0001` |

## 3. Top 15 Class II Directionally Biased Convergent Dials — Tier 1 (6 Calibration Genes)

| Gene | rCRS Pos | Human -> Alt | `k_ext` vs `k_ctl` | Entropy | Grantham | Cross-Hypoxia? | Supporting Extreme Species |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `COX1` | `57` | `I->V` | **8 vs 1** | 1.21 | 29.61 | YES | Pantholops hodgsonii (High), Ochotona curzoniae (High), Nannospalax galili (Divi), Ursus maritimus (Arct), Ovibos moschatus (Arct), Oryx leucoryx (Dese), Camelus dromedarius (Dese), Jaculus jaculus (Dese) |
| `CYTB` | `115` | `I->V` | **6 vs 0** | 1.126 | 29.61 | YES | Pseudopodoces humilis (High), Physeter macrocephalus (Divi), Mirounga angustirostris (Divi), Vulpes lagopus (Arct), Oryx leucoryx (Dese), Camelus dromedarius (Dese) |
| `ND1` | `313` | `S->G` | **5 vs 1** | 1.323 | 55.27 | YES | Pantholops hodgsonii (High), Ochotona curzoniae (High), Rhinopithecus bieti (High), Physeter macrocephalus (Divi), Camelus dromedarius (Dese) |
| `COX1` | `469` | `I->V` | **5 vs 2** | 1.309 | 29.61 | YES | Ochotona curzoniae (High), Rhinopithecus bieti (High), Leptonychotes weddellii (Divi), Nannospalax galili (Divi), Camelus dromedarius (Dese) |
| `ATP6` | `106` | `I->V` | **3 vs 0** | 0.777 | 29.61 | YES | Ochotona curzoniae (High), Physeter macrocephalus (Divi), Ziphius cavirostris (Divi) |
| `COX3` | `84` | `I->V` | **3 vs 1** | 0.823 | 29.61 | YES | Pseudopodoces humilis (High), Heterocephalus glaber (Divi), Jaculus jaculus (Dese) |
| `COX3` | `165` | `I->T` | **2 vs 0** | 0.689 | 89.28 | YES | Panthera uncia (High), Heterocephalus glaber (Divi) |
| `COX3` | `128` | `E->D` | **3 vs 0** | 0.505 | 44.6 | no | Pseudopodoces humilis (High), Vulpes lagopus (Arct), Jaculus jaculus (Dese) |
| `CYTB` | `181` | `F->L` | **3 vs 0** | 0.875 | 21.82 | no | Physeter macrocephalus (Divi), Mirounga angustirostris (Divi), Mammuthus primigenius (Arct) |
| `ND1` | `7` | `L->I` | **3 vs 0** | 0.875 | 4.86 | no | Ochotona curzoniae (High), Vulpes lagopus (Arct), Vulpes zerda (Dese) |
| `ATP6` | `195` | `I->V` | **3 vs 0** | 1.009 | 29.61 | no | Ursus maritimus (Arct), Oryx leucoryx (Dese), Camelus dromedarius (Dese) |
| `CYTB` | `95` | `F->Y` | **3 vs 1** | 0.991 | 21.61 | no | Physeter macrocephalus (Divi), Leptonychotes weddellii (Divi), Camelus dromedarius (Dese) |
| `ATP6` | `180` | `A->V` | **3 vs 1** | 1.034 | 64.43 | no | Ziphius cavirostris (Divi), Leptonychotes weddellii (Divi), Mirounga angustirostris (Divi) |
| `COX3` | `62` | `T->V` | **3 vs 1** | 1.118 | 29.61 | no | Nannospalax galili (Divi), Camelus dromedarius (Dese), Jaculus jaculus (Dese) |
| `COX1` | `31` | `T->M` | **2 vs 0** | 0.262 | 81.04 | no | Rhinopithecus bieti (High), Camelus dromedarius (Dese) |

## 4. Top 15 Class II Directionally Biased Convergent Dials — Tier 2 (7 Frozen Out-of-Sample Genes)

| Gene | rCRS Pos | Human -> Alt | `k_ext` vs `k_ctl` | Entropy | Grantham | Cross-Hypoxia? | Supporting Extreme Species |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `ND4` | `270` | `I->V` | **5 vs 0** | 1.034 | 29.61 | YES | Ochotona curzoniae (High), Panthera uncia (High), Heterocephalus glaber (Divi), Nannospalax galili (Divi), Jaculus jaculus (Dese) |
| `ND5` | `54` | `F->Y` | **4 vs 0** | 0.772 | 21.61 | YES | Rhinopithecus bieti (High), Physeter macrocephalus (Divi), Ziphius cavirostris (Divi), Jaculus jaculus (Dese) |
| `COX2` | `2` | `A->P` | **3 vs 0** | 0.353 | 26.87 | YES | Rhinopithecus bieti (High), Heterocephalus glaber (Divi), Camelus dromedarius (Dese) |
| `ND5` | `76` | `L->I` | **3 vs 0** | 0.624 | 4.86 | YES | Rhinopithecus bieti (High), Heterocephalus glaber (Divi), Camelus dromedarius (Dese) |
| `ND5` | `366` | `M->L` | **3 vs 0** | 0.971 | 14.3 | YES | Ochotona curzoniae (High), Physeter macrocephalus (Divi), Heterocephalus glaber (Divi) |
| `ND5` | `433` | `G->N` | **3 vs 0** | 1.161 | 79.35 | YES | Ochotona curzoniae (High), Physeter macrocephalus (Divi), Camelus dromedarius (Dese) |
| `ND4` | `455` | `T->M` | **3 vs 0** | 1.187 | 14.3 | YES | Rhinopithecus bieti (High), Heterocephalus glaber (Divi), Camelus dromedarius (Dese) |
| `ND4` | `173` | `S->T` | **4 vs 2** | 0.939 | 57.75 | YES | Ochotona curzoniae (High), Physeter macrocephalus (Divi), Ziphius cavirostris (Divi), Leptonychotes weddellii (Divi) |
| `ND4` | `407` | `S->T` | **3 vs 1** | 0.918 | 57.75 | YES | Ochotona curzoniae (High), Nannospalax galili (Divi), Jaculus jaculus (Dese) |
| `ND5` | `172` | `I->V` | **3 vs 1** | 1.159 | 29.61 | YES | Bos grunniens (High), Ziphius cavirostris (Divi), Jaculus jaculus (Dese) |
| `ND4L` | `5` | `Y->H` | **3 vs 1** | 1.3 | 83.33 | YES | Bos grunniens (High), Physeter macrocephalus (Divi), Oryx leucoryx (Dese) |
| `ND4L` | `11` | `A->T` | **2 vs 0** | 0.262 | 58.02 | YES | Pseudopodoces humilis (High), Nannospalax galili (Divi) |
| `COX2` | `47` | `T->S` | **2 vs 0** | 0.262 | 57.75 | YES | Ochotona curzoniae (High), Heterocephalus glaber (Divi) |
| `ND5` | `63` | `I->T` | **2 vs 0** | 1.126 | 89.28 | YES | Pseudopodoces humilis (High), Nannospalax galili (Divi) |
| `ND4` | `434` | `N->H` | **2 vs 0** | 1.17 | 68.35 | YES | Ochotona curzoniae (High), Physeter macrocephalus (Divi) |

## 5. Top 15 Class I Ultra-Conserved Radical Rule-Breakers (All 13 Genes; Tier 1 + Tier 2)

| Gene (Tier) | rCRS Pos | Cons -> Alt | Cons % | Entropy | Grantham | dCharge | Extreme Species & Stressor |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `COX1` (T1) | `31` | `T->M` | 95.6% | 0.262 | 81.04 | +0.0 | Rhinopithecus bieti (High_Altitude_Hypoxia), Camelus dromedarius (Desert_Heat) |
| `COX1` (T1) | `181` | `T->I` | 91.1% | 0.522 | 89.28 | +0.0 | Leptonychotes weddellii (Diving_Subterranean_Hypoxia), Mirounga angustirostris (Diving_Subterranean_Hypoxia) |
| `ND2` (T1) | `3` | `P->I` | 91.1% | 0.566 | 95.38 | +0.0 | Nannospalax galili (Diving_Subterranean_Hypoxia), Jaculus jaculus (Desert_Heat) |
| `COX3` (T1) | `135` | `S->A` | 86.7% | 0.567 | 99.13 | +0.0 | Heterocephalus glaber (Diving_Subterranean_Hypoxia), Nannospalax galili (Diving_Subterranean_Hypoxia) |
| `ND1` (T1) | `208` | `V->C` | 97.8% | 0.154 | 191.23 | +0.0 | Jaculus jaculus (Desert_Heat) |
| `ND4` (T2) | `334` | `Y->D` | 97.8% | 0.154 | 159.94 | -1.0 | Panthera uncia (High_Altitude_Hypoxia) |
| `CYTB` (T1) | `218` | `I->N` | 97.8% | 0.154 | 148.91 | +0.0 | Jaculus jaculus (Desert_Heat) |
| `COX1` (T1) | `498` | `C->S` | 97.8% | 0.154 | 111.67 | +0.0 | Rhinopithecus bieti (High_Altitude_Hypoxia) |
| `ND5` (T2) | `99` | `S->C` | 97.8% | 0.154 | 111.67 | +0.0 | Rhinopithecus bieti (High_Altitude_Hypoxia) |
| `ND4` (T2) | `347` | `G->V` | 97.8% | 0.154 | 108.79 | +0.0 | Panthera uncia (High_Altitude_Hypoxia) |
| `ND4` (T2) | `206` | `K->A` | 97.8% | 0.154 | 105.57 | -1.0 | Rangifer tarandus (Arctic_Cold) |
| `COX3` (T1) | `116` | `W->R` | 97.8% | 0.154 | 101.29 | +1.0 | Anser indicus (High_Altitude_Hypoxia) |
| `ND1` (T1) | `50` | `A->S` | 97.8% | 0.154 | 99.13 | +0.0 | Jaculus jaculus (Desert_Heat) |
| `COX1` (T1) | `89` | `A->S` | 97.8% | 0.154 | 99.13 | +0.0 | Ziphius cavirostris (Diving_Subterranean_Hypoxia) |
| `CYTB` (T1) | `152` | `A->S` | 97.8% | 0.154 | 99.13 | +0.0 | Panthera uncia (High_Altitude_Hypoxia) |
