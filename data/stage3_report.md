# Stage 3B Report — 3D Mitonuclear Compatibility, Live ClinVar/gnomAD Safety & DdCBE/TALED Base-Editor Gate

## 1. Answer-Key Calibration & Out-of-Sample Validation (5/5 Controls)

| Control Role | Gene & AA | Human rCRS Edit | Expected Gate | Observed Gate | gnomAD Hom AF | Base Editor | Verdict |
| :--- | :---: | :---: | :--- | :--- | :---: | :---: | :---: |
| `Tier1_Adaptive_Positive` | `ND1 30->H` | `m.3394T>C` | `Safe_Natural_Poly` | `Safe_Natural_Poly` | `0.00911024` | `TALED_AT_GC` | **PASS** |
| `Tier1_Pathogenic_Negative` | `ND1 52->T` | `m.3460G>A` | `Contraindicated_Pathogenic` | `Contraindicated_Pathogenic` | `0.0` | `DdCBE_DddA11_HC` | **PASS** |
| `Tier1_Pathogenic_Negative` | `ATP6 156->R` | `m.8993T>G` | `Contraindicated_Pathogenic` | `Contraindicated_Pathogenic` | `0.0` | `Not_Base_Editable_Transversion` | **PASS** |
| `Tier2_Frozen_Pathogenic_OOS` | `ND4 340->H` | `m.11778G>A` | `Contraindicated_Pathogenic` | `Contraindicated_Pathogenic` | `0.00019496` | `DdCBE_DddA11_HC` | **PASS** |
| `Tier2_Frozen_Pathogenic_OOS` | `ND6 64->V` | `m.14484T>C` | `Contraindicated_Pathogenic` | `Contraindicated_Pathogenic` | `0.00053166` | `TALED_AT_GC` | **PASS** |

## 2. Automated Funnel Summary Across All 131 3D-Mapped Candidates
- **Transfer Strategy Breakdown:** `{'Motif_Grafting_BaseEdit': 42, 'Allotopic_Nuclear_Express': 70, 'Co_Adapted_Cassette': 14, 'Incompatible': 5}`
- **Human Safety Gate (`ClinVar` + `gnomAD` v3):** `{'Neutral_Unreported': 105, 'Safe_Natural_Poly': 25, 'Epistatic_Rescue_Required': 1}`
- **Base-Editor Class (`NCBI Table 2` on `rCRS`):** `{'DdCBE_DddA11_HC': 14, 'TALED_AT_GC': 32, 'Not_Base_Editable_Transversion': 78, 'Dual_Transition': 7}`
- **Actionable `Motif_Grafting_BaseEdit` Finalists:** **`42`** (`12` with **0 nonsynonymous bystanders** in the ±4 nt editing window!)
- **Top Contacting Human Nuclear Partners (for AlphaGenome Audit):** `{'NDUFA13': 2, 'NDUFB6': 1, 'ATP5MJ': 1, 'COX7C': 1, 'COX8A': 1, 'UQCRFS1': 1, 'NDUFS5': 1}`

## 3. Top 20 Actionable `Motif_Grafting_BaseEdit` Candidates (Safe, Compatible, Base-Editable Today)

| Rank | Gene (T) | Site | rCRS Edit | Base Editor (Bystanders) | Safety (`gnomAD` Hom) | Clash / Compat | Nuclear Partners (<=4.5 A) | FBT Class & Axes |
| :---: | :--- | :---: | :---: | :--- | :--- | :---: | :--- | :--- |
| 1 | `ND4L` (T2) | `A11T` | `m.10500G>A` | `DdCBE_DddA11_HC` (bys=0) | `Neutral_Unreported` (0.0000) | `0.69` / **0.836** | None (core/lipid) | `POLARITY_GAIN_CHAN` (Divi/High) |
| 2 | `ND1` (T1) | `S313G` | `m.4243A>G` | `TALED_AT_GC` (bys=1) | `Safe_Natural_Poly` (0.0001) | `0.0` / **0.99** | NDUFA13 | `CONSERVATIVE_TRANS` (Dese/Divi/High) |
| 3 | `ND4L` (T2) | `Y5H` | `m.10482T>C` | `TALED_AT_GC` (bys=2) | `Neutral_Unreported` (0.0000) | `0.0` / **0.99** | None (core/lipid) | `CONSERVATIVE_TRANS` (Dese/Divi/High) |
| 4 | `COX3` (T1) | `I165T` | `m.9700T>C` | `TALED_AT_GC` (bys=1) | `Safe_Natural_Poly` (0.0001) | `0.0` / **0.99** | None (core/lipid) | `CONSERVATIVE_TRANS` (Divi/High) |
| 5 | `ND5` (T2) | `I63T` | `m.12524T>C` | `TALED_AT_GC` (bys=2) | `Neutral_Unreported` (0.0000) | `0.0` / **0.99** | NDUFB6 | `CONSERVATIVE_TRANS` (Divi/High) |
| 6 | `COX1` (T1) | `I57V` | `m.6072A>G` | `TALED_AT_GC` (bys=0) | `Safe_Natural_Poly` (0.0000) | `0.0` / **0.99** | None (core/lipid) | `CONSERVATIVE_TRANS` (Arct/Dese/Divi/High) |
| 7 | `ATP6` (T1) | `I106V` | `m.8842A>G` | `TALED_AT_GC` (bys=0) | `Safe_Natural_Poly` (0.0023) | `0.0` / **0.99** | ATP5MJ | `CONSERVATIVE_TRANS` (Divi/High) |
| 8 | `ND5` (T2) | `I172V` | `m.12850A>G` | `TALED_AT_GC` (bys=0) | `Safe_Natural_Poly` (0.0002) | `0.0` / **0.99** | None (core/lipid) | `CONSERVATIVE_TRANS` (Dese/Divi/High) |
| 9 | `CYTB` (T1) | `I115V` | `m.15089A>G` | `TALED_AT_GC` (bys=2) | `Neutral_Unreported` (0.0000) | `0.0` / **0.99** | None (core/lipid) | `CONSERVATIVE_TRANS` (Arct/Dese/Divi/High) |
| 10 | `ND4` (T2) | `I270V` | `m.11567A>G` | `TALED_AT_GC` (bys=2) | `Neutral_Unreported` (0.0000) | `0.0` / **0.99** | None (core/lipid) | `CONSERVATIVE_TRANS` (Dese/Divi/High) |
| 11 | `COX1` (T1) | `I469V` | `m.7308A>G` | `TALED_AT_GC` (bys=1) | `Safe_Natural_Poly` (0.0000) | `0.0` / **0.99** | COX7C, COX8A | `CONSERVATIVE_TRANS` (Dese/Divi/High) |
| 12 | `COX3` (T1) | `I84V` | `m.9456A>G` | `TALED_AT_GC` (bys=1) | `Safe_Natural_Poly` (0.0002) | `0.0` / **0.99** | None (core/lipid) | `CONSERVATIVE_TRANS` (Dese/Divi/High) |
| 13 | `ND6` (T2) | `S35G` | `m.14571T>C` | `TALED_AT_GC` (bys=0) | `Neutral_Unreported` (0.0000) | `0.0` / **0.99** | None (core/lipid) | `CONSERVATIVE_TRANS` (Dese/Divi) |
| 14 | `ND5` (T2) | `I443T` | `m.13664T>C` | `TALED_AT_GC` (bys=0) | `Neutral_Unreported` (0.0000) | `0.0` / **0.99** | None (core/lipid) | `CONSERVATIVE_TRANS` (Divi) |
| 15 | `ATP6` (T1) | `A180V` | `m.9065C>T` | `DdCBE_DddA11_HC` (bys=0) | `Safe_Natural_Poly` (0.0000) | `0.0` / **0.99** | None (core/lipid) | `HYDROPHOBIC_BULK_I` (Divi) |
| 16 | `ND1` (T1) | `I314T` | `m.4247T>C` | `TALED_AT_GC` (bys=0) | `Neutral_Unreported` (0.0000) | `0.0` / **0.99** | NDUFA13 | `CONSERVATIVE_TRANS` (Divi) |
| 17 | `COX3` (T1) | `W116R` | `m.9552T>C` | `TALED_AT_GC` (bys=0) | `Neutral_Unreported` (0.0000) | `0.75` / **0.821** | None (core/lipid) | `POLARITY_GAIN_CHAN` (High) |
| 18 | `CYTB` (T1) | `M53T` | `m.14904T>C` | `TALED_AT_GC` (bys=0) | `Neutral_Unreported` (0.0000) | `1.1` / **0.738** | UQCRFS1 | `POLARITY_GAIN_CHAN` (High) |
| 19 | `ND5` (T2) | `H264Y` | `m.13126C>T` | `DdCBE_DddA11_HC` (bys=2) | `Neutral_Unreported` (0.0000) | `0.22` / **0.948** | None (core/lipid) | `HYDROPHOBIC_BULK_I` (Arct/Dese/Divi) |
| 20 | `ND4` (T2) | `I100T` | `m.11058T>C` | `TALED_AT_GC` (bys=2) | `Neutral_Unreported` (0.0000) | `0.0` / **0.99** | None (core/lipid) | `CONSERVATIVE_TRANS` (Arct/Divi) |

## 4. Safety Traps Caught Automatically (Wild Extreme Adaptations That Hit Human Disease Hotspots!)

| Gene (T) | Site | Extreme Species (Axis) | Safety Gate | ClinVar / gnomAD Detail |
| :--- | :---: | :--- | :--- | :--- |
| `ND6` (T2) | `A74S` | Heterocephalus glaber | `Epistatic_Rescue_Required` | Same-site disease hotspot (1 ClinVar Pathogenic hit(s) at pos 74) |
