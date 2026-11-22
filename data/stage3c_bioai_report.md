# Stage 3C Report — AlphaFold DB v6 Structural Audit, Live AlphaGenome AVI Promoter Scans & Co-Folding Manifests

> *Research-Use Framing:* All AlphaGenome AVI scores reflect predicted molecular/regulatory functional impact (never clinical diagnosis or medical advice).

## 1. AlphaFold DB v6 Structural Audit of Stage-4 Compensatory Buffers (`pLDDT` + `PAE`)

| Buffer | UniProt | Length (AA) | Global pLDDT | Very Low (<50) | Very High (>90) | PAE Rigid Domain & Engineering Guardrail |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| `Ciona_AOX` | `H2XKF9` | 276 | **63.88** | 46.4% | 42.0% | 132-276 (145 AA rigid catalytic core; 1-131 flexible MTS/N-tail) |
| `Yeast_NDI1` | `P32340` | 513 | **90.31** | 8.2% | 81.3% | 46-513 (468 AA rigid FAD core; 1-45 MTS) |
| `Human_SOD2` | `P04179` | 222 | **93.19** | 5.0% | 88.3% | 24-222 (199 AA rigid core; 1-23 cleaved MTS) |
| `Human_UCP2` | `P55851` | 309 | **62.53** | 10.0% | 0.0% | 57-304 (248 AA dynamic transmembrane carrier) |
| `Human_PINK1` | `Q9BXM7` | 581 | **77.25** | 19.3% | 52.0% | 268-544 (277 AA rigid kinase domain; 1-267 MTS/TM sensor) |

## 2. Live AlphaGenome Variant Impact (`AVI`) Core Promoter Saturation Scans (`TSS ± 20 bp`, `123` SNVs/gene)

| Nuclear Gene | Role | Scanned Locus (hg38) | SNVs Phred>=10 | SNVs Phred>=15 | Max AVI Phred | Top Variant | Atlas Link |
| :--- | :--- | :--- | :---: | :---: | :---: | :--- | :--- |
| `NDUFA13` | Nuclear_OXPHOS_Partner | `chr19:19516204-19516244` | 123 | **38** | **32.61** | `chr19:19516241:G>T` | [Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=NDUFA13&m=entity) |
| `UQCRFS1` | Nuclear_OXPHOS_Partner | `chr19:29213131-29213171` | 123 | **39** | **22.78** | `chr19:29213160:C>T` | [Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=UQCRFS1&m=entity) |
| `ATP5MJ` | Nuclear_OXPHOS_Partner | `chr14:103921509-103921549` | 123 | **49** | **24.24** | `chr14:103921529:C>A` | [Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=ATP5MJ&m=entity) |
| `SOD2` | Endogenous_Human_Buffer | `chr6:159693221-159693261` | 123 | **0** | **11.54** | `chr6:159693241:T>G` | [Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=SOD2&m=entity) |
| `UCP2` | Endogenous_Human_Buffer | `chr11:73982823-73982863` | 123 | **6** | **19.53** | `chr11:73982830:C>G` | [Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=UCP2&m=entity) |

*(Updated `18` matching rows in SQLite `transferability_assessments.alphagenome_partner_summary`.)*

## 3. Exported AlphaFold 3 JSON & Boltz-2 YAML Co-Folding Manifests (Top 5 Actionable Finalists)

| Finalist | Mutation & rCRS Edit | Base Editor | Cofactors | AF3 JSON | Boltz-2 YAML |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `finalist1_ND4L_A11T` | `ND4L A11T (m.10500G>A)` | `DdCBE_DddA11_HC (0 bystanders)` | `None` | `data/cofolding_manifests/finalist1_ND4L_A11T_af3.json` | `data/cofolding_manifests/finalist1_ND4L_A11T_boltz2.yaml` |
| `finalist2_COX3_W116R` | `COX3 W116R (m.9552T>C)` | `TALED_AT_GC (0 bystanders)` | `HEA` | `data/cofolding_manifests/finalist2_COX3_W116R_af3.json` | `data/cofolding_manifests/finalist2_COX3_W116R_boltz2.yaml` |
| `finalist3_COX3_I165T` | `COX3 I165T (m.9700T>C)` | `TALED_AT_GC (1 bystander)` | `HEA` | `data/cofolding_manifests/finalist3_COX3_I165T_af3.json` | `data/cofolding_manifests/finalist3_COX3_I165T_boltz2.yaml` |
| `finalist4_ATP6_I106V` | `ATP6 I106V (m.8842A>G)` | `TALED_AT_GC (0 bystanders)` | `ATP, MG` | `data/cofolding_manifests/finalist4_ATP6_I106V_af3.json` | `data/cofolding_manifests/finalist4_ATP6_I106V_boltz2.yaml` |
| `finalist5_CYTB_M53T` | `CYTB M53T (m.14904T>C)` | `TALED_AT_GC (0 bystanders)` | `HEM, U10` | `data/cofolding_manifests/finalist5_CYTB_M53T_af3.json` | `data/cofolding_manifests/finalist5_CYTB_M53T_boltz2.yaml` |
