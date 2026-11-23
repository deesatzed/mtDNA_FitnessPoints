# A College Freshman's Guide to the MitoFitness Atlas Explorer

Welcome! This guide walks you through **every section, number, chart, and interactive control** in the **MitoFitness Atlas Explorer** (`app/mitofitness_explorer.html`). You do not need a PhD in biochemistry or statistics to read this—if you have taken Intro Biology and remember what DNA, proteins, mitochondria, and ATP are, you have everything you need.

---

## The Big Picture in 60 Seconds

Your cells are powered by tiny bean-shaped engines called **mitochondria**. Inside each mitochondrion is a mini-chromosome called **mitochondrial DNA (mtDNA)** that encodes **13 core engine parts** (proteins) inside four giant molecular machines called **OXPHOS Complexes I, III, IV, and V**. Together, these four machines burn oxygen to make **ATP** (cellular fuel).

Here is the core question behind the **MitoFitness Atlas**:

> *Wild extremophile animals—like bar-headed geese flying over Mount Everest, deep-diving whales, subterranean naked mole-rats, Arctic foxes, and desert camels—evolved tiny, single-letter changes in their mtDNA so their cellular engines keep running in brutal environments. Which of those wild upgrades could theoretically fit into a **human** mitochondrion, can actually be written with modern DNA base editors, won't cause human mitochondrial disease, and what happens if we turn them on during a heart-attack-like oxygen crash and sudden oxygen rush?*

The Explorer lets you inspect **every piece of real data** our 5-stage computational pipeline gathered to answer that question—plus the **6 formal stress-tests ("proofs")** we ran to try to break our own conclusions.

---

## Top Header Bar & Yellow Scope Banner (Always Visible)

At the very top of the page, you will see:

1. **Title & Subtitle:** Lists the six real public data sources feeding the Explorer:
   - **Comparative Vertebrate mtDNA Convergence** (44 wild animal mitogenomes + Human)
   - **2.5–2.9 Å Human OXPHOS Cryo-EM** (real 3D atomic maps of human mitochondrial machines from the Protein Data Bank)
   - **ClinVar / gnomAD** (real human medical genetics & healthy population variation databases)
   - **AlphaFold DB v6** (3D AI protein structure predictions for helper "buffer" proteins)
   - **AlphaGenome AVI** (AI predictions of how DNA changes affect gene "on/off" switches in the cell nucleus)
   - **Epistasis ODE** (a math simulation of oxygen crash + oxygen rush)
2. **Three Status Pills (Top Right):**
   - `● 6/6 Formal Proofs Verified` — All six automated stress-tests passed.
   - `Human rCRS NC_012920.1` — The official standard reference sequence for human mtDNA (`16,569` DNA letters long).
   - `100% Real Public Data` — Zero made-up sequences or fake structures.
3. **Yellow "Research-Use Scope" Banner:**
   - **Why it is here:** Everything in this tool is *in silico* (done inside a computer). It narrows down thousands of possibilities to the **top 2–5 best lab experiments** to run in human cells in a dish. It is **not** medical advice or a clinical treatment.

---

## Tab 1: `1. Overview & 6-Proof Scorecard`

Think of Tab 1 as the **executive dashboard and lie-detector test** for the whole project.

### Region 1A: The 6 KPI Cards (Top Row)

| Card on Screen | Exact Number | Plain-English Meaning |
| :--- | :---: | :--- |
| **Wild Mitogenomes** | `44 + Human` | We downloaded complete mitochondrial genomes from **44 wild vertebrates** paired into **23 matched "Extreme vs. Normal Cousin" pairs** (for example, *Bar-headed goose* vs. *lowland goose*; *Sperm whale* vs. *cow*) across 4 stress axes (High-Altitude Hypoxia, Diving/Underground Hypoxia, Arctic Cold, Desert Heat). |
| **3D-Mapped Candidates** | `131` | Out of `3,756` total protein positions mapped in 3D human mitochondria, **131 specific amino-acid changes** passed our strict evolutionary filter (showing clear, directional shifts in extreme animals rather than random evolutionary noise). All `131/131` (`100%`) were located on real 3D human atomic structures. |
| **Adj. Interface Enrichment** | `1.519×` (`p = 0.027`) | Extreme-animal adaptations are **~52% more likely** than random background sites to sit right at the **mitonuclear interface** (within `4.5 Å`, or roughly 2 atom-widths, of a partner protein encoded in the cell's nucleus). Even after statistically adjusting for burial depth, evolutionary speed, and gene identity, this result stays statistically significant (`p < 0.05`). |
| **Blind OOS ClinVar AUROC** | `0.888` (`5-CV: 0.890`) | **AUROC** is a `0.5` (coin flip) to `1.0` (perfect) grading score. When we train a simple model on **6 human genes (Tier 1)** to tell real disease mutations apart from harmless human variants—using only 3D structure and animal conservation—and then test it **completely blind** on the **other 7 genes (Tier 2)** it has never seen, it still scores **`0.888`** (an A-/B+ grade out-of-sample!). |
| **Ultra-Clean ±6 nt Edits** | `6` | Mitochondrial base editors act like molecular pencils, but sometimes accidentally smudge neighboring DNA letters within a few bases ("bystander edits"). Only **6 out of 131** wild adaptations stay completely smudge-free (`0` accidental protein changes) even across a wide **`±6`-letter** window—and our top two finalists (`COX3 W116R` and `ND4L A11T`) are among those 6. |
| **Top Epistasis Synergy** | `1.322×` (`93.5% LHS > 0`) | Pairing our #1 wild hypoxia upgrade (`COX3 W116R`) with a nuclear helper enzyme (`Yeast NDI1`) beats the sum of their individual effects by **32.2%** (`1.322×`), and that extra boost holds true in **93.5%** of `10,000` randomized parameter worlds. |

### Region 1B: Stage 5 Scientific Hardening & Formal Proof Scorecard

In science, finding a cool pattern is only half the job—the other half is trying hard to **disprove yourself** the way a tough peer reviewer would. This table shows all **6 formal proofs** we ran:

- **Proof 1 (Mitonuclear Interface Confounders):**
  - *Skeptical question:* "Maybe extreme adaptations just look enriched near nuclear partner proteins because surface spots mutate faster than buried core spots?"
  - *What the table proves:* When we control simultaneously for **3D burial**, **evolutionary speed (Shannon entropy)**, and **gene identity**, the enrichment actually gets *stronger* (`Adjusted OR = 1.519, p = 0.027`). Meanwhile, our **Matched Control-Branch Null** (doing the exact same analysis on the *normal cousin* animals instead of the extreme animals) is **not statistically significant** (`OR = 1.194, p = 0.136`).
- **Proof 2 (Proteome ClinVar Calibration & OOS):**
  - *Skeptical question:* "How do we know your 3D + conservation features actually capture real human biology across the whole mtDNA genome?"
  - *What the table proves:* Across `42` real ClinVar Pathogenic disease mutations vs. `935` harmless/common human variants, our features separate disease from benign with `AUROC = 0.890` (5-fold CV) and `0.888` blindly transferring from Tier 1 → Tier 2 genes (and `0.938` in reverse!).
- **Proof 3 (3D Cryo-EM Rotamer Steric Clashes):**
  - *Skeptical question:* "If you swap an animal amino acid into a human 3D protein, won't the atoms collide like putting a truck engine into a compact car?"
  - *What the table proves:* All **5 of 5** finalists fit into their real human 3D cryo-EM pockets with **zero atomic clashes (`< 2.2 Å`)**. In fact, Bar-headed goose **`COX3 W116R`** replaces a huge Tryptophan (`W`) with a smaller Arginine (`R`), opening `-46 Å³` of extra room (`13/81` clash-free shapes, closest atom `3.21 Å`).
- **Proof 4 (Base-Editor Window Stress Test `±2..±6 nt`):**
  - *Skeptical question:* "Real base editors aren't laser-sharp; what if the editing window spills out to `±5` or `±6` letters?"
  - *What the table proves:* Even at `±6` letters on either side, both **`COX3 W116R`** and **`ND4L A11T`** still produce **`0`** accidental amino-acid smudges.
- **Proof 5 (Supermodularity, Alternative Models & Adversarial Nulls):**
  - *Skeptical question:* "Is your Stage-4 synergy just a quirk of one computer equation?"
  - *What the table proves:* The synergy stays positive across **all 4 mathematical/stress models** (linear ROS, curved Hill ROS, severe 30-min ischemia, mild 10-min ischemia), has **`0`** physics-bound errors, gives **exact `0.0000`** when fed a fake **"Null Placebo"** buffer, and turns **strictly negative** (`-0.051` to `-0.156`) when fed a **"Toxic Uncoupler (FCCP)"**.
- **Proof 6 (Database Referential Integrity):**
  - *Skeptical question:* "Are your SQLite database tables actually lined up cleanly with zero broken IDs?"
  - *What the table proves:* All 8 tables pass foreign-key checks with `0` violations and `0` orphaned records.

### Region 1C: End-to-End Discovery Funnel (Bottom Left)

This horizontal bar chart shows how our strict filters narrow the entire human mitochondrial proteome down to actionable lab experiments:
1. **`3,756`** human OXPHOS protein positions resolved in 3D cryo-EM
2. → **`131`** wild extremophile candidates (directional, polarized evolutionary signal)
3. → **`42`** single-step base-editable with modern `TALED` or `DdCBE` mitochondrial editors (the other 78 require "transversions" that today's editors cannot make!)
4. → **`12`** zero-bystander clean edits inside a standard `±4 nt` editing window
5. → **`6`** ultra-clean zero-bystander edits even out to a harsh `±6 nt` window
6. → **`5`** top Stage-4 finalists tested in full ischemia-reperfusion synergy simulations

### Region 1D: Stage 3 Answer-Key Calibration Controls (Bottom Right)

Before trusting any automated safety filter on *new* wild mutations, you must test it on **5 famous human mutations whose answer we already know** ("the answer key"):
- **`ND1 Y30H` (`m.3394T>C`):** A famous **helpful/harmless** Tibetan high-altitude human variant (present in ~0.9% of humans in gnomAD) → correctly graded **`Safe_Natural_Poly`**.
- **`ND1 A52T`, `ATP6 L156R`, `ND4 R340H`, `ND6 M64V`:** Four classic **severe human mitochondrial disease mutations** (causing LHON blindness or Leigh syndrome, including `ND6 M64V` which sits on the tricky *reverse* DNA strand!) → all 4 correctly blocked as **`Contraindicated_Pathogenic`**.
- **Result:** **`5/5 PASS`**.

---

## Tab 2: `2. 131-Candidate Atlas`

This tab is your **interactive spreadsheet and magnifying glass** for every single one of the 131 wild extremophile adaptations mapped onto human mitochondria.

### Region 2A: Search & Filter Toolbar

You can slice and dice the 131 candidates in real time:
- **Search box:** Type a gene (`COX3`, `ND4L`, `ATP6`), a variant (`W116R`), a nuclear partner (`UQCRFS1`, `ATP5MJ`), or an rCRS coordinate (`m.9552`).
- **All Tiers:**
  - `Tier 1 Calibration Core` = The 6 genes we calibrated our rules on (`ND1, ND2, CYTB, COX1, COX3, ATP6`).
  - `Tier 2 Frozen OOS` = The other 7 genes (`ND3, ND4, ND4L, ND5, ND6, COX2, ATP8`) run under **frozen, untouched rules** so we could prove our pipeline wasn't overfit!
- **All Complexes:** Filter by molecular machine (`Complex I`, `III`, `IV`, or `V`).
- **All Strategies:**
  - `Motif_Grafting_BaseEdit` = Safe + fits in 3D + editable with one base-editor step (the engineering sweet spot!).
  - `Co_Adapted_Cassette` / `Allotopic_Nuclear_Express` = Interesting biology, but either requires a transversion editor we don't have yet or sits in a tight interface that needs partner co-engineering.
  - `Incompatible` = Blocked by human disease safety rules.
- **All Editors:** `TALED (A:T→G:C)` changes A↔G / T↔C; `DdCBE (5'-TC)` and `DdCBE DddA11 (HC)` change C↔T / G↔A; `Transversion` cannot be base-edited today.
- **Checkboxes:**
  - `Mitonuclear Interface (≤4.5 Å)` = Shows only candidates physically touching a human nuclear-encoded partner subunit.
  - `0 Bystanders (±4 nt)` = Shows only the **12 cleanest single-step edits** with zero accidental nearby protein changes.
- **Export CSV:** Downloads whatever filtered list you are currently looking at straight to Excel/Google Sheets.

### Region 2B: Understanding All 13 Columns in the Main Table

| Column | What It Means in Plain English |
| :--- | :--- |
| **Variant** | `Gene_HumanLetterPositionWildLetter`. Example: `COX3_W116R` means at position 116 of gene `COX3`, humans have Tryptophan (`W`), and the extreme animal (Bar-headed goose) swapped it to Arginine (`R`). |
| **Tier** | `T1` (Calibration Core gene) or `T2` (Frozen Out-of-Sample gene). |
| **Complex** | Which of the 4 OXPHOS engine complexes (`I`, `III`, `IV`, or `V`) this protein belongs to. |
| **Class** | `C1_RuleB` = **Class I Rule-Breaker** (an extreme animal breaks an amino acid that ~85–99% of all vertebrates keep frozen!); `C2_Conve` = **Class II Convergent Dial** (multiple independent extreme lineages evolved shifts at the same spot). |
| **Axes** | Which stress environments drove the hit (`High` = High-Altitude Hypoxia, `Diving` = Diving/Underground Hypoxia, `Arctic` = Cold, `Desert` = Heat). |
| **Cons%** | How conserved (unchanged) this spot is across all 45 vertebrate species. `97.8%` means almost every animal on Earth keeps the same amino acid here—so breaking it is a loud evolutionary signal! |
| **Grantham** | A chemistry score (`5` to `215`) measuring how physically different the new amino acid is from the old one. Small (`< 50`) = gentle tweak (like `I → V`); Big (`> 100`) = bold chemical change (like `W → R`, Grantham `101.3`). |
| **Nuc Dist** | Exact 3D distance in **Angstroms (Å)** (`1 Å = 0.1 nm`) to the nearest human nuclear-encoded partner protein in the cryo-EM map. Highlighted in **bold blue** when `≤ 4.5 Å` (direct physical contact!). |
| **Editor & rCRS** | Which mitochondrial base editor can write this change, plus the exact human mtDNA address (e.g., `m.9552T>C` = letter #9,552 on human mtDNA changes from `T` to `C`). |
| **Bys (±4)** | How many *unintended* amino-acid changes happen within `±4` DNA letters of the edit site. **`0` is ideal.** |
| **Safety** | `Safe Natural Poly` (already seen harmlessly in healthy humans in gnomAD/ClinVar!), `Neutral Unreported` (never linked to human disease), or `Epistatic Rescue Required` / `Contraindicated` (flagged by ClinVar!). |
| **Strategy** | Automated engineering verdict (`BaseEdit` is the ready-to-test lane). |
| **Compat** | Overall `0.0` to `1.0` 3D structural compatibility score (`> 0.70` is strong). |

### Region 2C: The Inline Deep-Dive Drawer (Click Any Row!)

Whenever you click a row (try clicking **`COX3_W116R`** or **`CYTB_M53T`**), a blue-bordered inspection card pops open right below the table showing four panels:
1. **Evolutionary Signal:** Exact Extreme vs. Control branch counts, consensus %, Shannon entropy, and stress axes.
2. **3D Cryo-EM Microenvironment:** Exact PDB structure ID & chain, nuclear partner names within `4.5 Å`, distance to the catalytic redox/ATP core, distance to membrane lipids/cardiolipin, local atom crowding (`8 Å burial`), and side-chain volume change (`ΔVol`).
3. **Human Safety & Base Editing:** ClinVar clinical classification, gnomAD homoplasmic allele frequency (`AF_hom`), exact `rCRS` edit, and `±2` to `±6 nt` bystander profile.
4. **AlphaGenome / Compatibility:** 3D compatibility score + live AlphaGenome regulatory constraint notes for any contacting nuclear partner gene.

---

## Tab 3: `3. 3D Cryo-EM, Rotamers & Bio-AI`

Tab 3 zooms down from DNA letters to **actual 3D atoms, AI-folded buffer proteins, and nuclear gene switches**.

### Region 3A: Proof 3 — 3D Experimental Cryo-EM Rotamer Steric-Clash Audit

An amino-acid side chain can rotate around its chemical bonds into several preferred 3D shapes called **rotamers** (like twisting a key in a lock). Here we tested all **5 Stage-3 Finalists** plus **2 real human disease controls** directly inside experimental human cryo-EM crystal/EM structures:
- **Look at `COX3_W116R` (`9I6F`):** People naturally worry that mutating Tryptophan (`W`) to positively charged Arginine (`R`) might jam the protein. The 3D atomic math proves the opposite: Tryptophan has a huge double-ring side chain (`228 Å³`), while Arginine is **`46 Å³` smaller** (`ΔVol = -46.0 Å³`). Out of `81` standard Arginine rotamers, **`13` fit completely clash-free**, with the nearest neighboring heavy atom a roomy **`3.21 Å`** away (`0` severe clashes `< 2.2 Å`)!
- **Contrast with pathogenic `ATP6_L156R` (`8H9S`, Leigh syndrome):** Here Arginine replaces a *smaller* Leucine (`ΔVol = +13.0 Å³`) inside Complex V's tight proton channel—only `4/81` rotamers squeeze in, pressed tightly up against a **`2.26 Å`** steric wall.
- **Exact 5 Å Contacts:** Notice the exact human nuclear residues touching our other finalists: **`COX3 I165T`** touches `COX6A1:PHE43`, **`ATP6 I106V`** touches `ATP5MJ:MET31 & GLY27`, and **`CYTB M53T`** touches `UQCRFS1:MET140` (the Rieske iron-sulfur protein!).

### Region 3B: AlphaFold DB v6 Nuclear Buffer Structural Audits (Bottom Left)

Because Stage 4 pairs wild mtDNA upgrades with **5 protective nuclear "buffer" proteins**, we checked every buffer's 3D confidence using **AlphaFold DB v6**:
- **`Mean pLDDT`** (`0–100`): AlphaFold's per-atom confidence score (`> 90` = rock-solid rigid fold; `< 50` = floppy/disordered noodle).
- **Key Engineering Discovery:** Look at **`Ciona_AOX` (`H2XKF9`)**—its global `pLDDT` is `63.88` because **`46.4%`** of the protein is a floppy, disordered N-terminal tail (residues `1–131`), while its catalytic core (`132–276`) is rock-solid! That tells wet-lab engineers to **trim off the floppy tail** and attach a clean human mitochondrial import signal (Cox8 MTS) directly to the rigid core `132–276`. Meanwhile, **`Yeast NDI1` (`90.31`)** and **`Human SOD2` (`93.19`)** are already ultra-ordered.

### Region 3C: Live AlphaGenome AVI Promoter Saturation Scans (Bottom Right)

When you engineer a mitochondrial protein, what about its **nuclear partner genes** or **endogenous human buffer genes**? We used Google DeepMind's **AlphaGenome** API to test every possible single-letter mutation (`123` SNVs per gene) right across each gene's core promoter ("on/off" switch):
- **`SNVs Phred ≥ 15`** = How many single-letter promoter mutations cause a strong predicted disruption (top ~3% impact across the genome).
- **The Cool Discovery:**
  - The 3 contacting **nuclear OXPHOS partners** (`NDUFA13`: `38/123`, max `32.61`; `UQCRFS1`: `39/123`, max `22.78`; `ATP5MJ`: `49/123`, max `24.24`) have **heavily constrained, fragile promoters**—meaning you should be very careful not to break their native expression.
  - By contrast, the endogenous human buffers behave completely differently: **`UCP2`** has a tight, localized 6-SNV tuning switch right at its promoter (`max 19.53`), while **`SOD2`** has `0` high-impact SNVs in its core TSS (consistent with known biology: `SOD2` is controlled by a distal stress enhancer rather than its core TSS!).
  - *Tip:* Click any gene name with `↗` (e.g., `NDUFA13 ↗`) to open that exact gene live in the public **AlphaGenome Atlas**.

### Region 3D: Proof 4 — Base-Editor Window Stress Test (`±2 nt` to `±6 nt`)

Shows the bystander count for all 5 finalists as we widen the hypothetical DNA base-editor smudge zone from `±2` letters up to `±6` letters. Notice **`ND4L_A11T`** and **`COX3_W116R`** stay straight **`0, 0, 0, 0, 0`** across every single column (`YES (0 bys)`), whereas `COX3_I165T` picks up `1` bystander at `±3 nt` and `CYTB_M53T` picks up `1` at `±5 nt`.

---

## Tab 4: `4. Epistasis & Heteroplasmy Simulator`

This is the **heart of the discovery**—showing *why* a single wild mutation isn't enough on its own, and *how* combining a wild mtDNA driver with a nuclear buffer creates a **"1 + 1 > 2" survival breakthrough** during acute ischemia-reperfusion (a 20-minute severe oxygen cutoff followed by a 40-minute oxygen rush, mimicking heart attack / stroke / transplant stress).

### Region 4A: 5×5 Driver × Nuclear Buffer Positive-Epistasis Matrix (Top Left)

- **What is Positive Epistasis?** In genetics, if Mutation A alone gives `+10%` survival and Helper B alone gives `+0%` survival, you'd expect A + B together to give `+10%`. If A + B together actually gives **`+41%`** survival, that extra `+31%` bonus is **positive epistasis** (`E_add > 0`, synergy `S > 1`).
- **Why does it happen here? (The "Two-Cliff" Trap):**
  1. **Unbuffered Wild Driver (`COX3 W116R` alone):** Awesome during the 20-minute oxygen crash—it keeps cellular ATP high (`0.80` vs. human WT crash to `0.10`)! **The catch:** when oxygen rushes back in at minute 20, its revved-up engine over-pressurizes the electron chain and triggers a massive **RET-ROS (oxidative stress) explosion** (`EP peak` jumps from `1.155` to `1.486`), holding overall survival down to just `11.4%`.
  2. **Nuclear Buffer Alone (`Yeast NDI1`, `Human SOD2`, `Ciona AOX`, `Human UCP2`):** Great at defusing the ROS explosion, **but** they don't help Complex IV grab scarce oxygen during the 20-minute crash—so the cell starves of ATP first (`~0.2%` survival).
  3. **Driver + Buffer Together (Dark Blue Top Row!):** `COX3 W116R` carries the cell safely across the *low-oxygen ATP cliff*, and then `Yeast NDI1` / `Human SOD2` / `Ciona AOX` defuses the *high-oxygen ROS cliff*! Survival jumps to **`41.1%`** (`Yeast NDI1`, `S = 1.32×`, positive in `93.5%` of 10,000 parameter draws) and **`37.0%`** (`Human SOD2`, `S = 1.39×`, positive in `95.1%` of draws).
- **Interactive feature:** Click any square in the 5×5 grid! The gray box below instantly updates with that pair's 90% uncertainty interval, EP peak drop, and how many extra kPa of oxygen range (`Cliff Expansion`) the pair unlocks.

### Region 4B: Interactive Wet-Lab `TALED` Heteroplasmy Slider (Top Right)

- **What is Heteroplasmy?** Each human cell doesn't just have 1 copy of mtDNA—it has **hundreds to thousands of copies**! If a base editor fixes 60% of the mtDNA copies in a cell while 40% stay wild-type human, that cell has **60% heteroplasmy**.
- **Try dragging the slider from `0%` → `100%`!**
  - Watch the gray dashed line (`COX3 W116R` alone) vs. the colored lines (`+ Yeast NDI1`, `+ Human SOD2`, `+ Ciona AOX`, `+ Human UCP2`).
  - **The Huge Practical Discovery:** Between `0%` and `65%` editing, the curves stay almost flat near zero because the unedited copies still let ATP drop below the survival threshold during severe ischemia. Right around **`80%–95%`**, the curves shoot sharply upward!
  - **Takeaway for lab biologists:** Don't test a messy `50%`-edited bulk culture—isolate a **clonal cell line at `≥ 85–95%` heteroplasmy** before running the ischemia-reperfusion assay, or you could get a false negative!

### Region 4C: Proof 5A — Robustness Across Model Forms & Adversarial Controls (Bottom Left)

Checks whether our synergy survives if you change the math rules or stress durations:
- Across **all 4 scenarios** (`Baseline 20m/40m`, `Hill ROS n=4`, `Severe 30m/60m`, `Mild 10m/40m`), `Yeast NDI1`, `Human SOD2`, and `Ciona AOX` stay **green and strongly positive** (`+0.062` to `+0.140`, `88%–100%` positive).
- **Adversarial Sanity Checks:**
  - `Null Placebo` (a dummy buffer that does nothing) scores **exactly `0.0000`** in all 4 scenarios (proving zero mathematical bias).
  - `Toxic FCCP` (a chemical poison that drains mitochondrial energy) scores **red and strongly negative** (`-0.051` to `-0.156`), proving the model penalizes bad buffers instead of rubber-stamping everything!

### Region 4D: Saltelli–Sobol Global Sensitivity (`10,240` Evals, Bottom Right)

Suppose you are heading into a real wet lab and can only measure a few numbers precisely—which biological variables drive the most uncertainty in our prediction?
- **Dark Blue (`S_T`, Total-Order Sensitivity):** Total impact of each parameter including all interactions with other parameters.
- **Light Blue (`S_1`, First-Order Sensitivity):** Its solo effect alone.
- **Top 2 Knobs:** **`f_load`** (how hard the cell is working / baseline oxygen demand, `S_T = 0.579`) and **`drv_scale`** (actual editing heteroplasmy & kinetic strength of the `COX3 W116R` edit, `S_T = 0.460`). That tells experimentalists *exactly* what to measure first (`Seahorse OCR` + deep NGS heteroplasmy).

---

## Tab 5: `5. Wet-Lab Cybrid Protocol`

Science isn't finished until a computer prediction tells a wet-lab scientist **exactly how to prove it wrong**. Tab 5 lays out the exact benchtop blueprint:

1. **4-Arm Factorial Design (Left Card):** Why 4 separate cell dishes?
   - **Arm 1 (WT):** Baseline human cells.
   - **Arm 2 (Driver Alone):** Clonal `COX3 W116R` (`≥ 90%` heteroplasmy, 0 bystanders).
   - **Arm 3 (Buffer Alone):** `Yeast NDI1`, `Human SOD2`, or trimmed `Ciona AOX` (`132–276`) alone.
   - **Arm 4 (Combined):** Driver + Buffer together.
   - *Why all 4?* Without Arms 2 and 3 alone, you can never mathematically prove that Arm 4 is true **positive epistasis (`1 + 1 > 2`)** rather than just one part doing all the work.
2. **Pre-Registered Falsification Cutoffs (Right Card):**
   - Explicitly writes down the exact numbers in an **Oroboros O2k / Seahorse XF** chamber (`20 min 0.5% O2` → `40 min 21% O2`) that would **falsify (disprove)** our model:
     - If Arm 2 does *not* keep `≥ 15%` higher hypoxic ATP or does *not* spike reperfusion `H2O2/O2` by `≥ 20%`, our "two-cliff trade-off" hypothesis is wrong.
     - If Arm 4 does *not* slash Arm 2's reperfusion `H2O2/O2` by `≥ 35%` while preserving `≥ 85%` of its hypoxic ATP gain (`E_add ≤ 0`), our positive-epistasis prediction is wrong.
3. **AlphaFold 3 & Boltz-2 Co-Folding Manifests (Bottom Card):**
   - Five pre-formatted `.json` / `.yaml` files sitting in `data/cofolding_manifests/` so anyone with a GPU or AlphaFold Server account can immediately co-fold our top 5 mutated human subunits together with their contacting nuclear partner proteins and hemes/lipids to double-check 3D interface energetics at atomic resolution.

---

## Quick Freshman Glossary

- **Angstrom (`Å`):** `0.1` nanometers—roughly the diameter of a single hydrogen atom. `≤ 4.5 Å` means two proteins are practically touching.
- **AUROC (`0.5`–`1.0`):** Area Under the ROC Curve—the probability that a classifier ranks a real disease mutation higher than a harmless variant (`0.5` = blind guess, `0.89` = strong discriminator, `1.0` = perfection).
- **Base Editor (`TALED` / `DdCBE`):** Molecular tools that chemically flip a single DNA letter inside mitochondria without cutting the DNA double helix.
- **Bystander Edit:** An accidental extra DNA letter flip inside the editor's active reach (`±2` to `±6` bases away).
- **Cybrid:** A standard lab human cell line engineered to carry a specific mtDNA sequence inside a uniform human nuclear background.
- **Epistasis (`E_add`, `S`):** When two genetic changes interact so their combined effect is *different* from simply adding their individual effects (`E_add > 0` and `S > 1` = **positive epistasis / synergy**).
- **Heteroplasmy (`0%`–`100%`):** The percentage of mtDNA copies inside a cell carrying a specific variant.
- **Ischemia-Reperfusion (I/R):** A temporary cutoff of oxygen/blood flow (*ischemia*), followed by sudden return of oxygen (*reperfusion*), which triggers a burst of reactive oxygen species (**RET-ROS**).
- **Odds Ratio (`OR`):** How many times higher the odds of an event are in one group versus another (`OR = 1.0` = no difference; `OR = 1.52` = 52% higher odds).
- **pLDDT (`0–100`):** AlphaFold's confidence score for how tightly folded a protein segment is (`> 90` very rigid, `< 50` floppy/disordered).
- **Rotamer:** One of the preferred 3D rotational shapes an amino-acid side chain can adopt inside a folded protein.
