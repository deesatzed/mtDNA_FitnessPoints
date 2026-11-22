-- =============================================================================
-- MitoFitness Atlas Relational Schema (SQLite / PostgreSQL-compatible DDL)
-- =============================================================================
-- Extends the 5-table FBT architecture from "Fitness Point mtDNA Adaptations.txt"
-- with paired-clade donor metadata, verified human rCRS reference proteins,
-- ClinVar/gnomAD/MITOMAP safety annotations, AlphaFold/AlphaGenome audit fields,
-- and DdCBE/TALED mitochondrial base-editor feasibility gates.
--
-- Note for PostgreSQL migration:
--   Replace "INTEGER PRIMARY KEY AUTOINCREMENT" with "SERIAL PRIMARY KEY"
--   and JSON-encoded TEXT array columns with TEXT[] / JSONB if desired.
-- =============================================================================

PRAGMA foreign_keys = ON;

-- 0A. HUMAN rCRS REFERENCE GENES (NC_012920.1, NCBI Translation Table 2)
CREATE TABLE IF NOT EXISTS human_rcrs_genes (
    gene_symbol VARCHAR(20) PRIMARY KEY,      -- e.g., 'ND1', 'CYTB', 'COX1', 'ATP6'
    oxphos_complex VARCHAR(20) NOT NULL,      -- 'Complex_I', 'Complex_III', 'Complex_IV', 'Complex_V'
    tier INT NOT NULL CHECK (tier IN (1, 2)), -- 1 = Calibration Core (6 genes), 2 = Frozen Sweep (7 genes)
    rcrs_start_nt INT NOT NULL,               -- 1-based rCRS start coordinate
    rcrs_end_nt INT NOT NULL,                 -- 1-based rCRS end coordinate
    strand CHAR(1) NOT NULL CHECK (strand IN ('+', '-')),
    uniprot_id VARCHAR(20),                   -- Canonical UniProt accession
    expected_aa_len INT NOT NULL,             -- Canonical mature polypeptide length
    cds_nt_seq TEXT NOT NULL,                 -- Coding DNA sequence (5'->3' on coding strand)
    protein_aa_seq TEXT NOT NULL,             -- Translated protein sequence (NCBI Table 2)
    pdb_holocomplex_id VARCHAR(10),           -- Human cryo-EM RCSB PDB ID (populated in Stage 3)
    pdb_chain_id VARCHAR(10)                  -- Matching mtDNA chain ID in PDB holocomplex
);

-- 0B. WILD DONOR & SISTER-CONTROL SPECIES PANEL
CREATE TABLE IF NOT EXISTS donor_species (
    species_id INTEGER PRIMARY KEY AUTOINCREMENT,
    scientific_name VARCHAR(150) NOT NULL UNIQUE,
    common_name VARCHAR(150) NOT NULL,
    taxonomic_clade VARCHAR(100) NOT NULL,    -- e.g., 'Artiodactyla', 'Anseriformes', 'Cetacea', 'Rodentia'
    mitogenome_accession VARCHAR(30) NOT NULL UNIQUE, -- NCBI RefSeq NC_* accession
    stressor_axis VARCHAR(60) NOT NULL,       -- 'High_Altitude_Hypoxia', 'Diving_Hypoxia', 'Subterranean_Hypoxia', 'Arctic_Cold', 'Desert_Heat', 'Control', 'Reference'
    role VARCHAR(20) NOT NULL CHECK (role IN ('Extreme', 'Control', 'Reference')),
    sister_pair_id VARCHAR(50),               -- Links each Extreme species to its closest Control taxon
    habitat_notes TEXT
);

-- 1. WORLD SPACE (W): Objective environmental stressors
CREATE TABLE IF NOT EXISTS environmental_stressors (
    stressor_id INTEGER PRIMARY KEY AUTOINCREMENT,
    stressor_name VARCHAR(100) NOT NULL UNIQUE, -- e.g., 'High_Altitude_Hypoxia', 'Ischemia_Reperfusion', 'Arctic_Cold', 'Desert_Heat'
    physical_variable VARCHAR(50) NOT NULL,     -- e.g., 'pO2_kPa', 'Temp_C', 'pH', 'Flux_f'
    unit VARCHAR(20) NOT NULL,
    notes TEXT
);

-- 2. GENOMIC INTERFACE (X): Convergent mtDNA alleles mapped to human rCRS
CREATE TABLE IF NOT EXISTS mt_alleles (
    allele_id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol VARCHAR(20) NOT NULL REFERENCES human_rcrs_genes(gene_symbol),
    source_species VARCHAR(150) NOT NULL,
    taxonomic_clade VARCHAR(100),
    stressor_axis VARCHAR(60) NOT NULL,
    human_rcrs_aa_pos INT NOT NULL,           -- 1-based orthologous AA position in human rCRS protein
    human_ref_aa CHAR(1) NOT NULL,            -- Human rCRS reference amino acid
    control_aa CHAR(1),                       -- Sister-control taxon amino acid
    alt_aa CHAR(1) NOT NULL,                  -- Extreme-adapted amino acid
    convergent_clade_count INT DEFAULT 1,     -- Independent non-sister extreme clades with same directional shift
    sister_pair_support INT DEFAULT 1,        -- Count of matched Extreme-vs-Control pairs showing shift
    dnds_ratio NUMERIC(6,3),                  -- Optional selection metric (omega)
    grantham_dist NUMERIC(6,2),               -- Physicochemical Grantham distance
    physicochem_shift_score NUMERIC(5,2),     -- Composite radicality score (charge, polarity, volume, hydropathy)
    column_gap_fraction NUMERIC(4,3),         -- Objective alignment QC metric
    column_entropy NUMERIC(5,3),              -- Objective evolutionary entropy
    qc_passed INT DEFAULT 1 CHECK (qc_passed IN (0, 1)),
    UNIQUE (gene_symbol, source_species, human_rcrs_aa_pos, alt_aa)
);

-- 2B. HUMAN CLINICAL & POPULATION SAFETY ANNOTATIONS (ClinVar / gnomAD / MITOMAP)
CREATE TABLE IF NOT EXISTS human_safety_annotations (
    safety_id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol VARCHAR(20) NOT NULL REFERENCES human_rcrs_genes(gene_symbol),
    human_rcrs_aa_pos INT NOT NULL,
    human_ref_aa CHAR(1) NOT NULL,
    alt_aa CHAR(1) NOT NULL,
    rcrs_nt_change VARCHAR(40),               -- e.g., 'm.3394T>C'
    gnomad_af_hom NUMERIC(10,8),              -- gnomAD mtDNA homoplasmic allele frequency
    gnomad_af_het NUMERIC(10,8),              -- gnomAD mtDNA heteroplasmic allele frequency
    clinvar_significance VARCHAR(80),         -- e.g., 'Pathogenic', 'Benign', 'VUS', 'Not_Reported'
    mitomap_status VARCHAR(80),               -- e.g., 'Cfrm_Pathogenic', 'Adaptive_Haplogroup', 'Polymorphism'
    safety_gate_status VARCHAR(40) NOT NULL CHECK (
        safety_gate_status IN ('Safe_Natural_Poly', 'Neutral_Unreported', 'Epistatic_Rescue_Required', 'Contraindicated_Pathogenic')
    ),
    evidence_pmid VARCHAR(40),
    UNIQUE (gene_symbol, human_rcrs_aa_pos, alt_aa)
);

-- 3. FBT PAYOFF ENGINE f(W, O, S, A): Equivalence classes & operating envelopes
CREATE TABLE IF NOT EXISTS fitness_payoff_profiles (
    profile_id INTEGER PRIMARY KEY AUTOINCREMENT,
    allele_id INT NOT NULL REFERENCES mt_alleles(allele_id) ON DELETE CASCADE,
    stressor_id INT REFERENCES environmental_stressors(stressor_id),

    -- FBT Equivalence Class (derived from Stage 2 structural/physicochemical clustering)
    payoff_equivalence_class VARCHAR(100) NOT NULL, -- e.g., 'O2_CHANNEL_HYPOXIA', 'PROTON_LEAK_THERMO', 'Q_SITE_ROS_CLAMP'
    cluster_id INT,
    dist_to_catalytic_center_A NUMERIC(6,2),  -- Distance (Angstroms) to heme, Fe-S, Q-site, or proton channel
    dist_to_nuclear_interface_A NUMERIC(6,2), -- Min distance (Angstroms) to nearest nuclear OXPHOS subunit

    -- Non-monotonic Operating Envelope parameters for f(w) (populated/anchored in Stage 2/4)
    optimal_w_peak NUMERIC(8,3),
    lower_survival_cliff NUMERIC(8,3),
    upper_survival_cliff NUMERIC(8,3),

    -- Thermodynamic Trade-Off Vector (anchored by literature PMIDs or LHS priors in Stage 4)
    atp_o_ratio_delta NUMERIC(5,3),
    ros_leak_delta NUMERIC(5,3),
    heat_dissipation_delta NUMERIC(5,3),

    evidence_tier VARCHAR(20) NOT NULL CHECK (evidence_tier IN ('In_Silico', 'Phylogenetic', 'In_Vitro', 'In_Vivo')),
    anchor_pmid VARCHAR(40)
);

-- 4. STRUCTURAL COMPATIBILITY, BASE-EDITOR FEASIBILITY & ALPHAGENOME PARTNER AUDIT
CREATE TABLE IF NOT EXISTS transferability_assessments (
    assessment_id INTEGER PRIMARY KEY AUTOINCREMENT,
    allele_id INT NOT NULL REFERENCES mt_alleles(allele_id) ON DELETE CASCADE,
    recipient_species VARCHAR(150) NOT NULL DEFAULT 'Homo sapiens',
    pdb_id VARCHAR(10),
    interacting_nuclear_subunits_json TEXT,   -- JSON array of nuclear partner gene symbols within 4.5 A
    min_nuclear_contact_A NUMERIC(6,2),
    cryoem_local_bfactor NUMERIC(6,2),
    afdb_partner_plddt_min NUMERIC(5,2),      -- AlphaFold DB v6 disorder guardrail on partner loop
    steric_clash_score NUMERIC(6,2),          -- Packing/clash perturbation score (kcal/mol proxy or clash index)

    -- Mitochondrial Base-Editor Feasibility (NCBI Table 2 codon check on human rCRS)
    rcrs_ref_codon CHAR(3),
    target_codon CHAR(3),
    required_nt_edits_json TEXT,              -- JSON list of exact rCRS nt substitutions needed
    base_editor_class VARCHAR(50) CHECK (
        base_editor_class IN ('DdCBE_TC', 'DdCBE_DddA11_HC', 'TALED_AT_GC', 'Dual_Transition', 'Not_Base_Editable_Transversion')
    ),
    bystander_nonsyn_count INT,               -- Count of unintended nonsynonymous bystander edits in window

    -- AlphaGenome Nuclear Partner Regulatory Context
    alphagenome_partner_summary TEXT,

    transfer_strategy VARCHAR(50) NOT NULL CHECK (
        transfer_strategy IN ('Motif_Grafting_BaseEdit', 'Allotopic_Nuclear_Express', 'Co_Adapted_Cassette', 'Incompatible')
    ),
    compatibility_score NUMERIC(4,3) CHECK (compatibility_score >= 0.0 AND compatibility_score <= 1.0),
    UNIQUE (allele_id, recipient_species)
);

-- 5. SYNERGY / EPISTASIS: Driver + Nuclear Buffer combinations (Stage 4)
CREATE TABLE IF NOT EXISTS epistatic_synergies (
    synergy_id INTEGER PRIMARY KEY AUTOINCREMENT,
    driver_allele_id INT NOT NULL REFERENCES mt_alleles(allele_id) ON DELETE CASCADE,
    buffer_gene_or_allele VARCHAR(100) NOT NULL, -- e.g., 'Ciona_AOX', 'Yeast_NDI1', 'Human_SOD2', 'Human_UCP2', 'Human_PINK1'
    buffer_afdb_uniprot VARCHAR(20),             -- AlphaFold DB v6 accession audited for buffer
    buffer_alphagenome_locus VARCHAR(80),        -- AlphaGenome regulatory locus for endogenous human buffers
    target_scenario VARCHAR(150) NOT NULL,       -- e.g., 'Human_Acute_Ischemia_Reperfusion', 'High_Flux_Hypoxia'

    -- Net FBT / Energy-Resistance Payoff Shift (median & 5th-95th LHS interval)
    synergy_index_median NUMERIC(6,3),           -- > 1.0 indicates positive epistasis
    synergy_index_p05 NUMERIC(6,3),
    synergy_index_p95 NUMERIC(6,3),
    expanded_cliff_margin NUMERIC(8,3),
    ep_peak_driver_only NUMERIC(6,3),
    ep_peak_driver_plus_buffer NUMERIC(6,3),
    mechanism_summary TEXT
);
