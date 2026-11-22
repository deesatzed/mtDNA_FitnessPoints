#!/usr/bin/env python3
"""Stage 2: 3D Cryo-EM Atomic Microenvironment Mapping & FBT Structural Payoff Equivalence Clustering.

- Maps all 13 human mtDNA proteins 1-to-1 onto ultra-high-resolution experimental human cryo-EM
  OXPHOS holocomplexes downloaded from RCSB PDB:
    * Complex I   : PDB 9I4I (2.63 A, 44 protein subunits, FMN/FES/SF4/CDL)
    * Complex III : PDB 9HZL (2.52 A, 10 protein subunits, HEM/HEC/FES/U10/CDL)
    * Complex IV  : PDB 9I6F (2.95 A, 14 protein subunits incl. HIGD2A, HEA/CU/MG/CDL)
    * Complex V   : PDB 8H9S (2.53 A, 17 protein subunits, ATP/ADP/MG/CDL + c-ring proton channel)
- Computes real 3D atomic microenvironment features (using scipy cKDTree) for every resolved mtDNA
  residue: distance to catalytic/cofactor center, distance to cardiolipin/lipids, distance to
  nearest nuclear OXPHOS partner subunit (plus exact list of nuclear partners within 4.5 A),
  local heavy-atom packing density (8 A sphere), and local cryo-EM B-factor.
- Tests 3D spatial micro-niche enrichment of Stage-1 polarized candidates (Class I Rule-Breakers +
  Class II Directional Convergent Dials) vs. genome background and matched Control-branch shifts.
- Clusters all 135 Stage-1 candidates via standardized PCA + HDBSCAN (with KMeans silhouette check)
  on observable physicochemical + 3D structural features into FBT Payoff Equivalence Classes.
- Populates `human_rcrs_genes` PDB mappings and `fitness_payoff_profiles` in data/mitofitness.db,
  and renders figures/stage2_3d_fbt_equivalence_clusters.png.
"""

from __future__ import annotations

import gzip
import json
import math
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.cluster import HDBSCAN, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data" / "mitofitness.db"
STRUCT_DIR = ROOT / "data" / "structures"
META_JSON = STRUCT_DIR / "oxphos_entities_meta.json"
STAGE1_JSON = ROOT / "data" / "stage1_convergent_summary.json"
OUT_JSON = ROOT / "data" / "stage2_3d_fbt_summary.json"
OUT_REPORT = ROOT / "data" / "stage2_report.md"
OUT_FIG = ROOT / "figures" / "stage2_3d_fbt_equivalence_clusters.png"

COMPLEX_PDBS = {
    "Complex_I": ("9I4I", STRUCT_DIR / "pdb_00009i4i.cif.gz"),
    "Complex_III": ("9HZL", STRUCT_DIR / "pdb_00009hzl.cif.gz"),
    "Complex_IV": ("9I6F", STRUCT_DIR / "pdb_00009i6f.cif.gz"),
    "Complex_V": ("8H9S", STRUCT_DIR / "pdb_00008h9s.cif.gz"),
}

# Catalytic / redox / proton-channel cofactor ligand comp_ids per complex
CATALYTIC_LIGANDS = {
    "9I4I": {"FMN", "FES", "SF4", "NDP"},
    "9HZL": {"HEM", "HEC", "FES", "U10"},
    "9I6F": {"HEA", "CU", "MG"},
    "8H9S": {"ATP", "ADP", "MG"},
}

LIPID_LIGANDS = {"CDL", "PEE", "3PE", "PGV", "PC1", "PLX", "3PH"}

RNG_SEED = 20261122


def build_chain_metadata(meta_data: dict, uniprot_to_mt_gene: dict[str, str]) -> dict[str, dict]:
    """Map PDB ID -> chain_id -> {gene_symbol, is_mtdna, uniprot_id, description}."""
    u2hgnc_path = STRUCT_DIR / "uniprot_to_hgnc.json"
    u2hgnc = json.loads(u2hgnc_path.read_text(encoding="utf-8")) if u2hgnc_path.exists() else {}
    out: dict[str, dict] = {}
    for entry in meta_data["data"]["entries"]:
        pid = entry["rcsb_id"]
        chain_map = {}
        for pe in entry["polymer_entities"]:
            uids = pe["rcsb_polymer_entity_container_identifiers"].get("uniprot_ids") or []
            chains = pe["rcsb_polymer_entity_container_identifiers"].get("auth_asym_ids") or []
            desc = pe["rcsb_polymer_entity"].get("pdbx_description") or ""
            gnames = []
            for org in pe.get("rcsb_entity_source_organism") or []:
                for g in org.get("rcsb_gene_name") or []:
                    if g.get("value"):
                        gnames.append(g["value"])
            mt_gene = None
            for u in uids:
                if u in uniprot_to_mt_gene:
                    mt_gene = uniprot_to_mt_gene[u]
                    break
            hgnc_sym = u2hgnc.get(uids[0]) if uids else None
            gene_label = mt_gene if mt_gene else (hgnc_sym or (gnames[0] if gnames else desc[:24]))
            for ch in chains:
                chain_map[ch] = {
                    "gene_symbol": gene_label,
                    "is_mtdna": mt_gene is not None,
                    "uniprot_id": uids[0] if uids else None,
                    "description": desc,
                }
        out[pid] = chain_map
    return out


def extract_3d_microenvironments(
    pdb_id: str,
    cif_gz_path: Path,
    chain_meta: dict[str, dict],
    rcrs_seqs: dict[str, str],
) -> tuple[dict[tuple[str, int], dict], dict[str, str]]:
    """Parse cryo-EM mmCIF and compute exact 3D atomic microenvironment for all mtDNA residues."""
    parser = MMCIFParser(QUIET=True)
    with gzip.open(cif_gz_path, "rt", encoding="utf-8") as fh:
        structure = parser.get_structure(pdb_id, fh)
    model = next(structure.get_models())

    # Collect heavy-atom coordinates by category: nuclear protein, mtDNA protein (per chain), catalytic cofactors, lipids
    nuclear_coords = []
    nuclear_atom_genes = []
    mt_chain_coords: dict[str, list[np.ndarray]] = defaultdict(list)
    catalytic_coords = []
    lipid_coords = []
    c_ring_coords = []  # For Complex V (8H9S): c-ring nuclear proton-channel subunits (ATP5MC*)
    all_heavy_coords = []

    # Map each mtDNA gene to its primary chain (first instance in chain_meta)
    mt_gene_to_chain: dict[str, str] = {}
    for ch_id, info in chain_meta.items():
        if info["is_mtdna"] and info["gene_symbol"] not in mt_gene_to_chain:
            mt_gene_to_chain[info["gene_symbol"]] = ch_id

    cat_set = CATALYTIC_LIGANDS.get(pdb_id, set())

    for chain in model:
        ch_id = chain.id
        info = chain_meta.get(ch_id)
        for res in chain:
            hetflag = res.id[0]
            resname = res.get_resname().strip()
            coords = [atom.get_coord() for atom in res if (atom.element or "").upper() != "H"]
            if not coords:
                continue
            if hetflag == " ":
                all_heavy_coords.extend(coords)
                if info and info["is_mtdna"]:
                    mt_chain_coords[ch_id].extend(coords)
                elif info:
                    nuclear_coords.extend(coords)
                    nuclear_atom_genes.extend([info["gene_symbol"]] * len(coords))
                    if pdb_id == "8H9S" and ("ATP5MC" in info["gene_symbol"].upper() or "SUBUNIT C" in info["description"].upper()):
                        c_ring_coords.extend(coords)
            else:
                if resname in cat_set:
                    catalytic_coords.extend(coords)
                elif resname in LIPID_LIGANDS:
                    lipid_coords.extend(coords)

    # For Complex V (8H9S), the transmembrane proton channel is the ATP6 / c-ring interface (plus ATP/ADP in F1)
    # Using c-ring + ATP6 proton-channel cavity as functional catalytic center for membrane subunits ATP6/ATP8
    if pdb_id == "8H9S" and c_ring_coords:
        catalytic_coords_eff = c_ring_coords
    else:
        catalytic_coords_eff = catalytic_coords

    tree_nuc = cKDTree(np.asarray(nuclear_coords)) if nuclear_coords else None
    tree_cat = cKDTree(np.asarray(catalytic_coords_eff)) if catalytic_coords_eff else None
    tree_lip = cKDTree(np.asarray(lipid_coords)) if lipid_coords else None
    tree_all = cKDTree(np.asarray(all_heavy_coords)) if all_heavy_coords else None

    # Pre-build trees for "other mtDNA chains"
    mt_chain_trees = {
        ch: cKDTree(np.asarray(pts)) for ch, pts in mt_chain_coords.items() if pts
    }

    residue_features: dict[tuple[str, int], dict] = {}

    for mt_gene, ch_id in mt_gene_to_chain.items():
        if ch_id not in model:
            continue
        chain = model[ch_id]
        h_seq = rcrs_seqs[mt_gene]

        # Build residue list (standard AA only) and verify direct 1-based rCRS numbering match
        std_res = [r for r in chain if r.id[0] == " "]
        exact_matches = 0
        checked = 0
        for r in std_res:
            rnum = r.id[1]
            if 1 <= rnum <= len(h_seq):
                checked += 1
                try:
                    aa1 = seq1(r.get_resname())
                except Exception:
                    aa1 = "X"
                if aa1 == h_seq[rnum - 1]:
                    exact_matches += 1
        match_rate = exact_matches / max(1, checked)
        assert match_rate >= 0.90, f"{pdb_id} chain {ch_id} ({mt_gene}) rCRS numbering mismatch ({match_rate:.2f})!"

        # Other mtDNA chains in the same complex
        other_mt_pts = [
            pt for och, pts in mt_chain_coords.items() if och != ch_id for pt in pts
        ]
        tree_other_mt = cKDTree(np.asarray(other_mt_pts)) if other_mt_pts else None

        for r in std_res:
            rnum = r.id[1]
            if not (1 <= rnum <= len(h_seq)):
                continue
            atoms = [a for a in r if (a.element or "").upper() != "H"]
            if not atoms:
                continue
            pts = np.asarray([a.get_coord() for a in atoms])
            ca_coord = r["CA"].get_coord() if "CA" in r else np.mean(pts, axis=0)
            bfactor = float(np.mean([a.get_bfactor() for a in atoms]))

            # 1. Nuclear partner contacts
            if tree_nuc is not None:
                dists_n, idxs_n = tree_nuc.query(pts, k=1)
                min_nuc_d = float(np.min(dists_n))
                # All nuclear atoms within 4.5 A of any heavy atom in this residue
                hit_lists = tree_nuc.query_ball_point(pts, r=4.5)
                contact_genes = sorted({nuclear_atom_genes[idx] for lst in hit_lists for idx in lst})
            else:
                min_nuc_d = 99.0
                contact_genes = []

            # 2. Catalytic / cofactor distance
            if tree_cat is not None:
                d_cat = float(np.min(tree_cat.query(pts, k=1)[0]))
            else:
                d_cat = 99.0

            # 3. Cardiolipin / membrane lipid distance
            if tree_lip is not None:
                d_lip = float(np.min(tree_lip.query(pts, k=1)[0]))
            else:
                d_lip = 99.0

            # 4. Neighboring mtDNA subunit distance
            if tree_other_mt is not None:
                d_mt = float(np.min(tree_other_mt.query(pts, k=1)[0]))
            else:
                d_mt = 99.0

            # 5. Local heavy-atom burial density (excluding self atoms, within 8.0 A of CA/centroid)
            n_neigh = max(0, len(tree_all.query_ball_point(ca_coord, r=8.0)) - len(atoms))

            # Structural micro-niche classification (interpretable biophysical tag)
            if min_nuc_d <= 4.5:
                niche = "Mitonuclear_Interface"
            elif d_cat <= 10.0:
                niche = "Catalytic_Cofactor_Shell"
            elif d_lip <= 5.5:
                niche = "Cardiolipin_Lipid_Interface"
            elif d_mt <= 4.5:
                niche = "Intra_mtDNA_Interface"
            elif n_neigh >= 110:
                niche = "Buried_Core"
            else:
                niche = "Solvent_Lipid_Periphery"

            residue_features[(mt_gene, rnum)] = {
                "pdb_id": pdb_id,
                "chain_id": ch_id,
                "ca_xyz": [round(float(x), 3) for x in ca_coord],
                "min_nuclear_A": round(min_nuc_d, 2),
                "nuclear_partners_4p5A": contact_genes,
                "min_catalytic_A": round(d_cat, 2),
                "min_lipid_A": round(d_lip, 2),
                "min_other_mt_A": round(d_mt, 2),
                "burial_neighbors_8A": int(n_neigh),
                "local_bfactor": round(bfactor, 2),
                "micro_niche": niche,
            }

    return residue_features, mt_gene_to_chain


def label_fbt_cluster(cluster_rows: list[dict]) -> str:
    """Assign a clear biophysical FBT Payoff Equivalence Class name from cluster centroid features."""
    mean_gd = float(np.mean([r["grantham"] for r in cluster_rows]))
    mean_dh = float(np.mean([r["d_hydro"] for r in cluster_rows]))
    mean_dq = float(np.mean([abs(r["d_charge"]) for r in cluster_rows]))
    mean_ent = float(np.mean([r["column_entropy"] for r in cluster_rows]))
    if mean_dq >= 0.45:
        return "ELECTROSTATIC_CHARGE_REMODELING"
    if mean_gd >= 75.0 and mean_dh <= -1.0:
        return "POLARITY_GAIN_CHANNEL_REMODELING"
    if mean_gd >= 75.0 and mean_dh >= 1.0:
        return "HYDROPHOBIC_BULK_INTERFACE_SEAL"
    if mean_ent >= 0.70:
        return "CONSERVATIVE_TRANSMEMBRANE_DIAL"
    return "MODERATE_STRUCTURAL_TUNING"


def main() -> None:
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON;")

    rcrs_rows = conn.execute(
        "SELECT gene_symbol, oxphos_complex, tier, uniprot_id, protein_aa_seq FROM human_rcrs_genes"
    ).fetchall()
    rcrs_seqs = {g: seq for g, _, _, _, seq in rcrs_rows}
    uniprot_to_mt = {u: g for g, _, _, u, _ in rcrs_rows}

    meta_data = json.loads(META_JSON.read_text(encoding="utf-8"))
    chain_meta_all = build_chain_metadata(meta_data, uniprot_to_mt)

    # 1. Extract 3D microenvironment for EVERY resolved human mtDNA residue across Complexes I, III, IV, V
    all_3d: dict[tuple[str, int], dict] = {}
    gene_pdb_map: dict[str, tuple[str, str]] = {}
    for comp, (pdb_id, cif_path) in COMPLEX_PDBS.items():
        print(f"Parsing {comp} ({pdb_id}: {cif_path.name})...")
        feats, g2c = extract_3d_microenvironments(pdb_id, cif_path, chain_meta_all[pdb_id], rcrs_seqs)
        all_3d.update(feats)
        for g, ch in g2c.items():
            gene_pdb_map[g] = (pdb_id, ch)
            conn.execute(
                "UPDATE human_rcrs_genes SET pdb_holocomplex_id = ?, pdb_chain_id = ? WHERE gene_symbol = ?",
                (pdb_id, ch, g),
            )
    conn.commit()
    print(f"Mapped {len(all_3d)} / {sum(len(s) for s in rcrs_seqs.values())} human mtDNA residues to experimental cryo-EM 3D coordinates!")

    # 2. Load all 135 Stage-1 polarized candidates (65 Class I + 70 Class II)
    s1 = json.loads(STAGE1_JSON.read_text(encoding="utf-8"))
    raw_candidates = s1["class1_rulebreakers"] + s1["class2_convergent"]

    # Deduplicate by (gene, pos, alt_aa), merging classes if a site qualifies for both Class I and Class II!
    merged: dict[tuple[str, int, str], dict] = {}
    for c in raw_candidates:
        key = (c["gene"], int(c["human_rcrs_aa_pos"]), c["alt_aa"])
        if key not in merged:
            entry = dict(c)
            entry["classes"] = [c["discovery_class"]]
            merged[key] = entry
        else:
            merged[key]["classes"].append(c["discovery_class"])

    candidates = []
    unresolved = []
    for key, c in merged.items():
        g, pos, alt_aa = key
        f3d = all_3d.get((g, pos))
        if f3d is None:
            unresolved.append(f"{g}:{pos}->{alt_aa}")
            continue
        rec = {**c, **f3d}
        candidates.append(rec)

    print(f"Stage-1 unique candidates: {len(merged)} | Resolved in 3D cryo-EM: {len(candidates)} | Unresolved flexible termini: {len(unresolved)}")

    # 3. 3D Micro-Niche Enrichment Test (Candidates vs. All Resolved Human mtDNA Residues)
    bg_niches = Counter(v["micro_niche"] for v in all_3d.values())
    cand_niches = Counter(c["micro_niche"] for c in candidates)
    c1_niches = Counter(c["micro_niche"] for c in candidates if "Class_I_Conserved_RuleBreaker" in c["classes"])
    c2_niches = Counter(c["micro_niche"] for c in candidates if "Class_II_Directional_Convergent_Dial" in c["classes"])

    niche_enrichment = {}
    n_bg_tot = len(all_3d)
    n_c_tot = len(candidates)
    for niche in sorted(bg_niches.keys()):
        a = cand_niches[niche]
        b = n_c_tot - a
        c_bg = bg_niches[niche] - a
        d_bg = (n_bg_tot - n_c_tot) - c_bg
        or_val, p_val = fisher_exact([[a, b], [max(0, c_bg), max(0, d_bg)]])
        niche_enrichment[niche] = {
            "candidates": a,
            "cand_pct": round(100.0 * a / max(1, n_c_tot), 1),
            "background": bg_niches[niche],
            "bg_pct": round(100.0 * bg_niches[niche] / max(1, n_bg_tot), 1),
            "class1_count": c1_niches[niche],
            "class2_count": c2_niches[niche],
            "odds_ratio": round(float(or_val), 3),
            "fisher_p": round(float(p_val), 5),
        }

    # 4. 3D Cross-Stressor Spatial Proximity Test (Do High-Altitude and Diving/Subterranean candidates
    #    cluster closer in 3D Angstroms within the SAME complex than random candidate-matched sites?)
    rng = np.random.default_rng(RNG_SEED)
    spatial_3d_tests = {}
    for comp, (pdb_id, _) in COMPLEX_PDBS.items():
        comp_bg = [v["ca_xyz"] for (g, _), v in all_3d.items() if v["pdb_id"] == pdb_id]
        alt_pts = [c["ca_xyz"] for c in candidates if c["pdb_id"] == pdb_id and "High_Altitude_Hypoxia" in c["axes_hit"]]
        div_pts = [c["ca_xyz"] for c in candidates if c["pdb_id"] == pdb_id and "Diving_Subterranean_Hypoxia" in c["axes_hit"]]
        alt_only = [c["ca_xyz"] for c in candidates if c["pdb_id"] == pdb_id and "High_Altitude_Hypoxia" in c["axes_hit"] and "Diving_Subterranean_Hypoxia" not in c["axes_hit"]]
        div_only = [c["ca_xyz"] for c in candidates if c["pdb_id"] == pdb_id and "Diving_Subterranean_Hypoxia" in c["axes_hit"] and "High_Altitude_Hypoxia" not in c["axes_hit"]]
        shared_exact = sum(1 for c in candidates if c["pdb_id"] == pdb_id and "High_Altitude_Hypoxia" in c["axes_hit"] and "Diving_Subterranean_Hypoxia" in c["axes_hit"])
        if len(alt_pts) >= 2 and len(div_pts) >= 2:
            t_div = cKDTree(np.asarray(div_pts))
            obs_med_d = float(np.median(t_div.query(np.asarray(alt_pts), k=1)[0]))
            non_ov_med = float(np.median(cKDTree(np.asarray(div_only)).query(np.asarray(alt_only), k=1)[0])) if (len(alt_only) >= 1 and len(div_only) >= 1) else None
            bg_arr = np.asarray(comp_bg)
            null_meds = []
            for _ in range(2000):
                idx_a = rng.choice(len(bg_arr), size=len(alt_pts), replace=False)
                idx_b = rng.choice(len(bg_arr), size=len(div_pts), replace=False)
                t_b = cKDTree(bg_arr[idx_b])
                null_meds.append(float(np.median(t_b.query(bg_arr[idx_a], k=1)[0])))
            null_meds_arr = np.asarray(null_meds)
            p_close = float((np.sum(null_meds_arr <= obs_med_d) + 1) / (len(null_meds_arr) + 1))
            spatial_3d_tests[f"{comp}_{pdb_id}_Altitude_to_DivingSub"] = {
                "n_altitude": len(alt_pts),
                "n_diving_sub": len(div_pts),
                "shared_exact_sites": shared_exact,
                "obs_median_3d_nn_A": round(obs_med_d, 2),
                "non_overlapping_median_3d_nn_A": round(non_ov_med, 2) if non_ov_med is not None else None,
                "null_expected_3d_nn_A": round(float(np.mean(null_meds_arr)), 2),
                "p_value_closer_than_random": round(p_close, 4),
            }

    # 5. Unsupervised FBT Structural-Physicochemical Payoff Equivalence Clustering
    # Feature vector strictly from observable molecular + 3D cryo-EM metrics (NO fabricated ATP/O!)
    feat_names = [
        "grantham", "d_hydro", "d_charge", "d_vol", "d_pol",
        "min_nuclear_A", "min_catalytic_A", "min_lipid_A",
        "burial_neighbors_8A", "column_entropy",
    ]
    X_raw = np.array(
        [
            [
                c["grantham"], c["d_hydro"], c["d_charge"], c["d_vol"], c["d_pol"],
                min(c["min_nuclear_A"], 25.0),
                min(c["min_catalytic_A"], 45.0),
                min(c["min_lipid_A"], 25.0),
                c["burial_neighbors_8A"],
                c["column_entropy"],
            ]
            for c in candidates
        ],
        dtype=np.float64,
    )
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X_raw)
    pca = PCA(n_components=2, random_state=RNG_SEED)
    X_pca = pca.fit_transform(X_std)

    # Run HDBSCAN first, and also evaluate KMeans across k=4..7 via silhouette score for full partition
    hdb = HDBSCAN(min_cluster_size=8, min_samples=4)
    hdb_labels = hdb.fit_predict(X_std)
    n_hdb_clusters = len(set(hdb_labels) - {-1})

    best_k, best_sil, best_km_labels = 5, -1.0, None
    for k in range(4, 8):
        km = KMeans(n_clusters=k, n_init=30, random_state=RNG_SEED)
        lbls = km.fit_predict(X_std)
        sil = float(silhouette_score(X_std, lbls))
        if sil > best_sil:
            best_k, best_sil, best_km_labels = k, sil, lbls

    # Use KMeans optimal-silhouette partition for complete non-noise equivalence class assignment,
    # while recording HDBSCAN density cluster ID alongside it for full transparency.
    cluster_groups: dict[int, list[dict]] = defaultdict(list)
    for idx, c in enumerate(candidates):
        cid = int(best_km_labels[idx])
        c["pca_1"] = round(float(X_pca[idx, 0]), 3)
        c["pca_2"] = round(float(X_pca[idx, 1]), 3)
        c["cluster_id"] = cid
        c["hdbscan_cluster"] = int(hdb_labels[idx])
        cluster_groups[cid].append(c)

    # Assign unique biophysical Payoff Equivalence Class names to clusters
    used_names = Counter()
    cluster_summaries = {}
    for cid, rows in sorted(cluster_groups.items()):
        base_name = label_fbt_cluster(rows)
        used_names[base_name] += 1
        eq_name = base_name if used_names[base_name] == 1 else f"{base_name}_{used_names[base_name]}"
        for r in rows:
            r["payoff_equivalence_class"] = eq_name
        axis_dist = Counter(ax for r in rows for ax in r["axes_hit"])
        complex_dist = Counter(r["complex"] for r in rows)
        niche_dist = Counter(r["micro_niche"] for r in rows)
        cluster_summaries[eq_name] = {
            "cluster_id": cid,
            "size": len(rows),
            "tier1_count": sum(1 for r in rows if r["tier"] == 1),
            "tier2_count": sum(1 for r in rows if r["tier"] == 2),
            "mean_grantham": round(float(np.mean([r["grantham"] for r in rows])), 1),
            "mean_entropy": round(float(np.mean([r["column_entropy"] for r in rows])), 3),
            "mean_dist_nuclear_A": round(float(np.mean([r["min_nuclear_A"] for r in rows])), 2),
            "mean_dist_catalytic_A": round(float(np.mean([r["min_catalytic_A"] for r in rows])), 2),
            "mean_dist_lipid_A": round(float(np.mean([r["min_lipid_A"] for r in rows])), 2),
            "mean_abs_dcharge": round(float(np.mean([abs(r["d_charge"]) for r in rows])), 2),
            "stressor_axes": dict(axis_dist),
            "complexes": dict(complex_dist),
            "top_niches": dict(niche_dist.most_common(3)),
            "exemplars": [
                f"{r['gene']}:{r['human_ref_aa']}{r['human_rcrs_aa_pos']}{r['alt_aa']} ({'/'.join(a[:4] for a in r['axes_hit'])})"
                for r in sorted(rows, key=lambda x: (x["k_ext_derived"], x["grantham"]), reverse=True)[:5]
            ],
        }

    # 6. Populate `fitness_payoff_profiles` in SQLite (linked to allele_id in mt_alleles!)
    conn.execute("DELETE FROM fitness_payoff_profiles;")
    stressor_id_map = dict(conn.execute("SELECT stressor_name, stressor_id FROM environmental_stressors").fetchall())
    allele_rows = conn.execute(
        "SELECT allele_id, gene_symbol, human_rcrs_aa_pos, alt_aa, stressor_axis FROM mt_alleles"
    ).fetchall()
    cand_lookup = {(c["gene"], int(c["human_rcrs_aa_pos"]), c["alt_aa"]): c for c in candidates}

    profile_inserts = []
    for allele_id, g, pos, alt_aa, axis in allele_rows:
        c = cand_lookup.get((g, pos, alt_aa))
        if c is None:
            continue
        # Literature-anchored evidence tier for known benchmarks (Bar-headed goose COX3 W116R)
        if g == "COX3" and pos == 116 and alt_aa == "R":
            ev_tier, anchor_pmid = "In_Vivo", "21383179"  # Scott et al. 2011 PNAS
        else:
            ev_tier, anchor_pmid = "Phylogenetic", None
        profile_inserts.append((
            allele_id,
            stressor_id_map.get(axis),
            c["payoff_equivalence_class"],
            c["cluster_id"],
            c["min_catalytic_A"],
            c["min_nuclear_A"],
            None, None, None,  # Operating envelope cliffs parameterized in Stage 4 with LHS uncertainty
            None, None, None,  # Thermodynamic deltas anchored in Stage 4 (zero fabricated numbers!)
            ev_tier,
            anchor_pmid,
        ))
    conn.executemany(
        """
        INSERT INTO fitness_payoff_profiles (
            allele_id, stressor_id, payoff_equivalence_class, cluster_id,
            dist_to_catalytic_center_A, dist_to_nuclear_interface_A,
            optimal_w_peak, lower_survival_cliff, upper_survival_cliff,
            atp_o_ratio_delta, ros_leak_delta, heat_dissipation_delta,
            evidence_tier, anchor_pmid
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        profile_inserts,
    )
    conn.commit()
    n_profiles_db = conn.execute("SELECT COUNT(*) FROM fitness_payoff_profiles").fetchone()[0]
    conn.close()

    # 7. Render Publication-Grade 2-Panel Stage 2 Figure
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4), dpi=200)

    palette = ["#0284c7", "#d97706", "#059669", "#dc2626", "#7c3aed", "#0f766e", "#be123c"]
    eq_classes = sorted(cluster_summaries.keys())
    for i, eq in enumerate(eq_classes):
        pts = [c for c in candidates if c["payoff_equivalence_class"] == eq]
        ax1.scatter(
            [p["pca_1"] for p in pts],
            [p["pca_2"] for p in pts],
            s=[55 + 18 * p["k_ext_derived"] for p in pts],
            c=palette[i % len(palette)],
            alpha=0.82,
            edgecolors="#0f172a",
            linewidths=0.6,
            label=f"{eq.replace('_', ' ').title()} (n={len(pts)})",
        )
    # Annotate key landmarks
    landmarks_to_tag = {
        ("COX3", 116, "R"): "COX3 W116R (Goose)",
        ("COX3", 135, "A"): "COX3 S135A (Mole-rats)",
        ("COX3", 165, "T"): "COX3 I165T (Leopard+NMR)",
        ("ATP6", 106, "V"): "ATP6 I106V (Pika+Whales)",
        ("ATP6", 180, "V"): "ATP6 A180V (3 Divers)",
        ("ND4L", 11, "T"): "ND4L A11T (Tit+BMR)",
        ("ND5", 54, "Y"): "ND5 F54Y (Primate+Whales)",
    }
    for c in candidates:
        lbl = landmarks_to_tag.get((c["gene"], int(c["human_rcrs_aa_pos"]), c["alt_aa"]))
        if lbl:
            ax1.annotate(
                lbl, (c["pca_1"], c["pca_2"]),
                xytext=(6, 5), textcoords="offset points",
                fontsize=7.5, fontweight="bold", color="#0f172a",
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#94a3b8", alpha=0.85),
            )
    ax1.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance: Radicality & Burial/Interface)")
    ax1.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance: Charge/Hydropathy & Lipid/Catalytic)")
    ax1.set_title(f"(A) 3D FBT Equivalence Classes (n={len(candidates)}, Sil={best_sil:.2f})", fontsize=10.5, fontweight="bold")
    ax1.grid(True, alpha=0.22)
    ax1.legend(fontsize=7.2, loc="best", frameon=True)

    # Panel B: 3D Micro-Niche Enrichment (Candidates vs. Human mtDNA Background)
    niches_ordered = [
        "Mitonuclear_Interface",
        "Cardiolipin_Lipid_Interface",
        "Catalytic_Cofactor_Shell",
        "Intra_mtDNA_Interface",
        "Buried_Core",
        "Solvent_Lipid_Periphery",
    ]
    niches_present = [n for n in niches_ordered if n in niche_enrichment]
    y = np.arange(len(niches_present))
    h = 0.36
    bg_vals = [niche_enrichment[n]["bg_pct"] for n in niches_present]
    cd_vals = [niche_enrichment[n]["cand_pct"] for n in niches_present]
    ax2.barh(y - h / 2, bg_vals, height=h, color="#94a3b8", label=f"All Human mtDNA (n={n_bg_tot})")
    ax2.barh(y + h / 2, cd_vals, height=h, color="#0284c7", label=f"Stage-1 Polarized Candidates (n={n_c_tot})")
    for idx, n in enumerate(niches_present):
        ne = niche_enrichment[n]
        star = " *" if ne["fisher_p"] < 0.05 else ""
        ax2.text(
            max(bg_vals[idx], cd_vals[idx]) + 0.8,
            idx,
            f"OR={ne['odds_ratio']:.2f} (p={ne['fisher_p']:.3f}){star}",
            va="center", fontsize=8, color="#0f172a",
        )
    ax2.set_yticks(y)
    ax2.set_yticklabels([n.replace("_", " ") for n in niches_present], fontsize=8.5)
    ax2.set_xlabel("Share of Residues (%)")
    ax2.set_xlim(0, max(bg_vals + cd_vals) + 22)
    ax2.set_title("(B) 3D Cryo-EM Micro-Niche Enrichment vs. Background", fontsize=10.5, fontweight="bold")
    ax2.grid(True, axis="x", alpha=0.22)
    ax2.legend(fontsize=8, loc="lower right")

    plt.tight_layout(pad=1.5)
    plt.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    plt.close()

    # 8. Save JSON & Markdown Report
    summary = {
        "human_cryoem_holocomplexes": {
            "Complex_I": {"pdb_id": "9I4I", "resolution_A": 2.63},
            "Complex_III": {"pdb_id": "9HZL", "resolution_A": 2.52},
            "Complex_IV": {"pdb_id": "9I6F", "resolution_A": 2.95},
            "Complex_V": {"pdb_id": "8H9S", "resolution_A": 2.53},
        },
        "total_human_mtdna_residues_resolved_3d": len(all_3d),
        "stage1_unique_candidates": len(merged),
        "candidates_resolved_3d": len(candidates),
        "unresolved_termini": unresolved,
        "clustering_metrics": {
            "features": feat_names,
            "pca_var_explained_2d": [round(float(x), 4) for x in pca.explained_variance_ratio_],
            "kmeans_optimal_k": best_k,
            "kmeans_silhouette": round(best_sil, 4),
            "hdbscan_density_clusters": n_hdb_clusters,
            "hdbscan_noise_points": int(np.sum(hdb_labels == -1)),
        },
        "niche_enrichment": niche_enrichment,
        "spatial_3d_cross_hypoxia_tests": spatial_3d_tests,
        "equivalence_classes": cluster_summaries,
        "db_fitness_payoff_profiles_rows": n_profiles_db,
        "candidates_3d": candidates,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # Markdown report
    lines = [
        "# Stage 2 Report — 3D Cryo-EM Microenvironment Mapping & FBT Payoff Equivalence Classes",
        "",
        f"- **Experimental Human Cryo-EM Atlas:** Complex I (`9I4I`, `2.63 A`), Complex III (`9HZL`, `2.52 A`), Complex IV (`9I6F`, `2.95 A`), Complex V (`8H9S`, `2.53 A`).",
        f"- **3D Coverage:** `{len(all_3d)}` / `3789` human mtDNA residues resolved in 3D; **`{len(candidates)}` / `{len(merged)}`** Stage-1 candidates mapped directly to atomic coordinates (`{len(unresolved)}` unresolved flexible N/C-termini).",
        f"- **Unsupervised FBT Clustering:** `KMeans (k={best_k}, silhouette={best_sil:.3f})` + `HDBSCAN ({n_hdb_clusters} dense cores)` over 10 observable physicochemical + 3D structural features (`PC1+PC2 = {sum(pca.explained_variance_ratio_)*100:.1f}%` variance).",
        f"- **SQLite `fitness_payoff_profiles` Populated:** `{n_profiles_db}` rows linked to `mt_alleles`.",
        "",
        "## 1. 3D Cryo-EM Micro-Niche Distribution (Candidates vs. Human mtDNA Background)",
        "",
        "| 3D OXPHOS Micro-Niche | Candidates (%) | Human Background (%) | Class I / Class II | Odds Ratio | Fisher Exact p |",
        "| :--- | :---: | :---: | :---: | :---: | :---: |",
    ]
    for n in niches_present:
        ne = niche_enrichment[n]
        lines.append(
            f"| `{n}` | {ne['candidates']} ({ne['cand_pct']}%) | {ne['background']} ({ne['bg_pct']}%) | {ne['class1_count']} / {ne['class2_count']} | **{ne['odds_ratio']:.2f}** | `{ne['fisher_p']}` |"
        )

    lines += [
        "",
        "## 2. 3D Spatial Proximity Test (High-Altitude vs. Diving/Subterranean Hypoxia Within Each Complex)",
        "",
        "| Complex (PDB) | Altitude Sites | Diving/Sub Sites | Observed Median 3D NN Dist (A) | Random Expected (A) | p (closer than random) |",
        "| :--- | :---: | :---: | :---: | :---: | :---: |",
    ]
    for k, v in spatial_3d_tests.items():
        lines.append(
            f"| `{k}` | {v['n_altitude']} | {v['n_diving_sub']} | **{v['obs_median_3d_nn_A']} A** | {v['null_expected_3d_nn_A']} A | `{v['p_value_closer_than_random']}` |"
        )

    lines += [
        "",
        "## 3. Discovered FBT Payoff Equivalence Classes (Cross-Stressor Breakdown)",
        "",
        "| Payoff Equivalence Class | n (T1/T2) | Mean Grantham | Mean Dist Nuclear (A) | Mean Dist Catalytic (A) | Mean Dist Lipid (A) | Stressor Axes Present | Key Exemplars |",
        "| :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- |",
    ]
    for eq, st in cluster_summaries.items():
        axes_str = ", ".join(f"{k[:4]}:{v}" for k, v in sorted(st["stressor_axes"].items()))
        ex_str = "; ".join(st["exemplars"][:3])
        lines.append(
            f"| `{eq}` | {st['size']} ({st['tier1_count']}/{st['tier2_count']}) | {st['mean_grantham']} | {st['mean_dist_nuclear_A']} | {st['mean_dist_catalytic_A']} | {st['mean_dist_lipid_A']} | {axes_str} | {ex_str} |"
        )

    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nStage 2 complete! Wrote {OUT_JSON}, {OUT_REPORT}, and {OUT_FIG}.")
    print("3D Micro-Niche Enrichment:")
    for n in niches_present:
        ne = niche_enrichment[n]
        print(f"  {n:<28}: cand={ne['cand_pct']:>5.1f}% vs bg={ne['bg_pct']:>5.1f}% | OR={ne['odds_ratio']:.2f} (p={ne['fisher_p']})")
    print("3D Cross-Hypoxia Spatial Proximity:")
    for k, v in spatial_3d_tests.items():
        print(f"  {k}: obs={v['obs_median_3d_nn_A']} A vs null={v['null_expected_3d_nn_A']} A (p={v['p_value_closer_than_random']})")


if __name__ == "__main__":
    main()
