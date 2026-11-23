#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "openmm",
#     "pdbfixer",
#     "biopython",
#     "numpy",
#     "scipy",
# ]
# ///
"""Stage 6C: Physics-Based OpenMM AMBER14 + OBC2 Implicit-Solvent Local Energy Relaxation.

Directly mutates and energy-minimizes complete interacting subunit assemblies extracted from
experimental human OXPHOS cryo-EM structures (9I6F, 9I4I, 9HZL, 8H9S) for all 5 Stage-3 Finalists
plus 2 Human Pathogenic Controls (ATP6 L156R, ND1 A52T), measuring:
  - Initial vs. Minimized AMBER14+OBC2 potential energy (WT vs. Mutant strain relief)
  - Post-relaxation 8 A pocket C-alpha & heavy-atom RMSD (structural perturbation)
  - Post-relaxation minimum non-bonded heavy-atom distance (0-clash verification with explicit H)
  - Newly formed polar H-bond / salt-bridge contacts (<= 3.5 A)

Saves minimized WT/Mutant PDBs to data/stage6_openmm/ and summary to data/stage6c_openmm_summary.json.
"""

from __future__ import annotations

import gzip
import io
import json
import math
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform, unit
from openmm.app import ForceField, Modeller, NoCutoff, PDBFile, Simulation
from pdbfixer import PDBFixer
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
STRUCT_DIR = ROOT / "data" / "structures"
OUT_DIR = ROOT / "data" / "stage6_openmm"
OUT_JSON = ROOT / "data" / "stage6c_openmm_summary.json"

AA1_TO_3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

STANDARD_AA3 = set(AA1_TO_3.values())

TARGETS = [
    ("COX3_W116R", "COX3", 116, "W", "R", "9I6F", STRUCT_DIR / "pdb_00009i6f.cif.gz", "Stage3_Finalist_1"),
    ("ND4L_A11T",  "ND4L",  11, "A", "T", "9I4I", STRUCT_DIR / "pdb_00009i4i.cif.gz", "Stage3_Finalist_2"),
    ("COX3_I165T", "COX3", 165, "I", "T", "9I6F", STRUCT_DIR / "pdb_00009i6f.cif.gz", "Stage3_Finalist_3"),
    ("ATP6_I106V", "ATP6", 106, "I", "V", "8H9S", STRUCT_DIR / "pdb_00008h9s.cif.gz", "Stage3_Finalist_4"),
    ("CYTB_M53T",  "CYTB",  53, "M", "T", "9HZL", STRUCT_DIR / "pdb_00009hzl.cif.gz", "Stage3_Finalist_5"),
    ("ATP6_L156R", "ATP6", 156, "L", "R", "8H9S", STRUCT_DIR / "pdb_00008h9s.cif.gz", "Pathogenic_Control_Leigh"),
    ("ND1_A52T",   "ND1",   52, "A", "T", "9I4I", STRUCT_DIR / "pdb_00009i4i.cif.gz", "Pathogenic_Control_LHON"),
]


def load_chain_gene_map() -> dict[str, dict[str, str]]:
    meta = json.loads((STRUCT_DIR / "oxphos_entities_meta.json").read_text(encoding="utf-8"))
    u2h = json.loads((STRUCT_DIR / "uniprot_to_hgnc.json").read_text(encoding="utf-8"))
    mt_map = {
        "MT-ND1": "ND1", "MT-ND2": "ND2", "MT-ND3": "ND3", "MT-ND4": "ND4",
        "MT-ND4L": "ND4L", "MT-ND5": "ND5", "MT-ND6": "ND6", "MT-CYB": "CYTB",
        "MT-CO1": "COX1", "MT-CO2": "COX2", "MT-CO3": "COX3", "MT-ATP6": "ATP6", "MT-ATP8": "ATP8",
    }
    out = {}
    for entry in meta["data"]["entries"]:
        pid = entry["rcsb_id"]
        cmap = {}
        for pe in entry["polymer_entities"]:
            uids = pe["rcsb_polymer_entity_container_identifiers"].get("uniprot_ids") or []
            chains = pe["rcsb_polymer_entity_container_identifiers"].get("auth_asym_ids") or []
            sym = u2h.get(uids[0], "") if uids else ""
            sym = mt_map.get(sym, sym)
            for c in chains:
                cmap[c] = sym or c
        out[pid] = cmap
    return out


class ChainSubsetSelect(Select):
    def __init__(self, allowed_chains: set[str]):
        self.allowed_chains = allowed_chains

    def accept_chain(self, chain):
        return chain.id in self.allowed_chains

    def accept_residue(self, residue):
        return residue.id[0] == " " and residue.get_resname().strip() in STANDARD_AA3


def extract_interacting_subassembly(
    cif_gz: Path,
    pdb_id: str,
    target_gene: str,
    pos1: int,
    chain_map: dict[str, str],
) -> tuple[str, str, int, list[str], np.ndarray]:
    """Extract target chain + every protein chain with an atom within 7.5 A of target residue."""
    parser = MMCIFParser(QUIET=True)
    with gzip.open(cif_gz, "rt", encoding="utf-8") as fh:
        model = next(parser.get_structure(pdb_id, fh).get_models())

    target_cid = next((c for c, g in chain_map.items() if g == target_gene and c in model), None)
    assert target_cid is not None, f"Chain for {target_gene} not found in {pdb_id}"

    # Find exact residue
    t_res = None
    for res in model[target_cid]:
        if res.id[0] == " " and res.id[1] == pos1:
            t_res = res
            break
    if t_res is None:
        std_r = [r for r in model[target_cid] if r.id[0] == " "]
        t_res = std_r[pos1 - 1]

    t_coords = np.array([a.coord for a in t_res.get_atoms() if a.element != "H"], dtype=float)
    t_ca = t_res["CA"].coord.astype(float)
    t_tree = cKDTree(t_coords)

    # Identify chains with heavy atoms within 7.5 A of target residue (cap at 3 nearest neighbors for speed)
    neighbor_dists = []
    for ch in model:
        if ch.id == target_cid:
            continue
        coords = [a.coord for r in ch if r.id[0] == " " and r.get_resname().strip() in STANDARD_AA3 for a in r.get_atoms() if a.element != "H"]
        if not coords:
            continue
        dmin = float(np.min(t_tree.query(np.array(coords, dtype=float))[0]))
        if dmin <= 7.5:
            neighbor_dists.append((dmin, ch.id))
    neighbor_dists.sort()
    selected_chains = {target_cid} | {cid for _, cid in neighbor_dists[:3]}

    # Map PDB chain IDs to clean single-character chain IDs if needed & serialize to PDB text
    # Note: Biopython PDBIO requires 1-char chain IDs (all our cryo-EM chain IDs here are 1-char)
    buf = io.StringIO()
    io_writer = PDBIO()
    io_writer.set_structure(model)
    io_writer.save(buf, ChainSubsetSelect(selected_chains))
    partner_syms = [f"{chain_map.get(cid, cid)}({cid})" for cid in sorted(selected_chains)]
    return buf.getvalue(), target_cid, int(t_res.id[1]), partner_syms, t_ca


def minimize_structure(
    pdb_text: str,
    target_chain: str,
    resseq: int,
    ref_aa: str,
    alt_aa: str | None,
    site_ca_nm: np.ndarray,
    out_pdb: Path,
) -> dict:
    """Fix missing atoms/H at pH 7.4 via PDBFixer, optionally mutate, and minimize in OpenMM AMBER14+OBC2."""
    fixer = PDBFixer(pdbfile=io.StringIO(pdb_text))
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    # Do not insert long unresolved cryo-EM terminal loops — keep experimental resolved core
    fixer.missingResidues = {}
    if alt_aa is not None:
        mut_str = f"{AA1_TO_3[ref_aa]}-{resseq}-{AA1_TO_3[alt_aa]}"
        fixer.applyMutations([mut_str], target_chain)
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
    system = ff.createSystem(fixer.topology, nonbondedMethod=NoCutoff, constraints=None)

    # Soft harmonic positional restraints on CA atoms > 10.0 A (1.0 nm) from mutated site
    # Everything within 10.0 A (including mutated side chain + surrounding pocket) relaxes 100% freely!
    restraint = CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint.addGlobalParameter("k", 500.0 * unit.kilojoules_per_mole / unit.nanometer**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    pos_nm = np.array([[p.x, p.y, p.z] for p in fixer.positions.value_in_unit(unit.nanometer)])
    atom_meta = []
    for atom in fixer.topology.atoms():
        idx = atom.index
        res = atom.residue
        ch_id = res.chain.id
        r_id = int(res.id)
        is_heavy = atom.element is not None and atom.element.symbol != "H"
        d_site = float(np.linalg.norm(pos_nm[idx] - site_ca_nm))
        atom_meta.append((idx, ch_id, r_id, res.name, atom.name, atom.element.symbol if atom.element else "", is_heavy, d_site))
        if atom.name == "CA" and d_site > 1.0:
            restraint.addParticle(idx, pos_nm[idx].tolist())
    system.addForce(restraint)

    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    try:
        platform = Platform.getPlatformByName("CPU")
    except Exception:
        platform = Platform.getPlatformByName("Reference")
    sim = Simulation(fixer.topology, system, integrator, platform)
    sim.context.setPositions(fixer.positions)

    e_init = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    sim.minimizeEnergy(tolerance=10.0 * unit.kilojoules_per_mole / unit.nanometer, maxIterations=800)
    state_fin = sim.context.getState(getEnergy=True, getPositions=True)
    e_fin = state_fin.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    pos_fin_nm = np.array([[p.x, p.y, p.z] for p in state_fin.getPositions().value_in_unit(unit.nanometer)])

    with out_pdb.open("w", encoding="utf-8") as fh:
        PDBFile.writeFile(fixer.topology, state_fin.getPositions(), fh)

    # Analyze relaxed geometry in Angstroms (1 nm = 10 A)
    pos_init_A = pos_nm * 10.0
    pos_fin_A = pos_fin_nm * 10.0

    # Target site CA displacement & 8 A pocket heavy-atom RMSD (excluding mutated sidechain atoms)
    site_ca_idx = next(m[0] for m in atom_meta if m[1] == target_chain and m[2] == resseq and m[4] == "CA")
    ca_shift_A = float(np.linalg.norm(pos_fin_A[site_ca_idx] - pos_init_A[site_ca_idx]))

    pocket_idxs = [m[0] for m in atom_meta if m[6] and m[7] <= 0.8 and not (m[1] == target_chain and m[2] == resseq and m[4] not in {"N", "CA", "C", "O"})]
    pocket_rmsd_A = float(np.sqrt(np.mean(np.sum((pos_fin_A[pocket_idxs] - pos_init_A[pocket_idxs]) ** 2, axis=1))))

    # Post-relaxation minimum non-bonded heavy-atom distance from target side chain to environment
    sc_idxs = [m[0] for m in atom_meta if m[1] == target_chain and m[2] == resseq and m[6] and m[4] not in {"N", "CA", "C", "O"}]
    env_idxs = [
        m[0] for m in atom_meta
        if m[6] and not (m[1] == target_chain and m[2] == resseq)
        and not (m[1] == target_chain and abs(m[2] - resseq) == 1 and m[4] in {"N", "CA", "C", "O"})
    ]
    if sc_idxs and env_idxs:
        env_tree = cKDTree(pos_fin_A[env_idxs])
        d_sc, _ = env_tree.query(pos_fin_A[sc_idxs])
        dmin_sc_A = float(np.min(d_sc))
    else:
        dmin_sc_A = 4.0

    # Polar H-bond / salt-bridge contacts (<= 3.5 A between target sidechain N/O and env N/O)
    sc_polar = [m for m in atom_meta if m[1] == target_chain and m[2] == resseq and m[5] in {"N", "O"} and m[4] not in {"N", "O"}]
    env_polar = [m for m in atom_meta if m[5] in {"N", "O"} and not (m[1] == target_chain and m[2] == resseq)]
    polar_contacts = []
    if sc_polar and env_polar:
        ep_coords = np.array([pos_fin_A[m[0]] for m in env_polar])
        ep_tree = cKDTree(ep_coords)
        for sp in sc_polar:
            hits = ep_tree.query_ball_point(pos_fin_A[sp[0]], r=3.5)
            for h in hits:
                ep = env_polar[h]
                d = float(np.linalg.norm(pos_fin_A[sp[0]] - pos_fin_A[ep[0]]))
                polar_contacts.append(f"{sp[4]}--{ep[1]}:{ep[3]}{ep[2]}.{ep[4]}({d:.2f}A)")

    return {
        "n_atoms": len(atom_meta),
        "E_init_kJ_mol": round(float(e_init), 1),
        "E_min_kJ_mol": round(float(e_fin), 1),
        "relaxation_drop_kJ_mol": round(float(e_fin - e_init), 1),
        "pocket_8A_heavy_rmsd_A": round(pocket_rmsd_A, 3),
        "site_ca_shift_A": round(ca_shift_A, 3),
        "post_relax_sc_dmin_A": round(dmin_sc_A, 2),
        "polar_contacts_3p5A": sorted(set(polar_contacts))[:8],
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    chain_maps = load_chain_gene_map()
    results = []

    for slug, gene, pos1, ref_aa, alt_aa, pid, cif_path, role in TARGETS:
        print(f"=== OpenMM AMBER14+OBC2 Relaxation: {slug} ({role}, PDB {pid}) ===", flush=True)
        pdb_text, t_cid, resseq, partners, t_ca_A = extract_interacting_subassembly(
            cif_path, pid, gene, pos1, chain_maps[pid]
        )
        t_ca_nm = t_ca_A / 10.0
        wt_res = minimize_structure(pdb_text, t_cid, resseq, ref_aa, None, t_ca_nm, OUT_DIR / f"{slug}_WT_min.pdb")
        mut_res = minimize_structure(pdb_text, t_cid, resseq, ref_aa, alt_aa, t_ca_nm, OUT_DIR / f"{slug}_MUT_min.pdb")

        delta_rmsd = round(mut_res["pocket_8A_heavy_rmsd_A"] - wt_res["pocket_8A_heavy_rmsd_A"], 3)
        delta_ca = round(mut_res["site_ca_shift_A"] - wt_res["site_ca_shift_A"], 3)
        new_polars = sorted(set(mut_res["polar_contacts_3p5A"]) - set(wt_res["polar_contacts_3p5A"]))
        rec = {
            "variant": slug,
            "role": role,
            "pdb_id": pid,
            "target_chain": t_cid,
            "resseq": resseq,
            "subassembly_chains": partners,
            "wt": wt_res,
            "mutant": mut_res,
            "delta_pocket_rmsd_A": delta_rmsd,
            "delta_site_ca_shift_A": delta_ca,
            "new_mutant_polar_contacts_3p5A": new_polars,
            "physics_verdict": "STABLE_LOW_PERTURBATION" if (mut_res["post_relax_sc_dmin_A"] >= 2.6 and abs(delta_rmsd) <= 0.15) else "MODERATE_LOCAL_REMODELING",
        }
        results.append(rec)
        print(
            f"  Chains={partners} | WT RMSD={wt_res['pocket_8A_heavy_rmsd_A']} A -> MUT RMSD={mut_res['pocket_8A_heavy_rmsd_A']} A "
            f"(Δ={delta_rmsd:+.3f} A) | MUT d_min={mut_res['post_relax_sc_dmin_A']} A | New H-bonds={len(new_polars)} {new_polars[:2]}",
            flush=True,
        )

    OUT_JSON.write_text(json.dumps({"forcefield": "AMBER14 + OBC2 implicit solvent", "results": results}, indent=2), encoding="utf-8")
    print(f"\nOpenMM relaxation complete -> {OUT_JSON}")


if __name__ == "__main__":
    main()
