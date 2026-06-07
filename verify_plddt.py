"""Verify pLDDT filtering hypothesis without P2Rank.

Hypothesis: filtering low-pLDDT (< 70) regions before pocket detection
reduces false-positive pockets while preserving true binding sites.

Test: For each Phase 1 protein (AF + PDB pair, not in training set):
  1. Parse AF structure → per-residue pLDDT (from B-factor)
  2. Parse PDB structure → ligand positions (ground truth binding sites)
  3. Align structures (via Bio.PDB.Superimposer)
  4. For each AF residue, check:
     a. Is it low-pLDDT (< 70)?
     b. Is the corresponding PDB residue near a ligand (< 6 Å)?
  5. Report: what fraction of binding-site residues would be lost by
     pLDDT filtering, and what fraction of non-binding residues would be
     filtered out.

Usage:
    uv run python verify_plddt.py
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.Polypeptide import is_aa

sys.path.insert(0, str(Path(__file__).parent / "src"))

# ── config ────────────────────────────────────────────────────────────
PHASE1_DIR = Path("data/01_raw/proteins_fresh")
PLDDT_THRESHOLD = 70.0
LIGAND_PROXIMITY_A = 6.0  # Residue is "near ligand" if any atom < this distance

# PDB IDs overlapping with training set (Phase 2 scraped structures)
TRAINING_OVERLAP = {"10DC", "10HY", "11OY", "21LN", "22RH", "24OK", "28MR", "29MO"}


def parse_cif_plddt(cif_path: Path) -> dict:
    """Parse gzipped AF mmCIF, return per-residue pLDDT + CA coords.

    Returns dict keyed by (chain_id, res_id) with (ca_xyz, plddt).
    """
    if cif_path.suffix == ".gz":
        with gzip.open(cif_path, "rt", encoding="utf-8") as f:
            content = f.read()
    else:
        content = cif_path.read_text(encoding="utf-8")

    import io
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("af", io.StringIO(content))

    result = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0].strip():  # HETATM — skip ligands/water
                    continue
                if not is_aa(residue, standard=True):
                    continue
                ca = _find_ca(residue)
                if ca is None:
                    continue
                key = (chain.id, residue.id[1])
                result[key] = {
                    "ca_xyz": np.array(ca.get_coord()),
                    "plddt": ca.get_bfactor(),
                }
    return result


def _find_ca(residue):
    """Safely get CA atom, handling Bio.PDB vagaries."""
    try:
        ca = residue["CA"]
        return ca
    except (KeyError, TypeError):
        pass
    for atom in residue:
        try:
            if atom.get_name() == "CA":
                return atom
        except (AttributeError, TypeError):
            pass
    return None


def parse_pdb_ligands(cif_path: Path) -> dict:
    """Parse gzipped PDB mmCIF, return ligand atom coordinates.

    Returns dict with 'ligand_atoms': (M,3) array, 'residues': dict of
    (chain_id, res_id) -> ca_xyz for alignment.
    """
    if cif_path.suffix == ".gz":
        with gzip.open(cif_path, "rt", encoding="utf-8") as f:
            content = f.read()
    else:
        content = cif_path.read_text(encoding="utf-8")

    import io
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("pdb", io.StringIO(content))

    ligand_atoms = []
    ca_residues = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0].strip()
                if hetflag and hetflag != "W":
                    # Ligand — collect heavy atoms
                    for atom in residue:
                        if atom.element != "H":
                            ligand_atoms.append(atom.get_coord())
                elif not hetflag and is_aa(residue, standard=True):
                    ca = _find_ca(residue)
                    if ca:
                        ca_residues[(chain.id, residue.id[1])] = np.array(ca.get_coord())

    return {
        "ligand_atoms": np.array(ligand_atoms) if ligand_atoms else np.empty((0, 3)),
        "ca_residues": ca_residues,
    }


def align_af_to_pdb(af_data: dict, pdb_data: dict) -> tuple:
    """Align AF CA atoms to PDB CA atoms. Returns rotation + translation."""
    # Find matching residues
    common_keys = sorted(set(af_data.keys()) & set(pdb_data["ca_residues"].keys()))
    if len(common_keys) < 10:
        raise ValueError(f"Too few overlapping residues: {len(common_keys)}")

    af_ca = np.array([af_data[k]["ca_xyz"] for k in common_keys])
    pdb_ca = np.array([pdb_data["ca_residues"][k] for k in common_keys])

    class _Atom:
        def __init__(self, c):
            self._c = c

        def get_coord(self):
            return self._c

    sup = Superimposer()
    sup.set_atoms([_Atom(v) for v in pdb_ca], [_Atom(v) for v in af_ca])
    rot, tran = sup.rotran

    return rot, tran, sup.rms, common_keys


def analyze_protein(protein_dir: Path) -> dict | None:
    """Analyze one protein: pLDDT vs binding-site proximity."""
    af_files = list(protein_dir.glob("AF-*.cif.gz"))
    pdb_files = list(protein_dir.glob("PDB-*.cif.gz"))

    if not af_files or not pdb_files:
        return None

    pdb_id = pdb_files[0].stem.replace("PDB-", "").split(".")[0]
    if pdb_id in TRAINING_OVERLAP:
        return None

    try:
        af_data = parse_cif_plddt(af_files[0])
        pdb_data = parse_pdb_ligands(pdb_files[0])
    except Exception as e:
        return {"protein": protein_dir.name, "pdb_id": pdb_id, "error": str(e)}

    if len(pdb_data["ligand_atoms"]) == 0:
        return {"protein": protein_dir.name, "pdb_id": pdb_id, "error": "no_ligands"}

    # Align
    try:
        rot, tran, rmsd, common_keys = align_af_to_pdb(af_data, pdb_data)
    except ValueError as e:
        return {"protein": protein_dir.name, "pdb_id": pdb_id, "error": str(e)}

    lig_coords = pdb_data["ligand_atoms"]

    # Classify each overlapping residue
    low_plddt_near_lig = 0  # Would be lost (false negative risk)
    low_plddt_far = 0  # Would be correctly filtered (true negative)
    high_plddt_near_lig = 0  # Would be kept (true positive)
    high_plddt_far = 0  # Would be kept (false positive risk)

    for key in common_keys:
        plddt = af_data[key]["plddt"]
        af_ca = af_data[key]["ca_xyz"]
        # Align AF CA to PDB frame
        af_ca_pdb = af_ca @ rot + tran

        # Distance to nearest ligand atom
        dist = float(np.linalg.norm(lig_coords - af_ca_pdb, axis=1).min())

        near_lig = dist <= LIGAND_PROXIMITY_A
        low = plddt < PLDDT_THRESHOLD

        if low and near_lig:
            low_plddt_near_lig += 1
        elif low and not near_lig:
            low_plddt_far += 1
        elif not low and near_lig:
            high_plddt_near_lig += 1
        else:
            high_plddt_far += 1

    total_low = low_plddt_near_lig + low_plddt_far
    total_near = low_plddt_near_lig + high_plddt_near_lig
    total = len(common_keys)

    return {
        "protein": protein_dir.name,
        "pdb_id": pdb_id,
        "rmsd": round(float(rmsd), 2),
        "n_residues": total,
        "n_near_ligand": total_near,
        "n_low_plddt": total_low,
        "n_low_near_lig": low_plddt_near_lig,
        "n_low_far": low_plddt_far,
        "n_high_near_lig": high_plddt_near_lig,
        "n_high_far": high_plddt_far,
        "pct_low": round(100 * total_low / total, 1) if total else 0,
        "pct_near_lost": round(100 * low_plddt_near_lig / total_near, 1) if total_near else 0,
        "pct_far_filtered": round(100 * low_plddt_far / (total - total_near), 1) if (total - total_near) else 0,
    }


def main():
    protein_dirs = sorted(
        d for d in PHASE1_DIR.iterdir() if d.is_dir() and not d.name.startswith(".")
    )
    results = []
    ok = fail = skip = 0

    for i, d in enumerate(protein_dirs):
        r = analyze_protein(d)
        if r is None:
            skip += 1
            continue
        if "error" in r:
            fail += 1
            print(f"  [{i+1}/{len(protein_dirs)}] {r['protein']}: {r['error']}")
            results.append(r)
        else:
            ok += 1
            print(f"  [{i+1}/{len(protein_dirs)}] {r['protein']}: "
                  f"RMSD={r['rmsd']:.1f}Å, "
                  f"low-pLDDT={r['n_low_plddt']}/{r['n_residues']} ({r['pct_low']}%), "
                  f"near-lig lost={r['n_low_near_lig']}/{r['n_near_ligand']} ({r['pct_near_lost']}%), "
                  f"far filtered={r['n_low_far']}/{r['n_residues']-r['n_near_ligand']} ({r['pct_far_filtered']}%)")
            results.append(r)

    # ── aggregate ──────────────────────────────────────────────────────
    ok_results = [r for r in results if "error" not in r]
    print(f"\n{'='*60}")
    print(f"SUMMARY: {ok} OK, {fail} failed, {skip} skipped ({len(protein_dirs)} total)")

    if ok_results:
        df = pd.DataFrame.from_records(ok_results)
        total_low = df["n_low_near_lig"].sum() + df["n_low_far"].sum()
        total_near = df["n_low_near_lig"].sum() + df["n_high_near_lig"].sum()
        total_res = df["n_residues"].sum()
        total_far = total_res - total_near

        print(f"\n  Total residues analyzed:      {total_res}")
        print(f"  Low pLDDT (< {PLDDT_THRESHOLD}):           {total_low} ({100*total_low/total_res:.1f}%)")
        print(f"  Near ligand (< {LIGAND_PROXIMITY_A}Å):          {total_near} ({100*total_near/total_res:.1f}%)")
        print(f"  ─────────────────────────────────────")
        print(f"  Low pLDDT, near ligand  (LOST):   {df['n_low_near_lig'].sum()} ({100*df['n_low_near_lig'].sum()/total_near:.1f}% of binding)")
        print(f"  Low pLDDT, far          (FILTERED): {df['n_low_far'].sum()} ({100*df['n_low_far'].sum()/total_far:.1f}% of non-binding)")
        print(f"  High pLDDT, near ligand (KEPT):    {df['n_high_near_lig'].sum()}")
        print(f"  High pLDDT, far         (KEPT):    {df['n_high_far'].sum()}")

        precision_gain = 100 * df["n_low_far"].sum() / max(total_low, 1)
        recall_loss = 100 * df["n_low_near_lig"].sum() / max(total_near, 1)
        print(f"\n  Filtering removes {precision_gain:.1f}% non-binding residues")
        print(f"  Filtering loses    {recall_loss:.1f}% binding-site residues")
        print(f"  Net benefit ratio:   {precision_gain/recall_loss:.1f}x" if recall_loss > 0 else "")

        # Save
        out_path = Path("data/08_reporting/plddt_verification.json")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
