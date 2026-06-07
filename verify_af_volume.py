"""Hypothesis: AlphaFold pockets are smaller than PDB pockets (no induced fit).

For each protein with both AF and PDB:
  1. Align AF → PDB
  2. Find ligand positions in PDB
  3. Count protein atoms within 8Å of each ligand in both structures
  4. Compare volumes

Usage:
    uv run python verify_af_volume.py
"""

from __future__ import annotations

import gzip, io, json, sys
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.Polypeptide import is_aa

sys.path.insert(0, str(Path(__file__).parent / "src"))

TEST_DIRS = [
    Path("data/01_raw/proteins"),        # Phase 1
    Path("data/01_raw/proteins_fresh"),  # Fresh test data
]
RADIUS = 8.0  # Å — pocket radius around ligand

SKIP_LIGANDS = {"HOH", "DOD", "WAT", "GOL", "EDO", "SO4", "PO4", "ACT", "PEG",
                "BME", "DMS", "FMT", "EPE", "CIT", "TRS", "MPD",
                "MG", "CL", "NA", "K", "CA", "ZN", "MN", "FE", "CO", "NI", "CD", "HG"}


def parse_cif(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as f:
            content = f.read()
    else:
        content = path.read_text(encoding="utf-8")
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure("s", io.StringIO(content))


def get_ca_atoms(structure):
    """Get (chain_id, res_id) → CA coord for all standard residues."""
    result = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0].strip():
                    continue
                if not is_aa(res, standard=True):
                    continue
                try:
                    ca = res["CA"]
                except KeyError:
                    continue
                result[(chain.id, res.id[1])] = np.array(ca.get_coord())
    return result


def get_ligand_centers(structure):
    """Get positions of non-solvent ligands."""
    ligands = []
    for model in structure:
        for chain in model:
            for res in chain:
                if not res.id[0].strip():
                    continue
                if res.resname.strip() in SKIP_LIGANDS:
                    continue
                coords = []
                for atom in res:
                    if atom.element != "H":
                        coords.append(atom.get_coord())
                if coords:
                    ligands.append({
                        "name": res.resname.strip(),
                        "center": np.array(coords).mean(axis=0),
                    })
    return ligands


def count_nearby_atoms(structure, center_xyz, radius=RADIUS):
    """Count protein atoms within radius of a point."""
    center = np.array(center_xyz)
    count = 0
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0].strip():  # HETATM — skip ligands
                    continue
                for atom in res:
                    if np.linalg.norm(np.array(atom.get_coord()) - center) <= radius:
                        count += 1
    return count


def analyze_protein_dir(protein_dir: Path) -> list[dict]:
    results = []
    af_files = list(protein_dir.glob("AF-*.cif.gz"))
    pdb_files = list(protein_dir.glob("PDB-*.cif.gz"))
    if not af_files or not pdb_files:
        return results

    try:
        af = parse_cif(af_files[0])
        pdb = parse_cif(pdb_files[0])
    except Exception:
        return results

    ligands = get_ligand_centers(pdb)
    if not ligands:
        return results

    # Align
    af_ca = get_ca_atoms(af)
    pdb_ca = get_ca_atoms(pdb)
    common = sorted(set(af_ca) & set(pdb_ca))
    if len(common) < 10:
        return results

    af_xyz = np.array([af_ca[k] for k in common])
    pdb_xyz = np.array([pdb_ca[k] for k in common])

    class A:
        def __init__(self, c): self._c = c
        def get_coord(self): return self._c
    sup = Superimposer()
    sup.set_atoms([A(v) for v in pdb_xyz], [A(v) for v in af_xyz])
    rot, tran = sup.rotran

    # For each ligand, count atoms in PDB and AF
    for lig in ligands:
        pdb_count = count_nearby_atoms(pdb, lig["center"])
        # Transform AF atoms? No — just check AF structure at the same position
        # (structures are aligned, coordinates are in PDB frame now)
        af_count = 0
        for model in af:
            for chain in model:
                for res in chain:
                    if res.id[0].strip():
                        continue
                    for atom in res:
                        af_coord = np.array(atom.get_coord()) @ rot + tran
                        if np.linalg.norm(af_coord - lig["center"]) <= RADIUS:
                            af_count += 1

        results.append({
            "protein": protein_dir.name,
            "ligand": lig["name"],
            "pdb_atoms": pdb_count,
            "af_atoms": af_count,
            "ratio": round(af_count / pdb_count, 3) if pdb_count > 0 else None,
        })

    return results


def main():
    all_results = []
    for test_dir in TEST_DIRS:
        if not test_dir.exists():
            continue
        for protein_dir in sorted(test_dir.iterdir()):
            if not protein_dir.is_dir():
                continue
            res = analyze_protein_dir(protein_dir)
            if res:
                all_results.extend(res)

    if not all_results:
        print("No results — no ligand-containing protein pairs found.")
        return

    ratios = [r["ratio"] for r in all_results if r["ratio"] is not None]
    pdb_vols = [r["pdb_atoms"] for r in all_results]
    af_vols = [r["af_atoms"] for r in all_results]

    print(f"  Ligand sites analyzed: {len(all_results)}")
    print(f"  Mean PDB volume (atoms):  {np.mean(pdb_vols):.1f}")
    print(f"  Mean AF volume (atoms):   {np.mean(af_vols):.1f}")
    print(f"  Mean AF/PDB ratio:        {np.mean(ratios):.3f}")
    print(f"  AF < PDB: {sum(1 for r in ratios if r < 1.0)}/{len(ratios)} ({100*sum(1 for r in ratios if r < 1.0)/len(ratios):.0f}%)")

    # Statistical test
    from scipy.stats import wilcoxon
    stat, pval = wilcoxon(af_vols, pdb_vols, alternative="less")
    print(f"\n  Wilcoxon signed-rank (AF < PDB): p = {pval:.4f}")

    if pval < 0.05:
        print(f"  ✓ Significant — AF pockets are smaller (p < 0.05)")
    else:
        print(f"  ✗ Not significant")

    # Save
    out = Path("data/08_reporting/af_volume_comparison.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump({
            "n_sites": int(len(all_results)),
            "mean_pdb_atoms": float(np.mean(pdb_vols)),
            "mean_af_atoms": float(np.mean(af_vols)),
            "mean_ratio": float(np.mean(ratios)),
            "af_smaller_pct": float(100 * sum(1 for r in ratios if r < 1.0) / len(ratios)),
            "wilcoxon_p": float(pval),
            "significant": bool(pval < 0.05),
            "per_ligand": all_results,
        }, f, indent=2, default=str)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
