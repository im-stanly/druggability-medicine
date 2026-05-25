"""mmCIF parser that extracts protein atoms and ligand atoms from PDB structures.

Each PDB mmCIF file contains entities of type:
  - polymer: protein chains (standard amino acids)
  - non-polymer: ligands, ions, cofactors
  - water: solvent molecules

We separate protein atoms (for surface generation) from ligand atoms
(for labeling binding pockets).
"""

import gzip
import io
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Structure import Structure

from .constants import STANDARD_AAS, IGNORE_RESIDUES


@dataclass
class ProteinAtoms:
    """Container for protein atom arrays extracted from mmCIF."""
    coords: np.ndarray       # (N, 3) -- xyz coordinates
    atom_names: list[str]    # (N,)   -- atom names
    residue_names: list[str] # (N,)   -- 3-letter residue codes
    residue_ids: list[int]   # (N,)   -- residue sequence numbers
    chain_ids: list[str]     # (N,)   -- chain identifiers
    elements: list[str]      # (N,)   -- element symbols
    sasa: np.ndarray | None = None  # (N,) -- per-atom SASA (computed later)

    def __len__(self) -> int:
        return len(self.coords)


@dataclass
class LigandInfo:
    """Container for ligand data extracted from mmCIF."""
    name: str                # Compound name
    coords: np.ndarray       # (M, 3) -- heavy atom coordinates
    elements: list[str]      # (M,)   -- element symbols
    center: np.ndarray       # (3,)   -- geometric center
    formula: str = ""        # Chemical formula if available
    n_atoms: int = 0

    def __post_init__(self):
        self.n_atoms = len(self.coords)
        if self.center is None:
            self.center = self.coords.mean(axis=0)


@dataclass
class ParsedStructure:
    """Complete parsed structure: protein atoms + list of ligands."""
    pdb_id: str
    protein: ProteinAtoms
    ligands: list[LigandInfo] = field(default_factory=list)
    raw_structure: Structure | None = None


def _is_ligand_residue(resname: str, hetfield: str) -> bool:
    """Check if a residue is a ligand (non-polymer, non-water, non-buffer).

    In mmCIF, hetero residues have hetfield like 'H_WPW' and standard
    amino acids have hetfield ' ' (blank).
    """
    resname = resname.strip()
    if resname in STANDARD_AAS:
        return False
    if resname in IGNORE_RESIDUES:
        return False
    if hetfield == "W":  # water
        return False
    return True


def _is_metal_ion(resname: str) -> bool:
    """Check if a residue is a single metal ion (not a complex ligand)."""
    resname = resname.strip()
    metals = {"MG", "CA", "ZN", "FE", "MN", "CU", "NA", "K", "CO", "NI", "CD", "HG"}
    return resname in metals


def parse_cif(filepath: str | Path) -> ParsedStructure:
    """Parse a gzipped mmCIF file and extract protein + ligand atoms.

    Args:
        filepath: Path to a .cif.gz file from the PDB.

    Returns:
        ParsedStructure with protein atoms and ligand info separated.
    """
    filepath = Path(filepath)

    # Decompress and parse
    if filepath.suffix == ".gz":
        with gzip.open(filepath, "rt", encoding="utf-8") as f:
            content = f.read()
    else:
        content = filepath.read_text(encoding="utf-8")

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("struct", io.StringIO(content))

    pdb_id = structure.id if structure.id != "struct" else filepath.stem.split(".")[0]

    # Collect atoms
    prot_coords, prot_names, prot_resnames = [], [], []
    prot_resids, prot_chains, prot_elements = [], [], []

    ligands: list[LigandInfo] = []

    model = structure[0]
    for chain in model:
        for residue in chain:
            resname = residue.resname.strip()
            hetfield = str(residue.id[0]).strip()

            if _is_ligand_residue(resname, hetfield):
                # Collect ligand atoms (heavy atoms only)
                atoms = [a for a in residue.get_atoms() if a.element != "H"]
                if not atoms:
                    continue
                coords = np.array([a.get_coord() for a in atoms])
                elements = [a.element for a in atoms]

                # Get compound name from the residue
                lig = LigandInfo(
                    name=resname,
                    coords=coords,
                    elements=elements,
                    center=coords.mean(axis=0),
                )
                if not _is_metal_ion(resname):
                    ligands.append(lig)

            elif resname in STANDARD_AAS:
                # Collect protein atoms (all atoms including H)
                for atom in residue:
                    prot_coords.append(atom.get_coord())
                    prot_names.append(atom.get_name())
                    prot_resnames.append(resname)
                    prot_resids.append(residue.id[1])
                    prot_chains.append(chain.id)
                    prot_elements.append(atom.element)

    if not prot_coords:
        raise ValueError(f"No protein atoms found in {filepath}")

    protein = ProteinAtoms(
        coords=np.array(prot_coords),
        atom_names=prot_names,
        residue_names=prot_resnames,
        residue_ids=prot_resids,
        chain_ids=prot_chains,
        elements=prot_elements,
    )

    return ParsedStructure(pdb_id=pdb_id, protein=protein, ligands=ligands,
                           raw_structure=structure)


def get_ca_atoms(protein: ProteinAtoms) -> np.ndarray:
    """Extract coordinates of alpha-carbon atoms only."""
    mask = np.array([name == "CA" for name in protein.atom_names])
    return protein.coords[mask]


def get_ligand_centers(parsed: ParsedStructure) -> np.ndarray:
    """Return (L, 3) array of ligand geometric centers."""
    if not parsed.ligands:
        return np.empty((0, 3))
    return np.array([lig.center for lig in parsed.ligands])


def get_ligand_atoms(parsed: ParsedStructure) -> np.ndarray:
    """Return (M, 3) array of all ligand heavy-atom coordinates."""
    all_coords = []
    for lig in parsed.ligands:
        all_coords.append(lig.coords)
    if not all_coords:
        return np.empty((0, 3))
    return np.vstack(all_coords)
