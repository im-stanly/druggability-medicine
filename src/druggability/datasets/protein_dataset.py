import gzip
import fsspec
from Bio.PDB import MMCIFParser, Structure, MMCIFIO, Select
from kedro.io import AbstractDataset
from typing import Any, Dict
import os

import io
from pathlib import Path
from typing import Any, Generator

from kedro.io import AbstractDataset


class ResidueSelect(Select):
    """Filter to keep only the specifically matched residues."""
    def __init__(self, keep_residues):
        self.keep_ids = {res.get_full_id() for res in keep_residues}

    def accept_residue(self, residue):
        if residue.get_full_id() in self.keep_ids:
            return 1
        return 0


class MmcifPairedDataset(AbstractDataset):
    def __init__(
        self,
        root_path: str,
        pdb_filename: str = "PDB.cif",
        af_filename: str = "AF.cif",
        credentials: dict[str, Any] | None = None,
    ) -> None:
        """
        Args:
            root_path:     Path to the directory containing one sub-folder per protein.
            pdb_filename:  Name of the PDB mmCIF file inside each protein folder.
            af_filename:   Name of the AlphaFold mmCIF file inside each protein folder.
            credentials:   Unused; kept for Kedro compatibility.
        """
        self._root = Path(root_path)
        self._pdb_filename = pdb_filename
        self._af_filename = af_filename
        self._parser = MMCIFParser(QUIET=True)

    def _load(self) -> Generator[dict[str, Any], None, None]:
        """Yield one dict per protein folder found under root_path."""
        protein_dirs = sorted(
            p for p in self._root.iterdir() if p.is_dir()
        )

        if not protein_dirs:
            raise FileNotFoundError(
                f"No protein sub-directories found in '{self._root}'"
            )

        for protein_dir in protein_dirs:
            pdb_path = protein_dir / self._pdb_filename
            af_path = protein_dir / self._af_filename

            self._assert_exists(pdb_path)
            self._assert_exists(af_path)

            protein_name = protein_dir.name

            yield {
                "name": protein_name,
                "pdb": self._parse(pdb_path, structure_id=f"{protein_name}_pdb"),
                "af": self._parse(af_path, structure_id=f"{protein_name}_af"),
            }

    def _save(self, data: list[dict]) -> None:
        os.makedirs(self._root, exist_ok=True)
        for protein_dict in data:
            dirname = protein_dict["name"]
            pdb_seq = protein_dict["pdb"]
            af_seq = protein_dict["af"]
            pdb_keep = protein_dict.get("keep_pdb")
            af_keep = protein_dict.get("keep_af")
            protein_root_path = self._root / dirname
            os.makedirs(protein_root_path, exist_ok=True)

            pdb_save_path = protein_root_path / self._pdb_filename
            af_save_path = protein_root_path / self._af_filename

            io = MMCIFIO()
            io.set_structure(pdb_seq)

            if pdb_keep:
                io.save(os.path.join(protein_root_path, self._pdb_filename), select=ResidueSelect(pdb_keep))
            else:
                io.save(os.path.join(protein_root_path, self._pdb_filename))

            io.set_structure(af_seq)
            if af_keep:
                io.save(os.path.join(protein_root_path, self._af_filename), select=ResidueSelect(af_keep))
            else:
                io.save(os.path.join(protein_root_path, self._af_filename))


    def _describe(self) -> dict[str, Any]:
        return {
            "root_path": str(self._root),
            "pdb_filename": self._pdb_filename,
            "af_filename": self._af_filename,
            "n_proteins": self._count_proteins(),
        }

    def _parse(self, path: Path, structure_id: str) -> Structure:
        return self._parser.get_structure(structure_id, str(path))

    def _assert_exists(self, path: Path) -> None:
        if not path.exists():
            raise FileNotFoundError(
                f"Expected mmCIF file not found: '{path}'. "
                f"Each protein folder must contain '{self._pdb_filename}' "
                f"and '{self._af_filename}'."
            )

    def _count_proteins(self) -> int:
        if not self._root.exists():
            return 0
        return sum(1 for p in self._root.iterdir() if p.is_dir())
    
   