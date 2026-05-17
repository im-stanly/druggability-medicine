import gzip
import io
from pathlib import Path
from typing import Any, BinaryIO, Generator, TextIO, cast

from Bio.PDB import MMCIFParser
from Bio.PDB.Structure import Structure
from kedro.io import AbstractDataset


class ProteinLigandPairsDataset(AbstractDataset):
    """Load protein-ligand pairs from a folder structure like:

    `root_path/<pdb_id>/PDB-<pdb_id>.cif[.gz]`
    `root_path/<pdb_id>/LIG-<ligand_id>.cif[.gz]`

    One item is yielded per `<pdb_id>` directory.

    Output shape (per yielded item):
        {
          "pdb_id": str,
          "protein": Bio.PDB.Structure.Structure,
          "ligands": {ligand_id: Bio.PDB.Structure.Structure, ...},
          "paths": {"protein": Path, "ligands": {ligand_id: Path, ...}},
        }

    Notes:
      * This dataset is read-only.
      * Supports both plain `.cif` and gzipped `.cif.gz` files.
      * By default, it will raise if a protein folder doesn't contain a protein
        structure or doesn't contain any ligands.
    """

    def __init__(
        self,
        root_path: str,
        protein_prefix: str = "PDB-",
        ligand_prefix: str = "LIG-",
        cif_ext: str = ".cif",
        allow_missing_ligands: bool = False,
        credentials: dict[str, Any] | None = None,
    ) -> None:
        self._root = Path(root_path)
        self._protein_prefix = protein_prefix
        self._ligand_prefix = ligand_prefix
        self._cif_ext = cif_ext
        self._allow_missing_ligands = allow_missing_ligands
        self._parser = MMCIFParser(QUIET=True)
        # credentials kept for Kedro compatibility

    def _load(self) -> Generator[dict[str, Any], None, None]:
        protein_dirs = sorted(p for p in self._root.iterdir() if p.is_dir())
        if not protein_dirs:
            raise FileNotFoundError(
                f"No protein-ligand sub-directories found in '{self._root}'"
            )

        for protein_dir in protein_dirs:
            pdb_id = protein_dir.name

            protein_path = self._find_single(
                protein_dir,
                stem=f"{self._protein_prefix}{pdb_id}",
                exts=(self._cif_ext, f"{self._cif_ext}.gz"),
            )
            if protein_path is None:
                raise FileNotFoundError(
                    f"Expected protein mmCIF file not found in '{protein_dir}'. "
                    f"Looked for '{self._protein_prefix}{pdb_id}{self._cif_ext}' "
                    f"or gzipped variant."
                )

            ligand_paths = self._find_many_with_prefix(
                protein_dir,
                prefix=self._ligand_prefix,
                exts=(self._cif_ext, f"{self._cif_ext}.gz"),
            )

            if not ligand_paths and not self._allow_missing_ligands:
                raise FileNotFoundError(
                    f"No ligand mmCIF files found in '{protein_dir}'. "
                    f"Expected files like '{self._ligand_prefix}<ID>{self._cif_ext}[.gz]'."
                )

            protein_structure = self._parse(protein_path, structure_id=f"{pdb_id}_protein")

            ligands: dict[str, dict[str, Any]] = {}
            for lig_id, lig_path in ligand_paths.items():
                complex_centroid = self._ligand_centroid_from_complex(
                    protein_structure, lig_id=lig_id
                )
                if complex_centroid is None:
                    raise ValueError(
                        f"Ligand '{lig_id}' has a LIG file ('{lig_path.name}') but was not found "
                        f"as a HET residue in the complex file '{protein_path.name}'."
                    )

                ligands[lig_id] = {
                    "id": lig_id,
                    "centroid": complex_centroid,
                    "source": "complex",
                    "complex_file": protein_path,
                    "ligand_file": lig_path,
                }

            yield {
                "pdb_id": pdb_id,
                "protein": protein_structure,
                "ligands": ligands,
                "paths": {"protein": protein_path, "ligands": ligand_paths},
            }

    def _save(self, data: Any) -> None:
        raise NotImplementedError(
            "ProteinLigandPairsDataset is read-only. "
            "Writing protein/ligand mmCIF files is not supported."
        )

    def _describe(self) -> dict[str, Any]:
        return {
            "root_path": str(self._root),
            "protein_prefix": self._protein_prefix,
            "ligand_prefix": self._ligand_prefix,
            "cif_ext": self._cif_ext,
            "allow_missing_ligands": self._allow_missing_ligands,
            "n_folders": self._count_folders(),
        }

    def _count_folders(self) -> int:
        if not self._root.exists():
            return 0
        return sum(1 for p in self._root.iterdir() if p.is_dir())

    def _open_binary(self, path: Path) -> BinaryIO:
        if path.suffix == ".gz":
            with gzip.open(path, "rb") as f:
                data = cast("bytes", f.read())
                return io.BytesIO(data)
        data = cast("bytes", path.read_bytes())
        return io.BytesIO(data)

    def _parse(self, path: Path, structure_id: str) -> Structure:
        # BioPython's parser accepts a file path or a file handle with `readline`.
        # For gzipped inputs we parse from an in-memory binary handle.
        if path.suffix == ".gz":
            handle = self._open_binary(path)
            wrapped = cast(io.BufferedIOBase, handle)
            text_handle: TextIO = io.TextIOWrapper(wrapped, encoding="utf-8")
            return self._parser.get_structure(structure_id, text_handle)
        return self._parser.get_structure(structure_id, str(path))

    def _ligand_centroid_from_complex(
        self, protein: Structure, lig_id: str
    ) -> tuple[float, float, float] | None:
        """Compute ligand centroid from its bound pose in the protein complex.

        It collects coordinates of all atoms belonging to residues with `resname == lig_id`
        and a non-empty hetflag.
        """
        coords: list[tuple[float, float, float]] = []
        for res in protein.get_residues():
            if not res.id[0].strip():
                continue
            if res.resname != lig_id:
                continue
            for atom in res.get_atoms():
                x, y, z = atom.coord
                coords.append((float(x), float(y), float(z)))

        if not coords:
            return None
        return self._centroid(coords)

    def _extract_ligand_coords_from_protein(
        self, protein: Structure, lig_id: str
    ) -> list[tuple[float, float, float]] | None:
        """Deprecated: use `_ligand_centroid_from_complex` instead."""
        return None

    def _extract_ligand_coords_from_ligand_cif(
        self, ligand_path: Path
    ) -> list[tuple[float, float, float]] | None:
        """Extract ligand coordinates directly from the ligand CIF.

        Supports:
          * coordinate mmCIF with `_atom_site.Cartn_*`
          * chemical component CIF with `_chem_comp_atom.model_Cartn_*`

        Returns:
            List of (x, y, z) tuples, or None if coordinates cannot be found.
        """
        text = self._read_text(ligand_path)

        # Chemical component CIFs can be either loop_ form (multiple atoms) or
        # single-record form (e.g. ions like MG) where fields are listed once.
        if "_chem_comp_atom.model_Cartn_x" in text:
            coords = self._extract_coords_from_chem_comp_atom(text)
            if coords is not None:
                return coords
            coords = self._extract_coords_from_single_chem_comp_atom(text)
            if coords is not None:
                return coords

        # Try parsing as mmCIF structure and reading atom coords
        try:
            s = self._parser.get_structure(ligand_path.stem, io.StringIO(text))
        except Exception:
            return None

        coords: list[tuple[float, float, float]] = []
        for atom in s.get_atoms():
            x, y, z = atom.coord
            coords.append((float(x), float(y), float(z)))
        return coords or None

    def _extract_coords_from_chem_comp_atom(self, cif_text: str) -> list[tuple[float, float, float]] | None:
        """Parse a `_chem_comp_atom` loop and return model_Cartn coordinates."""
        lines = cif_text.splitlines()

        chem_headers: list[str] = []
        header_end: int = -1
        for i, line in enumerate(lines):
            if line.strip() != "loop_":
                continue
            j = i + 1
            while j < len(lines) and lines[j].strip().startswith("_"):
                j += 1
            headers = [h.strip() for h in lines[i + 1 : j] if h.strip().startswith("_")]
            if any(h.startswith("_chem_comp_atom.") for h in headers):
                chem_headers = headers
                header_end = j
                break

        if header_end == -1:
            return None

        def idx(name: str) -> int | None:
            full = f"_chem_comp_atom.{name}"
            try:
                return chem_headers.index(full)
            except ValueError:
                return None

        x_i = idx("model_Cartn_x")
        y_i = idx("model_Cartn_y")
        z_i = idx("model_Cartn_z")
        if x_i is None or y_i is None or z_i is None:
            return None

        import shlex

        rows: list[str] = []
        k = header_end
        while k < len(lines):
            s = lines[k].strip()
            if not s or s.startswith("#"):
                if s.startswith("#"):
                    break
                k += 1
                continue
            if s == "loop_" or s.startswith("_"):
                break
            rows.append(lines[k])
            k += 1

        coords: list[tuple[float, float, float]] = []
        for row in rows:
            toks = shlex.split(row, posix=True)
            if len(toks) <= max(x_i, y_i, z_i):
                continue
            try:
                coords.append((float(toks[x_i]), float(toks[y_i]), float(toks[z_i])))
            except ValueError:
                continue

        return coords or None

    def _extract_coords_from_single_chem_comp_atom(
        self, cif_text: str
    ) -> list[tuple[float, float, float]] | None:
        """Extract coordinates from non-loop chem_comp_atom blocks.

        Example (ions):
            _chem_comp_atom.model_Cartn_x 0.000
            _chem_comp_atom.model_Cartn_y 0.000
            _chem_comp_atom.model_Cartn_z 0.000

        Returns a single coordinate triple.
        """
        x = y = z = None
        for line in cif_text.splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith("_chem_comp_atom.model_Cartn_x"):
                x = s.split()[-1]
            elif s.startswith("_chem_comp_atom.model_Cartn_y"):
                y = s.split()[-1]
            elif s.startswith("_chem_comp_atom.model_Cartn_z"):
                z = s.split()[-1]

        if x is None or y is None or z is None:
            return None
        try:
            return [(float(x), float(y), float(z))]
        except ValueError:
            return None

    def _read_text(self, path: Path) -> str:
        if path.suffix == ".gz":
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                return f.read()
        return path.read_text(encoding="utf-8", errors="replace")

    def _find_single(self, folder: Path, stem: str, exts: tuple[str, ...]) -> Path | None:
        for ext in exts:
            candidate = folder / f"{stem}{ext}"
            if candidate.exists():
                return candidate
        return None

    def _find_many_with_prefix(
        self, folder: Path, prefix: str, exts: tuple[str, ...]
    ) -> dict[str, Path]:
        out: dict[str, Path] = {}
        for p in folder.iterdir():
            if not p.is_file():
                continue

            matches_ext = any(str(p.name).endswith(ext) for ext in exts)
            if not matches_ext:
                continue

            name = p.name
            if not name.startswith(prefix):
                continue

            lig_id = name[len(prefix) :]
            for ext in exts:
                if lig_id.endswith(ext):
                    lig_id = lig_id[: -len(ext)]
                    break

            out[lig_id] = p

        return dict(sorted(out.items(), key=lambda kv: kv[0]))

    def _centroid(
        self, coords: list[tuple[float, float, float]]
    ) -> tuple[float, float, float]:
        """Return the mean coordinate triplet for a set of (x,y,z) points."""
        if not coords:
            raise ValueError("Cannot compute centroid of empty coordinate list")
        sx = sy = sz = 0.0
        for x, y, z in coords:
            sx += x
            sy += y
            sz += z
        n = float(len(coords))
        return (sx / n, sy / n, sz / n)
