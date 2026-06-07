"""Filter AlphaFold mmCIF by pLDDT — set occupancy=0 for low-confidence residues.

In AlphaFold mmCIF, pLDDT is stored in _atom_site.B_iso_or_equiv.
We set _atom_site.occupancy to 0.0 for atoms whose pLDDT < threshold.
P2Rank (BioJava) ignores zero-occupancy atoms, effectively hiding them
from pocket detection while keeping the structure intact.

Usage:
    from plddt_filter import filter_cif
    filter_cif("input.cif.gz", "output.cif", threshold=70)
"""

from __future__ import annotations

import gzip
import shlex
from pathlib import Path


def filter_cif(
    input_path: str | Path,
    output_path: str | Path,
    threshold: float = 70.0,
) -> dict:
    """Filter an AlphaFold mmCIF, setting occupancy=0 for pLDDT < threshold.

    Returns stats dict with n_atoms_before, n_atoms_after (kept), n_zeroed.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    if input_path.suffix == ".gz":
        with gzip.open(input_path, "rt", encoding="utf-8") as f:
            content = f.read()
    else:
        content = input_path.read_text(encoding="utf-8")

    lines = content.splitlines(keepends=True)

    # ── find _atom_site loop ─────────────────────────────────────────
    header_start: int | None = None
    data_start: int | None = None
    bfactor_col: int | None = None
    occupancy_col: int | None = None
    headers: list[str] = []

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped == "loop_":
            j = i + 1
            cand_headers = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                cand_headers.append(lines[j].strip())
                j += 1
            if any(h.startswith("_atom_site.") for h in cand_headers):
                header_start = i
                data_start = j
                headers = cand_headers
                for ci, h in enumerate(headers):
                    if h == "_atom_site.B_iso_or_equiv":
                        bfactor_col = ci
                    if h == "_atom_site.occupancy":
                        occupancy_col = ci
                break

    if bfactor_col is None or occupancy_col is None or data_start is None:
        raise ValueError(f"Required _atom_site columns not found in {input_path}")

    # ── parse atom data range ────────────────────────────────────────
    atom_end = data_start
    while atom_end < len(lines):
        s = lines[atom_end].strip()
        if not s or s.startswith("#"):
            atom_end += 1
            continue
        if s == "loop_" or s.startswith("_"):
            break
        atom_end += 1

    # ── rewrite with occupancy modifications ─────────────────────────
    n_zeroed = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        # Write header (unchanged)
        for i in range(data_start):
            f.write(lines[i])

        # Write atom data with modified occupancy
        for i in range(data_start, atom_end):
            line = lines[i].strip()
            if not line or line.startswith("#"):
                f.write(lines[i])
                continue

            try:
                tokens = shlex.split(line, posix=True)
            except ValueError:
                f.write(lines[i])
                continue

            if len(tokens) <= max(bfactor_col, occupancy_col):
                f.write(lines[i])
                continue

            try:
                bf = float(tokens[bfactor_col])
            except (ValueError, IndexError):
                f.write(lines[i])
                continue

            if bf < threshold:
                tokens[occupancy_col] = "0.00"
                n_zeroed += 1

            # Reconstruct line preserving whitespace as best we can
            f.write(" ".join(tokens) + "\n")

        # Write remaining lines
        for i in range(atom_end, len(lines)):
            f.write(lines[i])

    n_atoms = atom_end - data_start
    return {
        "n_atoms_total": n_atoms,
        "n_atoms_kept": n_atoms - n_zeroed,
        "n_zeroed": n_zeroed,
    }
