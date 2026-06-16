"""
Nodes for the `extract_pockets` pipeline.

This first iteration focuses on:
- running P2Rank on each protein complex
- parsing its output (pockets CSV)

Later iterations can add:
- pocket-to-ligand distance matching
- druggable/non-druggable labels
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import csv
import os
import shutil
import subprocess
import tempfile
import platform as sys_platform
import pandas as pd


@dataclass(frozen=True)
class P2RankConfig:
    p2rank_dir: str
    threads: int = 1
    # Optional: point directly to the prank launcher; otherwise derived from p2rank_dir
    prank_bin: str | None = None
    # Optional: set JAVA_HOME for prank
    java_home: str | None = None
    # Optional: choose a built-in profile like "alphafold" OR a path to a groovy file
    config: str | None = None
    # Disable heavy visualization generation by default (faster)
    visualizations: int = 0
    # pLDDT threshold: set occupancy=0 for atoms with pLDDT < this (0 = disabled)
    plddt_threshold: float = 0.0


def run_p2rank_and_parse_pockets(
    protein_ligand_ds: Iterable[dict[str, Any]],
    p2rank: dict[str, Any],
) -> pd.DataFrame:
    """Run P2Rank for each complex and parse pocket predictions.

    Notes on I/O contract:
      * Input is the output of `ProteinLigandPairsDataset`, i.e. an iterable of
        dicts with keys like `pdb_id`, `paths`, `ligands`.
      * Output is a *flat*, tabular dataset: one row per predicted pocket.

    Returns:
        A pandas DataFrame with (at minimum):
          - pdb_id
          - pocket_id
          - center_x, center_y, center_z
          - score, probability
    """
    cfg = P2RankConfig(**p2rank)

    p2rank_dir = Path(cfg.p2rank_dir)
    
    prank = Path(cfg.prank_bin) if cfg.prank_bin else (p2rank_dir / "distro" / "prank")
    if sys_platform.system() == "Windows" and prank.suffix != ".bat":
        prank = prank.with_name(f"{prank.name}.bat")
    if not prank.exists():
        raise FileNotFoundError(f"P2Rank launcher not found at: {prank}")

    records: list[dict[str, Any]] = []

    def _iter_with_progress(it: Iterable[dict[str, Any]]):
        """Add a tqdm progress bar if tqdm is installed."""
        try:
            from tqdm.auto import tqdm  # type: ignore
        except Exception:
            yield from it
            return

        # Some iterables (like generators) don't support len(). That's fine.
        total = None
        try:
            total = len(it)  # type: ignore[arg-type]
        except Exception:
            total = None

        yield from tqdm(it, total=total, desc="P2Rank", unit="protein")

    for item in _iter_with_progress(protein_ligand_ds):
        pdb_id = item["pdb_id"]
        protein_path = Path(item["paths"]["protein"])

        with tempfile.TemporaryDirectory(prefix=f"p2rank_{pdb_id}_") as tmp_dir:
            tmp_dir_p = Path(tmp_dir)
            in_path = tmp_dir_p / protein_path.name

            # ── pLDDT filtering ──────────────────────────────────────
            if cfg.plddt_threshold > 0:
                from plddt_filter import filter_cif

                filtered_path = tmp_dir_p / f"filtered_{protein_path.name}"
                filter_cif(protein_path, filtered_path, threshold=cfg.plddt_threshold)
                shutil.copy2(filtered_path, in_path)
            else:
                shutil.copy2(protein_path, in_path)
            # ─────────────────────────────────────────────────────────

            out_dir = tmp_dir_p / "out"
            out_dir.mkdir(parents=True, exist_ok=True)

            cmd = [
                str(prank),
                "predict",
                "-f",
                str(in_path),
                "-o",
                str(out_dir),
                "-threads",
                str(cfg.threads),
                "-visualizations",
                str(cfg.visualizations),
            ]
            if cfg.config:
                cmd.extend(["-c", cfg.config])

            env = os.environ.copy()
            if cfg.java_home:
                env["JAVA_HOME"] = cfg.java_home

            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    "P2Rank failed.\n"
                    f"PDB: {pdb_id}\n"
                    f"Command: {' '.join(cmd)}\n"
                    f"stdout:\n{e.stdout}\n\n"
                    f"stderr:\n{e.stderr}\n"
                ) from e

            predictions_csv = _find_predictions_csv(out_dir=out_dir, input_name=in_path.name)

            ligand_centroids = {
                lig_id: lig["centroid"] for lig_id, lig in item.get("ligands", {}).items()
            }

            pocket_records = parse_p2rank_predictions_csv(
                predictions_csv,
                pdb_id=pdb_id,
                ligand_centroids=ligand_centroids,
            )
            records.extend(pocket_records)

    df = pd.DataFrame.from_records(records)

    # Keep output stable even when empty.
    expected_cols = [
        "pdb_id",
        "pocket_id",
        "rank",
        "center_x",
        "center_y",
        "center_z",
        "score",
        "probability",
    ]
    for c in expected_cols:
        if c not in df.columns:
            df[c] = pd.Series(dtype="float64" if c not in {"pdb_id", "pocket_id"} else "object")

    return df[expected_cols + [c for c in df.columns if c not in expected_cols]]


def parse_p2rank_predictions_csv(
    predictions_csv: Path,
    *,
    pdb_id: str,
    ligand_centroids: dict[str, tuple[float, float, float]] | None = None,
) -> list[dict[str, Any]]:
    """Parse a P2Rank `<input>_predictions.csv` into flat records.

    This is separated for unit testing (no need to execute P2Rank).
    """
    ligand_centroids = ligand_centroids or {}

    with predictions_csv.open("r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    def norm_key(k: str) -> str:
        # P2Rank sometimes pads headers (e.g. "   center_x").
        return k.strip()

    def norm_row(row: dict[str, Any]) -> dict[str, Any]:
        return {norm_key(k): v for k, v in row.items()}

    records: list[dict[str, Any]] = []
    for i, row0 in enumerate(rows, start=1):
        row = norm_row(row0)

        pocket_id = (
            row.get("name")
            or row.get("pocket")
            or row.get("id")
            or row.get("pocket_id")
            or f"pocket_{i}"
        )

        rec: dict[str, Any] = {
            "pdb_id": pdb_id,
            "pocket_id": str(pocket_id).strip(),
            "rank": _to_int(row.get("rank") or row.get("order") or row.get("#")),
            "center_x": _to_float(row.get("center_x") or row.get("centerX") or row.get("center.x")),
            "center_y": _to_float(row.get("center_y") or row.get("centerY") or row.get("center.y")),
            "center_z": _to_float(row.get("center_z") or row.get("centerZ") or row.get("center.z")),
            "score": _to_float(row.get("score")),
            "probability": _to_float(row.get("probability") or row.get("prob")),
            # Keeping these around for the next step (distances & labels)
            "ligand_centroids": ligand_centroids,
            "raw": row0,
        }
        records.append(rec)

    return records


def _find_predictions_csv(*, out_dir: Path, input_name: str) -> Path:
    """Find the predictions CSV produced by P2Rank.

    P2Rank typically writes `<input_file_name>_predictions.csv` but may sometimes
    use the stem of the file.
    """
    candidates = [
        out_dir / f"{input_name}_predictions.csv",
        out_dir / f"{Path(input_name).stem}_predictions.csv",
    ]
    for c in candidates:
        if c.exists():
            return c

    # Fallback: pick the only *_predictions.csv in out_dir if unambiguous
    matches = sorted(out_dir.glob("*_predictions.csv"))
    if len(matches) == 1:
        return matches[0]

    raise FileNotFoundError(
        f"Expected P2Rank predictions CSV not found for input '{input_name}'. "
        f"Tried: {[str(c) for c in candidates]}. "
        f"out_dir contents: {[p.name for p in out_dir.iterdir()]}"
    )


def _to_float(v: str | None) -> float | None:
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _to_int(v: str | None) -> int | None:
    if v is None:
        return None
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None
