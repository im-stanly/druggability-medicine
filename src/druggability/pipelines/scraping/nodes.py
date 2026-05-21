"""Kedro node that runs the RCSB PDB scraper.

Wraps scrappingData.py (Stanly's scraper) so it runs as a Kedro node.
Output lands in data/scrapped/structures/ with one subdirectory per entry,
each containing PDB-*.cif.gz and optionally LIG-*.cif.gz files.
"""

import logging
import runpy
from pathlib import Path

logger = logging.getLogger(__name__)


def run_scraper(
    output_dir: str = "data/scrapped/structures",
    output_json: str = "data/scrapped/rcsb_hits.json",
    years: int = 4,
    resolution_min: float = 2.0,
    resolution_max: float = 3.0,
    rfree_min: float = 0.15,
    rfree_max: float = 0.25,
    limit: int = 500,
    polite_delay: float = 0.25,
    debug: bool = False,
) -> dict:
    """Run the RCSB scraper via runpy (triggers __main__ guard).

    Returns a summary dict with entry/file counts.
    """
    scraper_path = str(
        Path(__file__).resolve().parent.parent.parent.parent.parent
        / "scrappingData.py"
    )

    # Inject CLI args so scraper.main() sees them
    import sys
    old_argv = sys.argv
    sys.argv = [
        "scrappingData.py",
        "--years", str(years),
        "--resolution-min", str(resolution_min),
        "--resolution-max", str(resolution_max),
        "--rfree-min", str(rfree_min),
        "--rfree-max", str(rfree_max),
        "--limit", str(limit),
        "--output-dir", output_dir,
        "--output-json", output_json,
        "--polite-delay", str(polite_delay),
    ]
    if debug:
        sys.argv.append("--debug")

    try:
        logger.info("Scraping RCSB PDB (years=%d, limit=%d) → %s",
                    years, limit, output_dir)
        runpy.run_path(scraper_path, run_name="__main__")
    finally:
        sys.argv = old_argv

    # Count what we got
    out = Path(output_dir)
    n_entries = sum(1 for p in out.iterdir() if p.is_dir()) if out.exists() else 0
    n_pdb = sum(1 for p in out.rglob("PDB-*.cif.gz")) if out.exists() else 0

    logger.info("Scraper done: %d entries, %d PDB files", n_entries, n_pdb)
    return {"n_entries": n_entries, "n_pdb_files": n_pdb,
            "output_dir": output_dir, "output_json": output_json}
