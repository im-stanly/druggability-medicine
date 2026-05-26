# Druggability-medicine

[![Powered by Kedro](https://img.shields.io/badge/powered_by-kedro-ffc900?logo=kedro)](https://kedro.org)

## Overview

Project investigating whether AlphaFold-predicted protein structures
are accurate enough to assess **druggability** of binding pockets. For a curated set
of human proteins we collect both experimental ([RCSB PDB](https://www.rcsb.org))
and predicted ([AlphaFold DB](https://alphafold.ebi.ac.uk)) structures, detect
binding pockets, compare them, and train a classifier that predicts whether a pocket
is druggable from geometric and physicochemical descriptors. The project is built
on top of [Kedro](https://kedro.org) for pipeline orchestration.

## Prerequisites

- **Python 3.13** — see `.python-version`. `pyproject.toml` accepts `>=3.10`, but
  the lockfile is resolved against 3.13.
- **[`uv`](https://docs.astral.sh/uv/)** — the package manager used throughout.
  Install with `curl -LsSf https://astral.sh/uv/install.sh | sh` (macOS / Linux)
  or follow the [official installer instructions](https://docs.astral.sh/uv/getting-started/installation/).

The current pipelines are pure Python — no JVM, Docker, or system libraries
required. Pocket-detection tools (fpocket, P2Rank) will be needed once that step
is integrated; see [Roadmap](#roadmap).

## Setup

Clone the repository and install dependencies into a project-local virtual
environment:

```bash
make init
# equivalent to:
#   uv init && uv sync --all-groups
```

For an existing checkout, re-sync after pulling changes:

```bash
make sync          # uv sync --all-groups
```

Verify the install:

```bash
uv run kedro info
```

You should see the pipelines `unzip` and `compare` listed.

## Project layout

```
src/druggability/
  pipelines/
    protein_unzip/      # decompress raw .cif.gz files
    protein_compare/    # chain matching + sequence alignment + RMSD
    protein_parsing/    # placeholder for downstream parsing
  datasets/
    protein_dataset.py  # MmcifPairedDataset (paired PDB + AlphaFold mmCIF)
    path_dataset.py     # PathDataset (passes file paths through the catalog)
  pipeline_registry.py
conf/
  base/                 # catalog.yml, parameters_*.yml, logging.yml
  local/                # gitignored — local credentials / overrides
data/
  01_raw/proteins/      # committed sample data (PDB-*.cif.gz, LIG-*.cif.gz, AF-*.cif.gz)
  02_intermediate/      # decompressed CIFs (output of `unzip` pipeline)
  03_primary/           # matched PDB↔AlphaFold pairs
  08_reporting/         # final outputs (e.g. protein_alignment_results.csv)
notebooks/              # exploratory work (e.g. mmcif_parsing.ipynb)
docs/pocket-finding/    # fpocket vs P2Rank comparison + sample outputs
tests/                  # pytest suite
scrappingData.py        # standalone RCSB scraper (see below)
```

## Running the pipelines

Run everything registered in `pipeline_registry.py`:

```bash
uv run kedro run
```

Run a specific pipeline:

```bash
# Decompress data/01_raw/proteins/**/*.cif.gz into data/02_intermediate/proteins/
uv run kedro run --pipeline=unzip

# Match PDB↔AlphaFold chains, align them, write data/08_reporting/protein_alignment_results.csv
uv run kedro run --pipeline=compare
```

Alignment thresholds (chain identity cutoff, gap penalties) live in
`conf/base/parameters_protein_compare.yml`. Catalog entries — input paths,
intermediate locations, and the final CSV — are defined in `conf/base/catalog.yml`.

Visualize the DAG:

```bash
uv run kedro viz
```

For interactive exploration (`context`, `catalog`, `pipelines` are pre-loaded):

```bash
uv run kedro jupyter notebook
uv run kedro ipython
```

## Fetching new data (optional)

`scrappingData.py` queries RCSB for human X-ray structures with bound non-water
ligands and downloads the matching CIF files (plus a `rcsb_hits.json` manifest)
into `data/01_raw/proteins/<ENTRY_ID>/`.

> **Heads up:** the scraper imports `rcsbapi`, which is **not yet listed in
> `pyproject.toml`**. Install it once before running:
>
> ```bash
> uv add rcsbapi
> ```

Then:

```bash
uv run python scrappingData.py --help
uv run python scrappingData.py --output-dir data/01_raw/proteins
```

AlphaFold structures (`AF-<UniProt>-*.cif.gz`) are currently fetched manually
from [AlphaFold DB](https://alphafold.ebi.ac.uk) and dropped alongside the PDB
file in the same protein folder. Automating this is on the roadmap.

## Testing & formatting

```bash
make test      # uv run pytest
make format    # uv run ruff format
```

Coverage settings live under `[tool.coverage.report]` in `pyproject.toml`.

## Roadmap

**Implemented**

- RCSB scraper for human structures with bound ligands (`scrappingData.py`).
- `unzip` pipeline — decompress paired PDB / AlphaFold mmCIF inputs.
- `compare` pipeline — chain matching by sequence identity, structural
  alignment, and per-pair RMSD reporting.

**Planned**

- Pocket detection on both PDB and AlphaFold structures using **fpocket** and/or
  **P2Rank** — see `docs/pocket-finding/comparison.md`, `docs/pocket-finding/fpocket.md`,
  and `docs/pocket-finding/p2rank.md` for a tool-by-tool comparison and sample
  outputs.
- Pocket descriptor extraction (volume, hydrophobicity, polarity, compactness).
- Druggability classifier (logistic regression / random forest) trained on
  PDBbind / DrugBank ligand–protein pairs.
- pLDDT-based filtering of low-confidence AlphaFold regions before pocket
  detection, and analysis of how prediction confidence affects druggability
  agreement with experimental structures.

## References

- Kedro documentation — https://docs.kedro.org
- AlphaFold Protein Structure Database — https://alphafold.ebi.ac.uk
- RCSB Protein Data Bank — https://www.rcsb.org
- In-repo notes: `docs/pocket-finding/comparison.md`, `docs/pocket-finding/fpocket.md`,
  `docs/pocket-finding/p2rank.md`
