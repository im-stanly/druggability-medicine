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
- **[`P2Rank`](https://github.com/rdk/p2rank)** — for pocket detection (see below).


## P2Rank setup
The project uses P2Rank for pocket detection. In order to setup the P2Rank download the binary from [P2Rank releases](https://github.com/rdk/p2rank/releases).

Put the downloaded binary in the root folder of the project. If done correctly the `extract_pockets` pipeline should run without any issues.


## Training data setup

The pipeline needs a protein-ligand pair dataset to train the druggability classifier. 
The current implementation expects protein-ligand pairs in the `data/01_raw/protein_ligand_raw/` directory,
where each protein has its own subfolder containing the corresponding PDB and ligand mmcif files. 

In order to download the data use a script `scrappingData.py` that scrapes the RCSB PDB for human X-ray structures with bound non-water ligands. The script will download the matching CIF files (plus a `rcsb_hits.json` manifest) into `data/01_raw/protein_ligand_raw/<ENTRY_ID>/`.

You can control the output using various command-line arguments, e.g. to specify a different output directory and select only data that is not older than a certain release date:

```bash
uv run python scrappingData.py --output-dir data/01_raw/protein_ligand_raw --years 8 # only include structures released in the last 8 years
```

To see other options, run:

```bash
uv run python scrappingData.py --help
```

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

Run default pipeline from `pipeline_registry.py` (parse downloaded data and train the model):

```bash
uv run kedro run
```

Run a specific pipeline:

```bash
# Decompress data/01_raw/proteins/**/*.cif.gz into data/02_intermediate/proteins/
uv run kedro run --pipeline=unzip

# Match PDB↔AlphaFold chains, align them, write data/08_reporting/protein_alignment_results.csv
uv run kedro run --pipeline=compare

# Train the model on downloaded protein-ligand pairs 
uv run kedro run --pipeline=pocket_ml
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
