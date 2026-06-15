# Druggability-medicine

[![Powered by Kedro](https://img.shields.io/badge/powered_by-kedro-ffc900?logo=kedro)](https://kedro.org)

## Overview

Project investigating whether AlphaFold-predicted protein structures
are accurate enough to assess **druggability** of binding pockets. We collect
experimental ([RCSB PDB](https://www.rcsb.org)) and predicted ([AlphaFold DB](https://alphafold.ebi.ac.uk)) structures, detect
binding pockets, compare them, and train a classifier that predicts whether a pocket
is druggable from geometric and physicochemical descriptors. Built on [Kedro](https://kedro.org).

## Prerequisites

- **Python 3.13** — see `.python-version`.
- **[`uv`](https://docs.astral.sh/uv/)** — install per project docs.

## Setup

Clone the repository and install dependencies into a project-local virtual
environment:

```bash
make init
# equivalent to: uv init && uv sync --all-groups
```

For existing checkout, re-sync after pulling changes:

```bash
make sync
```

Verify the install:

```bash
uv run kedro info
```

## P2Rank setup
P2Rank is included as a git submodule under `libs/p2rank`. Do not add a built distro into the main repo.

Clone with submodules:

```bash
git clone --recurse-submodules <repo-url>
# or for an existing clone:
git submodule update --init --recursive
```

Build P2Rank (produces `distro` with binaries):

```bash
cd libs/p2rank
./make-distro.sh    # Linux / macOS
# or on Windows:
./gradlew.bat build
```

Paths used by project configuration:

- `p2rank_dir`: `./libs/p2rank/distro`
- `prank_bin` (Linux/macOS): `./libs/p2rank/distro/prank`
- `prank_bin` (Windows): `./libs/p2rank/distro/prank.bat`

To change these, edit `conf/base/parameters_extract_pockets.yml` and set `p2rank_dir` and `prank_bin` accordingly.

After building and configuring paths, run the extract_pockets pipeline.

## How to install dependencies

You should see the pipelines `unzip` and `compare` listed.

## Project layout

```
src/druggability/
  pipelines/
    protein_unzip/      # decompress raw .cif.gz files
    protein_compare/    # chain matching + sequence alignment + RMSD
    protein_parsing/    # downstream parsing utilities
    extract_pockets/    # P2Rank-based pocket detection (requires P2Rank distro)
  datasets/
    protein_dataset.py  # MmcifPairedDataset (paired PDB + AlphaFold mmCIF)
    path_dataset.py     # PathDataset (passes file paths through the catalog)
  pipeline_registry.py
conf/
  base/                 # catalog.yml, parameters_*.yml, logging.yml
  local/                # gitignored — local credentials / overrides
data/
  01_raw/proteins/      # committed sample data
  02_intermediate/      # decompressed CIFs
  03_primary/           # matched PDB↔AlphaFold pairs
  08_reporting/         # final outputs
notebooks/
  # exploratory work
docs/pocket-finding/    # fpocket vs P2Rank comparison + sample outputs
tests/                  # pytest suite
```

## Running the pipelines

Run everything registered in `pipeline_registry.py`:

```bash
uv run kedro run
```

Run a specific pipeline:

```bash
uv run kedro run --pipeline=unzip
uv run kedro run --pipeline=compare
uv run kedro run --pipeline=extract_pockets
```

## Testing & formatting

```bash
make test      # uv run pytest
make format    # uv run ruff format
```

## Roadmap

See `docs/pocket-finding/` for pocket-finding specifics.
