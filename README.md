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


**Java requirement:** P2Rank needs Java 17+. On macOS:
```bash
brew install openjdk@17
```
Set `java_home` in `conf/base/parameters_extract_pockets.yml` to your JDK path.

**xgboost requirement:** the `pocket_ml` pipeline needs xgboost + OpenMP:
```bash
uv add xgboost
brew install libomp
```

## Setup
Clone the repository using: 
```bash
git clone --recursive https://github.com/im-stanly/druggability-medicine.git
# or with ssh
git clone --recursive git@github.com:im-stanly/druggability-medicine.git
```
Install dependencies into a project-local virtual
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
git clone --recursive 
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

## Project layout

```
src/druggability/
  pipelines/
    compare_pockets/    # compare and match pockets between PDB and AlphaFold
    extract_pockets/    # P2Rank-based pocket detection (requires P2Rank distro)
    pocket_ml/          # pocket machine-learning pipeline
    protein_align/      # chain matching + sequence alignment + RMSD
    protein_unzip/      # decompress raw .cif.gz files
  datasets/
    protein_dataset.py  # MmcifPairedDataset (paired PDB + AlphaFold mmCIF)
    path_dataset.py     # PathDataset (passes file paths through the catalog)
  pipeline_registry.py
conf/
  base/                 # catalog.yml, parameters_*.yml, logging.yml
  local/                # gitignored — local credentials / overrides
data/
  01_raw/                # raw inputs (protein_ligand_raw, proteins)
  02_intermediate/       # decompressed CIFs and intermediate artifacts
  03_primary/            # matched PDB ↔ AlphaFold pairs
  04_feature/            # feature extraction outputs
  05_model_input/
  06_models/
  07_model_output/
  08_reporting/          # final outputs and reports
  scrapped/              # scraped RCSB/AlphaFold hits and intermediate JSON
docs/
  pocket-finding/        # fpocket vs P2Rank comparison + examples
libs/
  p2rank/                # p2rank submodule (build produces distro/)
notebooks/               # exploratory notebooks
tests/                   # pytest suite and pipeline tests
```

Notes:
- P2Rank is included as a git submodule under `libs/p2rank` and must be built
  to produce the `distro` used by the `extract_pockets` pipeline.
- Configuration is driven from `conf/base` with local overrides in `conf/local`.
- Pipelines are registered in `src/druggability/pipeline_registry.py` and
  executed via Kedro (invoked through `uv` in this project).
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
uv run kedro run --pipeline=align
uv run kedro run --pipeline=extract_pockets
uv run kedro run --pipeline=pocket_ml
```

# Run P2Rank pocket detection on all training structures
uv run kedro run --pipeline=extract_pockets

# Train XGBoost druggability classifier on P2Rank-detected pockets
uv run kedro run --pipeline=pocket_ml
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

### pLDDT filtering (optional)

Before running P2Rank, you can filter out low-confidence AlphaFold regions
by setting `plddt_threshold` in `conf/base/parameters_extract_pockets.yml`:

```yaml
p2rank:
  plddt_threshold: 70.0   # atoms with pLDDT < 70 → occupancy=0 → hidden from P2Rank
```

Set to `0.0` (default) to disable. The filter only affects AlphaFold structures
(where B-factor = pLDDT); PDB structures pass through unchanged.

```bash
uv run kedro run --pipeline=extract_pockets
```

### Hypothesis verification scripts

Three standalone verification scripts (no Kedro needed):

```bash
# H1: pLDDT filtering reduces false-positive pockets
uv run python verify_plddt.py

# H2: Simple pocket descriptors (volume, hydrophobicity, polarity)
#     classify druggable pockets with AUROC > 0.75
uv run python verify_simple_classifier.py

# H3: AlphaFold pockets are ~40% smaller than PDB pockets
#     (no induced fit → narrower binding sites)
uv run python verify_af_volume.py
```

Results are saved to `data/08_reporting/`:
- `plddt_verification.json`
- `simple_classifier_results.json`
- `af_volume_comparison.json`

## Fetching new data (optional)

`scrappingData.py` queries RCSB for human X-ray structures with bound non-water
ligands and downloads the matching CIF files (plus a `rcsb_hits.json` manifest)
into `data/01_raw/proteins/<ENTRY_ID>/`.

> **Heads up:** the scraper imports `rcsbapi`, which is **not yet listed in
> `pyproject.toml`**. Install it once before running:
>
> ```bash
> uv add rcsb-api
> ```

Then:

```bash
uv run python scrappingData.py --help
uv run python scrappingData.py --output-dir data/01_raw/protein_ligand_raw
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
- `protein_compare` pipeline — chain matching by sequence identity, structural
  alignment, and per-pair RMSD reporting.
- `extract_pockets` pipeline — P2Rank pocket detection on PDB structures with
  optional pLDDT-based filtering of low-confidence AlphaFold regions.
- `pocket_ml` pipeline — druggability classifier (XGBoost) trained on pocket
  descriptors and ECFP fingerprints.
- Hypothesis verification: pLDDT filtering, simple pocket classifier,
  AF vs PDB pocket volume comparison.

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
- fpocket — https://github.com/Discngine/fpocket
- P2Rank — https://github.com/rdk/p2rank
- In-repo notes: `docs/pocket-finding/comparison.md`, `docs/pocket-finding/fpocket.md`,
  `docs/pocket-finding/p2rank.md`
