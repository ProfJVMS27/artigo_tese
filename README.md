# Reproducibility package

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23387361.svg)](https://doi.org/10.5281/zenodo.23387361)

This package accompanies the manuscript *Uncertainty-aware surrogate-assisted multi-objective design of hybrid ballistic-radiological shields for RDD response*.

## Scope

The package contains the processed datasets used for model development, frozen split identifiers, trained surrogate artifacts, source code for post-processing, feature engineering, surrogate learning, constrained optimization, comparison studies, and high-fidelity re-evaluation summaries.  SHA-256 hashes are provided for every packaged file in `PROVENANCE/sha256_manifest.csv`.

The package deliberately excludes solver executables, raw LS-DYNA binary databases (`d3plot*`, `d3dump*`, `binout*`), MCNP runtime files (`runtpe*`), and other regenerated transient output.  These are not needed to inspect the reported tables or to rerun the Python analyses after installing the licensed solvers and reconstructing the original solver decks.  The included manifests preserve the associated simulation provenance.

## Directory map

- `DATA/FINAL_DATASETS`: final integrated MCNP6.3 DE/DF AP and LS-DYNA datasets and audit.
- `DATA/PHASE3`: fixed splits, calibration/test predictions, and cross-validation summaries.
- `MODELS`: trained operational surrogate bundles and metadata.
- `CODE`: scripts grouped by thesis phase.
- `RESULTS/PHASE4`: final NSGA-II output and 100-run optimization-comparison study.
- `RESULTS/PHASE5`: 30-candidate high-fidelity re-evaluation summaries and targeted numerical-sensitivity records.
- `MANUSCRIPT`: manuscript, supplementary information, bibliography, class/style files, and publication figures.
- `PROVENANCE`: package inventory, SHA-256 manifest, and release checklist.

## Environment

Create a Python 3.10 environment and install the pinned core dependencies in `requirements.txt`.  MCNP6.3 and LS-DYNA are proprietary licensed solvers and are not distributed with this archive.

## Reproduction order

1. Inspect the final datasets in `DATA/FINAL_DATASETS`.
2. Run the Phase 3 workflow in `CODE/PHASE3` to reproduce model-selection tables and model artifacts.
3. Run the Phase 4 workflow in `CODE/PHASE4` using the supplied operational models and datasets.
4. Inspect high-fidelity validation and numerical-sensitivity summaries in `RESULTS/PHASE5`.

## Important scope statement

The numerical-sensitivity materials document a targeted diagnostic, not a general mesh-convergence certification.  They must be interpreted with the qualification stated in the manuscript and supplementary information.

## Before a public release

Complete every item in `PROVENANCE/RELEASE_CHECKLIST.md`, particularly authorship, licensing, archival DOI, and any institutional review of material that may be publicly released.

## Como Citar

Se você utilizar este material, por favor cite:

> **MENDES DA SILVA, João Vitor**. *Uncertainty-aware surrogate-assisted multi-objective design of hybrid ballistic-radiological shields for RDD response*. Pacote de reprodutibilidade. Tese de Doutorado em Engenharia Nuclear, COPPE/UFRJ, Rio de Janeiro, 2026. DOI: [10.5281/zenodo.23387361](https://doi.org/10.5281/zenodo.23387361)
