# LDSA analysis of a 65-station observational network

Version 1.0.0 | Supplementary Software

This package provides the analysis code for a study of lung-deposited surface area (LDSA) at 65 monitoring stations: 64 in Europe and ROC_UB in the United States. It includes analysis scripts, configuration files, a processed particle number size distribution (PNSD) record from BCN_UB, and reference outputs for numerical verification.

## Quick start

Use Python 3.12 and run these commands from the extracted package directory:

```sh
python -m venv .venv
# macOS or Linux
. .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements-demo.txt
python scripts/run_demo.py
python -m unittest discover -s tests -v
python scripts/verify_release.py
```

The demonstration calculates LDSA on the native diameter grid for the head-airway (HA), tracheobronchial (TB) and alveolar (ALV) regions, together with total LDSA. It exports monthly means and a figure, compares the results with stored study outputs, and reproduces particle-size-coverage comparisons from station summaries. After dependency installation, it requires no additional data downloads, MPPD or R installation, GPU, or external input files.

## Package contents

| Directory or file | Description |
|---|---|
| `src/` | Integration, temporal preprocessing, coverage comparisons and R trend analysis |
| `scripts/` | Demonstration, numerical verification and source-path configuration |
| `examples/data/BCN_UB/` | 62,412 processed observations, 103 native channels, deposition fractions and metadata |
| `examples/data/site_statistics.csv` | Summary data for the 65-station coverage comparisons |
| `examples/expected/` | Reference numerical outputs and an example figure |
| `implementation/` | Source scripts for acquisition, quality control, network analyses, models, figures and tables |
| `config/` | Study parameters, record exclusions and workspace configuration |
| `docs/` | Installation, methods, data requirements, workflow and validation |
| `SHA256SUMS` | File integrity checksums |

## Running the analyses

The BCN_UB demonstration runs with the included data. The full-network scripts additionally require source observations, ancillary measurements, document templates and intermediate analysis tables. Input locations are configurable; see [WORKFLOW.md](docs/WORKFLOW.md). The single-station demonstration does not reproduce the full network on its own.

The 53- and 49-station sets are particle-size-coverage sensitivity subsets of the same 65-station network. Source adapters retain native diameter grids, units, deposition fractions and applicable provider quality flags.

Start with [INSTALL.md](docs/INSTALL.md), [BCN_EXAMPLE.md](docs/BCN_EXAMPLE.md) and [METHODS.md](docs/METHODS.md). Project-authored code is distributed under the MIT License. Third-party materials retain their applicable terms; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
