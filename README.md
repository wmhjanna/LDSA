# LDSA analysis of a 65-station observational network

Version 1.0.0 | Supplementary Software

[Repository](https://github.com/wmhjanna/LDSA) · [Download v1.0.0](https://github.com/wmhjanna/LDSA/releases/tag/v1.0.0)

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
| `autocorrelation_bootstrap/` | Moving-block trend uncertainty, paired trend differences, REML synthesis, multiplicity adjustment and BCN monthly examples |
| `implementation/` | Source scripts for acquisition, quality control, network analyses, models, figures and tables |
| `config/` | Study parameters, record exclusions and workspace configuration |
| `docs/` | Installation, methods, data requirements, workflow and validation |
| `SHA256SUMS` | File integrity checksums |

## Autocorrelation-aware trend inference

Use [autocorrelation_bootstrap/README.md](autocorrelation_bootstrap/README.md) for the moving-block bootstrap analyses of the 47 eligible full-record LDSA trends and 143 eligible exact-overlap site–pollutant pairs. This module includes prepared BCN_UB monthly observations, reference summaries and a command that reproduces and verifies its examples. It requires R with openair 2.19.0 and meta, Python, and a C++11 compiler. Cross-metric eligibility requires at least 60 common valid months and 60% calendar-month coverage.

The R commands in the integration demonstration and the general source scripts use their own trend settings; use the dedicated bootstrap module to reproduce the autocorrelation-aware uncertainty estimates and paired comparisons.

## Running the analyses

The BCN_UB demonstration runs with the included data. The full-network scripts additionally require source observations, ancillary measurements, document templates and intermediate analysis tables. Input locations are configurable; see [WORKFLOW.md](docs/WORKFLOW.md). The single-station demonstration does not reproduce the full network on its own.

The 53- and 49-station sets are particle-size-coverage sensitivity subsets of the same 65-station network. Source adapters retain native diameter grids, units, deposition fractions and applicable provider quality flags.

Start with [INSTALL.md](docs/INSTALL.md), [BCN_EXAMPLE.md](docs/BCN_EXAMPLE.md) and [METHODS.md](docs/METHODS.md). Project-authored code is distributed under the MIT License. Third-party materials retain their applicable terms; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
