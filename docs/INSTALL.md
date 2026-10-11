# Installation and system requirements

## Python demonstration

Use Python 3.12 and install the exact dependency versions in `requirements-demo.txt` in a virtual environment, as shown in the README. A network connection is required for dependency installation. The gzip-compressed CSV inputs are read directly without manual decompression.

The demonstration has been tested with Python 3.12.14 on macOS 26.6.2 (ARM64, 16 GB system memory). It runs on a standard CPU without specialised hardware. Allow approximately 1 GB of available memory and 300 MB of disk space for the extracted package and outputs, excluding Python and its dependencies. The same entry point is provided for Windows and Linux, which have not been tested for this version.

Installation typically takes an estimated 5–15 minutes when binary packages are available; compilation or a slow connection may require longer. This estimate is not a clean-machine benchmark.

```sh
python scripts/run_demo.py --output demo_output
python -m unittest discover -s tests -v
python scripts/verify_release.py
```

The numerical demonstration takes a few seconds on the tested system, excluding installation. Font cache initialisation and slower storage may increase the first-run time. Figure output uses Times New Roman when available and DejaVu Serif otherwise. Fonts are not distributed with the package.

## Full-network analysis

`requirements-analysis.txt` lists Python dependencies for the additional analyses. R 4.5.3 and the package versions listed in `requirements-R.txt` are used for Theil–Sen trends and random-effects meta-analysis. R packages must be installed separately using the R package manager.

LightGBM requires a compatible OpenMP runtime. On macOS, install a matching ARM64 or x86_64 `libomp`, for example with `brew install libomp`. A missing `libomp.dylib` prevents LightGBM from importing. LightGBM is not required for the Python demonstration.

Data acquisition requires access to the corresponding provider services. The Multiple-Path Particle Dosimetry (MPPD) model is external software; this package uses exported deposition fractions and does not include the MPPD application. The full analyses also require the datasets and resources described in [WORKFLOW.md](WORKFLOW.md).

Use `scripts/stage_sources.py` to create a copy of the source files with paths resolved from the workspace configuration. This command does not execute the analysis.

## Optional R demonstration

```sh
Rscript src/bcn_trends.R examples/data/BCN_UB/trend_input.csv demo_trends
python scripts/check_bcn_trends.py demo_trends
```

Use a new output directory for each fit. The single-station calculation is expected to take seconds to a few minutes on a standard CPU. See [BCN_EXAMPLE.md](BCN_EXAMPLE.md) for the input definitions and [VALIDATION.md](VALIDATION.md) for the scope of testing.

## Autocorrelation-aware bootstrap module

Follow `autocorrelation_bootstrap/README.md` for its Python and R installation commands. The module requires openair 2.19.0 because the wrapper adjusts its internal bootstrap settings, meta 8.2-1, and a C++11 compiler for joint paired resampling. Run its commands from `autocorrelation_bootstrap/`.
