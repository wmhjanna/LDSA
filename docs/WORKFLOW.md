# Analysis workflow

## Demonstration

Run `scripts/run_demo.py` from the package root. This entry point uses only the included observations, deposition fractions and reference outputs. It calculates station-level LDSA and network coverage comparisons.

## Autocorrelation-aware inference

The `autocorrelation_bootstrap/` module starts from prepared monthly observations. Its README provides the full-record and exact-overlap moving-block analyses, joint paired bootstrap, REML synthesis and BH/BY adjustment commands. Use this module for the 60%-coverage, 143-pair inference; the general source scripts below include other processing and comparison settings. The supplied BCN monthly inputs reproduce the station examples, and the reference summaries reproduce the 47-site synthesis and 143-pair multiplicity adjustment. Refitting every station requires the corresponding monthly input records.

## Full-network source scripts

The `implementation/` directory contains source scripts and templates for the study analyses. They require the complete observation files and intermediate outputs associated with each processing stage. They do not form a single-command pipeline from the example dataset.

Set the workspace locations in `config/workspace.example.json`, then run:

```sh
python scripts/stage_sources.py --config config/workspace.example.json --destination staged_sources
```

This command copies source files and resolves configurable path tokens. It does not download data, execute analyses or overwrite an existing destination. `project/` contains processing scripts, `external/` contains input adapters and plotting templates, and `delivery/` contains figure assets and composition scripts. The configured runtime directories must match these locations. Check each stage's input files and schemas before execution.

## Processing stages and inputs

Paths in the following table are relative to `implementation/project/work/` unless stated otherwise.

| Operation | Source location | Required inputs |
|---|---|---|
| Provider quality flags and source matching | `quality_control/`, `provider_flag_audit/` | Native PNSD observations and provider metadata |
| Workbook PNSD integration | `quality_control/workbook/`, `build_site_tables.py` | Native grids, PNSD, per-grid deposition fractions and metadata |
| NetCDF integration and station registry | `network_trends/`, `station_registry/` | NetCDF observations and channel configuration tables |
| Study-specific record exclusions | `record_exclusions/`, `gio_record_filter/` | Station series and exclusion configurations |
| Pollutants and meteorology | `german_pollutants/`, `netcdf_pollutants/`, `meteorology/`, `network_analysis/` | Station-matched ancillary observations |
| Trends and meta-analysis | `network_analysis/`, `record_exclusions/workbook/trends.R` | Monthly series, phase definitions and eligibility tables |
| LightGBM and SHAP | `network_analysis/model/model.py` | Weekly observations, phase ranges and eligibility tables |
| Spatial figures | `spatial_statistics/`, `geographic_inset/` | Station statistics, coordinates and geographic resources |
| Station diagnostics | `record_exclusions/layout/`, `record_exclusions/render_s6.py`, `gio_record_filter/` | Trend estimates, breakpoints and time series |
| Pollutant and model figures | `record_exclusions/pollutants/`, `network_analysis/render.py` | Exact-overlap comparisons and model outputs |
| Integration sensitivity | `integration_sensitivity/` | Comparable native-grid integrations |
| Coverage sensitivity | `coverage_sensitivity/`, package-level `src/coverage.py` | Station summaries and native-grid coverage |
| Supplementary tables | `supplementary_tables/build_tables.py` | Analysis tables and Word document template |

The processing sequence is acquisition and source matching, native-grid integration and quality control, configured record exclusions, temporal and spatial aggregation, statistical analysis, then figure and table generation. The directories contain shared functions and templates as well as executable scripts; running every source file in filename order is not a supported execution sequence. Some figure scripts load selected functions or code sections from another source template.

## External resources

Full monitoring datasets, ancillary observations, original workbook/document templates, all intermediate data tables and the MPPD application are not included. The paths in `source_path_inventory.csv` identify configurable source references; individual file-reading calls define their required schemas. Supply these resources at the configured paths before running the corresponding stage. Access conditions and dataset citations follow the manuscript's Data availability statement and provider records.

The European map panel uses the supplied raster asset. The northeastern United States inset uses Natural Earth geometry. Figure appearance depends on the fonts and source assets available on the system. The numerical demonstration does not validate the geographic content of the supplied raster.

Full-network run time depends on data acquisition, input size and model fitting. See [VALIDATION.md](VALIDATION.md) for the tested components.

## Workbook and directory conventions

The station diagnostic adapter reads `17-500LDSAbreakpoint.xlsx` and per-station workbooks from `clean17-500/station_series/`. A station workbook contains `Time` and `LDSA_17_500nm` columns in `LDSA_part_1` and/or `LDSA_part_2` worksheets, or in a single `Total` worksheet. These names define the input schema expected by the distributed adapter.

Other commonly used external resources are `station_locations.xlsx`, `spatial_reference_tables/station_spatial_statistics.xlsx`, `EBAS_station_inventory/`, and `provider_quality_flags/`. Input data must match the names and schemas used in the corresponding source scripts. Data-provider filenames and identifiers embedded in source metadata are retained for provenance.
