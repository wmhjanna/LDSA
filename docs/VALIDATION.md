# Numerical validation

The included checks cover:

- Integration of 62,412 BCN_UB observations on 103 native diameter channels, producing HA, TB, ALV and total LDSA.
- Agreement of all timestamps, numerical values and missing-value locations with stored study outputs.
- BCN_UB spatial statistics, breakpoint and paired-month availability.
- Coverage comparisons for the 65-, 53- and 49-station sets, including Holm-adjusted pairwise tests.
- Optional R trends, confidence limits and unavailable estimates.
- Analytic tests for surface-area units, channel selection, zero and missing values, input representation and negative-value rejection.

The Python and R demonstrations pass the supplied reference comparisons. The maximum absolute integration difference is approximately 1.14 × 10⁻¹³ μm² cm⁻³. Numerical comparisons use explicit floating-point tolerances; pixel-identical figures are not expected across fonts, plotting-library versions and operating systems. The After phase in the R demonstration has unavailable slope estimates, which are retained and checked as missing values.

The autocorrelation-aware module also passes numerical checks for five BCN full-record bootstrap settings, 16 exact-overlap series, 40 paired bootstrap fits, five REML estimates and 47 site weights, 715 pair-setting BH/BY adjustments, and 48 independent direct Theil–Sen fits on resampled observations. The maximum paired-estimator difference from reference values is approximately 7.11 × 10⁻¹⁵. Run `autocorrelation_bootstrap/scripts/run_examples.py` to reproduce these checks; its output directory contains the generated validation report and timing information.

These tests validate the bundled example and summary-data comparisons. They do not establish end-to-end validation of the full-network processing, external data acquisition, machine-learning retraining or all manuscript graphics. Those analyses require the complete inputs described in [WORKFLOW.md](WORKFLOW.md).

`demo_validation.json` records the tested Python environment and numerical checks. The R demonstration writes `R_session_info.txt` alongside its results.
