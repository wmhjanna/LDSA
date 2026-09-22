# BCN_UB example dataset

The demonstration uses the complete available processed BCN_UB PNSD record: 62,412 timestamped rows from 1 January 2013 to 14 February 2020, with 103 native channels. Of these rows, 48,663 have a finite primary total LDSA value. The example contains observed study data rather than simulated data.

`pnsd_bin_concentrations.csv.gz` stores cleaned **number concentration per size bin, Delta N, in cm⁻³**. Do not multiply these values by log-diameter widths again. `Time` retains the precision of the source timestamps. `Source_part` identifies the native configuration. Blank cells are missing.

Native channel centres span 12.2–478.3 nm. Selection by the 17–500 nm target retains centres from 17.5 to 478.3 nm. There is no extrapolation to 500 nm. The demo is not evidence of complete 17–500 nm measurement coverage.

`deposition_fractions.csv` contains the corresponding MPPD-derived HA, TB, ALV and Total fractions. The Total column is the exact tabulated coefficient used for total LDSA. Rounded regional coefficients may not sum exactly to the separately tabulated Total column; the demonstration preserves that input convention.

The demo writes:

- `BCN_UB_ldsa.csv.gz`: LDSA values for the HA, TB and ALV regions, together with total LDSA, for every input timestamp.
- `BCN_UB_monthly.csv`: direct arithmetic means by calendar month, without a minimum monthly data-completeness criterion. This output does not implement the separate monthly completeness rule used for paired comparisons or meta-analysis.
- `BCN_UB_spatial_statistics.csv`: 2015–2019 site statistics.
- `BCN_UB_demo.png` and `.pdf`: monthly total LDSA.
- coverage comparison tables and `validation.json`.

All timestamps and all four LDSA output columns are compared with stored study results, including the locations of missing values. The maximum observed absolute discrepancy was approximately 1.14e-13 μm² cm⁻³, attributable to floating-point arithmetic and CSV serialization.

The numerical reference outputs come from the study's stored calculations, not from the demonstration run that is being checked. Input and output SHA-256 digests are included in the release manifest.

The MIT licence applies to project code. Observation and third-party data rights remain with their providers; inclusion as supplementary demonstration data does not relicense the entire source monitoring dataset.

## Optional trend calculation in R

With R and openair installed, run from the package root:

```sh
Rscript src/bcn_trends.R examples/data/BCN_UB/trend_input.csv demo_trends
python scripts/check_bcn_trends.py demo_trends
```

Use a new output directory for a fresh bootstrap fit. This demonstration runs the study Theil–Sen routine on the actual prepared BCN_UB monthly trend inputs. Those inputs retain the distinct monthly aggregation/phase conventions of each analysis and are not interchangeable with the direct monthly CSV plotted by the Python demo. The fitted Before, Full and MetaFull results agree with the stored study outputs. The After phase has only five monthly input values and yields unavailable slope/interval estimates under the fitted routine; these missing values are included in the numerical verification. R session information is saved with the output.

The Python demo also reproduces the BCN_UB observation-level breakpoint (5 October 2019 at 03:30) and 75 months meeting the daily-data completeness criterion for paired comparisons. This is distinct from the 83 eligible MetaFull input months.
