# Autocorrelation-aware bootstrap for the 65-site LDSA study

Part of LDSA Supplementary Software v1.0.0. Analysis scripts with BCN_UB
example data and reference results from the 65-site study.

This package contains the uncertainty calculations used for the 47 eligible
full-record LDSA trends and the exact-overlap cross-metric comparisons. It starts
from **prepared monthly observations**. It does not perform raw PNSD processing,
MPPD deposition calculations, breakpoint testing, or figure production.

## Analysis settings

| Analysis | Script | Bootstrap replicates | Main setting |
|---|---|---:|---|
| Individual full-record LDSA trends | `01_theilsen_block_bootstrap.R` | 1,999 | openair moving blocks, length round(n^(1/3)) |
| Exact-overlap trend classification | `01_theilsen_block_bootstrap.R` | 599 | Same moving-block setting |
| Overall and station-class synthesis | `02_reml_synthesis.R` | Uses the preceding site CIs | REML; conventional normal-approximation CI |
| Paired LDSA-minus-pollutant slope differences | `03_paired_trend_difference.py` | 9,999 | Joint moving blocks, length round(n^(1/3)) |
| Multiple-comparison adjustment | `04_adjust_multiplicity.py` | Uses paired bootstrap results | BH and BY over all 143 eligible pairs, per block setting |

Here n is the number of retained monthly observations. All trends use the original
dates and are expressed as % yr^-1 relative to the fitted starting concentration.
The cross-metric inclusion rule is **at least 60 common valid months and at least
60% coverage** of calendar months from the first to the last common observation,
inclusive. The 143 pairs involve eight pollutant metrics; they are not 143 sites.
The four main metrics contain 12 PM2.5, 19 PM10, 7 BC, and 27 NO2 pairs.

## Installation

Validated with Python 3.12.14, NumPy 2.3.5, pandas 2.2.3, R 4.5.3,
openair 2.19.0, and meta 8.2-1 on macOS. A C++11 compiler is needed for the paired
bootstrap engine. No GPU or proprietary software is required. Linux and Windows
build commands are supplied but have not been tested in this release.

From the repository root, enter this module directory with
`cd autocorrelation_bootstrap`, then run:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
Rscript -e 'remotes::install_version("openair", version="2.19.0", repos="https://cloud.r-project.org", upgrade="never")'
Rscript -e 'remotes::install_version("meta", version="8.2-1", repos="https://cloud.r-project.org", upgrade="never")'
.venv/bin/python scripts/build_engine.py
```

On Windows, use `.venv/Scripts/python.exe` and a compatible C++ compiler such as
MinGW; pass its command with `scripts/build_engine.py --compiler g++` if needed.
The R installation commands require internet access. Runtime analysis uses local
files only. R package dependencies are installed by `remotes`.

## Reproduce the examples

```sh
# BCN full-record LDSA uncertainty: 83 eligible monthly observations.
Rscript scripts/01_theilsen_block_bootstrap.R data/BCN_LDSA_MetaFull.csv results/BCN_Auto_1999.csv 1999 auto

# BCN exact-overlap trends: eight pollutants and their matched LDSA series.
Rscript scripts/01_theilsen_block_bootstrap.R data/BCN_cross_metric_monthly.csv results/BCN_cross_599.csv 599 auto

# Synthesis from the supplied 47-site CI table.
Rscript scripts/02_reml_synthesis.R reference_results/site_trends_47.csv results/REML

# Joint paired bootstrap; add --all-settings for the block-length sensitivities.
.venv/bin/python scripts/03_paired_trend_difference.py data/BCN_paired_monthly.csv results/paired --all-settings

# Reproduce the full 143-pair multiplicity adjustment from reference summaries.
.venv/bin/python scripts/04_adjust_multiplicity.py reference_results/paired_differences_143_all_settings.csv results/multiplicity_143.csv --expected-family-size 143
```

The overall REML result should be -2.906241944 % yr^-1, with a 95% CI of
[-3.493120875, -2.319363014]. BCN paired point estimates and intervals are in
`reference_results/BCN_paired_differences.csv`. The example q values calculated
from BCN alone use **eight pairs**, not the paper's 143-pair family. Use the full
reference table to reproduce the paper's q values.

For the R block-length sensitivity analyses, replace `auto` by `3`, `6`, or `12`.
Use `iid` for the independent-resampling comparator. Keep 1,999 replicates for all
these full-record comparisons. The synthesis script accepts a final `HK` argument
for the Hartung-Knapp sensitivity analysis; `classic` is the main result.

## What is resampled and how the intervals are calculated

### Individual trends and classification

`openair::TheilSen` is called with `avg.time="month"`, `deseason=FALSE`,
`slope.percent=TRUE`, and `autocor=TRUE`. The input is already monthly; the wrapper
checks for one finite observation per month. Openair resamples contiguous blocks
of retained (date, concentration) records with replacement, without circular
wrapping, until n observations are obtained. Dates remain attached to their
observations. The fit is repeated for each bootstrap sample. This is a pairs
bootstrap, not a residual bootstrap. Openair's internal seed is 2 in this version.

The default block length is round(n^(1/3)), **not 20 months**. The number 20 in the
sampler's function signature is an unused default when the main regression
function supplies its selected block length.

The wrapper returns `slope.percent`, `lower.percent`, and `upper.percent` directly
from openair. In version 2.19.0 these percentage bounds are constructed from the
package's coefficient limits. They must not be described as the direct 2.5th and
97.5th percentiles of replicate-level percentage slopes. The latter construction
is used for paired differences below, not for these individual percentage CIs.

Openair does not expose `nboot` as a public `TheilSen` parameter in this version.
The wrapper temporarily changes the default of its internal `regci` function and,
for fixed-block sensitivity runs, replaces the internal sampler. It restores both
functions after every fit. It stops on other openair versions to avoid silently
changing this behavior. No openair source code is redistributed here.

Classification is `Declining` when a valid CI is entirely below zero, `Increasing`
when it is entirely above zero, and `Uncertain` otherwise. Nonfinite or reversed
intervals, zero-width intervals, and intervals that do not contain the estimate
are marked `CI unavailable`. This descriptive classification is not adjusted for
multiple comparisons. Agreement includes two uncertain classifications and is
not a test that two slopes are equal.

### Random-effects synthesis

For each site, SE = (upper CI - lower CI) / (2 x 1.959963985).
This is the **CI-width approximation used in the analysis**, not the empirical SD
of its bootstrap slopes. `meta::metagen` estimates tau-squared by REML, with a
conventional normal-approximation CI in the primary analysis. Each station class
is fitted separately. The supplied site weights refer to the overall 47-site
fit. This bootstrap step does not model dependence between different sites.

### Paired differences

For every exact-overlap site-pollutant pair, the same block indices jointly
resample LDSA and the pollutant. For each indicator and each replicate, the
Theil-Sen slope is the median of all pairwise slopes at distinct dates and the
intercept is median(y) - slope x median(t), with t in days since the first common
observation. The annual percentage slope is 36500 x slope / intercept.

For replicate b, delta[b] = beta_LDSA[b] - beta_pollutant[b]. Its CI is the 2.5th
and 97.5th percentiles of the joint bootstrap differences and its SE is their
sample SD. The subtraction therefore includes the covariance between the two
slopes. A positive delta means the LDSA relative trend is more positive (or less
negative). Delta has units of percentage points yr^-1.

Percentage inference is flagged unavailable if either original fitted baseline
is nonpositive, either bootstrap baseline's 2.5th percentile is nonpositive, or
any joint percentage draw is nonfinite. Numeric intervals are retained for audit
but should not be interpreted as reliable percentage inference in those cases.

The approximate two-sided tail probability is
min(1, 2 x min((1 + count(delta <= 0))/(B + 1),
(1 + count(delta >= 0))/(B + 1))). It is a percentile-bootstrap tail probability,
not an exact permutation p value or a null-centered bootstrap test. BH and BY are
applied across all 143 eligible site-pollutant pairs for each block setting;
unavailable pairs remain in the family with p = 1. The primary Auto results have
21 BH-supported pairs over all eight pollutants, including eight among the four
main metrics.

The C++ engine computes the same Theil-Sen estimates from observation
multiplicities, avoiding explicit expansion of repeated observations. Random
seeds are derived from site, pollutant, and block-setting identifiers with the
SHA-256 rule in the script. A repeated block length within a pair reuses its
indices. The default B is 9,999; reducing B is useful for a smoke run but will not
reproduce the reference intervals.

### Missing months

No missing values are imputed. The default and fixed-size moving blocks are
defined in **retained observation order** and may bridge calendar gaps. The
`Calendar6` paired sensitivity uses six-calendar-month windows instead; nonempty
windows are resampled while retaining their original dates. Coverage, missing
months, and gap transitions are reported for every pair. The term
"autocorrelation-aware" describes these dependence-sensitive resampling choices;
it does not establish that every series' dependence is fully captured.

## Data and reproducibility

`data/README.md` describes the BCN example and all input columns.
`reference_results/` contains the current study's 47-site CI/REML results and
60%-coverage paired comparison summaries. These allow checking the released code
against the current analysis. Only BCN monthly observations are bundled; refitting
all other sites requires the corresponding prepared monthly input data in the
same schema. Aggregate summaries cannot replace those inputs.

The numerical verification checks cover all five BCN 1,999-replicate settings, 16 exact-overlap
599-replicate series, 40 paired 9,999-replicate fits, all five REML estimates and
47 weights, all 715 pair-setting BH/BY adjustments, and independent direct
Theil-Sen calculations on example resamples. Run them again with:

```sh
.venv/bin/python scripts/run_examples.py --output results
```

This convenience command builds the engine, runs all examples and fixed-block
checks, and verifies the outputs. To check results that already exist, run
`.venv/bin/python scripts/verify_results.py --results results`.

`scripts/run_examples.py` writes `validation_report.json` and
`execution_times.json` to the selected output directory. The package-level
`SHA256SUMS` covers the distributed files. The scripts use configurable input
paths. The software is available in the [LDSA repository](https://github.com/wmhjanna/LDSA)
and the [v1.0.0 release](https://github.com/wmhjanna/LDSA/releases/tag/v1.0.0).

## License and software references

Study scripts are distributed under the MIT License; see [LICENSE](../LICENSE). Third-party
dependencies retain their own licenses. Example observations and reference tables
are supplied to reproduce this study's calculations; the code license does not
alter the underlying data providers' access or reuse conditions.

- [openair TheilSen documentation](https://openair-project.github.io/openair/reference/TheilSen.html)
- [openair source repository](https://github.com/davidcarslaw/openair)
- [meta package](https://cran.r-project.org/package=meta)
- [NumPy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
