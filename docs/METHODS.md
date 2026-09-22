# Computational methods and study conventions

## Terminology

LDSA denotes lung-deposited surface area; PNSD, particle number size distribution; and MPPD, Multiple-Path Particle Dosimetry. HA, TB and ALV denote the head-airway, tracheobronchial and alveolar regions. Total LDSA is reported separately. UB, TR, SUB and RB denote urban background, traffic, suburban background and regional background stations, respectively.

## Study population

The study consists of 65 stations, including 64 European stations and ROC_UB in the United States. Source-format adapters handle workbook and NetCDF inputs within this single observational network. Final station classes are UB, TR, SUB and RB; U maps to UB, and F/M/P map to RB. Configuration changes are handled on their native grids. PAR_SUB and PRA_SUB each contribute one station record.

## Number and surface-area integration

For a native channel i with diameter D_i in nm, bin number concentration Delta N_i in cm⁻³ and MPPD deposition fraction DF_(r,i),

```
LDSA_r = sum_i [Delta N_i * pi * (D_i / 1000)^2 * DF_(r,i)]
```

The output unit is μm² cm⁻³. If a source stores dN/dlog10(Dp), convert it to Delta N using its native log10 bin width first. The workbook conversion uses forward differences of log10 channel centres, repeating the final spacing for the last bin. Native NetCDF configuration weights are retained by the relevant adapter. A source already expressed as Delta N is not weighted a second time.

All stations use MPPD-derived deposition fractions. Preserve the diameter-to-fraction mapping supplied for each configuration. The included model profile describes the exported curve calculation; settings for other tabulated curves are configuration-specific. The BCN example uses the tabulated coefficients without refitting.

The target interval is 17–500 nm, with the documented LND_TR interval of 31.65–500 nm. Sum native channel centres within the inclusive interval without extrapolation. Missing channels are omitted, all-missing rows remain missing, and zero is retained. Partial sums are not renormalized to represent unmeasured bins. This can underestimate full-window LDSA when some bins are unavailable.

The early PAL_RB configuration is treated as Delta N under the study's adopted interpretation. This source-unit interpretation must remain explicit when describing uncertainty.

## Quality control and exceptions

Apply missing/sentinel and source-quality checks at the measurement or channel level. Nonfinite and negative values are invalid. Provider flags classified as invalid are excluded only after source/configuration correspondence is established; valid warnings are not automatically rejected. The provider's valid override flag is respected. A channel-level flag does not by itself invalidate every other channel at the same timestamp.

`config/record_exclusions.json` lists explicit weekly intervals removed at GIO_RB, KOS_RB and PAL_RB, separately from provider flagging. `config/study.json` additionally lists the four excluded GIO_RB hours in July 2011. These exclusions are not an automatic statistical outlier threshold.

HEL_TR SO2 is excluded. PAR_SUB and PRA_SUB are retained in the LDSA network but have separate eligibility for conventional-pollutant comparisons. IPR_RB temperature data are excluded where incompatible source units cannot be resolved. Missing ancillary variables are not imputed from incompatible sources.

## Spatial and coverage analyses

Spatial statistics use available 2015–2019 measurements; if none exist in that interval, use the nearest available year, choosing the later year in a tie. The 53-station upper-bound subset requires a measured upper bound of at least 500 nm. The 49-station full-coverage subset requires coverage of the full 17–500 nm interval. These are sensitivity subsets, not separate versions of the study.

The comparisons use station-level means, Kruskal–Wallis statistics, pairwise two-sided Mann–Whitney tests and Holm multiplicity correction. See `src/coverage.py`. Additional integration sensitivity is restricted to configurations that actually support the required intervals.

## Temporal analysis

The source workflow selects a single linear-kernel change point on the ordered observation series, subject to the recorded eligibility rules. The 20-month condition is an overall data-availability criterion, not a guarantee of 20 months on either side. Direct monthly aggregation for full and phase trends is distinct from the monthly averaging of equally weighted daily means with a minimum data-completeness criterion for exact-overlap comparisons. Meta-analysis uses the separately prepared eligible monthly series and usable uncertainty estimates.

The R engine applies openair TheilSen and inverse-variance random-effects synthesis with REML and classical random-effects confidence intervals. Percentage-change confidence intervals with a lower bound greater than the upper bound are excluded from the synthesis. Inspect the eligibility tables and R session information alongside reported estimates. See `src/trends_meta.R` and the analysis scripts.

## Predictive models

LightGBM uses weekly observations, class/phase grouping, random seed 42, and a random 60/20/20 train/validation/test split. The model uses a learning rate of 0.05, num_leaves = 31, max_depth = 6, feature_fraction = 0.8, and up to 1,000 boosting rounds, with early stopping after 50 rounds without improvement. Eligibility requires at least 52 weeks in each phase. SHAP contributions are evaluated on held-out rows. The reported R² values are based on the random held-out test set; temporal and leave-one-site-out validation were not performed. The final class groups used for this analysis are UB, TR and RB.
