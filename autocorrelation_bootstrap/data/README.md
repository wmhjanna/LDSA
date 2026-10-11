# Example inputs: BCN_UB

These are the prepared monthly observations used in the current 65-site LDSA
study. They have been extracted from the analysis inputs, without new aggregation
or imputation. BCN_UB is the Barcelona urban-background station.

- `BCN_LDSA_MetaFull.csv`: 83 monthly observations, January 2013 to February 2020,
  for the eligible full-record LDSA trend used in the 47-site synthesis.
- `BCN_cross_metric_monthly.csv`: 16 series, eight LDSA-pollutant exact-overlap
  pairs, represented in long format (1,158 rows). Each pair contains 69-74 common
  monthly observations. LDSA is repeated as a separate matched series for each
  pollutant because the common month sets can differ.
- `BCN_paired_monthly.csv`: the same eight pairs in wide format (579 rows),
  ready for joint resampling. Pollutants are BC, CO, NO, NO2, O3, PM10, PM2.5,
  and SO2. All eight BCN pairs pass the >=60-month and >=60%-coverage rule.

## Long-format input for the R trend script

| Column | Meaning |
|---|---|
| date | Original monthly timestamp in ISO format |
| value | Prepared monthly concentration |
| Site | Station identifier |
| Station_class | UB, TR, SUB, or RB; the alias Type is also accepted |
| Metric | LDSA or pollutant name |
| Phase | MetaFull, or Overlap_ followed by the pollutant identifier |
| group_id | Unique series identifier: Site\|Metric\|Phase |

## Wide-format input for the joint Python bootstrap

| Column | Meaning |
|---|---|
| Site | Station identifier |
| Station_class | UB, TR, SUB, or RB |
| Pollutant | Pollutant identifier |
| Date | Common monthly timestamp in ISO format |
| LDSA | Monthly LDSA concentration on that common month |
| Pollutant_value | Monthly pollutant concentration on that same month |

Both values must be finite, with one row per pair per month. The inputs retain
the concentration units of the prepared analysis series; LDSA is in um^2 cm^-3.
These scripts do not convert pollutant units. A constant positive unit conversion
within an entire series cancels in the percentage-trend normalization. Absolute
`Slope_Unit` outputs retain the input concentration unit per year and must be
interpreted with the parent dataset's unit metadata.

The full-record LDSA input and exact-overlap LDSA inputs have different inclusion
sets. Do not substitute one for the other. Neither dataset includes inferred or
filled values for missing months. The package starts after the study's source
quality control and monthly preparation; it does not repeat those upstream steps.
