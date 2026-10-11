#!/usr/bin/env python3
"""Compare run_examples.py outputs with the current study reference results."""
import argparse
import json
from pathlib import Path
import platform
import numpy as np
import pandas as pd

def read(path):
    return pd.read_csv(path, float_precision="round_trip")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("results"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    ref = root / "reference_results"
    out = args.results.resolve()
    checks = []
    def compare(label, observed, expected, keys, columns, atol=1e-9):
        a = observed.set_index(keys).sort_index()
        b = expected.set_index(keys).sort_index()
        assert a.index.equals(b.index), f"{label}: keys differ"
        x = a[columns].to_numpy(float); y = b[columns].to_numpy(float)
        np.testing.assert_allclose(x, y, rtol=1e-10, atol=atol, equal_nan=True, err_msg=label)
        checks.append(dict(check=label, rows=len(a), columns=columns,
            max_absolute_difference=float(np.nanmax(np.abs(x-y))), passed=True))
    trends = pd.concat([read(out / f"BCN_{s}_1999.csv") for s in ["Auto", "iid", "3", "6", "12"]])
    trend_cols = ["n_months", "nboot", "block_length", "Slope_Unit", "CI_Lower_Unit", "CI_Upper_Unit",
                  "Slope_Percent", "CI_Lower_Percent", "CI_Upper_Percent", "P_value"]
    compare("BCN full-record bootstrap, all settings", trends,
            read(ref / "BCN_trends_1999_all_settings.csv"), ["Site", "Method"], trend_cols)
    cross = read(out / "BCN_cross_599.csv")
    compare("BCN exact-overlap trend bootstrap", cross, read(ref / "BCN_cross_trends_599.csv"), ["group_id"], trend_cols)
    pair_ref = read(ref / "paired_trends_143.csv")
    bcn = pair_ref[pair_ref.Site.eq("BCN_UB")]
    by_id = cross.set_index("group_id")
    for p in bcn.itertuples():
        for metric, expected in [("LDSA", p.LDSA_class_block), (p.Pollutant, p.Pollutant_class_block)]:
            assert by_id.loc[f"{p.Site}|{metric}|Overlap_{p.Pollutant}", "Trend_class"] == expected
    checks.append(dict(check="BCN trend classifications", series=16, passed=True))
    compare("REML estimates", read(out / "REML/REML_summary.csv"), read(ref / "REML_47.csv"),
            ["Station_class"], ["k", "effect", "lower", "upper", "p", "tau2", "I2", "Q"])
    compare("REML site SEs and overall weights", read(out / "REML/REML_site_weights.csv"), read(ref / "site_trends_47.csv"),
            ["Site"], ["SE", "weight_percent"])
    pairs = read(out / "paired/paired_difference_results.csv")
    numeric_cols = ["N_months", "Coverage", "Block_length", "Replicates", "LDSA_slope", "Pollutant_slope",
        "Delta", "CI_lower", "CI_upper", "Delta_SE", "Approx_p", "Bootstrap_bias", "LDSA_baseline_p025",
        "Pollutant_baseline_p025", "Nonpositive_baseline_draws", "Finite_draws", "Bootstrap_slope_correlation", "Unpaired_SE_diagnostic"]
    expected_pairs = read(ref / "BCN_paired_differences.csv")
    compare("BCN joint bootstrap, all settings", pairs, expected_pairs, ["Site", "Pollutant", "Method"], numeric_cols)
    for name in ["Inference_available", "CI_excludes_zero", "Status"]:
        a = pairs.set_index(["Pollutant", "Method"])[name].sort_index()
        b = expected_pairs.set_index(["Pollutant", "Method"])[name].sort_index()
        assert a.equals(b), f"Paired diagnostic differs: {name}"
    checks.append(dict(check="BCN inference availability and interval flags", rows=len(pairs), passed=True))
    adjusted = read(out / "multiplicity_143.csv")
    compare("BH and BY over all 143 pairs", adjusted.rename(columns={"BH_q_recomputed":"BH", "BY_q_recomputed":"BY"}),
        read(ref / "paired_differences_143_all_settings.csv").rename(columns={"BH_q_143":"BH", "BY_q_143":"BY"}),
        ["Site", "Pollutant", "Method"], ["BH", "BY"], atol=1e-12)
    primary = adjusted[adjusted.Method.eq("Auto")]
    assert len(primary) == 143
    supported = primary.BH_q_recomputed.lt(.05) & primary.Inference_available
    main = primary.Pollutant.isin(["BC", "PM2.5", "PM10", "NO2"])
    assert int(supported.sum()) == 21 and int((supported & main).sum()) == 8 and int(main.sum()) == 65
    checks.append(dict(check="Current 60%-coverage multiplicity family", pairs=143, main_metric_pairs=65,
                       BH_supported_all=21, BH_supported_main_metrics=8, passed=True))
    # Direct, unweighted estimator on explicitly expanded records, independently
    # of the C++ multiplicity implementation.
    audit = json.loads((out / "paired/independent_reference_samples.json").read_text())
    errors = []
    for record in audit:
        t = np.asarray(record["t"], float); y = np.asarray(record["y"], float)
        i, j = np.triu_indices(len(t), 1)
        keep = t[j] != t[i]
        slope = np.median((y[j[keep]]-y[i[keep]])/(t[j[keep]]-t[i[keep]]))
        start = np.median(y) - slope*np.median(t)
        calculated = np.array([slope, start, 36500*slope/start])
        expected = np.array([record[k] for k in ["Slope_per_day", "Start", "Percent"]])
        np.testing.assert_allclose(calculated, expected, atol=1e-9, rtol=1e-10)
        errors.append(float(np.max(np.abs(calculated-expected))))
    assert len(audit) == 48
    checks.append(dict(check="Independent direct Theil-Sen estimator on resampled records", fits=len(audit),
                       max_absolute_difference=max(errors), passed=True))
    report = dict(status="PASS", python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, checks=checks)
    (out / "validation_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
