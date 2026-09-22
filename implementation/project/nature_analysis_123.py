#!/usr/bin/env python3
"""LDSA trends, breakpoint offsets and synchronized pollutant trajectories.

1. Exact-overlap site-level trend comparison (LDSA vs co-pollutants).
2. Metric-specific KernelCPD breakpoint offsets relative to LDSA.
3. Synchronized London and Stockholm monthly trajectories.

The script creates clean CSV inputs/results and publication-style PNG/PDF figures.
The site-level Theil-Sen trends and their confidence intervals are calculated by
the companion R script ``nature_analysis_trends.R`` using openair::TheilSen.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import ruptures as rpt
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr
from docx import Document


ROOT = Path("@PROJECT_ROOT@")
OUT = ROOT / "outputs" / "nature_analysis_1_2_3"
OUT.mkdir(parents=True, exist_ok=True)

LDSA_DIR = Path("@INPUT_ROOT@/data_sources/plotting/breakpoint_plots/clean17-500")
POLLUTANT_BOOK = Path(
    "@INPUT_ROOT@/data_sources/predictive_models/pollutant_data/all_sites/"
    "paired_date_value_measurements.xlsx"
)
RAW_POLLUTANT_DIR = Path("@INPUT_ROOT@/data_sources/predictive_models/pollutant_data/all_sites")
SP_DOCX = Path(
    "@ATTACHMENT_ROOT@/"
    "xwechat_files/wxid_u2wh4eh9sunh22_ea96/temp/drag/SP(1).docx"
)

EXTERNAL_METRICS = ["PM2.5", "PM10", "BC", "NO2", "SO2", "O3", "CO", "NO"]
FOCUS_METRICS = ["PM2.5", "PM10", "BC", "NO2", "SO2"]
MIN_VALID_DAYS_MONTH = 14
MIN_PRIMARY_MONTHS = 60
MIN_PRIMARY_COVERAGE = 0.60
RANDOM_SEED = 20260816

TYPE_COLORS = {"UB": "#2C7FB8", "SUB": "#6A51A3", "TR": "#E6550D", "RB": "#31A354"}
METRIC_COLORS = {
    "LDSA": "#111111",
    "PM2.5": "#D73027",
    "PM10": "#FC8D59",
    "BC": "#542788",
    "NO2": "#4575B4",
    "SO2": "#91BFDB",
    "O3": "#1A9850",
    "CO": "#FDAE61",
    "NO": "#74ADD1",
    "Ntotal": "#4D4D4D",
    "Nucleation": "#7B3294",
    "Aitken": "#008837",
    "Accumulation": "#C51B7D",
}


def set_plot_style() -> None:
    candidates = ["Times New Roman", "Times", "Liberation Serif", "DejaVu Serif"]
    available = {f.name for f in mpl.font_manager.fontManager.ttflist}
    chosen = next((x for x in candidates if x in available), "DejaVu Serif")
    mpl.rcParams.update(
        {
            "font.family": chosen,
            "font.size": 11,
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "axes.linewidth": 1.1,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def site_metadata() -> pd.DataFrame:
    doc = Document(SP_DOCX)
    rows = [[c.text.strip() for c in r.cells] for r in doc.tables[0].rows]
    meta = pd.DataFrame(rows[1:], columns=rows[0])
    meta = meta.rename(
        columns={
            "Site code": "site",
            "Country": "country",
            "Station type": "site_type",
            "Measurement period": "measurement_period",
        }
    )
    return meta[["site", "country", "site_type", "measurement_period"]].drop_duplicates("site")


def month_start(x: pd.Series) -> pd.Series:
    return pd.to_datetime(x, errors="coerce").dt.to_period("M").dt.to_timestamp()


def monthly_from_daily(date: pd.Series, value: pd.Series, min_days: int = MIN_VALID_DAYS_MONTH) -> pd.DataFrame:
    d = pd.DataFrame({"date": pd.to_datetime(date, errors="coerce"), "value": pd.to_numeric(value, errors="coerce")})
    d = d.dropna().sort_values("date")
    if d.empty:
        return pd.DataFrame(columns=["date", "value", "valid_days"])
    # Equal weighting of calendar days even if the source is hourly.
    daily = d.assign(day=d["date"].dt.floor("D")).groupby("day", as_index=False)["value"].mean()
    daily["date"] = daily["day"].dt.to_period("M").dt.to_timestamp()
    out = daily.groupby("date", as_index=False).agg(value=("value", "mean"), valid_days=("value", "count"))
    out.loc[out["valid_days"] < min_days, "value"] = np.nan
    return out.dropna(subset=["value"])


def read_ldsa_monthly(site: str) -> pd.DataFrame:
    path = LDSA_DIR / f"{site}_Integrated_LDSA.xlsx"
    cols = ["Time", "LDSA_17_500nm", "Nucleation_17_25nm", "Aitken_25_100nm", "Accumulation_100_500nm"]
    df = pd.read_excel(path, sheet_name="Total", usecols=cols)
    out = None
    rename = {
        "LDSA_17_500nm": "LDSA",
        "Nucleation_17_25nm": "LDSA_Nucleation",
        "Aitken_25_100nm": "LDSA_Aitken",
        "Accumulation_100_500nm": "LDSA_Accumulation",
    }
    for col, new in rename.items():
        m = monthly_from_daily(df["Time"], df[col]).rename(columns={"value": new, "valid_days": f"{new}_valid_days"})
        out = m if out is None else out.merge(m, on="date", how="outer")
    return out.sort_values("date")


def read_pollutant_monthly(site: str) -> dict[str, pd.DataFrame]:
    df = pd.read_excel(POLLUTANT_BOOK, sheet_name=site)
    result: dict[str, pd.DataFrame] = {}
    for i in range(0, len(df.columns) - 1, 2):
        date_col, metric = df.columns[i], str(df.columns[i + 1]).strip()
        if metric not in EXTERNAL_METRICS:
            continue
        result[metric] = monthly_from_daily(df[date_col], df.iloc[:, i + 1])
    return result


def months_between(a: pd.Timestamp, b: pd.Timestamp) -> int:
    return (b.year - a.year) * 12 + (b.month - a.month)


def max_month_gap(dates: pd.Series) -> int:
    ds = sorted(pd.to_datetime(dates).dropna().tolist())
    if len(ds) < 2:
        return 0
    return max(months_between(a, b) for a, b in zip(ds[:-1], ds[1:]))


def prepare_monthly_pairs(meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    sheets = pd.ExcelFile(POLLUTANT_BOOK).sheet_names
    meta_idx = meta.set_index("site")
    pair_rows: list[pd.DataFrame] = []
    coverage_rows: list[dict] = []
    ldsa_all: dict[str, pd.DataFrame] = {}
    for site in sheets:
        ldsa = read_ldsa_monthly(site)
        ldsa_all[site] = ldsa
        pollutants = read_pollutant_monthly(site)
        for metric, pol in pollutants.items():
            pair = ldsa[["date", "LDSA"]].merge(pol[["date", "value"]], on="date", how="inner").dropna()
            if pair.empty:
                continue
            pair = pair.sort_values("date")
            start, end = pair["date"].min(), pair["date"].max()
            n = len(pair)
            span = months_between(start, end) + 1
            gap = max_month_gap(pair["date"])
            pair_id = f"{site}__{metric.replace('.', '')}"
            info = {
                "pair_id": pair_id,
                "site": site,
                "country": meta_idx.loc[site, "country"] if site in meta_idx.index else "",
                "site_type": meta_idx.loc[site, "site_type"] if site in meta_idx.index else site.rsplit("_", 1)[-1],
                "pollutant": metric,
                "overlap_start": start,
                "overlap_end": end,
                "n_months": n,
                "span_months": span,
                "coverage_fraction": n / span,
                "max_gap_months": gap,
                "eligible_5y": n >= MIN_PRIMARY_MONTHS and span >= MIN_PRIMARY_MONTHS,
                "primary_included": n >= MIN_PRIMARY_MONTHS
                and span >= MIN_PRIMARY_MONTHS
                and n / span >= MIN_PRIMARY_COVERAGE
                and site != "PAR_SUB",
                "qc_note": "PAR_SUB documented discontinuity"
                if site == "PAR_SUB"
                else ("long gap >12 months" if gap > 12 else ("intermittent monthly gaps" if gap > 3 else "")),
            }
            coverage_rows.append(info)
            common = pair.rename(columns={"LDSA": "ldsa_value", "value": "pollutant_value"})
            for role, col in [("LDSA", "ldsa_value"), ("Pollutant", "pollutant_value")]:
                tmp = common[["date", col]].rename(columns={col: "value"})
                for k, v in info.items():
                    tmp[k] = v
                tmp["role"] = role
                pair_rows.append(tmp)
    pairs = pd.concat(pair_rows, ignore_index=True)
    coverage = pd.DataFrame(coverage_rows).sort_values(["pollutant", "site"])
    pairs.to_csv(OUT / "01_exact_overlap_monthly_pairs.csv", index=False)
    coverage.to_csv(OUT / "01_pair_coverage_qc.csv", index=False)
    return pairs, coverage, ldsa_all


def run_r_trends() -> pd.DataFrame:
    rscript = ROOT / "nature_analysis_trends.R"
    proc = subprocess.run(
        ["Rscript", str(rscript), str(OUT / "01_exact_overlap_monthly_pairs.csv"), str(OUT / "01_exact_overlap_trends.csv")],
        check=True,
        text=True,
        capture_output=True,
    )
    (OUT / "01_R_trend_log.txt").write_text(proc.stdout + "\n" + proc.stderr, encoding="utf-8")
    return pd.read_csv(OUT / "01_exact_overlap_trends.csv", parse_dates=["overlap_start", "overlap_end"])


def classify_trend(lower: float, upper: float) -> str:
    if pd.isna(lower) or pd.isna(upper):
        return "Unavailable"
    lo, hi = min(lower, upper), max(lower, upper)
    if hi < 0:
        return "Declining"
    if lo > 0:
        return "Increasing"
    return "Stable/uncertain"


def summarize_trend_pairs(trends: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    trends["qc_note"] = trends["qc_note"].fillna("")
    wide = trends.pivot_table(
        index=["pair_id", "site", "country", "site_type", "pollutant", "overlap_start", "overlap_end", "n_months", "span_months", "coverage_fraction", "max_gap_months", "eligible_5y", "primary_included", "qc_note"],
        columns="role",
        values=["slope_percent", "lower_percent", "upper_percent", "p_value", "slope_absolute", "lower_absolute", "upper_absolute"],
        aggfunc="first",
    ).reset_index()
    wide.columns = ["_".join([str(x) for x in c if str(x) != ""]).strip("_") if isinstance(c, tuple) else c for c in wide.columns]
    wide = wide.rename(
        columns={
            "slope_percent_LDSA": "ldsa_slope_pct_yr",
            "lower_percent_LDSA": "ldsa_ci_low_pct_yr",
            "upper_percent_LDSA": "ldsa_ci_high_pct_yr",
            "p_value_LDSA": "ldsa_p",
            "slope_absolute_LDSA": "ldsa_slope_abs_yr",
            "lower_absolute_LDSA": "ldsa_ci_low_abs_yr",
            "upper_absolute_LDSA": "ldsa_ci_high_abs_yr",
            "slope_percent_Pollutant": "pollutant_slope_pct_yr",
            "lower_percent_Pollutant": "pollutant_ci_low_pct_yr",
            "upper_percent_Pollutant": "pollutant_ci_high_pct_yr",
            "p_value_Pollutant": "pollutant_p",
            "slope_absolute_Pollutant": "pollutant_slope_abs_yr",
            "lower_absolute_Pollutant": "pollutant_ci_low_abs_yr",
            "upper_absolute_Pollutant": "pollutant_ci_high_abs_yr",
        }
    )
    wide["ldsa_trend"] = [classify_trend(a, b) for a, b in zip(wide["ldsa_ci_low_pct_yr"], wide["ldsa_ci_high_pct_yr"])]
    wide["pollutant_trend"] = [classify_trend(a, b) for a, b in zip(wide["pollutant_ci_low_pct_yr"], wide["pollutant_ci_high_pct_yr"])]
    wide["trend_pair"] = wide["pollutant_trend"] + " / " + wide["ldsa_trend"]
    wide["declining_pollutant_ldsa_mismatch"] = (wide["pollutant_trend"] == "Declining") & (wide["ldsa_trend"] != "Declining")
    wide["sensitivity_3y_included"] = (
        (wide["n_months"] >= 36)
        & (wide["span_months"] >= 36)
        & (wide["coverage_fraction"] >= MIN_PRIMARY_COVERAGE)
        & (wide["site"] != "PAR_SUB")
    )
    wide = wide.sort_values(["pollutant", "site"])
    wide.to_csv(OUT / "01_exact_overlap_site_results.csv", index=False)

    levels = ["Declining", "Stable/uncertain", "Increasing"]
    matrix_rows = []
    for pollutant in EXTERNAL_METRICS:
        d = wide[(wide["pollutant"] == pollutant) & (wide["primary_included"] == True)]
        for pstate in levels:
            for lstate in levels:
                n = int(((d["pollutant_trend"] == pstate) & (d["ldsa_trend"] == lstate)).sum())
                matrix_rows.append({"pollutant": pollutant, "pollutant_trend": pstate, "ldsa_trend": lstate, "n": n, "percent": 100 * n / len(d) if len(d) else np.nan, "N": len(d)})
    matrix = pd.DataFrame(matrix_rows)
    matrix.to_csv(OUT / "01_trend_classification_matrix.csv", index=False)

    summary_rows = []
    for pollutant in EXTERNAL_METRICS:
        d = wide[(wide["pollutant"] == pollutant) & (wide["primary_included"] == True)].dropna(subset=["ldsa_slope_pct_yr", "pollutant_slope_pct_yr"])
        r, p = (pearsonr(d["pollutant_slope_pct_yr"], d["ldsa_slope_pct_yr"]) if len(d) >= 3 else (np.nan, np.nan))
        dec = d[d["pollutant_trend"] == "Declining"]
        mismatch = dec[dec["ldsa_trend"] != "Declining"]
        summary_rows.append(
            {
                "pollutant": pollutant,
                "N_primary": len(d),
                "Pearson_r": r,
                "Pearson_P": p,
                "pollutant_declining_N": len(dec),
                "declining_pollutant_ldsa_mismatch_N": len(mismatch),
                "mismatch_percent_among_pollutant_declining": 100 * len(mismatch) / len(dec) if len(dec) else np.nan,
                "median_ldsa_slope_pct_yr": d["ldsa_slope_pct_yr"].median(),
                "median_pollutant_slope_pct_yr": d["pollutant_slope_pct_yr"].median(),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "01_trend_summary.csv", index=False)

    sensitivity_rows = []
    for pollutant in EXTERNAL_METRICS:
        d = wide[(wide["pollutant"] == pollutant) & (wide["sensitivity_3y_included"] == True)].dropna(subset=["ldsa_slope_pct_yr", "pollutant_slope_pct_yr"])
        r, p = (pearsonr(d["pollutant_slope_pct_yr"], d["ldsa_slope_pct_yr"]) if len(d) >= 3 else (np.nan, np.nan))
        dec = d[d["pollutant_trend"] == "Declining"]
        mismatch = dec[dec["ldsa_trend"] != "Declining"]
        sensitivity_rows.append(
            {
                "pollutant": pollutant,
                "N_sensitivity_3y": len(d),
                "Pearson_r": r,
                "Pearson_P": p,
                "pollutant_declining_N": len(dec),
                "declining_pollutant_ldsa_mismatch_N": len(mismatch),
                "mismatch_percent_among_pollutant_declining": 100 * len(mismatch) / len(dec) if len(dec) else np.nan,
            }
        )
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(OUT / "01_trend_summary_3y_sensitivity.csv", index=False)
    return wide, matrix, summary, sensitivity


def kernel_breakpoint(values: np.ndarray, min_size: int) -> int | None:
    y = np.asarray(values, dtype=float)
    if len(y) < 2 * min_size + 1 or np.nanstd(y) == 0:
        return None
    z = (y - np.nanmean(y)) / np.nanstd(y)
    bkps = rpt.KernelCPD(kernel="linear", min_size=min_size, jump=1).fit(z.reshape(-1, 1)).predict(n_bkps=1)
    return int(bkps[0])


def breakpoint_analysis(pairs: pd.DataFrame, coverage: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cov_idx = coverage.set_index("pair_id")
    rows = []
    for pair_id, g in pairs.groupby("pair_id"):
        date_value = {}
        for role in ["LDSA", "Pollutant"]:
            r = g[g["role"] == role].sort_values("date")
            date_value[role] = (pd.to_datetime(r["date"]).reset_index(drop=True), r["value"].to_numpy(float))
        dates = date_value["LDSA"][0]
        if not dates.equals(date_value["Pollutant"][0]):
            raise RuntimeError(f"Paired dates differ for {pair_id}")
        info = cov_idx.loc[pair_id].to_dict()
        out = {"pair_id": pair_id, **info}
        for min_size, tag in [(2, "implemented"), (12, "robust12")]:
            idx_l = kernel_breakpoint(date_value["LDSA"][1], min_size)
            idx_p = kernel_breakpoint(date_value["Pollutant"][1], min_size)
            if idx_l is not None and idx_p is not None:
                bp_l, bp_p = dates.iloc[idx_l], dates.iloc[idx_p]
                out[f"ldsa_breakpoint_{tag}"] = bp_l
                out[f"pollutant_breakpoint_{tag}"] = bp_p
                out[f"offset_months_{tag}"] = months_between(bp_l, bp_p)
                out[f"ldsa_left_months_{tag}"] = idx_l
                out[f"ldsa_right_months_{tag}"] = len(dates) - idx_l
                out[f"pollutant_left_months_{tag}"] = idx_p
                out[f"pollutant_right_months_{tag}"] = len(dates) - idx_p
            else:
                out[f"ldsa_breakpoint_{tag}"] = pd.NaT
                out[f"pollutant_breakpoint_{tag}"] = pd.NaT
                out[f"offset_months_{tag}"] = np.nan
        out["breakpoint_sensitivity_difference_months"] = abs(out.get("offset_months_implemented", np.nan) - out.get("offset_months_robust12", np.nan))
        out["stable_within_12_months"] = out["breakpoint_sensitivity_difference_months"] <= 12
        rows.append(out)
    result = pd.DataFrame(rows).sort_values(["pollutant", "site"])
    result.to_csv(OUT / "02_breakpoint_offsets_site_results.csv", index=False)

    sum_rows = []
    for metric in EXTERNAL_METRICS:
        d = result[(result["pollutant"] == metric) & (result["primary_included"] == True)].dropna(subset=["offset_months_robust12"])
        v = d["offset_months_robust12"].astype(float)
        sum_rows.append(
            {
                "pollutant": metric,
                "N_primary": len(v),
                "median_months": v.median() if len(v) else np.nan,
                "Q1_months": v.quantile(0.25) if len(v) else np.nan,
                "Q3_months": v.quantile(0.75) if len(v) else np.nan,
                "P10_months": v.quantile(0.10) if len(v) else np.nan,
                "P90_months": v.quantile(0.90) if len(v) else np.nan,
                "earlier_N": int((v < 0).sum()),
                "same_month_N": int((v == 0).sum()),
                "later_N": int((v > 0).sum()),
                "within_12_months_N": int((v.abs() <= 12).sum()),
                "within_12_months_percent": 100 * (v.abs() <= 12).mean() if len(v) else np.nan,
                "12_to_24_months_N": int(((v.abs() > 12) & (v.abs() <= 24)).sum()),
                "over_24_months_N": int((v.abs() > 24).sum()),
                "sensitivity_stable_N": int(d["stable_within_12_months"].fillna(False).sum()),
            }
        )
    summary = pd.DataFrame(sum_rows)
    summary.to_csv(OUT / "02_breakpoint_offsets_summary.csv", index=False)
    return result, summary


def monthly_raw_pnc(site: str) -> pd.DataFrame:
    path = RAW_POLLUTANT_DIR / f"{site}.xlsx"
    cols = ["TIME UTC, start", "Nucleation", "Aitken", "Accumulation", "Ntotal"]
    df = pd.read_excel(path, sheet_name="Data", usecols=cols)
    out = None
    for c in cols[1:]:
        m = monthly_from_daily(df[cols[0]], df[c]).rename(columns={"value": c, "valid_days": f"{c}_valid_days"})
        out = m if out is None else out.merge(m, on="date", how="outer")
    return out


def full_ldsa_breakpoint(ldsa: pd.DataFrame) -> pd.Timestamp:
    d = ldsa[["date", "LDSA"]].dropna().sort_values("date")
    idx = kernel_breakpoint(d["LDSA"].to_numpy(), 12)
    if idx is None:
        raise RuntimeError("Unable to estimate LDSA breakpoint")
    return pd.to_datetime(d.iloc[idx]["date"])


def prepare_city_series(site: str, ldsa: pd.DataFrame) -> tuple[pd.DataFrame, pd.Timestamp, pd.DataFrame]:
    p = read_pollutant_monthly(site)
    pnc = monthly_raw_pnc(site)
    merged = ldsa[["date", "LDSA"]].copy()
    for metric in ["PM2.5", "PM10", "BC", "NO2"]:
        if metric in p:
            merged = merged.merge(p[metric][["date", "value"]].rename(columns={"value": metric}), on="date", how="outer")
    merged = merged.merge(pnc[["date", "Ntotal", "Nucleation", "Aitken", "Accumulation"]], on="date", how="outer")
    merged = merged.sort_values("date")
    bp = full_ldsa_breakpoint(ldsa)
    base_start = bp - pd.DateOffset(months=24)
    base_end = bp - pd.DateOffset(months=1)
    norm = merged.copy()
    summary_rows = []
    for c in [x for x in merged.columns if x != "date"]:
        base = merged.loc[(merged["date"] >= base_start) & (merged["date"] <= base_end), c].mean()
        norm[c] = 100 * merged[c] / base if pd.notna(base) and base != 0 else np.nan
        pre = norm.loc[(norm["date"] >= bp - pd.DateOffset(months=36)) & (norm["date"] < bp), c]
        post = norm.loc[(norm["date"] >= bp) & (norm["date"] < bp + pd.DateOffset(months=36)), c]
        summary_rows.append(
            {
                "site": site,
                "metric": c,
                "ldsa_breakpoint": bp,
                "baseline_start": base_start,
                "baseline_end": base_end,
                "baseline_mean_raw": base,
                "pre_3y_mean_index": pre.mean(),
                "post_3y_mean_index": post.mean(),
                "post_vs_pre_percent": 100 * (post.mean() / pre.mean() - 1) if pre.notna().sum() and post.notna().sum() and pre.mean() != 0 else np.nan,
                "pre_months": int(pre.notna().sum()),
                "post_months": int(post.notna().sum()),
            }
        )
    return norm, bp, pd.DataFrame(summary_rows)


def policy_context() -> pd.DataFrame:
    rows = [
        ["LND_TR", "Total road vehicle-km", "2013–2016", -2.4, "%", "Westminster major roads"],
        ["LND_TR", "Diesel-car vehicle-km", "2013–2016", 6.8, "%", "Westminster major roads"],
        ["LND_TR", "Diesel LGV vehicle-km", "2013–2016", 4.3, "%", "Westminster major roads"],
        ["LND_TR", "Diesel passenger cars Euro 5/6", "2016", 62.4, "% of fleet", "Central London CCZ"],
        ["LND_TR", "Diesel passenger cars Euro 4 or older", "2016", 37.6, "% of fleet", "Central London CCZ"],
        ["STO_TR", "Annual average daily traffic", "2013–2014", -4.0, "%", "Stockholm context"],
        ["STO_TR", "Diesel particulate-filter penetration", "2013", 31.9, "% of fleet", "Stockholm context"],
        ["STO_TR", "Diesel particulate-filter penetration", "2018", 48.4, "% of fleet", "Stockholm context"],
        ["STO_TR", "Diesel passenger cars Euro 5/6", "2017", 86.4, "% of fleet", "Stockholm context"],
        ["STO_TR", "Diesel light commercial vehicles Euro 5/6", "2017", 77.0, "% of fleet", "Stockholm context"],
        ["STO_TR", "Trucks Euro V/VI", "2017", 85.4, "% of fleet", "Stockholm context"],
        ["STO_TR", "Buses Euro VI", "2017", 97.0, "% of fleet", "Stockholm context"],
        ["STO_TR", "BC emission factor", "2013–2019", -66.0, "%", "Stockholm context"],
        ["STO_TR", "Particle-number emission factor", "2013–2019", -55.0, "%", "Stockholm context"],
    ]
    return pd.DataFrame(rows, columns=["site", "indicator", "period", "value", "unit", "evidence_scale"])


def figure_trend_mismatch(site_results: pd.DataFrame) -> None:
    set_plot_style()
    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.4), constrained_layout=True)
    metrics = ["PM2.5", "PM10", "BC", "NO2", "SO2"]
    for ax, metric, letter in zip(axes.flat, metrics, "abcde"):
        d = site_results[(site_results["pollutant"] == metric) & (site_results["primary_included"] == True)].copy()
        for typ, g in d.groupby("site_type"):
            ax.scatter(g["pollutant_slope_pct_yr"], g["ldsa_slope_pct_yr"], s=44, color=TYPE_COLORS.get(typ, "grey"), marker={"UB": "o", "SUB": "^", "TR": "s", "RB": "D"}.get(typ, "o"), alpha=0.9, label=typ)
        for _, r in d[d["declining_pollutant_ldsa_mismatch"]].iterrows():
            ax.annotate(r["site"], (r["pollutant_slope_pct_yr"], r["ldsa_slope_pct_yr"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
        ax.axhline(0, color="#A6BDD8", lw=1)
        ax.axvline(0, color="#A6BDD8", lw=1)
        display = {"PM2.5": r"PM$_{2.5}$", "PM10": r"PM$_{10}$", "NO2": r"NO$_2$", "SO2": r"SO$_2$"}.get(metric, metric)
        ax.set_xlabel(f"{display} slope (% yr$^{{-1}}$)")
        ax.set_ylabel("LDSA slope (% yr$^{-1}$)")
        ax.set_title(f"{letter}  {display}", loc="left")
        ax.grid(True, ls="--", lw=0.5, alpha=0.35)
        dec = d[d["pollutant_trend"] == "Declining"]
        mm = d[d["declining_pollutant_ldsa_mismatch"]]
        ax.text(0.98, 0.96, f"N = {len(d)}\nMismatch = {len(mm)}/{len(dec) if len(dec) else 0}", transform=ax.transAxes, ha="right", va="top", fontsize=9)
    ax = axes.flat[-1]
    ax.axis("off")
    handles = [Line2D([0], [0], marker={"UB": "o", "SUB": "^", "TR": "s", "RB": "D"}[t], color="none", markerfacecolor=TYPE_COLORS[t], markeredgecolor=TYPE_COLORS[t], markersize=8, label=t) for t in ["UB", "SUB", "TR", "RB"]]
    ax.legend(handles=handles, loc="center", frameon=False, ncol=2, title="Station type")
    fig.savefig(OUT / "Figure_A1_exact_overlap_trend_mismatch.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "Figure_A1_exact_overlap_trend_mismatch.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_breakpoints(results: pd.DataFrame, summary: pd.DataFrame) -> None:
    set_plot_style()
    metrics = FOCUS_METRICS
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    rng = np.random.default_rng(RANDOM_SEED)
    for i, metric in enumerate(metrics):
        d = results[(results["pollutant"] == metric) & (results["primary_included"] == True)].dropna(subset=["offset_months_robust12"])
        jitter = rng.normal(0, 0.055, len(d))
        colors = [TYPE_COLORS.get(t, "grey") for t in d["site_type"]]
        ax.scatter(np.full(len(d), i) + jitter, d["offset_months_robust12"], c=colors, s=38, alpha=0.8, zorder=3)
        s = summary[summary["pollutant"] == metric]
        if not s.empty:
            med, q1, q3, p10, p90 = [float(s.iloc[0][x]) for x in ["median_months", "Q1_months", "Q3_months", "P10_months", "P90_months"]]
            ax.plot([i, i], [p10, p90], color="black", lw=1.4, zorder=4)
            ax.plot([i - 0.16, i + 0.16], [q1, q1], color="black", lw=4, zorder=4)
            ax.plot([i - 0.16, i + 0.16], [q3, q3], color="black", lw=4, zorder=4)
            ax.scatter([i], [med], color="white", edgecolor="black", s=62, zorder=5)
    ax.axhline(0, color="black", lw=1, ls="--")
    ax.set_xticks(range(len(metrics)), metrics)
    ax.set_ylabel("Breakpoint offset relative to LDSA (months)")
    ax.set_xlabel("Negative = metric earlier; positive = metric later")
    ax.grid(axis="y", ls="--", lw=0.5, alpha=0.4)
    handles = [Line2D([0], [0], marker={"UB": "o", "SUB": "^", "TR": "s", "RB": "D"}[t], color="none", markerfacecolor=TYPE_COLORS[t], markeredgecolor=TYPE_COLORS[t], markersize=7, label=t) for t in ["UB", "SUB", "TR", "RB"]]
    ax.legend(handles=handles, frameon=False, ncol=4, loc="upper center")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_A2_breakpoint_offsets.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "Figure_A2_breakpoint_offsets.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_city_series(city_data: dict[str, tuple[pd.DataFrame, pd.Timestamp]]) -> None:
    set_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 7.1), sharex="col", constrained_layout=True)
    for col, site in enumerate(["LND_TR", "STO_TR"]):
        d, bp = city_data[site]
        panels = [(["LDSA", "Ntotal", "Nucleation", "Aitken", "Accumulation"], "Particle metrics"), (["BC", "PM2.5", "PM10", "NO2"], "Conventional metrics")]
        for row, (metrics, subtitle) in enumerate(panels):
            ax = axes[row, col]
            for metric in metrics:
                if metric in d.columns and d[metric].notna().sum() >= 6:
                    ax.plot(d["date"], d[metric], lw=1.6, label=metric, color=METRIC_COLORS.get(metric))
            ax.axvline(bp, color="black", ls="--", lw=1.3)
            ax.axhline(100, color="#BDBDBD", ls=":", lw=1)
            ax.set_ylabel("Index (pre-breakpoint 24-month mean = 100)")
            ax.set_title(f"{site}: {subtitle}" if row == 0 else subtitle)
            ax.grid(True, ls="--", lw=0.4, alpha=0.35)
            ax.legend(frameon=False, ncol=2 if row == 0 else 4, fontsize=8, loc="best")
            if row == 0:
                ax.text(bp, ax.get_ylim()[1], f"  LDSA breakpoint: {bp:%Y-%m}", ha="left", va="top", fontsize=8)
    fig.savefig(OUT / "Figure_A3_London_Stockholm_synchronized_series.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "Figure_A3_London_Stockholm_synchronized_series.pdf", bbox_inches="tight")
    plt.close(fig)


def write_readme(summary: pd.DataFrame, sensitivity: pd.DataFrame, bp_summary: pd.DataFrame, city_summary: pd.DataFrame) -> None:
    pm = summary[summary["pollutant"] == "PM2.5"].iloc[0]
    pm3 = sensitivity[sensitivity["pollutant"] == "PM2.5"].iloc[0]
    focus = bp_summary[bp_summary["pollutant"].isin(FOCUS_METRICS)].copy()
    city_pivot = city_summary.pivot(index="metric", columns="site", values="post_vs_pre_percent")
    lines = [
        "# LDSA and co-pollutant analyses",
        "",
        "## 1. Exact-overlap LDSA–pollutant trends",
        "",
        f"Primary analysis required at least {MIN_PRIMARY_MONTHS} paired monthly observations and at least {MIN_PRIMARY_COVERAGE:.0%} coverage of the common span, and excluded PAR_SUB from the headline summary because of its documented record discontinuity.",
        f"For PM2.5, {int(pm['N_primary'])} site pairs qualified. Among {int(pm['pollutant_declining_N'])} sites with a significantly declining PM2.5 trend, {int(pm['declining_pollutant_ldsa_mismatch_N'])} ({pm['mismatch_percent_among_pollutant_declining']:.1f}%) did not have a significantly declining LDSA trend.",
        f"The exact-overlap Pearson correlation for PM2.5 was r = {pm['Pearson_r']:.3f} (P = {pm['Pearson_P']:.3g}). This correlation is descriptive; the 3×3 trend-classification matrix is the primary mismatch result.",
        f"A less restrictive sensitivity analysis requiring at least 36 paired months retained {int(pm3['N_sensitivity_3y'])} PM2.5 site pairs; {int(pm3['declining_pollutant_ldsa_mismatch_N'])}/{int(pm3['pollutant_declining_N'])} sites with declining PM2.5 lacked a significant LDSA decline.",
        "",
        "## 2. Metric-specific breakpoint timing",
        "",
        "Positive offsets mean the co-pollutant breakpoint occurred after the LDSA breakpoint; negative offsets mean it occurred earlier. The robust result uses KernelCPD with a minimum segment length of 12 months; the as-implemented two-month setting is retained as a sensitivity analysis.",
    ]
    for _, r in focus.iterrows():
        lines.append(f"- {r['pollutant']}: N={int(r['N_primary'])}, median {r['median_months']:.1f} months (IQR {r['Q1_months']:.1f} to {r['Q3_months']:.1f}); {r['within_12_months_percent']:.1f}% within ±12 months.")
    lines += ["", "## 3. London and Stockholm synchronized trajectories", ""]
    for site in ["LND_TR", "STO_TR"]:
        vals = []
        for metric in ["LDSA", "Ntotal", "BC", "PM2.5", "PM10", "NO2"]:
            if metric in city_pivot.index and site in city_pivot.columns and pd.notna(city_pivot.loc[metric, site]):
                vals.append(f"{metric} {city_pivot.loc[metric, site]:+.1f}%")
        lines.append(f"- {site}, mean index in the first 36 months after versus 36 months before the LDSA breakpoint: " + "; ".join(vals) + ".")
    lines += [
        "",
        "The policy/fleet indicators are provided as contextual annual or multi-year values and are not interpolated into monthly observations. Therefore, the city plots support temporal co-evolution but do not by themselves establish policy causality.",
        "",
        "## Sign convention and quality control",
        "",
        "- Trend classes: declining if the entire 95% CI is below zero; increasing if the entire 95% CI is above zero; otherwise stable/uncertain.",
        "- Monthly values require at least 14 valid calendar days.",
        "- The exact same paired months are used for the LDSA and co-pollutant trend at each site.",
        "- PAR_SUB is retained in the full tables but excluded from primary summaries because of its documented discontinuity.",
    ]
    (OUT / "README_analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    np.random.seed(RANDOM_SEED)
    meta = site_metadata()
    meta.to_csv(OUT / "00_site_metadata.csv", index=False)
    pairs, coverage, ldsa_all = prepare_monthly_pairs(meta)
    trends = run_r_trends()
    site_results, matrix, trend_summary, trend_sensitivity = summarize_trend_pairs(trends)
    bp_results, bp_summary = breakpoint_analysis(pairs, coverage)

    city_data = {}
    city_summaries = []
    for site in ["LND_TR", "STO_TR"]:
        norm, bp, csum = prepare_city_series(site, ldsa_all[site])
        norm.to_csv(OUT / f"03_{site}_normalized_monthly_series.csv", index=False)
        city_data[site] = (norm, bp)
        city_summaries.append(csum)
    city_summary = pd.concat(city_summaries, ignore_index=True)
    city_summary.to_csv(OUT / "03_city_pre_post_summary.csv", index=False)
    policy = policy_context()
    policy.to_csv(OUT / "03_policy_fleet_context.csv", index=False)

    figure_trend_mismatch(site_results)
    figure_breakpoints(bp_results, bp_summary)
    figure_city_series(city_data)
    write_readme(trend_summary, trend_sensitivity, bp_summary, city_summary)

    manifest = {
        "analysis_date": pd.Timestamp.today().date().isoformat(),
        "random_seed": RANDOM_SEED,
        "minimum_valid_days_per_month": MIN_VALID_DAYS_MONTH,
        "minimum_primary_paired_months": MIN_PRIMARY_MONTHS,
        "minimum_primary_common_span_coverage": MIN_PRIMARY_COVERAGE,
        "primary_exclusion": "PAR_SUB documented discontinuity",
        "trend_method": "openair::TheilSen, deseason=FALSE, slope.percent=TRUE, same paired months",
        "breakpoint_method": "ruptures KernelCPD linear kernel; min_size=12 primary, min_size=2 sensitivity",
        "output_files": sorted([p.name for p in OUT.iterdir() if p.is_file()]),
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
