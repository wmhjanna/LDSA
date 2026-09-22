#!/usr/bin/env python3
"""Redraw the LND_TR/STO_TR breakpoint-context figure from local source data.

The time-series panels use cubic B-spline smoothing for display only.
All three-year percentage changes are calculated from the unsmoothed monthly
means, with the established LDSA breakpoints fixed at March 2013 (LND_TR) and
October 2013 (STO_TR).
"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import nature_analysis_123 as source


ROOT = Path("@PROJECT_ROOT@")
OUTPUT_DIR = ROOT / "outputs" / "lnd_sto_breakpoint_context_smooth"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

BREAKPOINTS = {
    "LND_TR": pd.Timestamp("2013-03-01"),
    "STO_TR": pd.Timestamp("2013-10-01"),
}

BSPLINE_SMOOTHING_FRACTION = 0.045
BSPLINE_POINTS_PER_MONTH = 8
BSPLINE_MAX_GAP_DAYS = 100

SITE_CONFIG = {
    "LND_TR": {
        "city": "Marylebone Road",
        "panel_line": "a",
        "panel_change": "b",
        "line_metrics": ["LDSA", "Ntotal", "BC", "PM2.5", "NO2"],
        "change_metrics": ["LDSA", "Ntotal", "BC", "PM2.5", "PM10", "NO2"],
        "note": (
            "Road VKM −2.4%; diesel-car +6.8%; diesel-LGV +4.3%\n"
            "Diesel passenger cars Euro 5/6: 62.4% (2016)"
        ),
    },
    "STO_TR": {
        "city": "Hornsgatan",
        "panel_line": "c",
        "panel_change": "d",
        "line_metrics": ["LDSA", "Ntotal", "PM2.5", "PM10", "NO2"],
        "change_metrics": ["LDSA", "Ntotal", "PM2.5", "PM10", "NO2"],
        "note": "AADT −4.0%; DPF 31.9→48.4%\nBC EF −66%; PN EF −55% (2013–2019)",
    },
}

DISPLAY = {
    "LDSA": "LDSA",
    "Ntotal": "PNC total",
    "BC": "BC",
    "PM2.5": r"PM$_{2.5}$",
    "PM10": r"PM$_{10}$",
    "NO2": r"NO$_2$",
}

COLORS = {
    "LDSA": "#202020",
    "Ntotal": "#009E73",
    "BC": "#7B3294",
    "PM2.5": "#D55E00",
    "PM10": "#E69F00",
    "NO2": "#0072B2",
}


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman"],
            "mathtext.fontset": "stix",
            "font.size": 11.5,
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "axes.linewidth": 1.2,
            "xtick.major.width": 1.1,
            "ytick.major.width": 1.1,
            "xtick.major.size": 4.5,
            "ytick.major.size": 4.5,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def load_site_data(site: str) -> pd.DataFrame:
    ldsa = source.read_ldsa_monthly(site)[["date", "LDSA"]]
    pollutants = source.read_pollutant_monthly(site)
    pnc = source.monthly_raw_pnc(site)

    merged = ldsa.copy()
    for metric in ["PM2.5", "PM10", "BC", "NO2"]:
        if metric in pollutants:
            monthly = pollutants[metric][["date", "value"]].rename(columns={"value": metric})
            merged = merged.merge(monthly, on="date", how="outer")
    merged = merged.merge(pnc[["date", "Ntotal"]], on="date", how="outer")
    return merged.sort_values("date").reset_index(drop=True)


def normalize_to_pre_breakpoint_mean(data: pd.DataFrame, breakpoint: pd.Timestamp) -> pd.DataFrame:
    result = data.copy()
    start = breakpoint - pd.DateOffset(months=24)
    end = breakpoint - pd.DateOffset(months=1)
    for column in [c for c in result.columns if c != "date"]:
        baseline = result.loc[(result["date"] >= start) & (result["date"] <= end), column].mean()
        result[column] = 100.0 * result[column] / baseline if pd.notna(baseline) and baseline != 0 else np.nan
    return result


def three_year_changes(data: pd.DataFrame, breakpoint: pd.Timestamp, metrics: list[str]) -> dict[str, float]:
    values: dict[str, float] = {}
    for metric in metrics:
        pre = data.loc[
            (data["date"] >= breakpoint - pd.DateOffset(months=36))
            & (data["date"] < breakpoint),
            metric,
        ].dropna()
        post = data.loc[
            (data["date"] >= breakpoint)
            & (data["date"] < breakpoint + pd.DateOffset(months=36)),
            metric,
        ].dropna()
        values[metric] = 100.0 * (post.mean() / pre.mean() - 1.0)
    return values


def bspline_segments(dates: pd.Series, values: pd.Series):
    """Yield smooth cubic B-spline segments without bridging long data gaps."""
    valid = pd.DataFrame({"date": pd.to_datetime(dates), "value": values}).dropna().sort_values("date")
    if valid.empty:
        return

    gap_group = valid["date"].diff().dt.days.gt(BSPLINE_MAX_GAP_DAYS).cumsum()
    for _, segment in valid.groupby(gap_group):
        x = mdates.date2num(segment["date"].to_numpy())
        y = segment["value"].to_numpy(dtype=float)
        if len(segment) < 4 or np.allclose(np.nanvar(y), 0):
            yield segment["date"].to_numpy(), y
            continue

        smoothing = len(y) * np.nanvar(y) * BSPLINE_SMOOTHING_FRACTION
        spline = UnivariateSpline(x, y, k=3, s=smoothing)
        span_months = max(1.0, (x[-1] - x[0]) / 30.4375)
        n_dense = max(len(x), int(np.ceil(span_months * BSPLINE_POINTS_PER_MONTH)))
        x_dense = np.linspace(x[0], x[-1], n_dense)
        yield mdates.num2date(x_dense), spline(x_dense)


def style_axis(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#222222")
        spine.set_linewidth(1.2)
    ax.tick_params(labelsize=10.8)
    for label in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
        label.set_fontweight("bold")


def add_panel_letter(ax: plt.Axes, letter: str) -> None:
    ax.text(
        -0.14,
        1.025,
        letter,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=16,
        fontweight="bold",
        clip_on=False,
    )


def draw_time_series(ax: plt.Axes, site: str, data: pd.DataFrame) -> None:
    config = SITE_CONFIG[site]
    breakpoint = BREAKPOINTS[site]
    normalized = normalize_to_pre_breakpoint_mean(data, breakpoint)

    displayed_values: list[np.ndarray] = []
    for metric in config["line_metrics"]:
        # B-spline smoothing affects only the displayed curve. The lollipop
        for segment_number, (smooth_dates, smooth_values) in enumerate(
            bspline_segments(normalized["date"], normalized[metric])
        ):
            displayed_values.append(np.asarray(smooth_values, dtype=float))
            ax.plot(
                smooth_dates,
                smooth_values,
                color=COLORS[metric],
                linewidth=1.9,
                label=DISPLAY[metric] if segment_number == 0 else "_nolegend_",
                solid_capstyle="round",
                zorder=3,
            )

    ax.axhline(100, color="#9A9A9A", linestyle=(0, (2, 2)), linewidth=0.9, zorder=1)
    ax.axvline(breakpoint, color="#333333", linestyle="--", linewidth=1.25, zorder=2)
    breakpoint_label_y = 250.0 if site == "STO_TR" else 20.0
    ax.annotate(
        f"LDSA bp {breakpoint:%Y-%m}",
        xy=(breakpoint, breakpoint_label_y),
        xycoords="data",
        xytext=(5, 2),
        textcoords="offset points",
        rotation=0,
        ha="left",
        va="bottom",
        fontsize=10.8,
        fontweight="bold",
        color="#333333",
        zorder=5,
    )

    displayed = np.concatenate(displayed_values)
    low = np.nanpercentile(displayed, 1)
    high = np.nanpercentile(displayed, 99)
    padding = max(8.0, 0.08 * (high - low))
    display_max = np.nanmax(displayed)
    if site == "STO_TR":
        # for both the PM10 maximum and the vertical breakpoint annotation.
        ymax = np.ceil((display_max + 20.0) / 25.0) * 25.0
        ax.set_ylim(max(0, low - padding), ymax)
    else:
        # Extra headroom keeps the upper-right legend clear of the curves and
        # makes the breakpoint annotation legible without changing the data.
        ymax = max(160.0, np.ceil((display_max + 20.0) / 20.0) * 20.0)
        ax.set_ylim(max(0, low - padding), ymax)
    ax.set_xlim(normalized["date"].min() - pd.DateOffset(months=2), normalized["date"].max() + pd.DateOffset(months=2))
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.set_ylabel("Index (pre-breakpoint 24-month mean = 100)", fontsize=11.5, fontweight="bold")
    ax.set_title(f"{site} ({config['city']})", loc="left", fontsize=13.5, pad=8, fontweight="bold")
    ax.grid(True, color="#D8D8D8", linestyle=":", linewidth=0.75, alpha=0.75, zorder=0)
    ax.legend(
        loc="upper right",
        ncol=2,
        fontsize=9.2,
        handlelength=2.0,
        columnspacing=1.0,
        labelspacing=0.35,
        borderaxespad=0.4,
        prop={"family": "Times New Roman", "weight": "bold", "size": 9.2},
    )
    add_panel_letter(ax, config["panel_line"])
    style_axis(ax)


def draw_change_panel(ax: plt.Axes, site: str, data: pd.DataFrame) -> None:
    config = SITE_CONFIG[site]
    breakpoint = BREAKPOINTS[site]
    metrics = config["change_metrics"]
    values = three_year_changes(data, breakpoint, metrics)
    y = np.arange(len(metrics))

    for ypos, metric in zip(y, metrics):
        value = values[metric]
        ax.hlines(ypos, min(value, 0), max(value, 0), color="#B9B9B9", linewidth=1.5, zorder=1)
        ax.scatter(value, ypos, s=64, color=COLORS[metric], edgecolor="white", linewidth=0.8, zorder=3)
        offset = 5 if value <= 0 else 5
        ax.annotate(
            f"{value:+.1f}%" if value > 0 else f"{value:.1f}%",
            (value, ypos),
            xytext=(offset, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=9.5,
            fontweight="bold",
        )

    ax.axvline(0, color="#333333", linewidth=1.1, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([DISPLAY[m] for m in metrics], fontsize=10.8, fontweight="bold")
    # Reserve a dedicated blank band below the last metric for the contextual
    # note so that it cannot overlap the NO2 point or its horizontal line.
    ax.set_ylim(len(metrics) + 0.9, -0.5)
    all_values = np.array(list(values.values()), dtype=float)
    xmin = min(-45.0, np.nanmin(all_values) - 5.0)
    xmax = max(12.0, np.nanmax(all_values) + 8.0)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("3-y post vs pre change (%)", fontsize=11.5, fontweight="bold")
    ax.set_title("3-y changes around LDSA breakpoint", fontsize=13.5, pad=8, fontweight="bold")
    ax.grid(axis="x", color="#D8D8D8", linestyle=":", linewidth=0.75, alpha=0.75, zorder=0)
    ax.text(
        0.985,
        0.035,
        config["note"],
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.8,
        fontweight="bold",
        color="#555555",
        zorder=4,
    )
    add_panel_letter(ax, config["panel_change"])
    style_axis(ax)


def main() -> None:
    set_style()
    data = {site: load_site_data(site) for site in SITE_CONFIG}

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12.8, 8.2),
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )
    draw_time_series(axes[0, 0], "LND_TR", data["LND_TR"])
    draw_change_panel(axes[0, 1], "LND_TR", data["LND_TR"])
    draw_time_series(axes[1, 0], "STO_TR", data["STO_TR"])
    draw_change_panel(axes[1, 1], "STO_TR", data["STO_TR"])

    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.095, top=0.94, wspace=0.22, hspace=0.24)
    base = OUTPUT_DIR / "LND_STO_breakpoint_context_smooth"
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)

    rows = []
    for site, frame in data.items():
        for metric, value in three_year_changes(frame, BREAKPOINTS[site], SITE_CONFIG[site]["change_metrics"]).items():
            rows.append(
                {
                    "site": site,
                    "breakpoint": BREAKPOINTS[site].strftime("%Y-%m"),
                    "metric": DISPLAY[metric].replace("$_{", "").replace("}$", ""),
                    "post_vs_pre_change_percent": value,
                }
            )
    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "LND_STO_3y_changes_fixed_breakpoints.csv", index=False)
    print(base.with_suffix(".png"))


if __name__ == "__main__":
    main()
