#!/usr/bin/env python3
"""Render station-level LDSA time series, breakpoints and phase trends."""

from __future__ import annotations

import argparse
from pathlib import Path
import math

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D


TARGET_POLLUTANT = "LDSA_17_500nm"
SHEET_NAMES = ("LDSA_part_1", "LDSA_part_2")


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Render individual station breakpoint diagnostics")
    parser.add_argument("--root", type=Path, default=script_dir, help="Directory containing breakpoint and station time-series workbooks")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir / "breakpoint_all_sites_individual",
        help="Directory for station diagnostic figures",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution in dots per inch (default: 300)")
    parser.add_argument(
        "--legend-file",
        type=Path,
        default=None,
        help="Optional output path for a separate legend figure",
    )
    return parser.parse_args()


def configure_style() -> None:
    sns.set_theme(style="ticks", context="paper", font_scale=2.4)
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman"],
            "font.size": 19.2,
            "axes.labelsize": 19.2,
            "xtick.labelsize": 17.6,
            "ytick.labelsize": 17.6,
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def load_breakpoints(root: Path) -> pd.DataFrame:
    breakpoint_file = root / "17-500LDSAbreakpoint.xlsx"
    frame = pd.read_excel(breakpoint_file)
    frame = frame[frame["Pollutant"].eq(TARGET_POLLUTANT)].copy()
    if frame.empty:
        raise ValueError(f"No breakpoint records found for {TARGET_POLLUTANT}.")
    return frame


def load_site_weekly_data(site: str, site_data_dir: Path) -> pd.DataFrame:
    site_file = site_data_dir / f"{site}.xlsx"
    if not site_file.exists():
        raise FileNotFoundError(f"Station data file not found: {site_file}")

    frames: list[pd.DataFrame] = []
    with pd.ExcelFile(site_file) as workbook:
        available = [name for name in SHEET_NAMES if name in workbook.sheet_names]
        if not available and "Total" in workbook.sheet_names:
            available = ["Total"]
        if not available:
            raise ValueError(f"No supported LDSA worksheet found in {site_file}.")
        for sheet_name in available:
            frames.append(
                pd.read_excel(
                    workbook,
                    sheet_name=sheet_name,
                    usecols=["Time", TARGET_POLLUTANT],
                )
            )

    data = pd.concat(frames, ignore_index=True)
    data["Time"] = pd.to_datetime(data["Time"], errors="coerce")
    data[TARGET_POLLUTANT] = pd.to_numeric(data[TARGET_POLLUTANT], errors="coerce")
    data = data.dropna(subset=["Time", TARGET_POLLUTANT])

    weekly = (
        data.set_index("Time")[TARGET_POLLUTANT]
        .resample("W")
        .mean()
        .dropna()
        .reset_index()
        .rename(columns={"Time": "Date", TARGET_POLLUTANT: "Value"})
    )
    return weekly.sort_values("Date").reset_index(drop=True)


def first_row(frame: pd.DataFrame, period: str) -> pd.Series | None:
    rows = frame[frame["Period"].eq(period)]
    return None if rows.empty else rows.iloc[0]


def to_number(row: pd.Series | None, column: str) -> float:
    if row is None:
        return np.nan
    return float(pd.to_numeric(row.get(column), errors="coerce"))


def breakpoint_info(site_rows: pd.DataFrame) -> tuple[str | None, pd.Timestamp, float, float]:
    pre_row = first_row(site_rows, "Before")
    post_row = first_row(site_rows, "After")

    if pre_row is not None:
        breakpoint_text = str(pre_row["End_Time_Range"])
    elif post_row is not None:
        breakpoint_text = str(post_row["Start_Time_Range"])
    else:
        breakpoint_text = None

    breakpoint_date = pd.to_datetime(breakpoint_text, errors="coerce")
    if pd.notna(breakpoint_date):
        breakpoint_text = breakpoint_date.strftime("%Y/%m")
    pre_slope = to_number(pre_row, "Slope_Unit")
    post_slope = to_number(post_row, "Slope_Unit")
    return breakpoint_text, breakpoint_date, pre_slope, post_slope


def trend_values(segment: pd.DataFrame, slope: float) -> np.ndarray:
    numeric_dates = mdates.date2num(segment["Date"])
    trend = (slope / 365.25) * numeric_dates
    return trend + segment["Value"].mean() - trend.mean()


def format_slope(value: float) -> str:
    return f"{value:.2f}" if pd.notna(value) else "NA"


def draw_site(
    site: str,
    weekly: pd.DataFrame,
    site_rows: pd.DataFrame,
    output_file: Path,
    dpi: int,
) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    fig.subplots_adjust(left=0.125, right=0.985, bottom=0.095, top=0.86)

    ax.scatter(
        weekly["Date"], weekly["Value"],
        color="#95a5a6", alpha=0.50, s=18, label="Weekly Mean", zorder=2,
    )
    ax.plot(
        weekly["Date"], weekly["Value"],
        color="#95a5a6", alpha=0.20, linewidth=1.0, zorder=1,
    )

    bp_text, bp_date, pre_slope, post_slope = breakpoint_info(site_rows)
    if pd.notna(bp_date):
        ax.axvline(bp_date, color="black", linestyle="--", linewidth=1.8, zorder=3)

        if pd.notna(pre_slope):
            pre_segment = weekly[weekly["Date"] <= bp_date]
            if not pre_segment.empty:
                ax.plot(
                    pre_segment["Date"], trend_values(pre_segment, pre_slope),
                    color="#e74c3c", linewidth=3.0,
                    label="Pre-Breakpoint Trend", zorder=4,
                )

        if pd.notna(post_slope):
            post_segment = weekly[weekly["Date"] >= bp_date]
            if not post_segment.empty:
                ax.plot(
                    post_segment["Date"], trend_values(post_segment, post_slope),
                    color="#2980b9", linewidth=3.0,
                    label="Post-Breakpoint Trend", zorder=4,
                )
    annotation = (
        f"Breakpoint: {bp_text if bp_text is not None else 'NA'}\n"
        f"Pre-slope: {format_slope(pre_slope)}\n"
        f"Post-slope: {format_slope(post_slope)}"
    )
    ax.text(
        0.03, 0.97, annotation,
        transform=ax.transAxes, fontsize=20,
        fontfamily="Times New Roman", fontweight="bold",
        ha="left", va="top",
    )
    ax.set_ylabel(
        "LDSA (μm²/cm³)",
        fontfamily="Times New Roman",
        fontsize=22,
        fontweight="bold",
    )
    ax.set_xlabel("")
    ax.grid(True, linestyle=":", alpha=0.60)
    if not weekly.empty:
        span_years = max(
            1.0,
            (weekly["Date"].max() - weekly["Date"].min()).days / 365.25,
        )
        year_interval = max(1, math.ceil(span_years / 6.0))
        ax.xaxis.set_major_locator(mdates.YearLocator(base=year_interval))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.margins(x=0.02)
    for tick_label in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        tick_label.set_fontfamily("Times New Roman")
        tick_label.set_fontsize(18)
        tick_label.set_fontweight("bold")
    x_tick_labels = ax.get_xticklabels()
    if x_tick_labels:
        x_tick_labels[-1].set_ha("right")
    fig.text(
        0.555, 0.91, site,
        fontfamily="Times New Roman", fontsize=26, fontweight="bold",
        ha="center", va="center",
    )
    fig.savefig(output_file, dpi=dpi, facecolor="white")
    plt.close(fig)


def draw_legend(output_file: Path, dpi: int) -> None:
    """Save the shared legend for the station diagnostic figures."""
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.axis("off")
    handles = [
        Line2D(
            [0], [0], marker="o", linestyle="none",
            markerfacecolor="#95a5a6", markeredgecolor="#95a5a6",
            alpha=0.50, markersize=24, label="Weekly Mean",
        ),
        Line2D(
            [0], [0], color="#e74c3c", linewidth=10.0,
            label="Pre-Breakpoint Trend",
        ),
        Line2D(
            [0], [0], color="#2980b9", linewidth=10.0,
            label="Post-Breakpoint Trend",
        ),
    ]
    ax.legend(
        handles=handles, loc="center", frameon=False,
        prop={"family": "Times New Roman", "size": 44, "weight": "bold"},
        handlelength=1.6, handletextpad=0.5,
        borderaxespad=0,
    )
    fig.savefig(output_file, dpi=dpi, facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    root = args.root.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    site_data_dir = root / "clean17-500" / "station_series"

    configure_style()
    breakpoints = load_breakpoints(root)
    sites = sorted(breakpoints["Site"].dropna().astype(str).unique())

    outputs: list[Path] = []
    for site in sites:
        weekly = load_site_weekly_data(site, site_data_dir)
        site_rows = breakpoints[breakpoints["Site"].astype(str).eq(site)]
        output_file = output_dir / f"{site}_breakpoint.png"
        draw_site(site, weekly, site_rows, output_file, args.dpi)
        outputs.append(output_file)
        print(f"Saved figure: {output_file.name}")

    legend_file = (
        args.legend_file.expanduser().resolve()
        if args.legend_file is not None
        else output_dir / "breakpoint_legend.png"
    )
    legend_file.parent.mkdir(parents=True, exist_ok=True)
    draw_legend(legend_file, args.dpi)
    print(f"Saved legend: {legend_file.name}")

    print(f"Generated plots for {len(outputs)} stations.")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
