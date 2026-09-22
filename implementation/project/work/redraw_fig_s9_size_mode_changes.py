from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


LDSA_DIR = Path(
    "@INPUT_ROOT@/data_sources/plotting/breakpoint_plots/clean17-500"
)
RAW_PNSD_DIR = Path(
    "@INPUT_ROOT@/data_sources/predictive_models/pollutant_data/all_sites"
)
OUTPUT_DIR = Path(
    "@PROJECT_ROOT@/outputs/fig_s9_python_redraw"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "mathtext.fontset": "stix",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "axes.linewidth": 1.35,
        "xtick.major.width": 1.25,
        "ytick.major.width": 1.25,
        "xtick.major.size": 5.0,
        "ytick.major.size": 5.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def monthly_mean(date, value, min_valid_days=14):
    """Calendar-day-weighted monthly means used in the main analysis."""
    d = pd.DataFrame(
        {
            "date": pd.to_datetime(date, errors="coerce"),
            "value": pd.to_numeric(value, errors="coerce"),
        }
    ).dropna()
    d["day"] = d["date"].dt.floor("D")
    daily = d.groupby("day", as_index=False)["value"].mean()
    daily["date"] = daily["day"].dt.to_period("M").dt.to_timestamp()
    monthly = daily.groupby("date", as_index=False).agg(
        value=("value", "mean"), valid_days=("value", "count")
    )
    monthly.loc[monthly["valid_days"] < min_valid_days, "value"] = pd.NA
    return monthly.dropna(subset=["value"])


def calculate_site_changes(site, breakpoint):
    """Calculate 36-month post-vs-pre changes from the established local data."""
    ldsa_path = LDSA_DIR / f"{site}_Integrated_LDSA.xlsx"
    raw_path = RAW_PNSD_DIR / f"{site}.xlsx"

    ldsa_raw = pd.read_excel(
        ldsa_path,
        sheet_name="Total",
        usecols=["Time", "LDSA_17_500nm"],
    )
    ldsa = monthly_mean(ldsa_raw["Time"], ldsa_raw["LDSA_17_500nm"]).rename(
        columns={"value": "LDSA"}
    )

    pnsd_raw = pd.read_excel(
        raw_path,
        sheet_name="Data",
        usecols=["TIME UTC, start", "Nucleation", "Aitken", "Accumulation", "Ntotal"],
    )
    merged = ldsa[["date", "LDSA"]].copy()
    for column in ["Ntotal", "Nucleation", "Aitken", "Accumulation"]:
        monthly = monthly_mean(pnsd_raw["TIME UTC, start"], pnsd_raw[column]).rename(
            columns={"value": column}
        )
        merged = merged.merge(monthly[["date", column]], on="date", how="outer")

    bp = pd.Timestamp(breakpoint)
    columns = {
        "LDSA": "LDSA",
        "Total PNC": "Ntotal",
        "Nucleation PNC": "Nucleation",
        "Aitken PNC": "Aitken",
        "Accumulation PNC": "Accumulation",
    }
    values = {}
    for label, column in columns.items():
        pre = merged.loc[
            (merged["date"] >= bp - pd.DateOffset(months=36))
            & (merged["date"] < bp),
            column,
        ]
        post = merged.loc[
            (merged["date"] >= bp)
            & (merged["date"] < bp + pd.DateOffset(months=36)),
            column,
        ]
        values[label] = 100 * (post.mean() / pre.mean() - 1)
    return values


metric_order = [
    "LDSA",
    "Total PNC",
    "Nucleation PNC",
    "Aitken PNC",
    "Accumulation PNC",
]
site_config = [
    ("LND_TR", "London", "2013-03-01", "Mar 2013", "a"),
    ("STO_TR", "Stockholm", "2013-10-01", "Oct 2013", "b"),
]

fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.7), sharey=False)
point_color = "#2C7FB8"

for ax, (site, city, breakpoint, breakpoint_label, panel) in zip(axes, site_config):
    site_values = calculate_site_changes(site, breakpoint)
    values = [site_values[m] for m in metric_order]
    y = list(range(len(metric_order)))

    ax.scatter(
        values,
        y,
        s=155,
        color=point_color,
        edgecolors=point_color,
        linewidths=0.8,
        zorder=3,
    )
    for xval, ypos in zip(values, y):
        ax.annotate(
            f"{xval:.1f}%",
            (xval, ypos),
            xytext=(7, 0),
            textcoords="offset points",
            va="center",
            ha="left",
            fontsize=13,
            fontweight="bold",
            clip_on=False,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(metric_order, fontsize=13, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlim(-52, 3)
    ax.set_xticks([-50, -40, -30, -20, -10, 0])
    ax.tick_params(axis="x", labelsize=12.5, pad=5)
    for tick in ax.get_xticklabels():
        tick.set_fontweight("bold")

    ax.axvline(0, color="#6E6E6E", linewidth=1.4, zorder=1)
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.9, alpha=0.75, zorder=0)
    ax.set_xlabel(
        "Post- vs pre-transition mean change (%)",
        fontsize=14,
        fontweight="bold",
        labelpad=10,
    )
    ax.set_title(
        f"{site} ({city}; {breakpoint_label} LDSA transition)",
        fontsize=15,
        fontweight="bold",
        pad=11,
    )
    ax.text(
        -0.075,
        1.055,
        panel,
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#222222")
        spine.set_linewidth(1.35)

fig.subplots_adjust(left=0.115, right=0.985, bottom=0.17, top=0.90, wspace=0.46)

base = OUTPUT_DIR / "Fig_S9_size_mode_changes_python_redraw"
fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
plt.close(fig)

print(base.with_suffix(".png"))
print(base.with_suffix(".pdf"))
print(base.with_suffix(".svg"))
