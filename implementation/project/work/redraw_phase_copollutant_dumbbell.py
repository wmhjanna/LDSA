from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path("@PROJECT_ROOT@")
SOURCE = Path(
    "@INPUT_ROOT@/data_sources/plotting/phase_trends/"
    "phase_trend_summary.xlsx"
)
OUT_DIR = ROOT / "outputs" / "phase_copollutant_dumbbell_no_par_sub"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "font.size": 12,
        "axes.linewidth": 1.1,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
    }
)

# Read the phase-specific pollutant and LDSA slopes and apply the analysis-wide
# exclusion of PAR_SUB.
pollutants = pd.read_excel(SOURCE, sheet_name=0)
ldsa = pd.read_excel(SOURCE, sheet_name=1)
for df in (pollutants, ldsa):
    df.columns = df.columns.astype(str).str.strip()
    for col in ("Site", "Period", "Pollutant"):
        df[col] = df[col].astype(str).str.strip()
    df.drop(df.index[df["Site"].eq("PAR_SUB")], inplace=True)

median_slope = (
    pollutants.groupby(["Pollutant", "Period"])["Slope_Percent"]
    .median()
    .unstack()
)

ldsa_phase = (
    ldsa.loc[ldsa["Pollutant"].str.contains("LDSA", case=False, na=False),
             ["Site", "Period", "Slope_Percent"]]
    .rename(columns={"Slope_Percent": "Slope_LDSA"})
)
paired = pollutants.merge(ldsa_phase, on=["Site", "Period"], how="inner")
paired = paired.dropna(subset=["Slope_Percent", "Slope_LDSA"])
paired["Same_Sign"] = paired["Slope_Percent"] * paired["Slope_LDSA"] > 0
same_sign = paired.groupby(["Pollutant", "Period"])["Same_Sign"].mean().unstack()

pollutant_order = ["BC", "PM10", "PM2.5", "NO2", "NO", "SO2", "CO", "O3"]
labels = {
    "BC": "BC",
    "PM10": r"PM$_{10}$",
    "PM2.5": r"PM$_{2.5}$",
    "NO2": r"NO$_2$",
    "NO": "NO",
    "SO2": r"SO$_2$",
    "CO": "CO",
    "O3": r"O$_3$",
}

before_color = "#2c7fb8"
after_color = "#f47a2a"
connector_color = "#8f8f8f"


def draw_panel(ax, table, panel_letter, xlabel, reference):
    y = list(range(len(pollutant_order)))

    for yi, pollutant in zip(y, pollutant_order):
        if pollutant not in table.index:
            continue
        before = table.loc[pollutant, "Before"]
        after = table.loc[pollutant, "After"]
        if pd.notna(before) and pd.notna(after):
            y_before = yi - 0.10
            y_after = yi + 0.10
            ax.plot(
                [before, after], [y_before, y_after],
                color=connector_color, linewidth=1.25, zorder=1,
            )
            ax.scatter(before, y_before, s=48, color=before_color, edgecolor="none",
                       zorder=3, label="Before" if yi == 0 else None)
            ax.scatter(after, y_after, s=48, color=after_color, edgecolor="none",
                       zorder=3, label="After" if yi == 0 else None)

    ax.axvline(reference, color="#a6a6a6", linestyle="--", linewidth=1.0,
               zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[p] for p in pollutant_order],
                       fontsize=12, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel(xlabel, fontsize=13, fontweight="bold", labelpad=8)
    ax.tick_params(axis="both", labelsize=11, length=4)
    for tick in ax.get_xticklabels():
        tick.set_fontweight("bold")
    ax.grid(axis="x", color="#dedede", linewidth=0.8, alpha=0.7)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.1)
        spine.set_color("#333333")
    ax.text(-0.08, 1.015, panel_letter, transform=ax.transAxes,
            ha="left", va="bottom", fontsize=15, fontweight="bold")


fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4))
draw_panel(
    axes[0], median_slope, "a",
    r"Median co-pollutant slope (% yr$^{-1}$)", 0,
)
draw_panel(
    axes[1], same_sign, "b",
    "Fraction of sites with same-sign LDSA slope", 0.5,
)

axes[0].legend(
    loc="lower left", frameon=False, fontsize=11,
    prop={"family": "Times New Roman", "weight": "bold", "size": 11},
    handletextpad=0.6,
)
axes[1].set_xlim(-0.02, 1.02)

fig.subplots_adjust(left=0.105, right=0.985, bottom=0.16, top=0.96, wspace=0.10)

stem = "Phase_specific_copollutant_dumbbell_no_PAR_SUB_revised"
for suffix in ("png", "pdf", "svg"):
    path = OUT_DIR / f"{stem}.{suffix}"
    kwargs = {"bbox_inches": "tight", "facecolor": "white"}
    if suffix == "png":
        kwargs["dpi"] = 300
    fig.savefig(path, **kwargs)
    print(path)

plt.close(fig)
