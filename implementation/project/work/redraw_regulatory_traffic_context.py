from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


OUTDIR = Path(
    "@PROJECT_ROOT@/outputs/regulatory_traffic_context"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "mathtext.fontset": "stix",
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "axes.linewidth": 1.15,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

fig = plt.figure(figsize=(14.2, 8.3))
gs = fig.add_gridspec(
    2,
    2,
    height_ratios=[0.78, 1.0],
    width_ratios=[1.08, 1.0],
    hspace=0.36,
    wspace=0.26,
)
ax_a = fig.add_subplot(gs[0, :])
ax_b = fig.add_subplot(gs[1, 0])
ax_c = fig.add_subplot(gs[1, 1])


def panel_label(ax, label, x=-0.055):
    ax.text(
        x,
        1.025,
        label,
        transform=ax.transAxes,
        fontsize=15,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


# a. Regulatory timeline
ax_a.set_title(
    "Regulatory timing relative to the LND_TR breakpoint",
    loc="left",
    fontsize=14.5,
    pad=9,
)
panel_label(ax_a, "a", x=-0.025)

ax_a.axhline(0, color="#5B9BD5", linewidth=1.6, zorder=1)
events = [
    (2009.00, 0.0, "Low-sulfur road fuel\nJan 2009", "#1f77b4", 0.22, "center"),
    (2009.67, 0.0, "Euro 5 type approval\nSep 2009", "#ff7f0e", -0.30, "center"),
    (2011.00, 0.0, "Euro 5 new passenger vehicles\nJan 2011", "#2ca02c", 0.22, "center"),
    (2011.67, 0.0, "Euro 5 PN limit, new types\nSep 2011", "#ff7f0e", -0.30, "center"),
    (2013.00, 0.0, "Euro 5 PN limit, all new vehicles\nJan 2013", "#9467bd", 0.22, "center"),
]
for year, y, text, color, text_y, ha in events:
    ax_a.scatter(year, y, s=42, color=color, edgecolors="white", linewidth=0.6, zorder=3)
    ax_a.text(year, text_y, text, ha=ha, va="center", fontsize=10.2, fontweight="bold")

breakpoint = 2013 + 2 / 12
ax_a.axvline(breakpoint, color="#2585D9", linestyle="--", linewidth=1.8, zorder=2)
ax_a.text(
    breakpoint + 0.06,
    0.40,
    "LND_TR LDSA breakpoint\nMar 2013",
    ha="left",
    va="center",
    fontsize=10.5,
    fontweight="bold",
)

ax_a.axvline(2016.0, color="#2C6CA3", linestyle=":", linewidth=1.6, zorder=2)
ax_a.text(
    2016.03,
    0.25,
    "2016 CCZ ANPR\nfleet snapshot",
    ha="left",
    va="center",
    fontsize=10.0,
    fontweight="bold",
)

ax_a.set_xlim(2008.65, 2016.35)
ax_a.set_ylim(-0.70, 0.70)
ax_a.set_xticks(np.arange(2009, 2017))
ax_a.set_yticks([])
ax_a.set_xlabel("Year", fontsize=12)
for side in ["left", "right", "top"]:
    ax_a.spines[side].set_visible(False)
ax_a.spines["bottom"].set_color("#444444")
ax_a.tick_params(axis="x", labelsize=10.5)


# b. Westminster major-road traffic activity
ax_b.set_title("Westminster major-road traffic activity", loc="left", fontsize=14.0, pad=9)
panel_label(ax_b, "b")
years = np.array([2010, 2013, 2016])
series = {
    "Total": (np.array([108.3, 100.0, 97.6]), "#1f77b4"),
    "Diesel car": (np.array([96.0, 100.0, 106.8]), "#ff7f0e"),
    "Diesel LGV": (np.array([103.3, 100.0, 104.3]), "#2ca02c"),
}
for name, (values, color) in series.items():
    ax_b.plot(
        years,
        values,
        color=color,
        marker="o",
        markersize=6.5,
        linewidth=2.1,
        label=name,
        zorder=3,
    )
ax_b.axhline(100, color="#4D91C6", linestyle=":", linewidth=1.3)
ax_b.axvline(2013, color="#2585D9", linestyle="--", linewidth=1.4)
ax_b.set_xlim(2009.7, 2016.3)
ax_b.set_ylim(92.0, 112.0)
ax_b.set_xticks(years)
ax_b.set_xlabel("Year", fontsize=12)
ax_b.set_ylabel("Traffic activity index (2013 = 100)", fontsize=12)
ax_b.grid(axis="y", color="#D9D9D9", linewidth=0.85, alpha=0.8)
ax_b.legend(
    loc="upper right",
    ncol=3,
    frameon=False,
    fontsize=10.2,
    columnspacing=1.5,
    handlelength=1.8,
)
ax_b.text(
    2009.82,
    92.7,
    "2013–2016: Total −2.4%; Diesel car +6.8%; Diesel LGV +4.3%",
    fontsize=9.4,
    ha="left",
    va="bottom",
    fontweight="bold",
)


# c. 2016 Central London on-road fleet composition
ax_c.set_title(
    "2016 Central London on-road fleet composition",
    loc="left",
    fontsize=14.0,
    pad=9,
)
panel_label(ax_c, "c")

categories = ["Diesel M1 cars", "Diesel N1 LGVs", "Diesel N2/N3 HGVs"]
# Percentages reproduce the proportions and the Euro 5/6 combined shares shown in the source figure.
composition = np.array(
    [
        [9.8, 27.8, 52.0, 10.4],
        [12.5, 31.8, 53.4, 2.3],
        [8.8, 23.0, 48.4, 19.8],
    ]
)
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
labels = ["Euro 3/III or older", "Euro 4/IV", "Euro 5/V", "Euro 6/VI"]
y = np.arange(len(categories))
left = np.zeros(len(categories))
for j, (color, label) in enumerate(zip(colors, labels)):
    ax_c.barh(y, composition[:, j], left=left, height=0.58, color=color, label=label)
    left += composition[:, j]

ax_c.text(98.2, 0, "62.4% Euro 5/6", ha="right", va="center", fontsize=10.0, fontweight="bold")
ax_c.text(98.2, 1, "55.7% Euro 5/6", ha="right", va="center", fontsize=10.0, fontweight="bold")
ax_c.text(98.2, 2, "68.2% Euro V/VI", ha="right", va="center", fontsize=10.0, fontweight="bold")
ax_c.set_yticks(y)
ax_c.set_yticklabels(categories, fontsize=10.5)
ax_c.invert_yaxis()
ax_c.set_xlim(0, 100)
ax_c.set_xlabel("Share of observed diesel vehicles (%)", fontsize=12)
ax_c.grid(axis="x", color="#D9D9D9", linewidth=0.8, alpha=0.75)
ax_c.set_axisbelow(True)
ax_c.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.19),
    ncol=2,
    frameon=False,
    fontsize=10.0,
    columnspacing=1.7,
    handlelength=1.8,
)

for ax in (ax_b, ax_c):
    ax.tick_params(axis="both", labelsize=10.5)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight("bold")
    for spine in ax.spines.values():
        spine.set_color("#444444")
        spine.set_linewidth(1.05)

fig.subplots_adjust(left=0.075, right=0.985, top=0.94, bottom=0.12)

base = OUTDIR / "Regulatory_timing_traffic_fleet_context"
fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
plt.close(fig)

print(base.with_suffix(".png"))
print(base.with_suffix(".pdf"))
print(base.with_suffix(".svg"))
