import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False
df = pd.read_excel("17-500LDSAbreakpoint.xlsx")
df_ldsa = df[df['Pollutant'] == 'LDSA_17_500nm'].copy()
df_ldsa['Site_Type'] = df_ldsa['Site'].apply(lambda x: x.split('_')[-1] if '_' in x else 'Unknown')
def parse_date_to_decimal(date_str):
    if pd.isna(date_str):
        return np.nan
    parts = str(date_str).split('/')
    if len(parts) == 2:
        year, month = int(parts[0]), int(parts[1])
        return year + (month - 1) / 12.0
    return np.nan
df_before = df_ldsa[df_ldsa['Period'] == 'Before'].copy()
df_after = df_ldsa[df_ldsa['Period'] == 'After'].copy()
table_c = pd.merge(df_before, df_after, on='Site', suffixes=('_pre', '_post'), how='outer')
table_c['Breakpoint_Date'] = table_c['End_Time_Range_pre']
table_c['Breakpoint_Year_Decimal'] = table_c['Breakpoint_Date'].apply(parse_date_to_decimal)
table_c['Start_Time'] = table_c['Start_Time_Range_pre'].fillna(table_c['Start_Time_Range_post'])
table_c['End_Time'] = table_c['End_Time_Range_post'].fillna(table_c['End_Time_Range_pre'])
table_c['Record_Length'] = table_c['Start_Time'].astype(str) + " to " + table_c['End_Time'].astype(str)
table_c_final = pd.DataFrame({
    'Site': table_c['Site'],
    'Type': table_c['Site_Type_pre'].fillna(table_c['Site_Type_post']),
    'Record_Length': table_c['Record_Length'],
    'Breakpoint_Date': table_c['Breakpoint_Date'],
    'Breakpoint_Year_Decimal': table_c['Breakpoint_Year_Decimal'],
    'Pre_Slope_(Unit)': table_c['Slope_Unit_pre'].round(3),
    'Post_Slope_(Unit)': table_c['Slope_Unit_post'].round(3),
    'Slope_Change_(Unit)': (table_c['Slope_Unit_post'] - table_c['Slope_Unit_pre']).round(3),
    'Pre_Slope_(%)': table_c['Slope_Percent_pre'].round(2),
    'Post_Slope_(%)': table_c['Slope_Percent_post'].round(2),
    'Slope_Change_(%)': (table_c['Slope_Percent_post'] - table_c['Slope_Percent_pre']).round(2),
    'Pre_Significance': table_c['Significance_pre'].fillna('ns'),
    'Post_Significance': table_c['Significance_post'].fillna('ns')
})
plot_data = table_c_final.dropna(subset=['Breakpoint_Year_Decimal']).sort_values('Breakpoint_Year_Decimal')
plot_data.to_csv("Supplementary_Table_C.csv", index=False)
type_palette = {"UB": "#2ca02c", "TR": "#d62728", "SUB": "#1f77b4", "RB": "#ff7f0e"}
type_order = [site_type for site_type in type_palette if site_type in plot_data['Type'].unique()]

sns.set_theme(style="whitegrid")
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
fig, (ax_a, ax_b) = plt.subplots(
    1,
    2,
    figsize=(15, 7),
    sharey=True,
    gridspec_kw={'width_ratios': [2.4, 1.6], 'wspace': 0.08}
)
ax_a.scatter(
    plot_data['Site'],
    plot_data['Breakpoint_Year_Decimal'],
    c=plot_data['Type'].map(type_palette),
    s=85,
    edgecolors='none',
    linewidths=0,
    zorder=3
)

handles = [
    plt.Line2D(
        [0], [0], marker='o', linestyle='', markerfacecolor=type_palette[site_type],
        markeredgecolor='none', markersize=9, label=site_type
    )
    for site_type in type_order
]
ax_a.legend(handles=handles, title='Station Type', loc='upper left', frameon=False)
ax_a.set_title('Station-specific breakpoint year', fontsize=15, fontweight='bold')
ax_a.text(
    -0.055, 1.025, 'a', transform=ax_a.transAxes,
    fontsize=15, fontweight='bold', ha='left', va='bottom', clip_on=False
)
ax_a.set_xlabel('Station', fontsize=13, fontweight='bold')
ax_a.set_ylabel('Breakpoint year', fontsize=13, fontweight='bold')
ax_a.tick_params(
    axis='x', which='major', labelrotation=45, labelsize=10,
    bottom=True, top=False, length=5, width=1.0,
    direction='out', color='black'
)
plt.setp(ax_a.get_xticklabels(), ha='right')
ax_a.tick_params(
    axis='y', which='major', labelsize=11,
    left=True, right=False, length=6, width=1.2,
    direction='out', color='black'
)
sns.boxplot(
    ax=ax_b,
    data=plot_data,
    x='Type',
    y='Breakpoint_Year_Decimal',
    order=type_order,
    hue='Type',
    hue_order=type_order,
    palette=type_palette,
    width=0.5,
    boxprops={'alpha': 0.45},
    showfliers=False,
    legend=False
)
sns.swarmplot(
    ax=ax_b,
    data=plot_data,
    x='Type',
    y='Breakpoint_Year_Decimal',
    order=type_order,
    hue='Type',
    hue_order=type_order,
    palette=type_palette,
    size=6.5,
    edgecolor='none',
    linewidth=0,
    legend=False
)
ax_b.set_title('Distribution by Station Type', fontsize=15, fontweight='bold')
ax_b.text(
    -0.08, 1.025, 'b', transform=ax_b.transAxes,
    fontsize=15, fontweight='bold', ha='left', va='bottom', clip_on=False
)
ax_b.set_xlabel('Station Type', fontsize=13, fontweight='bold')
ax_b.set_ylabel('')
ax_b.tick_params(
    axis='x', which='major', labelsize=11,
    bottom=True, top=False, length=6, width=1.2,
    direction='out', color='black'
)
for ax in (ax_a, ax_b):
    ax.grid(axis='y', linestyle='-', linewidth=0.8, alpha=0.65)
    ax.grid(axis='x', visible=False)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color('#555555')
        spine.set_linewidth(0.75)
    for tick_label in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        tick_label.set_fontweight('bold')

fig.subplots_adjust(left=0.065, right=0.985, bottom=0.22, top=0.90)
fig.savefig('Plot_AB_Breakpoint_Combined.png', dpi=300, bbox_inches='tight')
plt.show()
