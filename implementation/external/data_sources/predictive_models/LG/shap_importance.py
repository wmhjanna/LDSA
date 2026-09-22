import pandas as pd
import matplotlib.pyplot as plt
import os
categories = ['UB', 'TR']
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['font.size'] = 12
plt.rcParams['mathtext.default'] = 'regular'
def format_feature_name(name):
    name = str(name)
    name = name.replace('PM2.5', 'PM$_{2.5}$')
    name = name.replace('O3', 'O$_3$')
    name = name.replace('NO2', 'NO$_2$')
    name = name.replace('SO2', 'SO$_2$')
    return name
# Load SHAP comparison tables.
dfs = {}
for cat in categories:
    excel_file = f'SHAP_Comparison_{cat}_Weekly.xlsx'

    if not os.path.exists(excel_file):
        print(f"Input file not found: {excel_file}. Skipping this station class.")
        continue

    df = pd.read_excel(excel_file)
    if 'SHAP_Before' not in df.columns or 'SHAP_After' not in df.columns:
        print(f"Required columns SHAP_Before and/or SHAP_After are missing from {excel_file}. Skipping this station class.")
        continue

    df['Feature'] = df['Feature'].apply(format_feature_name)
    dfs[cat] = df

if not dfs:
    print("No valid SHAP comparison data were loaded.")
    exit()
# Configure shared axes.
all_data = pd.concat(dfs.values())
max_x_val = all_data[['SHAP_Before', 'SHAP_After']].max().max() * 1.05
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 12), sharey=False)

for i, cat in enumerate(categories):
    df = dfs[cat].copy()
    ax_before = axes[i, 0]
    df_before = df.sort_values(by='SHAP_Before', ascending=True)

    ax_before.barh(df_before['Feature'], df_before['SHAP_Before'], color='#2a9d8f', edgecolor='black', linewidth=0.8)
    ax_before.set_title(f"{cat} stations: before the breakpoint", fontsize=16, fontweight='bold', pad=10)
    ax_before.set_xlabel("Mean |SHAP value|", fontsize=14, fontweight='bold')
    ax_before.set_ylabel("Environmental predictors", fontsize=16, fontweight='bold')
    ax_before.set_xlim(0, max_x_val)
    ax_before.grid(axis='x', linestyle='--', alpha=0.6)
    ax_before.tick_params(axis='both', which='major', labelsize=12)
    ax_after = axes[i, 1]
    df_after = df.sort_values(by='SHAP_After', ascending=True)

    ax_after.barh(df_after['Feature'], df_after['SHAP_After'], color='#e76f51', edgecolor='black', linewidth=0.8)
    ax_after.set_title(f"{cat} stations: after the breakpoint", fontsize=16, fontweight='bold', pad=10)
    ax_after.set_xlabel("Mean |SHAP value|", fontsize=14, fontweight='bold')
    ax_after.set_xlim(0, max_x_val)
    ax_after.grid(axis='x', linestyle='--', alpha=0.6)
    ax_after.yaxis.tick_right()
    ax_after.tick_params(axis='both', which='major', labelsize=12)
    ax_after.tick_params(axis='y', labelright=True, right=True)
    ax_after.set_ylabel("Environmental predictors", fontsize=16, fontweight='bold')
    ax_after.yaxis.set_label_position("right")
fig.suptitle("Importance of environmental predictors before and after the breakpoint", fontsize=20, fontweight='bold', y=1.02)
plt.tight_layout()
save_filename = 'SHAP_Absolute_Bar_2x2_Weekly_Sorted.png'
plt.savefig(save_filename, dpi=300, bbox_inches='tight')
plt.show()
plt.close()

print(f"Saved figure: {save_filename}")
