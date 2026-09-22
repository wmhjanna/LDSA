import pandas as pd
import numpy as np
import lightgbm as lgb
import shap
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import math
import re
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
# Configure fonts and predictor labels.
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 2.0
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['mathtext.default'] = 'regular'
def format_predictor_name(name):
    """Format predictor names and units for axis labels."""
    name_str = str(name).strip()
    mapping = {
        'BC': r'BC (ng/m$^3$)',
        'PM10': r'PM$_{10}$ (µg/m$^3$)',
        'PM2.5': r'PM$_{2.5}$ (µg/m$^3$)',
        'PM1': r'PM$_1$ (µg/m$^3$)',
        'SO2': r'SO$_2$ (µg/m$^3$)',
        'NO': r'NO (µg/m$^3$)',
        'NO2': r'NO$_2$ (µg/m$^3$)',
        'O3': r'O$_3$ (µg/m$^3$)',
        'CO': r'CO (mg/m$^3$)',
        'TEMP': r'TEMP ($^\circ$C)',
        'RH': r'RH (%)',
        'PRESS': r'PRESS (mbar)',
        'RAD': r'RAD (W/m$^2$)',
        'WS': r'WS (m/s)',
        'WD': r'WD ($^\circ$)',
        'RAIN': r'RAIN (mm)'
    }
    for k, v in mapping.items():
        if name_str.upper() == k.upper():
            return v
    m = re.match(r'^([A-Za-z]+)(\d+\.?\d*)$', name_str)
    if m:
        return f"{m.group(1)}$_{{{m.group(2)}}}$"
    return name_str


def format_predictor_name_no_unit(name):
    """Format predictor names without units."""
    name_str = str(name).strip()
    mapping = {
        'PM10': r'PM$_{10}$',
        'PM2.5': r'PM$_{2.5}$',
        'PM1': r'PM$_1$',
        'SO2': r'SO$_2$',
        'NO2': r'NO$_2$',
        'NO': 'NO',
        'O3': r'O$_3$',
        'CO': 'CO',
        'BC': 'BC'
    }
    for k, v in mapping.items():
        if name_str.upper() == k.upper():
            return v
    m = re.match(r'^([A-Za-z]+)(\d+\.?\d*)$', name_str)
    if m:
        return f"{m.group(1)}$_{{{m.group(2)}}}$"
    return name_str
# Input file paths.
pollutants_file = 'weekly_observations.xlsx'
ldsa_file = 'LDSA-Total.xlsx'
breakpoints_file = '@INPUT_ROOT@/data_sources/predictive_models/LG/17-500LDSAbreakpoint.xlsx'

TARGET_POLLUTANT = 'LDSA_17_500nm'
MIN_VALID_WEEKS = 52
MIN_RECORDS_PER_WEEK = 1

predictors_to_drop = ['Accumulation', 'Aitken', 'Nucleation', 'Ntotal', 'PM1', 'PM10']
categories = ['TR', 'UB']
def parse_and_align_sheet(df):
    current_date_col = None
    aligned_series = []

    for col in df.columns:
        col_str = str(col).lower().strip()
        if col_str.startswith('date') or col_str.startswith('time'):
            current_date_col = col
        else:
            if current_date_col is not None:
                temp = df[[current_date_col, col]].copy()
                clean_col_name = str(col).split('(')[0].split('（')[0].strip()
                temp.columns = ['date', clean_col_name]

                temp['date'] = pd.to_datetime(temp['date'], errors='coerce')
                temp[clean_col_name] = pd.to_numeric(temp[clean_col_name], errors='coerce')
                temp = temp.dropna(subset=['date', clean_col_name])

                if not temp.empty:
                    temp['date'] = temp['date'].dt.to_period('W').dt.to_timestamp()
                    grouped = temp.groupby('date')[clean_col_name].agg(['mean', 'count'])
                    grouped = grouped[grouped['count'] >= MIN_RECORDS_PER_WEEK]
                    temp_result = grouped[['mean']].rename(columns={'mean': clean_col_name})
                    aligned_series.append(temp_result)

    if not aligned_series: return pd.DataFrame()
    aligned_df = pd.concat(aligned_series, axis=1, sort=False)
    return aligned_df.rename_axis('date').reset_index()
print("Loading LDSA data...")
ldsa_raw = pd.read_excel(ldsa_file)
ldsa_list = []

for i in range(0, ldsa_raw.shape[1], 2):
    temp_df = ldsa_raw.iloc[:, [i, i + 1]].copy()
    site_name = str(ldsa_raw.columns[i + 1]).strip().upper()
    temp_df.columns = ['date', 'LDSA']
    temp_df['Site'] = site_name

    temp_df['date'] = pd.to_datetime(temp_df['date'], errors='coerce')
    temp_df['LDSA'] = pd.to_numeric(temp_df['LDSA'], errors='coerce')
    temp_df = temp_df.dropna(subset=['date', 'LDSA'])

    temp_df['date'] = temp_df['date'].dt.to_period('W').dt.to_timestamp()
    grouped_ldsa = temp_df.groupby(['date', 'Site'], as_index=False).agg(LDSA=('LDSA', 'mean'),
                                                                         valid_count=('LDSA', 'count'))
    grouped_ldsa = grouped_ldsa[grouped_ldsa['valid_count'] >= MIN_RECORDS_PER_WEEK].drop(columns=['valid_count'])
    ldsa_list.append(grouped_ldsa)

ldsa_combined = pd.concat(ldsa_list)
print(f"Preparing predictor data for {TARGET_POLLUTANT}...")
bp_df = pd.read_excel(breakpoints_file)
bp_df = bp_df[bp_df['Pollutant'] == TARGET_POLLUTANT]
bp_df['Site'] = bp_df['Site'].astype(str).str.strip().str.upper()

bp_df['Start_Time_Range'] = pd.to_datetime(bp_df['Start_Time_Range'], errors='coerce').dt.to_period(
    'W').dt.to_timestamp()
bp_df['End_Time_Range'] = pd.to_datetime(bp_df['End_Time_Range'], errors='coerce').dt.to_period('W').dt.to_timestamp()
print("Matching station observations and predictor data...")
pollutants_dict = pd.read_excel(pollutants_file, sheet_name=None)

before_dataset_list, after_dataset_list, passed_sites = [], [], []

for site, df in pollutants_dict.items():
    site = str(site).strip().upper()
    site_bp = bp_df[bp_df['Site'] == site]
    if site_bp.empty: continue

    before_info = site_bp[site_bp['Period'] == 'Before']
    after_info = site_bp[site_bp['Period'] == 'After']
    if before_info.empty or after_info.empty: continue

    start_before, end_before = before_info['Start_Time_Range'].values[0], before_info['End_Time_Range'].values[0]
    start_after, end_after = after_info['Start_Time_Range'].values[0], after_info['End_Time_Range'].values[0]

    aligned_features = parse_and_align_sheet(df)
    if aligned_features.empty: continue

    site_ldsa = ldsa_combined[ldsa_combined['Site'] == site]
    if site_ldsa.empty: continue

    merged = pd.merge(site_ldsa, aligned_features, on='date', how='inner')
    merged_before = merged[(merged['date'] >= start_before) & (merged['date'] <= end_before)].dropna(subset=['LDSA'])
    merged_after = merged[(merged['date'] >= start_after) & (merged['date'] <= end_after)].dropna(subset=['LDSA'])

    if len(merged_before) >= MIN_VALID_WEEKS and len(merged_after) >= MIN_VALID_WEEKS:
        before_dataset_list.append(merged_before)
        after_dataset_list.append(merged_after)
        passed_sites.append(site)
        print(f"Skipping station with unavailable inputs: {site}")

df_before_final = pd.concat(before_dataset_list).drop(columns=predictors_to_drop, errors='ignore')
df_after_final = pd.concat(after_dataset_list).drop(columns=predictors_to_drop, errors='ignore')


def clean_inf(df):
    numeric_cols = [col for col in df.columns if col not in ['LDSA', 'Site', 'date']]
    df[numeric_cols] = df[numeric_cols].replace([np.inf, -np.inf], np.nan)
    return df


df_before_final = clean_inf(df_before_final)
df_after_final = clean_inf(df_after_final)
periods = {'Before': df_before_final, 'After': df_after_final}


def train_and_get_data(cat_data, cat_name, period_name):
    print(f"Training model for station class {cat_name}, phase {period_name}.")

    if 'BC' in cat_data.columns:
        original_len = len(cat_data)
        cat_data = cat_data[cat_data['BC'] != 0].copy()
        filtered_len = len(cat_data)
        if original_len != filtered_len:
            print(f"Excluded {original_len - filtered_len} rows with BC = 0.")

    if len(cat_data) < 10:
        return None, None

    X = cat_data.drop(columns=['LDSA', 'Site', 'date'], errors='ignore')
    y = cat_data['LDSA']

    X_train_full, X_test, y_train_full, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    X_train, X_val, y_train, y_val = train_test_split(X_train_full, y_train_full, test_size=0.25, random_state=42)

    train_data = lgb.Dataset(X_train, label=y_train)
    val_data = lgb.Dataset(X_val, label=y_val, reference=train_data)

    params = {'objective': 'regression', 'metric': 'rmse', 'boosting_type': 'gbdt',
              'learning_rate': 0.05, 'num_leaves': 31, 'max_depth': 6,
              'feature_fraction': 0.8, 'verbosity': -1, 'seed': 42}

    model = lgb.train(params, train_data, num_boost_round=1000, valid_sets=[train_data, val_data],
                      callbacks=[lgb.early_stopping(stopping_rounds=50, verbose=False)])

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)

    mean_abs_shap = np.abs(shap_values).mean(axis=0)
    shap_importance_df = pd.DataFrame({'Predictor': X_test.columns, 'Importance': mean_abs_shap}).sort_values(
        by='Importance', ascending=False)

    pdp_list = []
    dropped_outliers_count = 0

    for col in X_test.columns:
        col_idx = X_test.columns.get_loc(col)
        temp_df = pd.DataFrame({
            'Predictor': col,
            'Predictor_Value': X_test[col].values,
            'SHAP_Value': shap_values[:, col_idx],
            'Category': cat_name,
            'Period': period_name
        })

        orig_len = len(temp_df)
        lower_bound = temp_df['SHAP_Value'].quantile(0.00)
        upper_bound = temp_df['SHAP_Value'].quantile(1)
        temp_df = temp_df[(temp_df['SHAP_Value'] >= lower_bound) & (temp_df['SHAP_Value'] <= upper_bound)]
        dropped_outliers_count += (orig_len - len(temp_df))

        pdp_list.append(temp_df)

    if dropped_outliers_count > 0:
        print(f"Excluded {dropped_outliers_count} rows outside the configured SHAP quantile bounds.")

    pdp_df = pd.concat(pdp_list, ignore_index=True)
    pdp_df = pdp_df.dropna(subset=['Predictor_Value'])
    return shap_importance_df, pdp_df


all_pdp_list = []
shap_imp_dict = {cat: {} for cat in categories}

for cat in categories:
    for period, df_period in periods.items():
        cat_period_data = df_period[df_period['Site'].str.endswith(f'_{cat}', na=False)].copy()
        if len(cat_period_data) == 0: continue

        imp_df, pdp_df = train_and_get_data(cat_period_data, cat, period)
        if imp_df is not None and pdp_df is not None:
            all_pdp_list.append(pdp_df)
            shap_imp_dict[cat][period] = imp_df
if all_pdp_list:
    final_pdp_df = pd.concat(all_pdp_list, ignore_index=True)
    final_pdp_df.to_excel('SHAP_dependence_data.xlsx', index=False)
    print("Saved SHAP dependence data: SHAP_dependence_data.xlsx")
    print("Rendering SHAP dependence plots...")

    unique_preds = final_pdp_df['Predictor'].unique()
    n_preds = len(unique_preds)

    cols = 4
    rows = math.ceil((n_preds + 1) / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(22, 5 * rows))
    axes = np.array(axes).flatten()

    style_dict = {
        ('UB', 'Before'): ('#1f77b4', 'UB - Before'),
        ('UB', 'After'): ('#ff7f0e', 'UB - After'),
        ('TR', 'Before'): ('#2ca02c', 'TR - Before'),
        ('TR', 'After'): ('#d62728', 'TR - After')
    }

    for idx, pred in enumerate(unique_preds):
        ax = axes[idx]
        pred_fmt_no_unit = format_predictor_name_no_unit(pred)
        pred_fmt_with_unit = format_predictor_name(pred)

        for (cat, period), (color, label_str) in style_dict.items():
            subset = final_pdp_df[(final_pdp_df['Predictor'] == pred) &
                                  (final_pdp_df['Category'] == cat) &
                                  (final_pdp_df['Period'] == period)]
            if not subset.empty:
                sns.scatterplot(x='Predictor_Value', y='SHAP_Value', data=subset,
                                color=color, edgecolor='w', s=50, alpha=0.3, ax=ax, legend=False)
                sns.regplot(x='Predictor_Value', y='SHAP_Value', data=subset,
                            scatter=False, order=3, color=color, line_kws={'linewidth': 3.0},
                            ax=ax)

        ax.axhline(0, color='gray', linestyle='--', linewidth=1.5, alpha=0.8)

        ax.set_title("")
        ax.set_xlabel(f"{pred_fmt_with_unit}", fontsize=14, fontweight='bold')
        ax.set_ylabel('SHAP value', fontsize=14, fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.6, linewidth=1.5)
        subset_for_axis = final_pdp_df[final_pdp_df['Predictor'] == pred]
        true_y_min = subset_for_axis['SHAP_Value'].min()
        true_y_max = subset_for_axis['SHAP_Value'].max()
        val_range = true_y_max - true_y_min
        ax.set_ylim(true_y_min - val_range * 0.10, true_y_max + val_range * 0.35)

        ax.text(0.04, 0.96, pred_fmt_no_unit, transform=ax.transAxes,
                fontsize=18, fontweight='bold', va='top', ha='left')
    custom_handles = [
        mlines.Line2D([], [], color='#1f77b4', lw=4, label='Trend: UB - Before'),
        mlines.Line2D([], [], marker='o', color='w', markerfacecolor='#1f77b4', markersize=14, alpha=0.5,
                      label='Scatter: UB - Before'),
        mlines.Line2D([], [], color='#ff7f0e', lw=4, label='Trend: UB - After'),
        mlines.Line2D([], [], marker='o', color='w', markerfacecolor='#ff7f0e', markersize=14, alpha=0.5,
                      label='Scatter: UB - After'),
        mlines.Line2D([], [], color='#2ca02c', lw=4, label='Trend: TR - Before'),
        mlines.Line2D([], [], marker='o', color='w', markerfacecolor='#2ca02c', markersize=14, alpha=0.5,
                      label='Scatter: TR - Before'),
        mlines.Line2D([], [], color='#d62728', lw=4, label='Trend: TR - After'),
        mlines.Line2D([], [], marker='o', color='w', markerfacecolor='#d62728', markersize=14, alpha=0.5,
                      label='Scatter: TR - After')
    ]
    custom_labels = [h.get_label() for h in custom_handles]

    legend_ax = axes[n_preds]
    legend_ax.axis('off')

    legend_ax.legend(custom_handles, custom_labels, loc='center', ncol=2,
                     frameon=False, prop={'family': 'Times New Roman', 'weight': 'bold', 'size': 16})

    for idx in range(n_preds + 1, len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle("SHAP dependence plots for UB and TR stations", fontsize=26, fontweight='bold', y=1.02)
    plt.tight_layout(pad=3.0)
    plt.savefig('Unified_Full_PDP_Combined.png', bbox_inches='tight', dpi=300)
    plt.close()
    # Predictor importance comparison.
    print("Rendering changes in predictor importance...")
    fig_bar, axes_bar = plt.subplots(1, 2, figsize=(18, 8))

    for i, cat in enumerate(['UB', 'TR']):
        ax = axes_bar[i]
        if 'Before' in shap_imp_dict[cat] and 'After' in shap_imp_dict[cat]:
            df_b = shap_imp_dict[cat]['Before'].rename(columns={'Importance': 'SHAP_Before'})
            df_a = shap_imp_dict[cat]['After'].rename(columns={'Importance': 'SHAP_After'})

            merged_shap = pd.merge(df_b, df_a, on='Predictor', how='outer').fillna(0)
            merged_shap = merged_shap.sort_values(by='SHAP_After', ascending=False)

            merged_shap['Predictor_Display'] = merged_shap['Predictor'].apply(format_predictor_name)

            merged_shap.set_index('Predictor_Display')[['SHAP_Before', 'SHAP_After']].plot(
                kind='bar', ax=ax, color=['#1f77b4', '#ff7f0e'], width=0.75
            )

            ax.set_title(f"Change in predictor importance: {cat} stations", fontsize=16, fontweight='bold', pad=15)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontweight='bold', fontsize=13)
            ax.set_xlabel("Predictors", fontsize=14, fontweight='bold')
            ax.set_ylabel("Mean |SHAP value|", fontsize=14, fontweight='bold')
            ax.grid(axis='y', linestyle='--', alpha=0.5, linewidth=1.5)

            y_min, y_max = ax.get_ylim()
            ax.set_ylim(y_min, y_max + abs(y_max - y_min) * 0.35)

            legend = ax.legend(['Before Breakpoint', 'After Breakpoint'], loc='upper left', frameon=False,
                               prop={'family': 'Times New Roman', 'weight': 'bold', 'size': 12})

    plt.tight_layout(pad=2.0)
    plt.savefig('Unified_Predictor_Importance_Shift.png', bbox_inches='tight', dpi=300)
    plt.close()

print("\n" + "=" * 50)
print("SHAP analysis outputs saved.")
print("=" * 50)
