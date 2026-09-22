from pathlib import Path
import sys,json,shutil
R=Path('@PROJECT_ROOT@');sys.path.insert(0,str(R/'work/ldsa_support/deps'))
import pandas as pd,numpy as np
W=Path(__file__).parent;A=R/'outputs/pal_channel_conversion/analysis';P=R/'outputs/network_trends/analysis'
O=Path('@DELIVERY_ROOT@/PAL_native_bin_analysis');O.mkdir(exist_ok=True)
t=pd.read_csv(A/'trends.csv');old=pd.read_csv(P/'trends.csv');keys=['Site','Metric','Phase']
pd.testing.assert_frame_equal(t[t.Site.ne('PAL_RB')].set_index(keys).sort_index(),old[old.Site.ne('PAL_RB')].set_index(keys).sort_index(),check_dtype=False,rtol=1e-12,atol=1e-12)
ni=pd.read_csv(W/'trend_input.csv');oi=pd.read_csv(R/'work/network_trends/trend_input.csv')
pd.testing.assert_frame_equal(ni[ni.Site.ne('PAL_RB')].sort_values(['group_id','date']).reset_index(drop=True),oi[oi.Site.ne('PAL_RB')].sort_values(['group_id','date']).reset_index(drop=True),check_dtype=False,rtol=1e-12,atol=1e-12)
b=pd.read_csv(A/'breakpoints.csv');bo=pd.read_csv(P/'breakpoints.csv')
pd.testing.assert_frame_equal(b[b.Site.ne('PAL_RB')].reset_index(drop=True),bo[bo.Site.ne('PAL_RB')].reset_index(drop=True))
ms=pd.read_csv(A/'meta_sites.csv');prior=pd.read_csv(P/'meta_sites.csv');assert set(ms.Site)==set(prior.Site);assert np.isclose(ms.weight_percent.sum(),100)
e=pd.read_csv(A/'site_eligibility.csv');assert len(e)==65;assert (ms.n_months>=60).all();assert len(t)==len(old)
audit=json.loads((W/'conversion_audit.json').read_text())
comparison=old[old.Site.eq('PAL_RB')].merge(t[t.Site.eq('PAL_RB')],on=keys,suffixes=('_before','_deltaN'));comparison.to_csv(O/'PAL_trend_comparison.csv',index=False,encoding='utf-8-sig')
m=pd.read_csv(A/'meta_summary.csv');mo=pd.read_csv(P/'meta_summary.csv');q=pd.read_csv(A/'meta_subgroup_test.csv').iloc[0]
meta_compare=mo.merge(m,on='Type',suffixes=('_before','_deltaN'));meta_compare.to_csv(O/'meta_analysis_comparison.csv',index=False,encoding='utf-8-sig')
pal=t[t.Site.eq('PAL_RB')&t.Phase.eq('MetaFull')].iloc[0];palo=old[old.Site.eq('PAL_RB')&old.Phase.eq('MetaFull')].iloc[0]
overall=m[m.Type.eq('Overall')].iloc[0];overold=mo[mo.Type.eq('Overall')].iloc[0]
report=f'''# PAL_RB native-bin integration check

The early PAL_RB configuration is interpreted as number concentration per size bin (ΔN). LDSA is calculated as the sum of ΔN_i × π × (D_i / 1000)² × DF_i, in μm² cm⁻³. Native diameters and deposition fractions are retained, without applying an additional logarithmic bin-width factor.

Valid observations: {audit['valid_updated_rows']}. Mean LDSA under the ΔN interpretation: {audit['after_mean']:.6f} μm² cm⁻³.

Breakpoint: {audit['new_breakpoint']}.
PAL_RB trend: {pal.Slope_Percent:.6f} [{pal.CI_Lower_Percent:.6f}, {pal.CI_Upper_Percent:.6f}] % yr⁻¹.
Overall meta-analysis estimate: {overall.effect:.6f} [{overall.lower:.6f}, {overall.upper:.6f}] % yr⁻¹.
The analysis includes {len(ms)} eligible meta-analysis records. Subgroup test: Q={q.Q:.6f}, df={q.df:g}, P={q.p:.8g}.

Reference comparisons and numerical verification are provided in the accompanying CSV and JSON outputs. The source-unit interpretation is an explicit study assumption.
'''
(O/'output_notes.md').write_text(report)
for name in ['trends.csv','breakpoints.csv','meta_sites.csv','meta_summary.csv','meta_subgroup_test.csv','site_eligibility.csv']:shutil.copy2(A/name,O/name)
shutil.copy2(W/'conversion_audit.json',O/'conversion_audit.json')
verified=dict(other64_inputs_unchanged=True,other64_trends_unchanged=True,other64_breakpoints_unchanged=True,meta_site_set_unchanged=True,meta_n=len(ms),meta_weight_sum=float(ms.weight_percent.sum()),PAL=pal.to_dict(),overall=overall.to_dict(),Q=q.to_dict())
(O/'verification.json').write_text(json.dumps(verified,ensure_ascii=False,indent=2));print(json.dumps(verified,ensure_ascii=False,indent=2))
