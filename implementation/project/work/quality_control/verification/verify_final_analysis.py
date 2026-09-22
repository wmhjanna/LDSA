from pathlib import Path
import json,sys
import pandas as pd,numpy as np
R=Path('@PROJECT_ROOT@');W=R/'work/quality_control';V=W/'verification';O=R/'outputs/quality_control';A=O/'analysis';OA=O/'workbook/analysis';B=Path('@DELIVERY_ROOT@');checks={};warnings=[]
assert (W/'workbook_done').exists()
def no_hel_so2(d):
 bad=d.Site.eq('HEL_TR')&d.Metric.eq('SO2')
 if 'Phase' in d:bad|=d.Site.eq('HEL_TR')&d.Phase.eq('Overlap_SO2')
 return not bad.any()
for filename in ['monthly_series.csv','full_monthly_series.csv','trends.csv']:
 d=pd.read_csv(OA/filename);assert no_hel_so2(d),filename
model=pd.read_csv(OA/'weekly_model_input.csv');assert model.loc[model.Site.eq('HEL_TR'),'SO2'].isna().all();checks['HEL_TR_SO2_excluded_trends_monthlies_model']=True
# Recompute IPR monthly values independently from quality-cleaned native inputs.
ms=pd.read_csv(OA/'monthly_series.csv',parse_dates=['date']);full=pd.read_csv(OA/'full_monthly_series.csv',parse_dates=['date']);ipr=pd.read_pickle(V/'pollutants/IPR_RB.pkl');iprchecks=[]
for m,g in ipr.groupby('metric'):
 d=g[['date','value']].dropna().copy();daily=d.groupby(d.date.dt.floor('D')).value.mean();expected=daily.resample('MS').agg(['mean','count']);expected=expected.loc[expected['count']>=14,'mean'];current=ms[ms.Site.eq('IPR_RB')&ms.Metric.eq(m)].set_index('date').value.sort_index();assert expected.index.equals(current.index),(m,'monthlyindex');assert np.allclose(expected,current,rtol=1e-12,atol=1e-9),(m,'monthly')
 d=d.drop_duplicates('date',keep='first');expectedfull=d.groupby(d.date.dt.to_period('M').dt.to_timestamp()).value.mean();currfull=full[full.Site.eq('IPR_RB')&full.Metric.eq(m)].set_index('date').value.sort_index();assert expectedfull.index.equals(currfull.index),(m,'fullindex');assert np.allclose(expectedfull,currfull,rtol=1e-12,atol=1e-9),(m,'fullmonthly')
 iprchecks.append(dict(metric=m,valid_native=int(g.value.notna().sum()),paired_months=len(current),full_months=len(currfull)))
checks['IPR_monthlies_independently_verified']=iprchecks
german=B/'german_pollutant_inputs/analysis_inputs';gc=[]
for p in sorted(german.glob('*.pkl')):
 old=pd.read_pickle(p);new=pd.read_pickle(V/'pollutants'/p.name);pd.testing.assert_frame_equal(old,new);gc.append(dict(Site=p.stem,rows=len(old),reference_rows=int(old.source.eq('other_databases').sum())))
checks['German_all17_canonical_exactly_preserved']=True;checks['German_retained_reference_rows']=sum(x['reference_rows'] for x in gc)
assert model.LDSA.notna().all();assert set(model.Type)<=set(['UB','TR']);perf=pd.read_csv(OA/'model_performance.csv');assert len(perf)==4;assert np.isfinite(perf[['R2','RMSE','MAE']]).all().all();assert (perf.RMSE>=0).all() and (perf.MAE>=0).all();checks['model_four_groups_finite_metrics']=True
# Central network checks, after the merge and pooled meta calculation.
e=pd.read_csv(A/'site_eligibility.csv');assert len(e)==65 and e.Site.nunique()==65;assert not e.Site.str.match(r'PRA_\d').any();assert e.Site.eq('PRA_SUB').sum()==1 and e.Site.eq('PAR_SUB').sum()==1;checks['network_65_unique_sites']=True
inp=pd.read_csv(W/'trend_input.csv');assert not inp.duplicated(['group_id','date']).any();assert np.isfinite(inp.value).all();checks['trend_input_unique_group_date']=True
tr=pd.read_csv(A/'trends.csv');assert not tr.duplicated(['Site','Metric','Phase']).any();assert no_hel_so2(tr)
for unit in ['Unit','Percent']:
 lo=tr['CI_Lower_'+unit];hi=tr['CI_Upper_'+unit];fi=lo.notna()&hi.notna();assert (lo[fi]<=hi[fi]).all(),unit
 pvalid=tr.P_value.dropna();assert pvalid.between(0,1).all()
checks['finite_CI_order_and_pvalue_ranges_valid']=True
meta=pd.read_csv(A/'meta_sites.csv');assert meta.n_months.ge(60).all();assert meta.SE.gt(0).all() and np.isfinite(meta.SE).all();assert meta.weight_percent.ge(0).all() and np.isclose(meta.weight_percent.sum(),100,atol=1e-8);assert not meta.Site.duplicated().any();checks['meta_sites']=len(meta);checks['meta_min_months']=int(meta.n_months.min());checks['meta_weights_sum']=float(meta.weight_percent.sum())
mom=pd.read_csv(A/'meta_monthly_series.csv');count=mom.dropna(subset=['value']).groupby('Site').size();assert np.array_equal(meta.n_months.to_numpy(),meta.Site.map(count).to_numpy());checks['meta_month_counts_equal_saved_series']=True
bp=pd.read_csv(A/'breakpoints.csv');assert bp.n_months.ge(20).all();checks['breakpoint_month_gate_verified']=True
summ=pd.read_csv(A/'meta_summary.csv');assert np.isfinite(summ[['effect','lower','upper']]).all().all();assert (summ.lower<=summ.upper).all();checks['overall_meta']=summ[summ.Type.eq('Overall')].to_dict('records')
checks['warnings']=warnings;(V/'final_analysis_verification.json').write_text(json.dumps(checks,indent=2,ensure_ascii=False));pd.DataFrame(gc).to_csv(V/'german_reference_retention_verified.csv',index=False);print(json.dumps(checks,indent=2,ensure_ascii=False))
