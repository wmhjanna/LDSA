from pathlib import Path
import sys,json,hashlib
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;A=R/'outputs/network_analysis/analysis';sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import pandas as pd,numpy as np
from scipy.stats import pearsonr,spearmanr
import lightgbm as lgb
from sklearn.metrics import r2_score
tr=pd.read_csv(A/'trends.csv');t=tr[~tr.Site.isin(['PAR_SUB','PRA_SUB'])];assert not t.duplicated(['Site','Metric','Phase']).any()
phase=t[t.Phase.isin(['Before','After'])];phase.groupby(['Metric','Phase']).Slope_Percent.agg(['count','mean','median','min','max']).reset_index().to_csv(A/'phase_trend_summary.csv',index=False)
tr[tr.Metric.eq('LDSA')&tr.Phase.isin(['Before','After'])].groupby(['Type','Phase']).Slope_Percent.agg(['count','mean','median','min','max']).reset_index().to_csv(A/'network_LDSA_phase_summary.csv',index=False)
phase[phase.n_months.lt(12)].to_csv(A/'short_phase_estimates_under_12_months.csv',index=False)
paired=phase[phase.Metric.ne('LDSA')].merge(phase[phase.Metric.eq('LDSA')][['Site','Phase','Slope_Percent']],on=['Site','Phase'],suffixes=('','_LDSA')).dropna(subset=['Slope_Percent','Slope_Percent_LDSA']);paired['Same_sign']=paired.Slope_Percent*paired.Slope_Percent_LDSA>0
paired.groupby(['Metric','Phase']).Same_sign.mean().unstack().to_csv(A/'phase_direction_fraction.csv');paired.to_csv(A/'phase_direction_pairs.csv',index=False)
paired.groupby(['Metric','Phase']).Same_sign.agg(['count','sum','mean']).reset_index().to_csv(A/'S8_direction_counts.csv',index=False)
full=t[t.Phase.eq('Full')].pivot(index=['Site','Type'],columns='Metric',values='Slope_Percent').reset_index();full.to_csv(A/'full_trend_pairs.csv',index=False);cor=[];sensitivity=[]
el=pd.read_csv(A/'site_eligibility.csv');long=set(el.loc[el.meta_eligible,'Site'])
for metric in [c for c in full if c not in ['Site','Type','LDSA']]:
 for label,z in [('All sites',full),('LDSA >=60 valid months',full[full.Site.isin(long)])]:
  d=z[['Site','LDSA',metric]].replace([np.inf,-np.inf],np.nan).dropna()
  if len(d)<3 or d.LDSA.std()==0 or d[metric].std()==0:continue
  p=pearsonr(d.LDSA,d[metric]);s=spearmanr(d.LDSA,d[metric]);loo=[pearsonr(d.drop(i).LDSA,d.drop(i)[metric]).statistic for i in d.index] if len(d)>3 else [np.nan]
  sensitivity.append(dict(Metric=metric,Subset=label,n=len(d),Pearson_r=p.statistic,Pearson_p=p.pvalue,Spearman_rho=s.statistic,Spearman_p=s.pvalue,leave_one_out_r_min=min(loo),leave_one_out_r_max=max(loo)))
  if label=='All sites':cor.append(dict(Metric=metric,n=len(d),r=p.statistic,p=p.pvalue,same_sign=int((d.LDSA*d[metric]>0).sum()),opposite_sign=int((d.LDSA*d[metric]<0).sum())))
pd.DataFrame(cor).to_csv(A/'slope_correlations.csv',index=False);pd.DataFrame(sensitivity).to_csv(A/'correlation_sensitivity.csv',index=False)
checks={};c=pd.read_csv(A/'correlation_temporal_sensitivity.csv');ov=pd.read_csv(A/'exact_overlap_pairs.csv');checks['S7_correlations_independent']=[]
for row in c.itertuples():
 if row.Strategy=='Full available record':d=full[['LDSA',row.Metric]].dropna();x=d.LDSA;y=d[row.Metric]
 else:
  col='primary_included' if '60' in row.Strategy else 'sensitivity_included';d=ov[(ov.Pollutant==row.Metric)&ov[col]];x=d['LDSA slope (% yr^-1) [Y]'];y=d['Pollutant slope (% yr^-1) [X]']
 assert len(d)==row.n and np.isclose(pearsonr(x,y).statistic,row.Pearson_r)
 checks['S7_correlations_independent'].append(dict(metric=row.Metric,strategy=row.Strategy,n=row.n))
assert not (t.Site.eq('HEL_TR')&t.Metric.eq('SO2')).any()
features=['BC','CO','NO','NO2','O3','PM2.5','SO2','TEMP','RH','PRESS','RAD','WS','WD','RAIN'];models=[];perf=pd.read_csv(A/'model_performance.csv');imp=pd.read_csv(A/'SHAP_importance.csv');mi=pd.read_csv(A/'weekly_model_input.csv',parse_dates=['date']);wk=pd.read_pickle(W/'weekly.pkl')
for r in perf.itertuples():
 h=pd.read_csv(A/f'heldout_SHAP_{r.Type}_{r.Phase}.csv',float_precision='round_trip');model=lgb.Booster(model_file=str(A/f'LightGBM_{r.Type}_{r.Phase}.txt'));pred=model.predict(h[features]);sc=model.predict(h[features],pred_contrib=True)
 assert np.allclose(pred,h.Prediction,atol=1e-8) and np.allclose(sc.sum(1),pred,atol=1e-8)
 assert np.allclose(sc[:,:-1],h[[f+'_SHAP' for f in features]],atol=1e-8)
 vals=imp[(imp.Type==r.Type)&(imp.Phase==r.Phase)].set_index('Feature').reindex(features).mean_abs_SHAP
 assert np.allclose(np.abs(sc[:,:-1]).mean(0),vals) and np.isclose(r2_score(h.LDSA,pred),r.R2)
 models.append(dict(Type=r.Type,Phase=r.Phase,sites=r.sites,test=r.test,R2=r.R2,SHAP_sum_verified=True))
checks['models']=models
joined=mi.merge(wk[['Site','date','LDSA']],on=['Site','date'],suffixes=('_model','_source'),validate='many_to_one');assert len(joined)==len(mi) and np.allclose(joined.LDSA_model,joined.LDSA_source)
checks['model_LDSA_matches_current_weekly']=len(joined)
assert mi[mi.Site.eq('HEL_TR')].SO2.isna().all() and mi[mi.Site.eq('IPR_RB')].TEMP.isna().all()
rb=mi[mi.Type.eq('RB')];assert rb.RAD.dropna().ge(0).all() and rb.TEMP.dropna().lt(100).all()
checks['RB_sites']=sorted(rb.Site.unique().tolist())
e=pd.read_csv(A/'model_site_eligibility.csv');assert ((e[e.included].before_weeks>=52)&(e[e.included].after_weeks>=52)).all()
s=pd.read_pickle(W/'meteorology/selected.pkl');assert s[s.official_invalid].value.isna().all();assert not s.duplicated(['Site','metric','date']).any();checks['met_provider_invalid_never_used']=True
# S9 two city cases use the same current LDSA/size-mode data, independent of newly downloaded ancillary data.
ms=pd.read_csv(A/'monthly_series.csv');ms.date=pd.to_datetime(ms.date,format='mixed');city=pd.read_csv(A/'city_36month_changes.csv');bps=pd.read_csv(A/'breakpoints.csv');verified=0
for row in city.itertuples():
 g=ms[ms.Site.eq(row.Site)&ms.Metric.eq(row.Metric)];b=pd.Timestamp(row.Breakpoint);pre=g[g.date.ge(b-pd.DateOffset(months=36))&g.date.lt(b)].value;post=g[g.date.ge(b)&g.date.lt(b+pd.DateOffset(months=36))].value
 expected=100*(post.mean()/pre.mean()-1);assert np.isclose(expected,row.change_percent,equal_nan=True);verified+=1
checks['S9_change_rows_recomputed']=verified
(W/'validation.json').write_text(json.dumps(checks,indent=2,ensure_ascii=False));print(json.dumps(checks,indent=2,ensure_ascii=False))
