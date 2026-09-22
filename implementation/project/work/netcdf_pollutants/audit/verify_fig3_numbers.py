from pathlib import Path
import json
import pandas as pd,numpy as np
import sys
sys.path.insert(0,"@PROJECT_ROOT@/work/ldsa_support/deps")
from scipy.stats import pearsonr
R=Path('@PROJECT_ROOT@');W=R/'work/netcdf_pollutants';A=W/'audit';out=R/'outputs/netcdf_pollutants/analysis';old=R/'outputs/quality_control/analysis';pp=pd.read_csv(W/'exact_overlap_pairs.csv');prior=pd.read_csv(old/'exact_overlap_pairs.csv');cl=pd.read_csv(W/'classification.csv');off=pd.read_csv(W/'exact_overlap_breakpoints_min12.csv');ms=pd.read_csv(out/'monthly_series.csv');ms['date']=pd.to_datetime(ms.date,format='mixed');tr=pd.read_csv(out/'trends.csv')
k=['Site','Pollutant'];merged=prior.merge(pp,on=k,suffixes=('_old','_new'),validate='one_to_one');assert len(merged)==len(prior)
for c in [x for x in prior if x not in k]:
 a=merged[c+'_old'];b=merged[c+'_new'];same=np.isclose(a,b,equal_nan=True,rtol=1e-12,atol=1e-10) if pd.api.types.is_numeric_dtype(a) else ((a==b)|(a.isna()&b.isna()));assert same.all(),c
report=[]
for metric,g in pp[pp.primary_included].groupby('Pollutant'):
 assert not g.Site.isin(['PAR_SUB','PRA_SUB']).any();assert not ((g.Site=='HEL_TR')&(metric=='SO2')).any();x=g['Pollutant slope (% yr^-1) [X]'];y=g['LDSA slope (% yr^-1) [Y]'];r=pearsonr(x,y).statistic if len(g)>=3 else None;report.append(dict(Metric=metric,n=len(g),r=r))
 for z in g.itertuples(index=False):
  site=z.Site;gr=tr[tr.Site.eq(site)&tr.Phase.eq('Overlap_'+metric)].set_index('Metric')
  def state(row):
   lo,hi=row.CI_Lower_Percent,row.CI_Upper_Percent
   if not np.isfinite([lo,hi]).all():return 'Stable/uncertain'
   return 'Declining' if hi<0 else 'Increasing' if lo>0 else 'Stable/uncertain'
  row=g[g.Site.eq(site)].iloc[0];assert row.LDSA_class==state(gr.loc['LDSA']) and row.Pollutant_class==state(gr.loc[metric]);assert np.isclose(row['Pollutant slope (% yr^-1) [X]'],gr.loc[metric].Slope_Percent) and np.isclose(row['LDSA slope (% yr^-1) [Y]'],gr.loc['LDSA'].Slope_Percent)
 if metric in set(cl.Pollutant):
  cnt={'Same classification (%)':0,'LDSA decline only (%)':0,'Pollutant decline only (%)':0,'Other mismatch (%)':0}
  for row in g.itertuples(index=False):
   l,p=row.LDSA_class,row.Pollutant_class;key='Same classification (%)' if l==p else 'LDSA decline only (%)' if l=='Declining' else 'Pollutant decline only (%)' if p=='Declining' else 'Other mismatch (%)';cnt[key]+=100/len(g)
  cr=cl[cl.Pollutant.eq(metric)].iloc[0];assert cr.n==len(g);assert all(np.isclose(cr[c],v) for c,v in cnt.items());assert np.isclose(sum(cnt.values()),100)
# Independently reproduce linear-kernel single CP by max between-segment mean-square contrast.
cpchecks=[]
for row in off.to_dict('records'):
 q=ms[ms.Site.eq(row['Site'])&ms.Metric.isin(['LDSA',row['Pollutant']])].pivot(index='date',columns='Metric',values='value').dropna().sort_index();n=len(q);k=np.arange(12,n-12+1);d={}
 for metric in ['LDSA',row['Pollutant']]:
  z=q[metric].to_numpy();z=(z-z.mean())/z.std();s=np.cumsum(z);j=int(k[np.argmax(s[k-1]**2/k+(s[-1]-s[k-1])**2/(n-k))]);d[metric]=q.index[j]
 delta=(d[row['Pollutant']].year-d['LDSA'].year)*12+d[row['Pollutant']].month-d['LDSA'].month;assert delta==row['Offset Δt (months) [Y]'],(row,delta);assert d['LDSA']==pd.Timestamp(row['LDSA_breakpoint']) and d[row['Pollutant']]==pd.Timestamp(row['Metric_breakpoint']);cpchecks.append(dict(Site=row['Site'],Metric=row['Pollutant'],offset=delta))
res=dict(prior_pair_rows_unchanged=len(prior),all_fig3_xy_and_classifications_verified=True,all_e_percentages_verified=True,all_f_breakpoints_independently_reproduced=len(cpchecks),correlations=report);(A/'fig3_numerical_verification.json').write_text(json.dumps(res,indent=2));pd.DataFrame(report).to_csv(A/'fig3_independent_correlations.csv',index=False);print(json.dumps(res,indent=2))
