import sys
sys.path.insert(0,'@PROJECT_ROOT@/work/ldsa_support/deps')
from pathlib import Path
import pandas as pd,numpy as np,json,shutil
from scipy.stats import pearsonr,spearmanr
R=Path('@PROJECT_ROOT@');O=R/'outputs/record_exclusions';A=O/'analysis';W=R/'work/record_exclusions'
sp=pd.read_csv(A/'site_spatial_statistics.csv');el=pd.read_csv(A/'site_eligibility.csv');sp.merge(el[['Site','valid_months','meta_eligible','breakpoint']],on='Site',validate='one_to_one').to_csv(A/'Table_S1_network.csv',index=False)
sp.groupby('Type').Mean.agg(['count','mean','median','std']).reset_index().to_csv(A/'spatial_group_summary.csv',index=False)
tr=pd.read_csv(A/'trends.csv');ldsa=tr[tr.Metric.eq('LDSA')&tr.Phase.isin(['Before','After'])];ldsa.groupby(['Type','Phase']).Slope_Percent.agg(['count','mean','median','min','max']).reset_index().to_csv(A/'network_LDSA_phase_summary.csv',index=False)
# Conventional summaries explicitly use the same established comparison cohort as Fig3/S7/S8.
t=tr[~tr.Site.isin(['PAR_SUB','PRA_SUB']) & tr.Site.isin(pd.read_csv(O/'workbook/analysis/site_eligibility.csv').Site)]
for name in ['phase_trend_summary.csv','phase_direction_fraction.csv','full_trend_pairs.csv','slope_correlations.csv','correlation_sensitivity.csv']:
 shutil.copy2(O/'workbook/analysis'/name,A/('workbook_all_sites_'+name))
phase=t[t.Phase.isin(['Before','After'])];phase.groupby(['Metric','Phase']).Slope_Percent.agg(['count','mean','median','min','max']).reset_index().to_csv(A/'phase_trend_summary.csv',index=False)
paired=phase[phase.Metric.ne('LDSA')].merge(phase[phase.Metric.eq('LDSA')][['Site','Phase','Slope_Percent']],on=['Site','Phase'],suffixes=('','_LDSA')).dropna(subset=['Slope_Percent','Slope_Percent_LDSA']);paired['Same_sign']=paired.Slope_Percent*paired.Slope_Percent_LDSA>0;paired.groupby(['Metric','Phase']).Same_sign.mean().unstack().to_csv(A/'phase_direction_fraction.csv')
full=t[t.Phase.eq('Full')].pivot(index=['Site','Type'],columns='Metric',values='Slope_Percent').reset_index();full.to_csv(A/'full_trend_pairs.csv',index=False)
rows=[];sensitivity=[];long=set(el.loc[el.meta_eligible,'Site'])
for metric in [c for c in full if c not in ['Site','Type','LDSA']]:
 for label,df in [('All established comparison sites',full),('LDSA >=60 valid months',full[full.Site.isin(long)])]:
  d=df[['Site','LDSA',metric]].dropna()
  if len(d)<3 or d.LDSA.std()==0 or d[metric].std()==0:continue
  c=pearsonr(d.LDSA,d[metric]);rho=spearmanr(d.LDSA,d[metric]);loo=[pearsonr(d.drop(i).LDSA,d.drop(i)[metric]).statistic for i in d.index] if len(d)>3 else [np.nan]
  sensitivity.append(dict(Metric=metric,Subset=label,n=len(d),Pearson_r=c.statistic,Pearson_p=c.pvalue,Spearman_rho=rho.statistic,Spearman_p=rho.pvalue,leave_one_out_r_min=np.min(loo),leave_one_out_r_max=np.max(loo)))
  if label.startswith('All'):rows.append(dict(Metric=metric,n=len(d),r=c.statistic,p=c.pvalue,same_sign=int((d.LDSA*d[metric]>0).sum()),opposite_sign=int((d.LDSA*d[metric]<0).sum())))
pd.DataFrame(rows).to_csv(A/'slope_correlations.csv',index=False);pd.DataFrame(sensitivity).to_csv(A/'correlation_sensitivity.csv',index=False)
(A/'README.md').write_text('For the study protocol and reproduction instructions, see docs/METHODS.md and docs/WORKFLOW.md.\n')
print('SCOPE_LABELS_AND_SUMMARIES_VERIFIED')
