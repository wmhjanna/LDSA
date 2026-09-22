from pathlib import Path
import sys,pickle
ROOT=Path('@PROJECT_ROOT@');sys.path.insert(0,str(ROOT/'work/ldsa_support/deps'))
import numpy as np,pandas as pd
from scipy.stats import pearsonr,spearmanr
W=ROOT/'work/integration_support';A=ROOT/'outputs/integration_support/analysis';OLD=ROOT/'outputs/pnsd_cleaning/analysis'
tr=pd.read_csv(A/'trends.csv');ms=pd.read_csv(A/'monthly_series.csv',parse_dates=['date']);bp=pd.read_csv(A/'breakpoints.csv',parse_dates=['Breakpoint']);sp=pd.read_csv(A/'site_spatial_statistics.csv');el=pd.read_csv(A/'site_eligibility.csv')
sp.groupby('Type').Mean.agg(['count','mean','median','std']).reset_index().to_csv(A/'spatial_group_summary.csv',index=False)
phase=tr[tr.Phase.isin(['Before','After'])];phase.groupby(['Metric','Phase']).Slope_Percent.agg(['count','mean','median','min','max']).reset_index().to_csv(A/'phase_trend_summary.csv',index=False)
phase[phase.n_months.lt(12)].to_csv(A/'short_phase_estimates_under_12_months.csv',index=False)
paired=phase[phase.Metric.ne('LDSA')].merge(phase[phase.Metric.eq('LDSA')][['Site','Phase','Slope_Percent']],on=['Site','Phase'],suffixes=('','_LDSA')).dropna(subset=['Slope_Percent','Slope_Percent_LDSA']);paired['Same_sign']=paired.Slope_Percent*paired.Slope_Percent_LDSA>0
paired.groupby(['Metric','Phase']).Same_sign.mean().unstack().to_csv(A/'phase_direction_fraction.csv')
full=tr[tr.Phase.eq('Full')].pivot(index=['Site','Type'],columns='Metric',values='Slope_Percent').reset_index();full.to_csv(A/'full_trend_pairs.csv',index=False);sensitivity=[];cor=[]
long=set(el.loc[el.meta_eligible,'Site'])
for metric in [c for c in full if c not in ['Site','Type','LDSA']]:
    for label,df in [('All sites',full),('LDSA >=60 valid months',full[full.Site.isin(long)])]:
        d=df[['Site','LDSA',metric]].dropna()
        if len(d)<3 or d.LDSA.std()==0 or d[metric].std()==0:continue
        r=pearsonr(d.LDSA,d[metric]);rho=spearmanr(d.LDSA,d[metric]);loo=[pearsonr(d.drop(i).LDSA,d.drop(i)[metric]).statistic for i in d.index] if len(d)>3 else [np.nan]
        sensitivity.append(dict(Metric=metric,Subset=label,n=len(d),Pearson_r=r.statistic,Pearson_p=r.pvalue,Spearman_rho=rho.statistic,Spearman_p=rho.pvalue,leave_one_out_r_min=np.min(loo),leave_one_out_r_max=np.max(loo)))
        if label=='All sites':cor.append(dict(Metric=metric,n=len(d),r=r.statistic,p=r.pvalue,same_sign=int((d.LDSA*d[metric]>0).sum()),opposite_sign=int((d.LDSA*d[metric]<0).sum())))
pd.DataFrame(cor).to_csv(A/'slope_correlations.csv',index=False);pd.DataFrame(sensitivity).to_csv(A/'correlation_sensitivity.csv',index=False)
changes=[]
for site in ['LND_TR','STO_TR']:
    d=ms[ms.Site.eq(site)].pivot(index='date',columns='Metric',values='value').reindex(columns=sorted(ms.Metric.unique()));b=bp.loc[bp.Site.eq(site)&bp.Metric.eq('LDSA'),'Breakpoint'].iloc[0];norm=d.copy()
    for metric in d:
        pre=d.loc[(d.index>=b-pd.DateOffset(months=36))&(d.index<b),metric];post=d.loc[(d.index>=b)&(d.index<b+pd.DateOffset(months=36)),metric]
        before=d.loc[(d.index>=b-pd.DateOffset(months=24))&(d.index<b),metric].mean();norm[metric]=100*d[metric]/before if before!=0 else np.nan
        changes.append(dict(Site=site,Metric=metric,Breakpoint=b,pre_months=pre.notna().sum(),post_months=post.notna().sum(),pre_mean=pre.mean(),post_mean=post.mean(),change_percent=100*(post.mean()/pre.mean()-1) if pd.notna(pre.mean()) and pre.mean()!=0 else np.nan))
    norm.to_csv(A/(site+'_normalized_monthly.csv'))
pd.DataFrame(changes).to_csv(A/'city_36month_changes.csv',index=False)
oldbp=pd.read_csv(OLD/'breakpoints.csv');bp[bp.Metric.eq('LDSA')].merge(oldbp[oldbp.Metric.eq('LDSA')],on='Site',suffixes=('_new','_old')).to_csv(A/'window_change_breakpoints.csv',index=False)
sp.merge(pd.read_csv(OLD/'site_spatial_statistics.csv'),on='Site',suffixes=('_new','_old')).to_csv(A/'window_change_spatial.csv',index=False)
oldtr=pd.read_csv(OLD/'trends.csv');tr[tr.Metric.eq('LDSA')].merge(oldtr[oldtr.Metric.eq('LDSA')],on=['Site','Phase'],suffixes=('_new','_old')).to_csv(A/'window_change_trends.csv',index=False)
# Native configuration ranges were already verified in prepare.py.
coords=pd.read_csv(ROOT/'outputs/german_RB_revision/site_map_coordinates.csv').rename(columns={'Site':'Physical_site'})
sp[['Site','Physical_site','Type','Period','n','lower_nm','upper_nm']].merge(el[['Site','valid_months','meta_eligible','breakpoint']],on='Site').merge(coords[['Physical_site','Country','City','Latitude','Longitude']],on='Physical_site',validate='many_to_one').to_csv(A/'Table_S1_network.csv',index=False)
print('TABLES',sp.Physical_site.nunique(),'physical sites',len(sp),'independent records',len(tr),'trend fits')
