from pathlib import Path
import sys
ROOT=Path('@PROJECT_ROOT@');STAGE=ROOT/'work/ldsa_support';sys.path.insert(0,str(STAGE/'deps'))
import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr
ANA=ROOT/'outputs/ldsa_cleaned_reanalysis/analysis';rows=[]
oldroot=Path('@INPUT_ROOT@/data_sources/plotting/breakpoint_plots/clean17-500')
REV={'SCH_RB':'SCH_SUB','WAL_RB':'WAL_SUB','ZUG_RB':'ZUG_SUB'}
for p in sorted(STAGE.glob('*_series.pkl')):
    site=p.name.removesuffix('_series.pkl');s=pd.read_pickle(p);oldpath=oldroot/f'{REV.get(site,site)}_Integrated_LDSA.xlsx'
    r=dict(Site=site,new_start=s.Time.min(),new_end=s.Time.max(),new_rows=len(s),new_valid=s.LDSA.notna().sum())
    if oldpath.exists():
        old=pd.read_excel(oldpath,sheet_name='Total',engine='calamine',usecols=['Time','LDSA_17_500nm']);t=pd.to_datetime(old.Time,errors='coerce');r.update(old_start=t.min(),old_end=t.max(),old_rows=len(old),old_valid=old.LDSA_17_500nm.notna().sum())
    rows.append(r)
pd.DataFrame(rows).to_csv(ANA/'old_new_time_coverage.csv',index=False)
wide=pd.read_csv(ANA/'full_trend_pairs.csv');elig=pd.read_csv(ANA/'site_eligibility.csv');long=set(elig[elig.meta_eligible].Site);checks=[]
for metric in wide.columns[2:]:
    if metric=='LDSA':continue
    for label,source in [('All sites',wide),('Exclude PRA_SUB',wide[wide.Site.ne('PRA_SUB')]),('LDSA at least 60 valid months',wide[wide.Site.isin(long)])]:
        d=source[['Site','LDSA',metric]].dropna();x=d.LDSA;y=d[metric]
        if len(d)<3 or x.std()==0 or y.std()==0:continue
        r,p=pearsonr(x,y);rho,sp=spearmanr(x,y);loo=[]
        if len(d)>3:
            for i in range(len(d)):
                z=d.drop(d.index[i]);loo.append(pearsonr(z.LDSA,z[metric]).statistic)
        checks.append(dict(Metric=metric,Subset=label,n=len(d),Pearson_r=r,Pearson_p=p,Spearman_rho=rho,Spearman_p=sp,leave_one_out_r_min=min(loo) if loo else np.nan,leave_one_out_r_max=max(loo) if loo else np.nan))
pd.DataFrame(checks).to_csv(ANA/'correlation_sensitivity.csv',index=False)
tr=pd.read_csv(ANA/'trends.csv');tr[tr.Phase.isin(['Before','After'])&tr.n_months.lt(12)].to_csv(ANA/'short_phase_estimates_under_12_months.csv',index=False)
