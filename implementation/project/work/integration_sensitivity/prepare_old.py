from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;A=R/'outputs/network_analysis/analysis';sp=pd.read_csv(A/'site_spatial_statistics.csv').set_index('Site');old=R/'work/quality_control/workbook/series';rows=[];audit=[]
for p in sorted(old.glob('*_series.pkl')):
 site=p.name.removesuffix('_series.pkl');s=pd.read_pickle(p);t=pd.to_datetime(s.Time);yr=sp.loc[site,'Period'];years=list(map(int,str(yr).split('-')));sel=t.dt.year.between(years[0],years[-1]);names=['primary','10_500nm','10_800nm'];sel&=s[['Total_LDSA_'+n for n in names]].notna().all(axis=1)
 for label,key in zip(['Primary','10-500','10-800'],names):
  for reg in ['Total','HA','TB','ALV']:
   v=s.loc[sel,reg+'_LDSA_'+key].dropna();rows.append(dict(Site=site,Type=site.rsplit('_',1)[1],Period=yr,Range=label,primary_low_nm=31.65 if site=='LND_TR' else 17,Region=reg,n=len(v),Mean=v.mean(),Median=v.median(),P25=v.quantile(.25),P75=v.quantile(.75),P95=v.quantile(.95)))
 audit.append(dict(Site=site,source=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),n=int(sel.sum()),Period=yr,mean=s.loc[sel,'Total_LDSA_primary'].mean(),current_spatial_mean=sp.loc[site,'Mean']))
 assert np.isclose(audit[-1]['mean'],sp.loc[site,'Mean'],rtol=1e-9),site
pd.DataFrame(rows).to_csv(W/'workbook_sensitivity.csv',index=False);pd.DataFrame(audit).to_csv(W/'workbook_sources.csv',index=False)
prev=pd.read_csv(A/'size_range_sensitivity.csv');now=pd.DataFrame(rows);j=prev.merge(now,on=['Site','Region','Range'],suffixes=('_prior','_now'));delta=(j.Mean_now-j.Mean_prior).abs();print('stations',now.Site.nunique(),'max old table change',delta.max());print(j.loc[delta>1e-8,['Site','Region','Range','Mean_prior','Mean_now']].to_string(index=False))
