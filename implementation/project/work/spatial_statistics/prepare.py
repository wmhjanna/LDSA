from pathlib import Path
import json,sys,hashlib
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;sys.path.insert(0,str(R/'work/ldsa_support/deps'))
import pandas as pd,numpy as np
A=R/'outputs/network_analysis/analysis';O=R/'outputs/spatial_statistics';O.mkdir(exist_ok=True);prev=pd.read_csv(A/'site_spatial_statistics.csv');codes={r['new']:r['old'] for r in json.loads((R/'work/station_registry/mapping.json').read_text())};rows=[];details=[];sources=[]
for r in prev.itertuples(index=False):
 site=r.Site
 if site in codes:
  p=R/'outputs/network_analysis/netcdf/hourly_ldsa'/f'{codes[site]}_LDSA.csv.gz';cs=[f'{x}_LDSA_um2_cm-3' for x in ['Head','TB','P','Total']];d=pd.read_csv(p,usecols=['Time_UTC']+cs);t=pd.to_datetime(d.Time_UTC)
 else:
  p=R/'work/record_exclusions/workbook/series'/f'{site}_series.pkl';cs=[f'{x}_LDSA_primary' for x in ['HA','TB','ALV','Total']];d=pd.read_pickle(p);t=pd.to_datetime(d.Time)
 x=d[cs].replace([np.inf,-np.inf],np.nan);valid=t.notna()&x.iloc[:,-1].notna();mask=valid&t.dt.year.between(2015,2019);period='2015-2019'
 if not mask.any():
  years=sorted(t[valid].dt.year.unique());year=min(years,key=lambda y:(max(2015-y,y-2019,0),-y));mask=valid&t.dt.year.eq(year);period=str(year)
 assert str(r.Period)==period,(site,r.Period,period)
 z=r._asdict();z.update(Period=period,n=int(mask.sum()),Record_Start=str(t[valid].min()),Record_End=str(t[valid].max()))
 for region,col in zip(['HA','TB','ALV','Total'],cs):
  v=x.loc[mask,col].dropna();q=v.quantile([.25,.5,.75,.95]);z[region]=float(v.mean());details.append(dict(Site=site,Type=r.Type,Period=period,Region=region,n=len(v),Mean=v.mean(),Median=q.loc[.5],P25=q.loc[.25],P75=q.loc[.75],P95=q.loc[.95]))
  if region=='Total':z.update(Mean=float(v.mean()),Median=float(q.loc[.5]),P25=float(q.loc[.25]),P75=float(q.loc[.75]),P95=float(q.loc[.95]))
 rows.append(z);sources.append(dict(Site=site,source=str(p),resolved_source=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),Period=period,actual_start=str(t[mask].min()),actual_end=str(t[mask].max()),records=int(mask.sum()),months=int(t[mask].dt.to_period('M').nunique()),source_primary_range='31.65-500' if site=='LND_TR' else '17-500'))
 print(site,period,int(mask.sum()),z['Mean'],flush=True)
now=pd.DataFrame(rows);assert len(now)==65 and now.Site.nunique()==65;now.to_csv(O/'site_spatial_statistics.csv',index=False);pd.DataFrame(details).to_csv(O/'regional_statistics.csv',index=False);pd.DataFrame(sources).to_csv(O/'spatial_provenance.csv',index=False);now.groupby('Type').Mean.agg(['count','mean','median','std']).reset_index().to_csv(O/'spatial_group_summary.csv',index=False)
a=prev.set_index('Site');b=now.set_index('Site');changes=[]
for site in a.index:
 for field in ['n','Mean','Median','P25','P75','P95','HA','TB','ALV']:
  changes.append(dict(Site=site,Field=field,previous=a.loc[site,field],current=b.loc[site,field],delta=b.loc[site,field]-a.loc[site,field],changed=not np.isclose(a.loc[site,field],b.loc[site,field],rtol=1e-10,atol=1e-9)))
c=pd.DataFrame(changes);c.to_csv(O/'numerical_check.csv',index=False);print('CHANGED',c[c.changed].to_string(index=False));print('COUNTS',now.Type.value_counts().to_dict())
