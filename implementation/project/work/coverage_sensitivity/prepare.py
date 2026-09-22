from pathlib import Path
import sys,json,hashlib,pickle
R=Path('@PROJECT_ROOT@');sys.path.insert(0,str(R/'work/ldsa_support/deps'))
import pandas as pd,numpy as np
W=Path(__file__).parent;O=R/'outputs/coverage_sensitivity';O.mkdir(exist_ok=True)
D=Path('@DELIVERY_ROOT@/coverage_sensitivity');D.mkdir(exist_ok=True)
reference=pd.read_csv(R/'outputs/spatial_statistics/site_spatial_statistics.csv').set_index('Site')
codes={r['new']:r['old'] for r in json.loads((R/'work/station_registry/mapping.json').read_text())}
registry_path=R/'outputs/network_analysis/netcdf/calculation_parameters/source_qc_registry.csv'
if not registry_path.exists():registry_path=R/'outputs/record_exclusions/netcdf/calculation_parameters/source_qc_registry.csv'
registry=pd.read_csv(registry_path).set_index('Source_ID')
rows=[];config=[];coverage=[];provenance=[]
for site,r in reference.iterrows():
 if site in codes:
  p=R/'outputs/network_analysis/netcdf/hourly_ldsa'/f'{codes[site]}_LDSA.csv.gz';z=pd.read_csv(p,usecols=['Time_UTC','Source_ID','Total_LDSA_um2_cm-3']);t=pd.to_datetime(z.Time_UTC);x=z['Total_LDSA_um2_cm-3'].replace([np.inf,-np.inf],np.nan)
 else:
  p=R/'work/record_exclusions/workbook/series'/f'{site}_series.pkl';z=pd.read_pickle(p);t=pd.to_datetime(z.Time);x=z.Total_LDSA_primary.replace([np.inf,-np.inf],np.nan)
 valid=t.notna()&x.notna();mask=valid&t.dt.year.between(2015,2019);period='2015-2019'
 if not mask.any():
  year=min(t[valid].dt.year.unique(),key=lambda y:(max(2015-y,y-2019,0),-y));period=str(year);mask=valid&t.dt.year.eq(year)
 v=x[mask];q=v.quantile([.25,.5,.75,.95]);assert period==str(r.Period)
 assert len(v)==r['n'] and np.isclose(v.mean(),r.Mean,rtol=1e-10,atol=1e-9)
 grids=[]
 if site in codes:
  for sid,n in z.loc[mask].Source_ID.value_counts().items():
   reg=registry.loc[sid];assert reg.code==codes[site];grids.append((reg.native_low_nm,reg.native_high_nm,int(n),f'Source_ID={sid}'))
 else:
  with (R/'work/record_exclusions/workbook/pnsd'/f'{site}.pkl').open('rb') as f:raw=pickle.load(f)
  assert len(raw['structural'])==len(z)
  masks,counts=np.unique(raw['structural'][mask],axis=0,return_counts=True)
  for m,n in zip(masks,counts):
   dp=np.asarray(raw['diameters'])[~m];assert len(dp)>0;grids.append((float(dp.min()),float(dp.max()),int(n),'native structural grid'))
  del raw
 low_bad=any(lo>17+1e-9 for lo,hi,n,s in grids);high_bad=any(hi<500-1e-9 for lo,hi,n,s in grids);keep=not(low_bad or high_bad or site=='LND_TR')
 for lo,hi,n,s in grids:config.append(dict(Site=site,code=codes.get(site,site),Period=period,source=s,native_low_nm=lo,native_high_nm=hi,valid_records=n,configuration_covers_17_500=lo<=17+1e-9 and hi>=500-1e-9))
 reason=[]
 if low_bad:reason.append('Native lower diameter exceeds 17 nm')
 if high_bad:reason.append('Native upper diameter is below 500 nm')
 if site=='LND_TR':reason.append('LND_TR integration lower bound is 31.65 nm')
 ranges='; '.join(f'{lo:g}–{hi:g}' for lo,hi in sorted(set((lo,hi) for lo,hi,n,s in grids)))
 coverage.append(dict(Site=site,Country=r.Country,code=codes.get(site,site),Type=r.Type,Period=period,native_ranges_nm=ranges,keep=keep,upper_complete=not high_bad,reason='; '.join(reason) if reason else 'Native configurations cover 17–500 nm'))
 rows.append(dict(Site=site,Country=r.Country,Type=r.Type,Period=period,n=len(v),Mean=v.mean(),Median=q.loc[.5],P25=q.loc[.25],P75=q.loc[.75],P95=q.loc[.95],actual_start=str(t[mask].min()),actual_end=str(t[mask].max())))
 provenance.append(dict(Site=site,source=str(p),resolved=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),records=len(v)))
 print(site,len(v),'strict',keep,flush=True)
d=pd.DataFrame(rows);cov=pd.DataFrame(coverage);cfg=pd.DataFrame(config)
assert len(d)==65 and cov.Site.is_unique and cfg.groupby('Site').valid_records.sum().sort_index().equals(d.set_index('Site')['n'].sort_index())
for name,table in [('site_spatial_statistics.csv',d),('coverage_selection.csv',cov),('configuration_coverage.csv',cfg)]:table.to_csv(O/name,index=False,encoding='utf-8-sig')
oldpath=Path('@DELIVERY_ROOT@/record_exclusion_analysis/diameter_coverage/coverage_selection.csv');old=pd.read_csv(oldpath).set_index('Site');assert set(cov.loc[cov.keep,'Site'])==set(old.index[old.keep])
(O/'coverage_validation.json').write_text(json.dumps(dict(current_analysis=str(R/'outputs/current_analysis.json'),reference=str(R/'outputs/spatial_statistics/site_spatial_statistics.csv'),records=provenance,native_registry=str(registry_path),native_registry_sha256=hashlib.sha256(registry_path.read_bytes()).hexdigest(),all_current_site_means_match_fig1=True,strict_membership_matches_reference=True),ensure_ascii=False,indent=2))
print('COUNTS',len(d),int(cov.upper_complete.sum()),int(cov.keep.sum()),flush=True)
