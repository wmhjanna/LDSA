from pathlib import Path
import sys,json,datetime,gc,hashlib
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;sys.path[:0]=[str(R/'work/ipr_download/deps'),str(R/'work/ldsa_support/deps')]
import pandas as pd,numpy as np,netCDF4
A=R/'outputs/network_analysis/analysis';sp=pd.read_csv(A/'site_spatial_statistics.csv').set_index('Site');N=R/'outputs/quality_control/netcdf';reg=pd.read_csv(N/'calculation_parameters/source_qc_registry.csv').set_index('Source_ID');coef=pd.read_csv(N/'calculation_parameters/diameter_deposition_weights.csv');mapping=json.loads((R/'work/station_registry/mapping.json').read_text());B=Path('@DELIVERY_ROOT@/EBAS_observations');man={r['file']:r for r in json.loads((B/'download_manifest.json').read_text())};audit=json.loads((R/'work/network_trends/recalculate/source_audit.json').read_text());paths={x['record']['file']:B/x['audit']['relative_path'] for x in audit};flags=json.loads((R/'outputs/pnsd_cleaning/PNSD/QC_protocol.json').read_text())['sentinel_values']+[999999.99999];ft=pd.read_csv(R/'work/provider_flag_audit/ebas_flags.csv',header=None,quotechar="'",names=['flag','description','status']);invalid=ft[ft.status.isin(['I','H','M'])].flag.to_numpy();rows=[];cover=[];verify=[]
# Native spectral configurations must cover both endpoints; no extrapolation or zero-fill of unmeasured sizes.
for m in mapping:
 site=m['new'];code=m['old'];p=R/'outputs/network_analysis/netcdf/hourly_ldsa'/f'{code}_LDSA.csv.gz';h=pd.read_csv(p);h.Time_UTC=pd.to_datetime(h.Time_UTC);ys=list(map(int,str(sp.loc[site,'Period']).split('-')));h=h[h.Time_UTC.dt.year.between(ys[0],ys[-1])&h['Total_LDSA_um2_cm-3'].notna()].copy();chunks=[]
 for sid,g in h.groupby('Source_ID'):
  z=reg.loc[sid];c=coef[coef.Source_ID.eq(sid)].sort_values('channel');dp=c.diameter_nm.to_numpy();good500=dp.min()<=10+1e-8 and dp.max()>=500-1e-8;good800=dp.min()<=10+1e-8 and dp.max()>=800-1e-8
  cover.append(dict(Site=site,Source_ID=int(sid),native_low_nm=dp.min(),native_high_nm=dp.max(),n=len(g),Period=sp.loc[site,'Period'],covers_10_500=good500,covers_10_800=good800,source=str(paths[z.file]),sha256=z.sha256))
  # Calculate also uncovered ranges as available-bin diagnostics, never plotted as full-range sensitivity.
  with netCDF4.Dataset(paths[z.file]) as ds:
   tv=ds['time'];fac={'days':86400,'hours':3600,'minutes':60,'seconds':1}[tv.units.split()[0]];ep=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units,calendar=getattr(tv,'calendar','standard'));t=pd.to_datetime(np.rint((np.asarray(ds['time_bnds'][:,0])-ep)*fac).astype('int64'),unit='s').round('h');ii=pd.Index(t).get_indexer(g.Time_UTC);assert (ii>=0).all()
   v=ds[z.variable];dim=v.dimensions[0];native=np.asarray(ds[dim][:]);jj=np.abs(native[:,None]-dp).argmin(axis=0);assert np.allclose(native[jj],dp,rtol=0,atol=1e-8)
   raw=np.asarray(np.ma.filled(v[:],np.nan),float).T[ii][:,jj];q=ds[z.variable+'_qc'];qa=np.asarray(np.ma.filled(q[:],0));qa=np.transpose(qa,[q.dimensions.index('time'),q.dimensions.index(dim)]+[i for i,d in enumerate(q.dimensions) if d not in ['time',dim]])[ii][:,jj];bad=np.isin(qa,invalid).any(2)&~(qa==100).any(2);bad|=~np.isfinite(raw)|(raw<0)|np.isin(raw,flags);clean=np.where(bad,np.nan,raw)
  dw=np.ones(len(dp)) if code=='FI0096G' and sid==13 else c.delta_log10.to_numpy();coeff=dw[:,None]*np.pi*(dp[:,None]/1000)**2*c[['DF_head','DF_TB','DF_ALV','DF_total']].to_numpy();out=pd.DataFrame({'Time':g.Time_UTC.to_numpy(),'Source_ID':sid})
  for name,low,high in [('Primary',17,500),('10-500',10,500),('10-800',10,800)]:
   mask=(dp>=low)&(dp<=high);v=np.nan_to_num(clean[:,mask],nan=0)@coeff[mask];v[np.isnan(clean[:,mask]).all(1)]=np.nan
   for k,r in enumerate(['HA','TB','ALV','Total']):out[r+'_'+name]=v[:,k]
  assert np.allclose(out.Total_Primary,g['Total_LDSA_um2_cm-3'],rtol=2e-10,atol=1e-8),site
  assert np.allclose(out[['HA_Primary','TB_Primary','ALV_Primary']].sum(axis=1),out.Total_Primary)
  # Independent row-wise scalar integration for each range and each region.
  for k in np.linspace(0,len(out)-1,min(7,len(out)),dtype=int):
   for name,low,high in [('Primary',17,500),('10-500',10,500),('10-800',10,800)]:
    inds=np.flatnonzero((dp>=low)&(dp<=high)&np.isfinite(clean[k]));expected=[sum(float(clean[k,j])*float(coeff[j,cidx]) for j in inds) for cidx in range(4)] if len(inds) else [np.nan]*4
    assert np.allclose(out.loc[k,[r+'_'+name for r in ['HA','TB','ALV','Total']]].to_numpy(float),expected,equal_nan=True,rtol=1e-12,atol=1e-8)
  chunks.append(out);del raw,qa,clean;gc.collect()
 allv=pd.concat(chunks,ignore_index=True);allv=allv[allv[['Total_'+n for n in ['Primary','10-500','10-800']]].notna().all(axis=1)];assert not allv.Time.duplicated().any()
 for name in ['Primary','10-500','10-800']:
  for region in ['Total','HA','TB','ALV']:
   v=allv[region+'_'+name];rows.append(dict(Site=site,Type=m['type'],Period=sp.loc[site,'Period'],Range=name,primary_low_nm=17,Region=region,n=len(v),Mean=v.mean(),Median=v.median(),P25=v.quantile(.25),P75=v.quantile(.75),P95=v.quantile(.95)))
 verify.append(dict(Site=site,n=len(allv),source_primary_verified=True,scalar_formula_verified=True,max_abs_primary_difference=float(abs(allv.Total_Primary.mean()-h['Total_LDSA_um2_cm-3'].mean()))));print(site,len(allv),'coverage',all(x['covers_10_500'] for x in cover if x['Site']==site),all(x['covers_10_800'] for x in cover if x['Site']==site),flush=True)
 pd.DataFrame(rows).to_csv(W/'netcdf_sensitivity.csv',index=False);pd.DataFrame(cover).to_csv(W/'netcdf_coverage.csv',index=False);pd.DataFrame(verify).to_csv(W/'netcdf_verification.csv',index=False)
print('DONE',len(verify),flush=True)
