from pathlib import Path
import sys,json,pickle,datetime
R=Path('@PROJECT_ROOT@');sys.path[:0]=[str(R/'work/ipr_download/deps'),str(R/'work/ldsa_support/deps')]
import pandas as pd,numpy as np,netCDF4
W=Path(__file__).parent;B=Path('@DELIVERY_ROOT@');D=B/'EBAS_observations/source_netcdf'
a=pd.read_csv(B/'provider_quality_flags/workbook_source_flags.csv');a=a[a.scope.eq('workbook_source_inventory')&a.invalid_flag_hours.gt(0)]
flags=pd.read_csv(R/'work/provider_flag_audit/ebas_flags.csv',header=None,quotechar="'",names=['flag','description','status']);invalid=flags[flags.status.isin(['I','H','M'])].flag.to_numpy()
rows=[];masks={}
for site,group in a.groupby('Site'):
 s=pickle.load(open(R/'work/workbook_processing/pnsd'/f'{site}.pkl','rb'));dp=s['diameters'];bins=s['bins'];ts=pd.DatetimeIndex(s['aux'].Time);width=np.r_[np.diff(np.log10(dp)),np.diff(np.log10(dp))[-1]];orig=bins/width;badout=np.zeros(bins.shape,bool)
 for r in group.itertuples():
  with netCDF4.Dataset(D/r.code/'PNSD'/r.file) as ds:
   v=ds[r.variable];da=next(x for x in v.dimensions if x!='time');ndp=np.asarray(ds[da][:]);q=ds[r.variable+'_qc'];qa=np.asarray(np.ma.filled(q[:],0));qa=np.transpose(qa,[q.dimensions.index('time'),q.dimensions.index(da)]+[i for i,x in enumerate(q.dimensions) if x not in ['time',da]]);bad=np.isin(qa,invalid).any(2)&~(qa==100).any(2)
   raw=np.moveaxis(np.asarray(np.ma.filled(v[:],np.nan)),v.dimensions.index('time'),0);tv=ds['time'];ep=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units);fac={'days':86400,'hours':3600,'seconds':1,'minutes':60}[tv.units.split()[0]];tt=pd.to_datetime(np.rint((np.asarray(ds['time_bnds'][:,0])-ep)*fac).astype('int64'),unit='s').round('h')
   ix=np.abs(dp[:,None]-ndp).argmin(1);matchdp=np.isclose(dp,ndp[ix],rtol=0,atol=.001)
   for offset in [0,-30,30]:
    ot=(ts+pd.Timedelta(minutes=offset)).round('min').floor('h');mi=pd.Index(tt).get_indexer(ot);ok=mi>=0;idx=np.flatnonzero(ok);rr=raw[mi[ok]][:,ix];xx=orig[ok];finite=np.isfinite(xx)&np.isfinite(rr)&matchdp;agree=finite&np.isclose(xx,rr,rtol=1e-5,atol=.001);n=finite.sum(1);ratio=agree.sum(1)/np.maximum(n,1);exact=(n>=10)&(ratio>=.95);badmatch=bad[mi[ok]][:,ix]&matchdp&np.isfinite(xx)&exact[:,None];badout[idx]|=badmatch
    rows.append(dict(Site=site,file=r.file,offset_min=offset,matched_diameters=int(matchdp.sum()),overlap_rows=len(idx),matching_spectra=int(exact.sum()),invalid_matched_rows=int(badmatch.any(1).sum()),invalid_matched_cells=int(badmatch.sum()),max_ratio=float(ratio.max(initial=0))))
 masks[site]=badout
 print(site,'invalid matched',badout.any(1).sum(),'cells',badout.sum(),flush=True)
pd.DataFrame(rows).to_csv(W/'old_source_match.csv',index=False);pickle.dump(masks,open(W/'old_qc_masks.pkl','wb'),protocol=5)
