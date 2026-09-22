from pathlib import Path
import sys,json,re,hashlib
R=Path('@PROJECT_ROOT@');W=R/'work/netcdf_pollutants';A=W/'audit';sys.path.insert(0,str(R/'work/ipr_download/deps'));import netCDF4,numpy as np,pandas as pd
raw=pd.read_pickle(W/'selected.pkl');ss=json.load(open(W/'sources.json'));sources={s['source_id']:s for s in ss};ff=pd.read_csv(R/'work/provider_flag_audit/ebas_flags.csv',header=None,quotechar="'",names=['code','desc','status']);invalid=set(ff.loc[ff.status.isin(['I','H','M']),'code']);exclude={e['file'] for e in json.load(open(A/'identity_exclusions.json'))};reports=[];checks={}
assert not raw.duplicated(['Site','metric','date']).any();assert raw.date.min()>=pd.Timestamp('1996-01-01') and raw.date.max()<pd.Timestamp('2024-01-01');assert (raw.end>raw.date).all();checks['unique_native_times_and_positive_duration']=True
bad=raw.qc.map(lambda q:(lambda c:bool(c&invalid) and 100 not in c)({int(x) for x in str(q).split(';') if x}));assert np.array_equal(bad,raw.official_invalid);assert raw.loc[bad,'value'].isna().all();assert raw.loc[raw.raw_value.lt(0),'value'].isna().all();finite=raw.value.notna();assert np.allclose(raw.loc[finite,'value'],raw.loc[finite,'raw_value']*raw.loc[finite,'conversion_factor'],rtol=1e-12,atol=1e-10);checks['all_official_invalid_and_negative_values_missing']=True
intervals=[]
for (site,m),g in raw.groupby(['Site','metric']):
 q=g.sort_values('date');start=q.date.to_numpy('datetime64[ns]').astype('int64');end=q.end.to_numpy('datetime64[ns]').astype('int64');overlap=int((start[1:]<np.maximum.accumulate(end)[:-1]-10**9).sum());intervals.append(dict(Site=site,metric=m,rows=len(g),valid=int(g.value.notna().sum()),interval_overlaps_over1second=overlap));assert overlap==0,(site,m,overlap)
checks['no_overlapping_native_intervals_over1second']=True
for i,(sid,g) in enumerate(raw.groupby('source_id')):
 s=sources[sid];assert s['file'] not in exclude;assert 'NRT' not in s['file'];assert not any('volatile' in str(m.get('matrix','')).lower() for m in s['metadata'])
 with netCDF4.Dataset(s['path']) as n:
  k=s['variable'];v=n[k];assert getattr(v,'cell_methods','') in ['time: mean','mean'];sel={};coords=s['coordinates']
  for dim,c in coords.items():
   vals=np.asarray(n[dim][:],float);idx=int(np.argmin(abs(vals-c['value'])));assert np.isclose(vals[idx],c['value']);sel[dim]=idx
   if 'wavelength' in dim.lower():assert idx==int(np.argmin(abs(vals-880)))
   elif 'height' in dim.lower():assert idx==int(np.argmin(vals))
  rawv=np.asarray(np.ma.filled(v[tuple(sel.get(d,slice(None)) for d in v.dimensions)],np.nan),float);tb=np.asarray(n['time_bnds'][:]);tt=pd.to_datetime(np.asarray(netCDF4.num2date(tb[:,0],n['time'].units,only_use_cftime_datetimes=False)).astype('datetime64[us]')).as_unit('ns');lookup=pd.Series(np.arange(len(tt)),index=tt);lookup=lookup[~lookup.index.duplicated(keep='first')];ix=lookup.reindex(g.date).to_numpy();assert np.isfinite(ix).all(),sid;ix=ix.astype(int);assert np.allclose(rawv[ix],g.raw_value,equal_nan=True,rtol=1e-12,atol=1e-10),(s['file'],k,'raw native mismatch')
  q=n[k+'_qc'];qd=[d for d in q.dimensions if d not in sel];qa=np.asarray(np.ma.filled(q[tuple(sel.get(d,slice(None)) for d in q.dimensions)],0));qa=np.moveaxis(qa,qd.index('time'),-1).reshape(-1,len(tt));qb=np.isin(qa,list(invalid)).any(0)&~(qa==100).any(0);assert np.array_equal(qb[ix],g.official_invalid)
  report=dict(Site=s['Site'],metric=s['metric'],source_id=sid,file=s['file'],variable=k,rows=len(g),valid=int(g.value.notna().sum()),invalid=int(g.official_invalid.sum()),native_data_exact=True);reports.append(report)
 if i%100==0:print(i,flush=True)
pd.DataFrame(reports).to_csv(A/'selected_source_verification.csv',index=False);pd.DataFrame(intervals).to_csv(A/'selected_interval_verification.csv',index=False);checks.update(selected_rows=len(raw),selected_sources=len(reports),station_count=raw.Site.nunique(),metric_pairs=raw.groupby(['Site','metric']).ngroups,finite_invalid_left=int((bad&finite).sum()),official_invalid_selected=int(bad.sum()),negative_selected=int(raw.raw_value.lt(0).sum()),valid_zero_retained=int((raw.raw_value.eq(0)&~bad).sum()),native_all_sources_independently_verified=True,MAC_conflicting_file_excluded=True,NRT_excluded=True);(A/'selected_verification.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
