from pathlib import Path
import sys,json,re,itertools,datetime,collections,math,hashlib
ROOT=Path('@PROJECT_ROOT@')
sys.path[:0]=[str(ROOT/'work/ipr_download/deps'),str(ROOT/'work/ldsa_support/deps')]
import netCDF4,numpy as np,pandas as pd
W=ROOT/'work/meteorology'
BASE=Path('@DELIVERY_ROOT@');O=BASE/'netcdf_meteorology'
for sub in ['station_meteorology','source_metadata']: (O/sub).mkdir(exist_ok=True)
overview=pd.read_csv(BASE/'netcdf_ldsa_inputs/netcdf_station_summary.csv').set_index('code')
mapping={r['old']:r for r in json.loads((ROOT/'work/station_registry/mapping.json').read_text())}
stations={s['code']:s for s in json.loads((ROOT/'work/deposition/mppd/inputs/stations.json').read_text()) if s['scope']=='europe31'}
manifest=json.loads((O/'download_manifest.json').read_text());manifest=sorted(manifest,key=lambda r:r['file'])
metrics={'temperature':'TEMP','relative_humidity':'RH','pressure':'PRESS','wind_speed':'WS','wind_direction':'WD','precipitation_amount':'RAIN','downward_solar_radiation_flux_density':'RAD'}
units={'TEMP':'deg C','RH':'%','PRESS':'hPa','WS':'m/s','WD':'deg','RAIN':'mm','RAD':'W/m2'}
def get(md,key):
 canon=lambda s:re.sub('[^a-z0-9]','',s.lower()).removeprefix('ebas')
 return next((v for k,v in md.items() if canon(k)==canon(key)),'')
def num(x):
 try:return float(str(x).split()[0])
 except:return None
def dates(v,unit,cal):
 factor={'days':86400,'day':86400,'hours':3600,'hour':3600,'minutes':60,'seconds':1}[unit.split()[0]]
 ep=netCDF4.date2num(datetime.datetime(1970,1,1),unit,calendar=cal)
 return pd.to_datetime(np.rint((np.asarray(v)-ep)*factor).astype('int64'),unit='s')
def distance(a,b,c,d):
 if any(v is None for v in [a,b,c,d]):return None
 return 6371*math.hypot(math.radians(c-a),math.radians(d-b)*math.cos(math.radians((a+c)/2)))
series=[];excluded=[];summary=[];files=[];allmd={}
for code,st in stations.items():
 label=mapping[code]['new'];ob=overview.loc[code];start=pd.Timestamp(ob.full_start);end=pd.Timestamp(ob.full_end)
 parts=[]
 for fi,r in enumerate(manifest,1):
  if r['station_code']!=code or r['status']!='downloaded':continue
  filepath=O/r['relative_path'];fileused=False
  with netCDF4.Dataset(filepath) as ds:
   tv=ds['time'];cal=getattr(tv,'calendar','standard');assert cal in ['standard','gregorian','proleptic_gregorian']
   tb=np.asarray(ds['time_bnds'][:]) if 'time_bnds' in ds.variables else np.column_stack([tv[:],tv[:]])
   ts=dates(tb[:,0],tv.units,cal);te=dates(tb[:,1],tv.units,cal)
   timeok=(ts<=end)&(te>start)
   for vn,v in ds.variables.items():
    mn=vn+'_ebasmetadata'
    if 'time' not in v.dimensions or mn not in ds.variables:continue
    mv=ds[mn];extra=[(dim,len(ds.dimensions[dim])) for dim in v.dimensions if dim!='time']
    combos=list(itertools.product(*[range(n) for dim,n in extra])) or [()]
    for combo in combos:
     ix=dict(zip([dim for dim,n in extra],combo));mix=tuple(ix.get(dim,0) for dim in mv.dimensions)
     md=json.loads(str(mv[mix]));component=get(md,'Component');matrix=get(md,'Matrix');metric=metrics.get(component)
     if metric is None:continue
     if matrix!='met' and not (metric=='RAIN' and matrix=='precip'):
      excluded.append(dict(file=r['file'],variable=vn,reason='Not confirmed ambient meteorology',component=component,matrix=matrix));continue
     # Do not merge min/max/gust with the arithmetic-mean meteorological records.
     statistic=get(md,'Statistics')
     if statistic and statistic.lower() not in ['arithmetic mean','sum','integral','vector mean']:
      excluded.append(dict(file=r['file'],variable=vn,reason='Non-primary statistic',statistic=statistic));continue
     raw=np.ma.filled(v[tuple(slice(None) if dim=='time' else ix[dim] for dim in v.dimensions)],np.nan).astype(float)
     good=timeok&np.isfinite(raw)
     if not good.any():continue
     qname=vn+'_qc'
     qc=np.full(len(raw),'',dtype=object)
     if qname in ds.variables:
      qv=ds[qname];qv.set_auto_mask(False)
      qi=tuple(ix[dim] if dim in ix else slice(None) for dim in qv.dimensions)
      qa=np.asarray(qv[qi]);qd=[dim for dim in qv.dimensions if dim not in ix]
      qa=np.moveaxis(qa,qd.index('time'),0).reshape(len(raw),-1)
      qc=np.array([';'.join(map(str,row)) for row in qa],dtype=object)
     unit=get(md,'Unit') or getattr(v,'units','');assert unit==units[metric],(r['file'],metric,unit)
     sid=len(series)+1
     mdlist=[]
     for mt in range(len(ds.dimensions.get('metadata_time',[])) or 1):
      mi=tuple(mt if dim=='metadata_time' else ix.get(dim,0) for dim in mv.dimensions);mdlist.append(json.loads(str(mv[mi])))
     allmd[str(sid)]=mdlist
     dur=np.asarray((te-ts).total_seconds());lat=num(get(md,'Measurement latitude')) or num(get(md,'Station latitude'));lon=num(get(md,'Measurement longitude'));lon=lon if lon is not None else num(get(md,'Station longitude'))
     height=get(md,'Measurement height') or '; '.join(f'{dim}={np.asarray(ds[dim][:])[ix[dim]]}' for dim,n in extra if dim in ds.variables)
     desc=dict(Series_ID=sid,Source_ID=fi,code=code,Site=label,station_name=st['name'],Metric=metric,Unit=unit,Component=component,Matrix=matrix,Statistics=statistic,native_resolution=r['file'].split('.')[7],valid_records=int(good.sum()),first_UTC=str(ts[good].min()),last_UTC=str(te[good].max()),min_duration_seconds=float(dur[good].min()),max_duration_seconds=float(dur[good].max()),height=height,measurement_latitude=lat,measurement_longitude=lon,pnsd_latitude=st['lat'],pnsd_longitude=st['lon'],distance_from_pnsd_km=distance(st['lat'],st['lon'],lat,lon),data_level=get(md,'Data level'),instrument=get(md,'Instrument name'),variable=vn,extra_dimension_indices=json.dumps(ix),file=r['file'],url=r['url'],sha256=r['sha256'])
     series.append(desc);fileused=True
     parts.append(pd.DataFrame({'Time_UTC':ts[good],'End_UTC':te[good],'Metric':metric,'Value':raw[good],'Unit':unit,'QC_flags':qc[good],'Source_ID':fi,'Series_ID':sid}))
  files.append(dict(Source_ID=fi,**r,contains_requested_valid_meteorology=fileused))
 long=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame(columns=['Time_UTC','End_UTC','Metric','Value','Unit','QC_flags','Source_ID','Series_ID'])
 long=long.sort_values(['Time_UTC','Metric','Series_ID']);long.to_csv(O/'station_meteorology'/f'{label}_{code}_meteorology.csv.gz',index=False,float_format='%.10g')
 per=dict(code=code,Site=label,station_name=st['name'],latitude=st['lat'],longitude=st['lon'],PNSD_start=str(start),PNSD_end=str(end),downloaded_files=sum(r['station_code']==code and r['status']=='downloaded' for r in manifest),ambient_series=len({x['Series_ID'] for x in series if x['code']==code}))
 for metric in metrics.values():
  z=long[long.Metric.eq(metric)];q=z[['Time_UTC','End_UTC']].drop_duplicates();hourly=q[((q.End_UTC-q.Time_UTC).dt.total_seconds()<=3601)&((q.End_UTC-q.Time_UTC).dt.total_seconds()>0)] if len(q) else q
  per[metric+'_unique_intervals']=len(q);per[metric+'_hourly_or_finer_intervals']=len(hourly);per[metric+'_first']=str(z.Time_UTC.min()) if len(z) else '';per[metric+'_last']=str(z.End_UTC.max()) if len(z) else ''
 summary.append(per);print(label,'files',per['downloaded_files'],'series',per['ambient_series'],'records',len(long),flush=True)
pd.DataFrame(series).to_csv(O/'source_metadata/meteorology_missingness.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(summary).to_csv(O/'netcdf_meteorology_summary.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(files).to_csv(O/'source_metadata/meteorology_input_records.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(excluded).to_csv(O/'source_metadata/excluded_meteorology_inputs.csv',index=False,encoding='utf-8-sig')
(O/'source_metadata/meteorology_provider_metadata.json').write_text(json.dumps(allmd,ensure_ascii=False,indent=2))
(W/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
print('DONE series',len(series),'files',len(manifest),'ambient sites',sum(s['ambient_series']>0 for s in summary))
