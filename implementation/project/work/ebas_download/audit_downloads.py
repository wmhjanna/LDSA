from pathlib import Path
import sys,json,time,re,datetime,traceback
sys.path.insert(0,'work/ipr_download/deps')
import netCDF4,numpy as np
W=Path('work/ebas_download');O=Path('@DELIVERY_ROOT@/EBAS_observations')
TARGET={'particle_number_size_distribution':'PNSD','pm25_mass':'PM2.5','pm10_mass':'PM10','pm1_mass':'PM1','nitrogen_dioxide':'NO2','nitrogen_monoxide':'NO','sulphur_dioxide':'SO2','ozone':'O3','carbon_monoxide':'CO','equivalent_black_carbon':'BC','black_carbon':'BC'}
def actualvar(n,v):
 if any(t in n for t in ['_qc','_ebasmetadata']):return None
 comp=getattr(v,'ebas_component','')
 for k,label in TARGET.items():
  if comp==k or n==k or n.startswith(k+'_'):return label
 return None
def mainvars(n):return not any(s in n.lower() for s in ['_max','_min','_stddev','_perc','_prec','_median','_variance'])
def audit(p):
 result={'file':p.name,'relative_path':str(p.relative_to(O)),'station_code':p.name.split('.')[0],'bytes':p.stat().st_size}
 with netCDF4.Dataset(p) as d:
  attr={k:str(d.getncattr(k)) for k in d.ncattrs()}
  result['metadata']={k:attr.get(k,'') for k in ['title','ebas_station_code','ebas_station_name','ebas_data_level','ebas_revision_date','time_coverage_start','time_coverage_end','time_coverage_resolution','ebas_resolution_code','ebas_instrument_type','ebas_instrument_name','ebas_matrix','ebas_volume_std_temperature','ebas_volume_std_pressure','doi','citation','license','ebas_export_filter','geospatial_lat_min','geospatial_lon_min']}
  assert attr.get('ebas_station_code',result['station_code'])==result['station_code']
  vtime=d['time'];tb=d['time_bnds'][:,0] if 'time_bnds' in d.variables else vtime[:]
  if 'time_bnds' in d.variables:
   end=d['time_bnds'][:,1];dt=(end-tb)*86400 if vtime.units.startswith('days') else end-tb
   result['median_interval_seconds']=float(np.ma.median(dt))
  low=netCDF4.date2num(datetime.datetime(1996,1,1),vtime.units,calendar=getattr(vtime,'calendar','standard'));high=netCDF4.date2num(datetime.datetime(2024,1,1),vtime.units,calendar=getattr(vtime,'calendar','standard'))
  mask=np.asarray((tb>=low)&(tb<high));result['rows_in_1996_2023']=int(mask.sum())
  variables=[];summary={}
  for name,v in d.variables.items():
   label=actualvar(name,v)
   if not label:continue
   vinfo={'name':name,'metric':label,'shape':list(v.shape),'dimensions':list(v.dimensions),'units':getattr(v,'units',''),'ebas_unit':getattr(v,'ebas_unit',''),'statistics':getattr(v,'ebas_statistics','')}
   if 'time' in v.dimensions and mainvars(name):
    axis=v.dimensions.index('time');good=np.zeros(len(tb),dtype=bool);complete=np.zeros(len(tb),dtype=bool)
    for st in range(0,len(tb),4096):
     ix=[slice(None)]*v.ndim;ix[axis]=slice(st,st+4096);a=v[tuple(ix)]
     valid=np.isfinite(np.ma.filled(a,np.nan));valid=np.moveaxis(valid,axis,-1).reshape((-1,valid.shape[axis]))
     good[st:st+valid.shape[1]]=valid.any(axis=0);complete[st:st+valid.shape[1]]=valid.all(axis=0)
    good &=mask;complete &=mask
    vinfo['finite_time_rows_in_1996_2023']=int(good.sum());vinfo['all_channels_finite_rows_in_1996_2023']=int(complete.sum())
    if good.any():
     dates=netCDF4.num2date(tb[good][[0,-1]],vtime.units,calendar=getattr(vtime,'calendar','standard'),only_use_cftime_datetimes=False)
     vinfo['first_valid_UTC']=str(dates[0]);vinfo['last_valid_UTC']=str(dates[1])
     summary.setdefault(label,[]).append({'name':name,'start':str(dates[0]),'end':str(dates[1]),'valid_rows':int(good.sum())})
    if label=='PNSD':
     axes=[dim for dim in v.dimensions if dim not in ['time','Location']]
     vinfo['channel_dimensions']={dim:len(d.dimensions[dim]) for dim in axes}
     daxis=next((dim for dim in axes if dim.startswith('D') and dim in d.variables),None)
     if daxis:
      ds=np.asarray(d[daxis][:]);vinfo['native_D_values']=ds.tolist();vinfo['diameter_coordinate']=daxis;vinfo['diameter_coordinate_units']=getattr(d[daxis],'units','')
      mdname=name+'_ebasmetadata'
      if mdname in d.variables:
       md=d[mdname];ix=tuple(0 for _ in md.dimensions)
       try:
        m=json.loads(str(md[ix]));vinfo['diameter_metadata_example']={k:str(val) for k,val in m.items() if 'characteristic' in k or 'diameter' in k}
       except Exception:pass
   variables.append(vinfo)
  result.update(status='verified',variables=variables,available_metrics=sorted(summary),valid_coverage=summary)
 return result
def run():
 log=O/'verified_metadata.jsonl';done={}
 if log.exists():
  for l in log.read_text().splitlines():
   a=json.loads(l);done[a['file']]=a
 while True:
  files=sorted((O/'source_netcdf').rglob('*.nc'))
  for p in files:
   if p.name in done:continue
   try:a=audit(p)
   except Exception as e:a={'file':p.name,'relative_path':str(p.relative_to(O)),'status':'audit_failed','error':str(e),'traceback':traceback.format_exc()}
   with log.open('a') as f:f.write(json.dumps(a,ensure_ascii=False)+'\n')
   done[p.name]=a;print(len(done),a['status'],p.name[:65],flush=True)
  if (O/'download_manifest.json').exists():break
  time.sleep(10)
 (O/'verified_metadata.json').write_text(json.dumps(list(done.values()),ensure_ascii=False,indent=2))
 print('AUDIT_DONE',len(done),flush=True)
if __name__=='__main__':run()
