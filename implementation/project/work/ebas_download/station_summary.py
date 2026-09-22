from pathlib import Path
import sys,json,datetime,collections,hashlib,re
import numpy as np
sys.path.insert(0,'work/ipr_download/deps')
import netCDF4
from audit_downloads import actualvar,mainvars

W=Path('work/ebas_download')
D=Path('@DELIVERY_ROOT@/EBAS_observations')
O=D.parent/'EBAS_station_inventory';O.mkdir(exist_ok=True)
stations=json.loads((W/'stations.json').read_text())
old={s['code']:s for s in json.loads((W/'local_comparison.json').read_text())}
audits=json.loads((D/'verified_metadata.json').read_text())
manifest={r['file']:r for r in json.loads((D/'download_manifest.json').read_text())}
METRICS=['PNSD','PM2.5','PM10','PM1','NO2','NO','SO2','O3','CO','BC']
TYPE={'U':'Urban','UB':'Urban background (UB)','RB':'Regional background','M':'Provider class M','F':'Provider class F','S':'Provider class S','P':'Provider class P'}
months=collections.defaultdict(set);files=[];file_metric_months={}
def typedesc(t):return TYPE[t.split(',')[0]]+(', coastal' if 'Coast' in t else '')
for ni,a in enumerate(audits,1):
 if not a['available_metrics']:continue
 code=a['station_code']
 with netCDF4.Dataset(D/a['relative_path']) as ds:
  tv=ds['time'];cal=getattr(tv,'calendar','standard');assert cal in ['standard','gregorian','proleptic_gregorian'],cal
  tb=ds['time_bnds'][:,0] if 'time_bnds' in ds.variables else tv[:]
  lo=netCDF4.date2num(datetime.datetime(1996,1,1),tv.units,calendar=cal);hi=netCDF4.date2num(datetime.datetime(2024,1,1),tv.units,calendar=cal)
  mask=np.asarray((tb>=lo)&(tb<hi));unit=tv.units.split()[0]
  factor={'days':86400,'day':86400,'hours':3600,'hour':3600,'minutes':60,'seconds':1}[unit]
  ep=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units,calendar=cal)
  ms=np.rint((np.asarray(tb)-ep)*factor*1000).astype('timedelta64[ms]')+np.datetime64('1970-01-01')
  mo=ms.astype('datetime64[M]').astype(int)
  for vi in a['variables']:
   name=vi['name'];metric=vi['metric']
   if not vi.get('finite_time_rows_in_1996_2023',0):continue
   v=ds[name];ta=v.dimensions.index('time');good=np.zeros(len(tb),bool)
   da=next((d for d in v.dimensions if d.startswith('D') and d in ds.variables),None) if metric=='PNSD' else None
   dgood=np.zeros(len(ds.dimensions[da]),bool) if da else None
   for start in range(0,len(tb),4096):
    ix=[slice(None)]*v.ndim;ix[ta]=slice(start,start+4096)
    z=np.isfinite(np.ma.filled(v[tuple(ix)],np.nan));tm=np.moveaxis(z,ta,0);n=len(tm)
    good[start:start+n]=tm.reshape(n,-1).any(axis=1)
    if da:
     dz=np.moveaxis(z,(v.dimensions.index(da),ta),(0,1))
     if dz.ndim>2:dz=dz.reshape(dz.shape[0],dz.shape[1],-1).any(axis=2)
     dgood|=dz[:,mask[start:start+n]].any(axis=1)
   good&=mask
   assert int(good.sum())==vi['finite_time_rows_in_1996_2023'],(a['file'],name)
   vm=set(map(int,np.unique(mo[good])));months[(code,metric)]|=vm
   file_metric_months.setdefault((a['file'],metric),set()).update(vm)
   if metric!='PNSD':continue
   assert da is not None,(a['file'],name)
   dp=np.asarray(ds[da][:],float);md={};mdlast={}
   mn=name+'_ebasmetadata'
   if mn in ds.variables:
    mv=ds[mn];md=json.loads(str(mv[tuple(0 for _ in mv.dimensions)]))
    lastix=[0 for _ in mv.dimensions]
    for j,dim in enumerate(mv.dimensions):
     if dim.startswith('D'):lastix[j]=-1
    mdlast=json.loads(str(mv[tuple(lastix)]))
   diameter_note=md.get('ebas_characteristic D',md.get('D',''))
   lower_note=md.get('ebas_characteristic Dmin',md.get('Dmin',''))
   upper_note=mdlast.get('ebas_characteristic Dmax',mdlast.get('Dmax',''))
   coord_unit=getattr(ds[da],'units','')
   assert 'nm' in diameter_note+lower_note+upper_note or coord_unit in ['nm','nanometer','nanometers'],(a['file'],diameter_note,coord_unit)
   if diameter_note:
    match=re.search(r'([-+0-9.eE]+)\s*nm',diameter_note)
    if match:assert abs(float(match[1])-dp[0])<max(0.001,abs(dp[0])*0.0001),(a['file'],diameter_note,dp[0])
   inst=a['metadata']['ebas_instrument_name'] or md.get('Instrument name','') or md.get('ebas_instrument_name','')
   gridkey=hashlib.sha256(np.asarray(dp,dtype='<f8').tobytes()).hexdigest()[:12]
   files.append({'code':code,'file':a['file'],'variable':name,'instrument_type':a['metadata']['ebas_instrument_type'] or a['file'].split('.')[3],'instrument':inst,'grid':gridkey,'diameters_nm':dp.tolist(),'native_low':float(dp.min()),'native_high':float(dp.max()),'published_lower_edge_nm':float(lower_note.split()[0]) if lower_note else None,'published_upper_edge_nm':float(upper_note.split()[0]) if upper_note else None,'channels':len(dp),'observed_low':float(dp[dgood].min()),'observed_high':float(dp[dgood].max()),'observed_channels':int(dgood.sum()),'start':vi['first_valid_UTC'],'end':vi['last_valid_UTC'],'months':sorted(vm),'months_count':len(vm),'valid_rows':int(good.sum()),'doi':a['metadata']['doi'],'url':manifest[a['file']]['url']})
 if ni%100==0:print('SCANNED',ni,'/',len(audits),flush=True)

config_groups=collections.defaultdict(list)
for f in files:config_groups[(f['code'],f['grid'],f['instrument_type'],f['instrument'])].append(f)
multi_months=collections.defaultdict(set)
for f in files:
 if f['channels']>=3:multi_months[f['code']].update(f['months'])
configs=[]
for key,ff in sorted(config_groups.items()):
 code,grid,it,ins=key
 configs.append({'code':code,'grid':grid,'instrument_type':it,'instrument':ins,'native_low':ff[0]['native_low'],'native_high':ff[0]['native_high'],'channels':ff[0]['channels'],'observed_low':min(r['observed_low'] for r in ff),'observed_high':max(r['observed_high'] for r in ff),'start':min(r['start'] for r in ff),'end':max(r['end'] for r in ff),'months_count':len(set().union(*(set(r['months']) for r in ff))),'file_count':len({r['file'] for r in ff})})

station_rows=[];pollutant_rows=[];type_counts=collections.Counter();metric_counts=collections.Counter();stats=[]
for s in sorted(stations,key=lambda s:s['name'].casefold()):
 code=s['code'];cc=[r for r in configs if r['code']==code and r['channels']>=3];ff=[r for r in files if r['code']==code and r['channels']>=3];fragments=[r for r in files if r['code']==code and r['channels']<3]
 assert ff,code
 ranges=sorted({(r['native_low'],r['native_high']) for r in cc});range_text='; '.join(f'{lo:g}–{hi:g}' for lo,hi in ranges)
 pub={m:[a for a in audits if a['station_code']==code and m in a['available_metrics']] for m in METRICS}
 found=[m for m in METRICS[1:] if months[(code,m)]]
 start=min(r['start'] for r in ff);end=max(r['end'] for r in ff)
 fragment_text='; '.join(sorted({f"{r['native_low']:g}" if r['channels']==1 else f"{r['native_low']:g}, {r['native_high']:g} (2 channels)" for r in fragments}))
 station_rows.append([code,s['name'],s['type'],typedesc(s['type']),range_text,start,end,len(multi_months[code]),', '.join(found) or 'Not available',len(cc),len(ranges),'; '.join(old[code]['current_analysis_sites']),old[code]['local_status'],code[:2],s['lat'],s['lon'],s['alt'],len(fragments),fragment_text])
 type_counts[s['type'].split(',')[0]]+=1
 for m in METRICS:
  mm=months[(code,m)]
  if mm:metric_counts[m]+=1
  if m=='PNSD':continue
  aa=pub[m];spans=[r for a in aa for r in a['valid_coverage'][m]]
  units=sorted({v.get('ebas_unit') or v.get('units','') for a in aa for v in a['variables'] if v['metric']==m and v.get('finite_time_rows_in_1996_2023',0)})
  pollutant_rows.append([code,s['name'],m,'Available' if mm else 'Not available',min((r['start'] for r in spans),default=None),max((r['end'] for r in spans),default=None),len(mm),len(mm & multi_months[code]),'; '.join(units),len(aa)])
 stats.append({'code':code,'name':s['name'],'type':s['type'],'ranges':ranges,'start':start,'end':end,'observed_months':len(multi_months[code]),'pollutants':found,'local_status':old[code]['local_status']})

names={s['code']:s['name'] for s in stations};seq=collections.Counter();config_rows=[]
for c in sorted(configs,key=lambda c:(names[c['code']].casefold(),c['start'],c['native_low'])):
 seq[c['code']]+=1;c['config_id']=f"{c['code']}_C{seq[c['code']]:02d}"
 config_rows.append([c['code'],names[c['code']],c['config_id'],c['instrument_type'],c['native_low'],c['native_high'],c['channels'],c['observed_low'],c['observed_high'],c['start'],c['end'],c['months_count'],c['file_count'],c['instrument'],c['grid'],'Multi-channel spectrum' if c['channels']>=3 else f"Partial spectrum: {c['channels']} channels"])
config_ids={(c['code'],c['grid'],c['instrument_type'],c['instrument']):c['config_id'] for c in configs}
file_rows=[]
for f in sorted(files,key=lambda f:(names[f['code']].casefold(),f['start'])):
 cid=config_ids[(f['code'],f['grid'],f['instrument_type'],f['instrument'])]
 file_rows.append([f['code'],cid,f['file'],f['variable'],f['native_low'],f['native_high'],f['channels'],f['observed_low'],f['observed_high'],f['observed_channels'],f['start'],f['end'],f['months_count'],f['valid_rows'],f['published_lower_edge_nm'],f['published_upper_edge_nm'],f['doi'],f['url']])
summary={'stations':len(stations),'configs':len(configs),'multi_channel_configs':sum(c['channels']>=3 for c in configs),'pnsd_files':len({f['file'] for f in files}),'multi_channel_files':len({f['file'] for f in files if f['channels']>=3}),'fragment_files':len({f['file'] for f in files if f['channels']<3}),'fragment_only_files':len({f['file'] for f in files if f['channels']<3}-{f['file'] for f in files if f['channels']>=3}),'mixed_files':len({f['file'] for f in files if f['channels']<3}&{f['file'] for f in files if f['channels']>=3}),'variable_records':len(files),'multi_channel_variables':sum(f['channels']>=3 for f in files),'fragment_variables':sum(f['channels']<3 for f in files),'type_counts':dict(type_counts),'metric_station_counts':dict(metric_counts),'global_native_low_nm':min(f['native_low'] for f in files if f['channels']>=3),'global_native_high_nm':max(f['native_high'] for f in files if f['channels']>=3),'multiple_range_sites':sum(s[10]>1 for s in station_rows),'multiple_grid_sites':sum(s[9]>1 for s in station_rows)}
payload={'stations':station_rows,'configs':config_rows,'pollutants':pollutant_rows,'files':file_rows,'summary':summary,'new_stations':[s for s in stats if s['local_status']=='PNSD_not_available_locally']}
(W/'station_summary_data.json').write_text(json.dumps(payload,ensure_ascii=False))
(O/'station_inventory_details.json').write_text(json.dumps({'summary':summary,'files':files,'configs':configs,'observed_months':{f'{k[0]}|{k[1]}':sorted(v) for k,v in months.items()},'multi_channel_PNSD_months':{k:sorted(v) for k,v in multi_months.items()}},ensure_ascii=False,indent=2))
print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
