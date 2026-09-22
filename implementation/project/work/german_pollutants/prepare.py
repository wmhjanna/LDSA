from pathlib import Path
import sys,json,re,hashlib
import numpy as np,pandas as pd
R=Path('@PROJECT_ROOT@');W=R/'work/german_pollutants'
sys.path.insert(0,str(R/'work/ipr_download/deps'));import netCDF4
from fetch import O,target
METRICS=list(dict.fromkeys(target.values()));UNITS={m:('ng m-3' if m=='BC' else 'mg m-3' if m=='CO' else 'ug m-3') for m in METRICS}
def norm(k):return k.lower().removeprefix('ebas_').replace(' ','_').replace('.','').replace('__','_')
def metadata(n,k):
 v=n[k];base={norm(a):v.getncattr(a) for a in v.ncattrs()}
 if k+'_ebasmetadata' not in n.variables:return [base],None
 vals=[{**base,**{norm(a):b for a,b in json.loads(str(x)).items()}} for x in n[k+'_ebasmetadata'][:].ravel()]
 return vals, np.asarray(n['metadata_time_bnds'][:]) if 'metadata_time_bnds' in n.variables else None
def factor(m,metric):
 u=str(m.get('unit',m.get('units',''))).replace('µ','u').replace('μ','u')
 if u=='ug N/m3':return ({'NO2':46.0055,'NO':30.0061}[metric]/14.0067),u,'elemental N to molecular mass'
 if u=='ug S/m3':return 64.066/32.065,u,'elemental S to SO2 mass'
 if u in ['ug/m3','ng/m3','mg/m3']:
  toug={'ug/m3':1,'ng/m3':.001,'mg/m3':1000}[u]
  return toug*(1000 if metric=='BC' else .001 if metric=='CO' else 1),u,'mass unit conversion'
 if u in ['nmol/mol','umol/mol']:
  masses={'NO2':46.0055,'NO':30.0061,'SO2':64.066,'CO':28.0101,'O3':47.9982}
  if metric not in masses:raise ValueError((metric,u))
  # Source standard state is used when documented; otherwise explicit project reporting reference.
  t=float(str(m.get('volume_std_temperature','293.15 K')).split()[0]);p=float(str(m.get('volume_std_pressure','1013.25 hPa')).split()[0])*100
  f=masses[metric]*p/(8.314462618*t)*.001*(1000 if u=='umol/mol' else 1)*(.001 if metric=='CO' else 1)
  return f,u,f'molar ratio at {t} K, {p} Pa; '+('source reference' if 'volume_std_temperature' in m else 'reporting reference assumed, not ambient')
 raise ValueError((metric,u))
def dt(a,units):return pd.to_datetime(np.asarray(netCDF4.num2date(a,units,only_use_cftime_datetimes=False)).astype('datetime64[us]'))
def extract():
 parts=[];audit=[];sources=[]
 manifest=json.load(open(O/'download_manifest.json'));sites=json.load(open(O/'station_mapping.json'))
 for row in manifest:
  if row['status']!='ok':continue
  f=Path(row['path']);code=f.name.split('.')[0]
  with netCDF4.Dataset(f) as n:
   tb=np.asarray(n['time_bnds'][:]);dates=dt(tb[:,0],n['time'].units);ends=dt(tb[:,1],n['time'].units)
   lo=min(pd.Timestamp(s['start']).floor('D') for s in sites if s['code']==code);hi=max(pd.Timestamp(s['end']).ceil('D') for s in sites if s['code']==code)
   timeok=(dates>=lo)&(dates<hi)
   if not timeok.any():continue
   candidates=[]
   for k,v in n.variables.items():
    if v.dimensions!=('time',) or getattr(v,'cell_methods','') not in ['time: mean','mean']:continue
    mm,bounds=metadata(n,k);comp=mm[0].get('component','');metric=target.get(comp)
    if not metric:continue
    if any(str(m.get('station_code',code))!=code for m in mm):raise ValueError('station mismatch '+f.name)
    if metric in ['NO','NO2','SO2','O3','CO'] and mm[0].get('matrix','air')!='air':continue
    # Prefer provider mass-concentration variable over an equivalent molar-ratio copy.
    priority=1 if 'mol/mol' in str(mm[0].get('unit',mm[0].get('units',''))) else 0
    candidates.append((metric,priority,k,mm,bounds))
   seen=set()
   for metric,priority,k,mm,bounds in sorted(candidates,key=lambda z:(z[0],z[1],z[2])):
    if metric in seen:continue
    seen.add(metric);raw=np.ma.filled(n[k][:],np.nan).astype(float);valid=timeok&np.isfinite(raw)&(raw>=0)
    if not valid.any():continue
    indices=np.zeros(len(raw),dtype=int)
    if bounds is not None:
     for j,(a,b) in enumerate(bounds):indices[(tb[:,0]>=a)&(tb[:,0]<b)]=j
    factors=np.ones(len(raw));rules=[]
    for j,m in enumerate(mm):
     fac,u,rule=factor(m,metric);factors[indices==j]=fac;rules.append(dict(index=j,unit=u,factor=fac,rule=rule))
    sid=hashlib.sha256((f.name+'|'+k).encode()).hexdigest()[:14]
    qc=np.asarray(n[k+'_qc'][:].filled(0)) if k+'_qc' in n.variables else np.zeros((1,len(raw)),int)
    qc=np.moveaxis(qc,n[k+'_qc'].dimensions.index('time'),-1).reshape(-1,len(raw)) if k+'_qc' in n.variables else qc
    qctext=[';'.join(str(x) for x in qc[:,i] if x) for i in np.flatnonzero(valid)]
    d=pd.DataFrame(dict(code=code,metric=metric,date=dates[valid],end=ends[valid],value=raw[valid]*factors[valid],unit=UNITS[metric],raw_value=raw[valid],conversion_factor=factors[valid],qc=qctext,source_id=sid))
    d['revision']=f.name.split('.')[2];d['level']=int(mm[0].get('data_level','0'));d['duration_s']=(d.end-d.date).dt.total_seconds();parts.append(d)
    sources.append(dict(source_id=sid,file=f.name,url=row['url'],sha256=row['sha256'],variable=k,metric=metric,code=code,metadata=mm,conversions=rules,rows=len(d)))
    audit.append(dict(code=code,metric=metric,source_id=sid,rows_in_window=int(timeok.sum()),retained=len(d),negative_removed=int((timeok&(raw<0)).sum()),first=str(d.date.min()),last=str(d.date.max())))
 (O/'sources.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2,default=str));pd.DataFrame(audit).to_csv(O/'source_audit.csv',index=False)
 allraw=pd.concat(parts,ignore_index=True);allraw.to_pickle(O/'EBAS_all_candidates.pkl')
 return allraw
def no_overlap(lower,higher):
 if lower.empty or higher.empty:return lower
 h=higher.sort_values('date')
 ns=lambda x:x.to_numpy(dtype='datetime64[ns]').astype('int64')
 starts=ns(h.date);ends=np.maximum.accumulate(ns(h.end));ld=ns(lower.date);le=ns(lower.end)
 i=np.searchsorted(starts,le,side='left')-1
 return lower[(i<0)|(ends[np.maximum(i,0)]<=ld)]
def merge_sources(raw):
 chunks=[]
 for (code,metric),g in raw.groupby(['code','metric']):
  accepted=g.iloc[:0];ranks=g.groupby('source_id').agg(level=('level','max'),duration_s=('duration_s','median'),revision=('revision','max')).sort_values(['level','duration_s','revision'],ascending=[False,True,False])
  for sid in ranks.index:
   d=g[g.source_id.eq(sid)].sort_values('date').drop_duplicates('date');d=no_overlap(d,accepted);accepted=pd.concat([accepted,d],ignore_index=True)
  chunks.append(accepted)
 return pd.concat(chunks,ignore_index=True)
if __name__=='__main__':
 merged=merge_sources(extract());merged.to_pickle(O/'EBAS_selected.pkl')
 print('EXTRACTED',len(merged),'observations',merged.groupby(['code','metric']).size().to_string(),flush=True)
