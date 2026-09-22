"""Read archived station EBAS pollutant means, with native QC and provenance."""
from pathlib import Path
import sys,json,re,hashlib
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent
sys.path[:0]=[str(R/'work/ipr_download/deps')]
import netCDF4
import numpy as np,pandas as pd
B=Path('@DELIVERY_ROOT@/EBAS_observations')
mapping={r['old']:r['new'] for r in json.load(open(R/'work/station_registry/mapping.json'))}
target={'pm25_mass':'PM2.5','pm10_mass':'PM10','pm1_mass':'PM1','nitrogen_dioxide':'NO2','nitrogen_monoxide':'NO','sulphur_dioxide':'SO2','ozone':'O3','carbon_monoxide':'CO','equivalent_black_carbon':'BC','black_carbon':'BC'}
flags=pd.read_csv(R/'work/provider_flag_audit/ebas_flags.csv',header=None,quotechar="'",names=['flag','description','status'])
invalid=set(flags.loc[flags.status.isin(['I','H','M']),'flag']);known=set(flags.flag)
def norm(k):return k.lower().removeprefix('ebas_').replace(' ','_').replace('.','').replace('__','_')
def selection(n,v):
 out={};coords={}
 for dim in v.dimensions:
  if dim=='time':continue
  vals=np.asarray(n[dim][:],float)
  if 'wavelength' in dim.lower():i=int(np.argmin(abs(vals-880)))
  elif 'height' in dim.lower():i=int(np.nanargmin(vals))
  elif len(vals)==1:i=0
  else:raise ValueError(f'Unsupported dimension {dim}: {vals}')
  out[dim]=i;coords[dim]=dict(value=float(vals[i]),unit=getattr(n[dim],'units',''))
 return out,coords
def sliced(v,sel):
 idx=tuple(sel[d] if d in sel else slice(None) for d in v.dimensions)
 return v[idx],tuple(d for d in v.dimensions if d not in sel)
def metadata(n,k,sel):
 base={norm(a):n.getncattr(a) for a in n.ncattrs()}
 base.update({norm(a):n[k].getncattr(a) for a in n[k].ncattrs()})
 if k+'_ebasmetadata' not in n.variables:return [base],None
 arr,dims=sliced(n[k+'_ebasmetadata'],sel)
 vals=[{**base,**{norm(a):b for a,b in json.loads(str(x)).items()}} for x in np.asarray(arr).ravel()]
 if len(dims)>1:raise ValueError(('metadata dimensions',dims))
 return vals,np.asarray(n['metadata_time_bnds'][:]) if 'metadata_time_bnds' in n.variables else None
def factor(m,metric):
 u=str(m.get('unit',m.get('units',''))).replace('µ','u').replace('μ','u')
 if u=='ug N/m3':return {'NO2':46.0055,'NO':30.0061}[metric]/14.0067,u,'elemental N to molecular mass'
 if u=='ug S/m3':return 64.066/32.065,u,'elemental S to SO2 mass'
 if u in ['ug/m3','ng/m3','mg/m3']:
  return {'ug/m3':1,'ng/m3':.001,'mg/m3':1000}[u]*(1000 if metric=='BC' else .001 if metric=='CO' else 1),u,'mass unit conversion'
 if u in ['nmol/mol','umol/mol']:
  mass={'NO2':46.0055,'NO':30.0061,'SO2':64.066,'CO':28.0101,'O3':47.9982}[metric]
  ts=str(m.get('volume_std_temperature','293.15 K'));ps=str(m.get('volume_std_pressure','1013.25 hPa'))
  t=float(ts.split()[0]);p=float(ps.split()[0]);p*=100 if 'hpa' in ps.lower() or 'mbar' in ps.lower() else 1000 if 'kpa' in ps.lower() else 1
  assert 200<t<330 and 50000<p<120000,(ts,ps)
  f=mass*p/(8.314462618*t)*.001*(1000 if u=='umol/mol' else 1)*(.001 if metric=='CO' else 1)
  return f,u,f'{t} K, {p} Pa; '+('source reference' if 'volume_std_temperature' in m and 'volume_std_pressure' in m else 'reporting reference assumed, not ambient')
 raise ValueError((metric,u))
def dt(a,units):return pd.to_datetime(np.asarray(netCDF4.num2date(a,units,only_use_cftime_datetimes=False)).astype('datetime64[us]')).as_unit('ns')
def no_overlap(lower,higher):
 if lower.empty or higher.empty:return lower
 h=higher.sort_values('date');ns=lambda x:x.to_numpy(dtype='datetime64[ns]').astype('int64')
 starts=ns(h.date);ends=np.maximum.accumulate(ns(h.end));ld=ns(lower.date);le=ns(lower.end)
 i=np.searchsorted(starts,le,side='left')-1
 return lower[(i<0)|(ends[np.maximum(i,0)]<=ld)]
def main():
 W.mkdir(exist_ok=True);(W/'native').mkdir(exist_ok=True)
 parts=[];sources=[];errors=[];audit=[];excluded_levels=[]
 manifest=json.load(open(B/'download_manifest.json'))
 # Identity exclusions are independently audited against each station's PNSD location.
 identity_path=W/'audit/identity_exclusions.json'
 identity=json.load(open(identity_path)) if identity_path.exists() else []
 excluded={x['file'] for x in identity}
 for row in manifest:
  code=row['station_code']
  if code not in mapping or row['category']=='PNSD' or row['status'] not in ['downloaded','existing','ok']:continue
  path=B/row['relative_path']
  if path.name in excluded:continue
  with netCDF4.Dataset(path) as n:
   tb=np.asarray(n['time_bnds'][:]);dates=dt(tb[:,0],n['time'].units);ends=dt(tb[:,1],n['time'].units)
   timeok=(dates>=pd.Timestamp('1996-01-01'))&(dates<pd.Timestamp('2024-01-01'))
   for k,v in n.variables.items():
    if 'time' not in v.dimensions or getattr(v,'cell_methods','') not in ['time: mean','mean']:continue
    # Exact component matching excludes volatile/nonvolatile parts and uncertainty products.
    basecomp=getattr(v,'ebas_component','')
    if basecomp not in target and not any(k==c or k.startswith(c+'_') for c in target):continue
    try:
     sel,coords=selection(n,v);mm,bounds=metadata(n,k,sel)
     metric=target.get(mm[0].get('component',basecomp))
     if metric is None:continue
     if metric in ['PM2.5','PM10','PM1'] and any('volatile' in str(m.get('matrix','')).lower() for m in mm):continue
     lev_match=re.search(r'\.lev([^/]+)\.nc$',path.name)
     level=str(mm[0].get('data_level',lev_match[1] if lev_match else 'unspecified archival'))
     if level not in ['2','3','3a','unspecified archival']:
      excluded_levels.append(dict(file=path.name,variable=k,level=level,reason='not final/manual-QA data level'));continue
     if any(m.get('component',basecomp) not in target or target[m.get('component',basecomp)]!=metric for m in mm):raise ValueError('changing component')
     if any(str(m.get('station_code',code))!=code for m in mm):raise ValueError('station code mismatch')
     if metric in ['NO','NO2','SO2','O3','CO'] and mm[0].get('matrix','air')!='air':continue
     rawarr,dims=sliced(v,sel);raw=np.ma.filled(rawarr,np.nan).astype(float);assert dims==('time',),dims
     indices=np.zeros(len(raw),int)
     if bounds is not None:
      for j,(a,b) in enumerate(bounds):indices[(tb[:,0]>=a)&(tb[:,0]<b)]=j
     factors=np.ones(len(raw));rules=[]
     for j,m in enumerate(mm):
      fac,u,rule=factor(m,metric);factors[indices==j]=fac;rules.append(dict(index=j,unit=u,factor=fac,rule=rule))
     if k+'_qc' in n.variables:
      qa,qd=sliced(n[k+'_qc'],sel);qc=np.ma.filled(qa,0).astype(int);qc=np.moveaxis(qc,qd.index('time'),-1).reshape(-1,len(raw))
     else:qc=np.zeros((1,len(raw)),int)
     bad=np.isin(qc,list(invalid)).any(0)&~(qc==100).any(0)
     missing=~np.isfinite(raw)|(raw<0);value=raw*factors;value[missing|bad]=np.nan
     ix=np.flatnonzero(timeok)
     if not np.isfinite(raw[ix]).any():continue
     sid=hashlib.sha256((path.name+'|'+k+'|'+json.dumps(coords,sort_keys=True)).encode()).hexdigest()[:14]
     qctext=[';'.join(str(x) for x in qc[:,i] if x) for i in ix]
     d=pd.DataFrame(dict(code=code,Site=mapping[code],metric=metric,date=dates[ix],end=ends[ix],value=value[ix],raw_value=raw[ix],conversion_factor=factors[ix],qc=qctext,official_invalid=bad[ix],source_id=sid))
     d['unit']='ng m-3' if metric=='BC' else 'mg m-3' if metric=='CO' else 'ug m-3'
     d['level']=3 if level.startswith('3') else 2 if level=='2' else 1.9;d['revision']=path.name.split('.')[2];d['duration_s']=(d.end-d.date).dt.total_seconds()
     d['molar_priority']=int('mol/mol' in rules[0]['unit']);d['campaign_priority']=int('_IMP_' in str(mm[0].get('method_ref','')))
     parts.append(d)
     audit.append(dict(Site=mapping[code],metric=metric,source_id=sid,rows=len(ix),finite_before=int(np.isfinite(raw[ix]).sum()),invalid_finite=int((bad[ix]&np.isfinite(raw[ix])).sum()),finite_after=int(np.isfinite(value[ix]).sum()),unknown_flags=';'.join(map(str,sorted(set(qc[:,ix].ravel())-known-{0})))))
     sources.append(dict(source_id=sid,file=path.name,path=str(path),url=row['url'],sha256=row['sha256'],variable=k,metric=metric,code=code,Site=mapping[code],coordinates=coords,metadata=mm,conversions=rules))
    except Exception as ex:errors.append(dict(file=path.name,variable=k,error=str(ex)))
 (W/'sources.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2,default=str));pd.DataFrame(audit).to_csv(W/'source_qc_audit.csv',index=False)
 (W/'extraction_errors.json').write_text(json.dumps(errors,ensure_ascii=False,indent=2))
 (W/'excluded_unfinalized_sources.json').write_text(json.dumps(excluded_levels,ensure_ascii=False,indent=2))
 if errors:raise RuntimeError(f'{len(errors)} extraction errors; resolve before merging')
 raw=pd.concat(parts,ignore_index=True);selected=[]
 for (site,metric),g in raw.groupby(['Site','metric']):
  ranks=g.groupby('source_id').agg(level=('level','max'),duration=('duration_s','median'),revision=('revision','max'),molar=('molar_priority','min'),campaign=('campaign_priority','min')).sort_values(['level','duration','revision','molar','campaign'],ascending=[False,True,False,True,True],kind='stable')
  accepted=g.iloc[:0]
  for sid in ranks.index:
   d=g[g.source_id.eq(sid)].sort_values('date').drop_duplicates('date');d=no_overlap(d,accepted);accepted=pd.concat([accepted,d],ignore_index=True)
  selected.append(accepted)
 merged=pd.concat(selected,ignore_index=True).sort_values(['Site','metric','date']);merged.to_pickle(W/'selected.pkl')
 for site,g in merged.groupby('Site'):g.to_pickle(W/'native'/f'{site}.pkl');g.to_csv(W/'native'/f'{site}.csv.gz',index=False)
 summary=merged.groupby(['Site','code','metric','unit']).agg(rows=('value','size'),valid=('value','count'),start=('date','min'),end=('end','max'),sources=('source_id','nunique'),official_invalid=('official_invalid','sum')).reset_index();summary.to_csv(W/'selected_summary.csv',index=False)
 print(summary.to_string(index=False),flush=True);print('COMPLETE',len(merged),len(summary),merged.Site.nunique(),flush=True)
if __name__=='__main__':main()
