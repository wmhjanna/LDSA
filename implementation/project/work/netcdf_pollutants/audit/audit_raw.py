from pathlib import Path
import sys,json,re,itertools,math,datetime,collections
R=Path('@PROJECT_ROOT@');sys.path.insert(0,str(R/'work/ipr_download/deps'));import netCDF4,numpy as np,pandas as pd
A=R/'work/netcdf_pollutants/audit';B=Path('@DELIVERY_ROOT@/EBAS_observations');M={r['old']:r for r in json.load(open(R/'work/station_registry/mapping.json'))};vv=[r for r in json.load(open(B/'verified_metadata.json')) if r['station_code'] in M and '/PNSD/' not in r['relative_path']];ff=pd.read_csv(R/'work/provider_flag_audit/ebas_flags.csv',header=None,quotechar="'",names=['code','description','status']);invalid=ff[ff.status.isin(['I','H','M'])].code.to_numpy();known=set(ff.code)|{0}
def norm(k):return re.sub('[^a-z0-9]','',k.lower()).removeprefix('ebas')
def normalized(d):return {norm(k):str(v) for k,v in d.items()}
def coords(d):
 def num(x):
  try:return float(x)
  except:return None
 return (num(d.get('measurementlatitude',d.get('stationlatitude'))),num(d.get('measurementlongitude',d.get('stationlongitude'))))
def distance(a,b):
 if None in a or None in b:return None
 p1,p2=map(math.radians,[a[0],b[0]]);dl=math.radians(a[1]-b[1]);z=math.sin((p1-p2)/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2;return 6371000*2*math.asin(min(1,math.sqrt(z)))
ref={c:[] for c in M};details=json.load(open(R/'work/network_trends/recalculate/metadata_details.json'));reg=pd.read_csv(R/'outputs/quality_control/netcdf/calculation_parameters/source_qc_registry.csv');selected=set(reg[reg.selected_rows.gt(0)].file)
for r in details:
 if r['code'] in M and r['file'] in selected:
  d=normalized({**r.get('varatts',{}),**r['metadata']});co=coords(d)
  if co not in [x['coords'] for x in ref[r['code']]]:ref[r['code']].append(dict(coords=co,name=d.get('stationname',''),platform=d.get('platformcode',''),height=d.get('measurementheight',''),file=r['file']))
rows=[];locations=[];errors=[];dups=[];coverage=collections.defaultdict(set)
for ri,r in enumerate(vv):
 try:
  code=r['station_code']
  with netCDF4.Dataset(B/r['relative_path']) as n:
   tv=n['time'];fac={'days':86400,'hours':3600,'minutes':60,'seconds':1}[tv.units.split()[0]];epoch=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units);tb=np.asarray(n['time_bnds'][:]);tt=pd.to_datetime(np.rint((tb[:,0]-epoch)*fac).astype('int64'),unit='s');inperiod=(tt>='1996-01-01')&(tt<'2024-01-01');dur=np.rint((tb[:,1]-tb[:,0])*fac);globalmd={k:n.getncattr(k) for k in n.ncattrs()};kept=[]
   for x in r['variables']:
    k=x['name'];v=n[k]
    if getattr(v,'cell_methods','') not in ['time: mean','mean'] or 'time' not in v.dimensions:continue
    base=normalized({**globalmd,**{a:v.getncattr(a) for a in v.ncattrs()}});ndims=[d for d in v.dimensions if d!='time'];idxs=list(itertools.product(*[range(len(n.dimensions[d])) for d in ndims])) if ndims else [()]
    for index in idxs:
     ix=dict(zip(ndims,index));key=tuple(ix[d] if d in ix else slice(None) for d in v.dimensions);raw=np.asarray(np.ma.filled(v[key],np.nan),float).reshape(-1);qn=k+'_qc';unknown=[]
     if qn in n.variables:
      q=n[qn];qkey=tuple(ix[d] if d in ix else slice(None) for d in q.dimensions);dims=[d for d in q.dimensions if d not in ix];qa=np.asarray(np.ma.filled(q[qkey],0));qa=np.moveaxis(qa,dims.index('time'),-1).reshape(-1,len(tt));bad=np.isin(qa,invalid).any(0)&~(qa==100).any(0);unknown=sorted(set(np.unique(qa))-known)
     else:bad=np.zeros(len(tt),bool)
     mds=[base];mn=k+'_ebasmetadata'
     if mn in n.variables:
      mdv=n[mn];mkey=tuple(ix[d] if d in ix else slice(None) for d in mdv.dimensions);mds=[{**base,**normalized(json.loads(str(z)))} for z in np.asarray(mdv[mkey]).ravel()]
     finite=inperiod&np.isfinite(raw)&(raw>=0);valid=finite&~bad;units=sorted({m.get('unit',m.get('units','')) for m in mds});matrices=sorted({m.get('matrix','') for m in mds});levels=sorted({m.get('datalevel',(re.search(r'\.lev(.+)\.nc$',r['file'])[1] if re.search(r'\.lev(.+)\.nc$',r['file']) else 'unknown')) for m in mds});configs={d:float(n[d][i]) for d,i in ix.items()};preferred=True
     for d,i in ix.items():
      vals=np.asarray(n[d][:],float)
      if 'wavelength' in d.lower():preferred&=i==int(np.abs(vals-880).argmin())
      elif 'height' in d.lower():preferred&=i==int(np.argmin(vals))
      else:preferred=False
     totalpm=not any('volatile' in ma or 'non_volatile' in ma for ma in matrices);levok=all((float(le)>=2 if re.fullmatch(r'\d+(\.\d+)?',str(le)) else le in ['3a','unknown']) for le in levels);nrt=('NRT' in r['file']) or any('NRT' in m.get('instrumentname','') for m in mds);identity_block=r['file'] in {z['file'] for z in json.load(open(A/'identity_exclusions.json'))};measureco=set(coords(m) for m in mds);dist=[min([distance(co,z['coords']) for z in ref[code] if distance(co,z['coords']) is not None],default=None) for co in measureco];maxdist=max([z for z in dist if z is not None],default=None)
     for co,distm in zip(measureco,dist):locations.append(dict(Site=M[code]['new'],code=code,file=r['file'],variable=k,configuration=json.dumps(configs),latitude=co[0],longitude=co[1],min_distance_to_PNSD_m=distm,station_names=';'.join(sorted({m.get('stationname','') for m in mds})),platforms=';'.join(sorted({m.get('platformcode','') for m in mds})),PNSD_reference=json.dumps(ref[code],ensure_ascii=False)))
     row=dict(Site=M[code]['new'],code=code,file=r['file'],metric=x['metric'],variable=k,configuration=json.dumps(configs),preferred_channel=bool(preferred),total_PM=totalpm,level_eligible=levok,NRT=nrt,identity_quarantined=identity_block,levels=';'.join(levels),matrix=';'.join(matrices),units_metadata=';'.join(units),units_variable=getattr(v,'units',''),time_rows=int(inperiod.sum()),finite_nonnegative=int(finite.sum()),official_invalid_finite=int((finite&bad).sum()),retained=int(valid.sum()),zero_retained=int((valid&(raw==0)).sum()),negative=int((inperiod&np.isfinite(raw)&(raw<0)).sum()),unknown_flags=';'.join(map(str,unknown)),first=str(tt[valid].min()) if valid.any() else '',last=str(tt[valid].max()) if valid.any() else '',median_duration_s=float(np.median(dur[inperiod])) if inperiod.any() else None,months=int(pd.DatetimeIndex(tt[valid]).to_period('M').nunique()),max_min_distance_to_PNSD_m=maxdist,methods=';'.join(sorted({m.get('methodref','') for m in mds})))
     rows.append(row);kept.append((row,raw,valid))
     if preferred and totalpm and levok and not nrt and not identity_block and valid.any():coverage[code,x['metric']].update(tt[valid].floor('D').strftime('%Y-%m-%d'))
   for a,b in itertools.combinations(kept,2):
    aa,ra,va=a;bb,rb,vb=b
    if aa['metric']==bb['metric'] and aa['units_metadata']==bb['units_metadata'] and aa['configuration']==bb['configuration']:
     both=va&vb;dups.append(dict(Site=M[code]['new'],file=r['file'],metric=aa['metric'],variable_a=aa['variable'],variable_b=bb['variable'],configuration=aa['configuration'],matrix_a=aa['matrix'],matrix_b=bb['matrix'],valid_a=int(va.sum()),valid_b=int(vb.sum()),both_valid=int(both.sum()),equal_both=int(np.isclose(ra[both],rb[both],rtol=1e-10,atol=1e-10).sum()),max_abs_diff=float(np.max(np.abs(ra[both]-rb[both]),initial=0))))
 except Exception as e:errors.append(dict(file=r['file'],error=repr(e)))
 if ri%50==0:print(ri,len(rows),flush=True)
pd.DataFrame(rows).to_csv(A/'all_mean_variables_native_QC.csv',index=False);pd.DataFrame(locations).to_csv(A/'pollutant_PNSD_locations.csv',index=False);pd.DataFrame(dups).to_csv(A/'within_file_candidate_overlap.csv',index=False);(A/'errors.json').write_text(json.dumps(errors,indent=2));(A/'PNSD_reference_locations.json').write_text(json.dumps(ref,indent=2,ensure_ascii=False))
summary=[]
for c in M:
 for metric in ['PM2.5','PM10','PM1','NO2','NO','SO2','O3','CO','BC']:
  days=sorted(coverage.get((c,metric),set()));months=pd.Series(pd.DatetimeIndex(days).to_period('M')).value_counts() if days else pd.Series(dtype=int);summary.append(dict(Site=M[c]['new'],code=c,metric=metric,valid_days_union=len(days),months_with_any=int(len(months)),months_ge14days_union=int((months>=14).sum()),first=days[0] if days else '',last=days[-1] if days else ''))
pd.DataFrame(summary).to_csv(A/'station_metric_coverage_candidate_union.csv',index=False);print('DONE',len(vv),len(rows),'errors',errors)
