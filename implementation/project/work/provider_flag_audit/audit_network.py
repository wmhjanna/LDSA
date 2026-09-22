from pathlib import Path
import sys,json,datetime,collections
R=Path('@PROJECT_ROOT@');sys.path[:0]=[str(R/'work/ipr_download/deps'),str(R/'work/ldsa_support/deps')]
import pandas as pd,numpy as np,netCDF4
W=Path(__file__).parent;O=Path('@DELIVERY_ROOT@/provider_quality_flags');O.mkdir(exist_ok=True)
D=O.parent/'EBAS_observations/source_netcdf';P=R/'outputs/network_trends/corrected_LDSA'
table=pd.read_csv(W/'ebas_flags.csv',header=None,quotechar="'",names=['flag','description','status']);classes=table.set_index('flag').status.to_dict();desc=table.set_index('flag').description.to_dict();invalid=table[table.status.isin(['I','H','M'])].flag.to_numpy();warnings=table[table.status.eq('V')&~table.flag.isin([100,110,111,120,210,211])].flag.to_numpy()
sentinels=json.load(open(R/'outputs/pnsd_cleaning/PNSD/QC_protocol.json'))['sentinel_values']+[999999.99999]
reg=pd.read_csv(P/'calculation_parameters/source_qc_registry.csv');coef=pd.read_csv(P/'calculation_parameters/diameter_deposition_weights.csv');mapping=json.load(open(R/'work/station_registry/mapping.json'));newnames={x['old']:x['new'] for x in mapping}
newdata={};jobs=[]
for code,rr in reg.groupby('code'):
 f=P/'hourly_ldsa'/f'{code}_LDSA.csv.gz'
 if code=='FI0096G':f=R/'outputs/pal_channel_conversion/FI0096G_LDSA_early_deltaN.csv.gz'
 newdata[code]=pd.read_csv(f,usecols=['Time_UTC','Source_ID','Total_LDSA_um2_cm-3'],parse_dates=['Time_UTC'])
 for r in rr.itertuples():
  c=coef[coef.Source_ID.eq(r.Source_ID)];dp=c[c.selected_17_500].diameter_nm.to_numpy();times=newdata[code].loc[newdata[code].Source_ID.eq(r.Source_ID),'Time_UTC']
  if len(times):jobs.append(dict(scope='current_netcdf',code=code,site=newnames[code],file=D/code/'PNSD'/r.file,variable=r.variable,sid=r.Source_ID,dp=dp,times=times))
oldmap=[r for r in json.load(open(R/'work/ebas_download/local_comparison.json')) if r['current_analysis_sites']]
for r in oldmap:
 spans=[]
 for site in r['current_analysis_sites']:
  f=R/'work/workbook_processing/series'/f'{site}_series.pkl'
  if f.exists():
   d=pd.read_pickle(f);spans.append((d.Time.min(),d.Time.max()))
 if not spans:continue
 lo=min(t[0] for t in spans);hi=max(t[1] for t in spans)
 for f in sorted((D/r['code']/'PNSD').glob('*.nc')):
  if '.lev1.' in f.name:continue
  jobs.append(dict(scope='workbook_source_inventory',code=r['code'],site=';'.join(r['current_analysis_sites']),file=f,lo=lo,hi=hi))
rows=[];flagrows=[];affected=[];errors=[]
def attrs(n,v):
 md={k:str(v.getncattr(k)) for k in v.ncattrs()}
 if v.name+'_ebasmetadata' in n.variables:
  a=n[v.name+'_ebasmetadata'];md.update(json.loads(str(a[tuple(0 for _ in a.dimensions)])))
 canon=lambda s:''.join(c for c in s.lower() if c.isalnum()).removeprefix('ebas')
 return {canon(k):v for k,v in md.items()}
for jj,job in enumerate(jobs):
 try:
  with netCDF4.Dataset(job['file']) as n:
   vv=[job['variable']] if 'variable' in job else [k for k,v in n.variables.items() if getattr(v,'cell_methods','') in ['time: mean','mean'] and (getattr(v,'ebas_component','')=='particle_number_size_distribution' or k.startswith('particle_number_size_distribution')) and 'time' in v.dimensions]
   tv=n['time'];tb=np.asarray(n['time_bnds'][:,0]);ep=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units);fac={'days':86400,'hours':3600,'seconds':1,'minutes':60}[tv.units.split()[0]];tt=pd.to_datetime(np.rint((tb-ep)*fac).astype('int64'),unit='s').round('h')
   timeok=tt.isin(job['times']) if 'times' in job else (tt>=job['lo'])&(tt<=job['hi'])
   if not timeok.any():continue
   for k in vv:
    v=n[k];md=attrs(n,v)
    if job['scope']=='workbook_source_inventory' and (md.get('datalevel')!='2' or 'non' in md.get('matrix','')):continue
    dims=[x for x in v.dimensions if x!='time'];da=next((x for x in dims if x.startswith('D') and x in n.variables),None)
    if da is None:continue
    dp=np.asarray(n[da][:]);sel=np.isin(dp,job['dp']) if 'dp' in job else (dp>= (31.65 if 'LND_TR' in job['site'] else 17))&(dp<=500)
    if not sel.any():continue
    if v.ndim!=2:errors.append(dict(file=job['file'].name,variable=k,error='additional data dimensions'));continue
    finite_n=inv_n=warn_n=all_invalid_n=inv_cells=warn_cells=qc_n=0;codes=collections.Counter();icode=collections.Counter();wcode=collections.Counter();unknown=set();dates_bad=[]
    hasqc=k+'_qc' in n.variables
    for st in range(0,len(tt),4096):
     stop=min(st+4096,len(tt));use=timeok[st:stop]
     if not use.any():continue
     ix=[slice(None)]*v.ndim;ix[v.dimensions.index('time')]=slice(st,stop);raw=np.moveaxis(np.asarray(np.ma.filled(v[tuple(ix)],np.nan),float),v.dimensions.index('time'),0)[:,sel][use]
     finite=np.isfinite(raw)&(raw>=0)&~np.isin(raw,sentinels);fh=finite.any(axis=1);finite_n+=int(fh.sum())
     if not hasqc:continue
     q=n[k+'_qc'];ix=[slice(None)]*q.ndim;ix[q.dimensions.index('time')]=slice(st,stop);qa=np.asarray(q[tuple(ix)].filled(0));qa=np.transpose(qa,[q.dimensions.index('time'),q.dimensions.index(da)]+[a for a,d in enumerate(q.dimensions) if d not in ['time',da]])[:,sel,:][use]
     override=(qa==100).any(axis=2);bad=np.isin(qa,invalid).any(axis=2)&~override;warn=np.isin(qa,warnings).any(axis=2)&~bad
     ih=(bad&finite).any(axis=1);wh=(warn&finite).any(axis=1);dates_bad.extend(tt[st:stop][use][ih].astype(str).tolist());inv_n+=int(ih.sum());warn_n+=int(wh.sum());all_invalid_n+=int((fh&~(finite&~bad).any(axis=1)).sum());inv_cells+=int((bad&finite).sum());warn_cells+=int((warn&finite).sum())
     for code in np.unique(qa):
      if not code:continue
      hit=(qa==code).any(axis=2)&finite;cnt=int(hit.sum())
      if cnt:codes[int(code)]+=cnt
      if code not in classes:unknown.add(int(code))
    row=dict(scope=job['scope'],Site=job['site'],code=job['code'],file=job['file'].name,variable=k,Source_ID=job.get('sid',''),has_qc=hasqc,finite_hours=finite_n,invalid_flag_hours=inv_n,all_bins_invalid_hours=all_invalid_n,warning_flag_hours=warn_n,invalid_cells=inv_cells,warning_cells=warn_cells,unknown_flags=';'.join(map(str,sorted(unknown))),nonzero_flags=';'.join(map(str,sorted(codes))))
    rows.append(row)
    for code,cnt in codes.items():flagrows.append(dict(scope=job['scope'],Site=job['site'],file=job['file'].name,flag=code,status=classes.get(code,'unknown'),description=desc.get(code,''),finite_cell_count=cnt))
    for date in dates_bad:affected.append(dict(scope=job['scope'],Site=job['site'],code=job['code'],Source_ID=job.get('sid',''),date=date,file=job['file'].name))
 except Exception as e:errors.append(dict(file=job['file'].name,error=repr(e)))
 if jj%10==0:print('AUDIT',jj,'/',len(jobs),'rows',len(rows),'errors',len(errors),flush=True)
pd.DataFrame(rows).to_csv(O/'workbook_source_flags.csv',index=False,encoding='utf-8-sig');pd.DataFrame(flagrows).to_csv(O/'quality_flag_counts.csv',index=False,encoding='utf-8-sig');pd.DataFrame(affected).to_csv(O/'flagged_observations.csv.gz',index=False);(O/'source_processing_errors.json').write_text(json.dumps(errors,ensure_ascii=False,indent=2))
df=pd.DataFrame(rows);agg=df.groupby(['scope','Site','code']).agg(source_variables=('file','size'),finite_hour_records=('finite_hours','sum'),invalid_hour_records=('invalid_flag_hours','sum'),all_bins_invalid=('all_bins_invalid_hours','sum'),warning_hour_records=('warning_flag_hours','sum'),invalid_cells=('invalid_cells','sum')).reset_index();agg.to_csv(O/'site_quality_flag_summary.csv',index=False,encoding='utf-8-sig');table.to_csv(O/'EBAS_source_quality_flags.csv',index=False,encoding='utf-8-sig')
print(agg[agg.invalid_hour_records.gt(0)].sort_values(['scope','invalid_hour_records'],ascending=[True,False]).to_string(index=False));print('ERRORS',errors)
