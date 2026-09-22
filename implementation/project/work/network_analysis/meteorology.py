from pathlib import Path
import sys,json,re,datetime,importlib.util
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent/'meteorology';W.mkdir(exist_ok=True);sys.path[:0]=[str(R/'work/ipr_download/deps'),str(R/'work/ldsa_support/deps')]
import netCDF4,numpy as np,pandas as pd
B=Path('@DELIVERY_ROOT@/netcdf_meteorology');reg=pd.read_csv(B/'source_metadata/meteorology_missingness.csv');files=pd.read_csv(B/'source_metadata/meteorology_input_records.csv').set_index('Source_ID');overview=pd.read_csv(B/'netcdf_meteorology_summary.csv').set_index('Site')
spec=importlib.util.spec_from_file_location('poll_helpers',R/'work/netcdf_pollutants/extract.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
# Choose the lowest available inlet for each native file/metric, retaining its exact coordinates.
def height(r):
 z=re.findall(r'[-+]?\d*\.?\d+',str(r.height));return float(z[-1]) if z else 1e9
reg['height_rank']=reg.apply(height,axis=1);reg=reg.sort_values(['file','Metric','height_rank','Series_ID']).drop_duplicates(['file','Metric']);audit=[];parts=[]
for r in reg.itertuples():
 if pd.notna(r.data_level) and r.data_level not in [2,3]:continue
 if pd.notna(r.distance_from_pnsd_km) and r.distance_from_pnsd_km>5:raise ValueError(('location',r.Site))
 path=B/files.loc[r.Source_ID,'relative_path'];sel=json.loads(r.extra_dimension_indices)
 with netCDF4.Dataset(path) as ds:
  v=ds[r.variable];rawarr,dims=h.sliced(v,sel);assert dims==('time',);raw=np.asarray(np.ma.filled(rawarr,np.nan),float);tb=np.asarray(ds['time_bnds'][:]);dates=h.dt(tb[:,0],ds['time'].units);ends=h.dt(tb[:,1],ds['time'].units);dur=(ends-dates).total_seconds().to_numpy()
  qname=r.variable+'_qc'
  if qname in ds.variables:
   qv=ds[qname];qv.set_auto_mask(False);qa,qd=h.sliced(qv,sel);qc=np.asarray(qa,int);qc=np.moveaxis(qc,qd.index('time'),0).reshape(len(raw),-1)
  else:qc=np.zeros((len(raw),1),int)
  invalid=np.isin(qc,list(h.invalid)).any(1)&~(qc==100).any(1);value=raw.copy();value[invalid|~np.isfinite(raw)]=np.nan
  # Keep native hourly meteorology; long accumulation precipitation is not an hourly rainfall predictor.
  timeok=(dates>=pd.Timestamp(overview.loc[r.Site,'PNSD_start']))&(dates<=pd.Timestamp(overview.loc[r.Site,'PNSD_end']))&(dur>0)&(dur<=3601)
  ix=np.flatnonzero(timeok)
  if not len(ix):audit.append(dict(Site=r.Site,Metric=r.Metric,Series_ID=r.Series_ID,used=False,reason='no native hourly observations overlapping LDSA'));continue
  d=pd.DataFrame(dict(date=dates[ix],end=ends[ix],value=value[ix],raw_value=raw[ix],official_invalid=invalid[ix],series_id=r.Series_ID,metric=r.Metric,Site=r.Site,duration_s=dur[ix],level=2 if pd.notna(r.data_level) else 1.9,revision=r.file.split('.')[2]))
  parts.append(d);audit.append(dict(Site=r.Site,Metric=r.Metric,Series_ID=r.Series_ID,used=True,reason='',observations=len(d),invalid_finite=int((invalid[ix]&np.isfinite(raw[ix])).sum()),height=r.height,file=r.file))
selected=[]
for (site,metric),g in pd.concat(parts,ignore_index=True).groupby(['Site','metric']):
 ranks=g.groupby('series_id').agg(level=('level','max'),duration=('duration_s','median'),revision=('revision','max')).sort_values(['level','duration','revision'],ascending=[False,True,False],kind='stable');accepted=g.iloc[:0]
 for sid in ranks.index:
  q=g[g.series_id.eq(sid)].sort_values('date').drop_duplicates('date');accepted=pd.concat([accepted,h.no_overlap(q,accepted)],ignore_index=True)
 selected.append(accepted)
selected=pd.concat(selected,ignore_index=True).sort_values(['Site','metric','date']);selected.to_pickle(W/'selected.pkl');pd.DataFrame(audit).to_csv(W/'source_selection_and_QC.csv',index=False);reg.to_csv(W/'selected_native_configurations.csv',index=False)
weekly=[]
for (site,metric),g in selected.groupby(['Site','metric']):
 z=g.groupby(g.date.dt.to_period('W').dt.start_time).value.mean().rename_axis('date').reset_index().assign(Site=site,Metric=metric);weekly.append(z)
pd.concat(weekly,ignore_index=True).to_pickle(W/'weekly.pkl');summary=selected.groupby(['Site','metric']).agg(rows=('value','size'),valid=('value','count'),invalid=('official_invalid','sum'),start=('date','min'),end=('end','max'));summary.to_csv(W/'summary.csv');print(summary.to_string(),flush=True);print('DONE',len(selected),flush=True)
