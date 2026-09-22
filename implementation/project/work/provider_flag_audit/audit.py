from pathlib import Path
import sys,json,datetime,ast,collections
R=Path('@PROJECT_ROOT@');sys.path[:0]=[str(R/'work/ipr_download/deps'),str(R/'work/ldsa_support/deps'),str(R/'work')]
import numpy as np,pandas as pd,netCDF4
from analyze_cleaned_ldsa import monthly
W=Path(__file__).parent;O=Path('@DELIVERY_ROOT@/provider_quality_flags');O.mkdir(exist_ok=True)
P=R/'outputs/network_trends/corrected_LDSA';D=O.parent/'EBAS_observations/source_netcdf'
reg=pd.read_csv(P/'calculation_parameters/source_qc_registry.csv');cf=pd.read_csv(P/'calculation_parameters/diameter_deposition_weights.csv')
flags=json.load(open(R/'outputs/pnsd_cleaning/PNSD/QC_protocol.json'))['sentinel_values']+[999999.99999]
cases={'GIO_RB':('MT0001R',['2013-02-18','2013-03-25','2013-04-01','2011-07-18']),'KOS_RB':('CZ0003R',['2019-11-18','2009-07-13']),'PAL_RB':('FI0096G',['2014-01-13'])}
for node in ast.parse((R/'work/workbook_processing/restore_original_sequence.py').read_text()).body:
 if isinstance(node,ast.FunctionDef) and node.name=='index_cp':exec(compile(ast.Module(body=[node],type_ignores=[]),'cp','exec'))
results=[];detail=[];profiles=[];sources=[];bp=[]
for site,(code,weeks) in cases.items():
 path=P/'hourly_ldsa'/f'{code}_LDSA.csv.gz'
 if site=='PAL_RB':path=R/'outputs/pal_channel_conversion/FI0096G_LDSA_early_deltaN.csv.gz'
 data=pd.read_csv(path,parse_dates=['Time_UTC']);data['week']=data.Time_UTC.dt.to_period('W').dt.start_time
 for week in weeks:
  t=pd.Timestamp(week);q=data[data.week.eq(t)].copy();v=q['Total_LDSA_um2_cm-3'];finite=v.notna();qvalid=q[finite].sort_values('Total_LDSA_um2_cm-3',ascending=False)
  row=dict(Site=site,week=week,n_hours=int(finite.sum()),n_days=qvalid.Time_UTC.dt.normalize().nunique(),week_mean=float(v.mean()),week_median=float(v.median()),max_hour=float(v.max()),max_date=str(qvalid.Time_UTC.iloc[0]),max_share_percent=float(v.max()/v.sum()*100),top5_share_percent=float(v.nlargest(5).sum()/v.sum()*100))
  local=data[data.Time_UTC.between(t-pd.Timedelta(days=28),t+pd.Timedelta(days=35))&~data.week.eq(t)]['Total_LDSA_um2_cm-3'];row['nearby_median']=float(local.median());row['nearby_valid_hours']=int(local.notna().sum())
  alldetail=[]
  for sid,g in qvalid.groupby('Source_ID'):
   rr=reg[reg.Source_ID.eq(sid)].iloc[0];c=cf[cf.Source_ID.eq(sid)].sort_values('channel');dp=c.diameter_nm.to_numpy();sel=(dp>=17)&(dp<=500);weight=c.delta_log10.to_numpy()*np.pi*(dp/1000)**2*c.DF_total.to_numpy()
   with netCDF4.Dataset(D/code/'PNSD'/rr.file) as n:
    tb=np.asarray(n['time_bnds'][:,0]);tv=n['time'];ep=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units);times=pd.to_datetime(np.rint((tb-ep)*86400).astype('int64'),unit='s').round('h');idx=pd.Index(times).get_indexer(g.Time_UTC);assert (idx>=0).all()
    raw=np.asarray(np.ma.filled(n[rr.variable][:,idx],np.nan),float).T;z=np.where(~np.isfinite(raw)|(raw<0)|np.isin(raw,flags),np.nan,raw);contrib=z*weight;calc=np.nansum(contrib[:,sel],axis=1)
    assert np.allclose(calc,g['Total_LDSA_um2_cm-3'],rtol=1e-10)
    qcvar=n[rr.variable+'_qc'];qc=np.asarray(qcvar[:,:,idx].filled(0));qcs=np.moveaxis(qc,2,0)
    md=n[rr.variable+'_ebasmetadata'];meta=json.loads(str(md[0,0]));sources.append(dict(site=site,sid=int(sid),file=rr.file,url=rr.url,metadata=meta))
    for j,(_,a) in enumerate(g.iterrows()):
     use=sel&np.isfinite(z[j]);qcvals=qcs[j,use,:].ravel();qs=sorted(set(map(int,qcvals))-{0});validvals=z[j,sel]
     order=np.argsort(np.nan_to_num(contrib[j],nan=-1))[::-1];order=[k for k in order if sel[k]][:5]
     entry=dict(Site=site,week=week,date=str(a.Time_UTC),LDSA=float(a['Total_LDSA_um2_cm-3']),PNC=float(a['PNC_cm-3']),valid_bins=int(use.sum()),window_bins=int(sel.sum()),max_raw=float(np.nanmax(validvals)),flags=';'.join(map(str,qs)),top_channel_nm=float(dp[order[0]]),top_channel_share=float(contrib[j,order[0]]/calc[j]*100),top5_channel_share=float(np.nansum(contrib[j,order])/calc[j]*100),Source_ID=int(sid))
     alldetail.append(entry)
     if a.Time_UTC in qvalid.head(5).Time_UTC.values:
      for k in range(len(dp)):profiles.append(dict(Site=site,week=week,date=str(a.Time_UTC),diameter_nm=dp[k],raw=raw[j,k],selected=bool(sel[k]),ldsa_contribution=contrib[j,k],flags=';'.join(map(str,sorted(set(map(int,qcs[j,k]))-{0})))))
  row['qc_flags']=sorted({f for z in alldetail for f in z['flags'].split(';') if f});detail+=alldetail;results.append(row)
 # Diagnostic only: remove entire specified plotted weeks, preserving production files.
 for label,exclude in [('current',[]),('circled_weeks',weeks[:3] if site=='GIO_RB' else weeks[:1]),('all_listed_high_weeks',weeks)]:
  d=data[~data.week.isin(pd.to_datetime(exclude))][['Time_UTC','Total_LDSA_um2_cm-3']].rename(columns={'Time_UTC':'date','Total_LDSA_um2_cm-3':'value'}).dropna().reset_index(drop=True);m=monthly(d.date,d.value);k=index_cp(d,len(m));bp.append(dict(Site=site,scenario=label,excluded_weeks=';'.join(exclude),n=len(d),breakpoint=str(d.date.iloc[k]) if k is not None else '',overall_mean=d.value.mean()))
pd.DataFrame(results).to_csv(O/'quality_flag_summary.csv',index=False,encoding='utf-8-sig');pd.DataFrame(detail).to_csv(O/'quality_flag_details.csv',index=False,encoding='utf-8-sig');pd.DataFrame(profiles).to_csv(O/'quality_flag_profiles.csv',index=False,encoding='utf-8-sig');pd.DataFrame(bp).to_csv(O/'breakpoint_quality_flags.csv',index=False,encoding='utf-8-sig');(O/'quality_flag_sources.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2))
print(pd.DataFrame(results).to_string(index=False));print(pd.DataFrame(bp).to_string(index=False));print('HOURLY TOP');print(pd.DataFrame(detail).sort_values('LDSA',ascending=False).head(18).to_string(index=False))
