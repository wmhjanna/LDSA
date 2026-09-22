from pathlib import Path
import sys,pickle,json,hashlib,re,gc
from zoneinfo import ZoneInfo
from datetime import timezone
ROOT=Path('@PROJECT_ROOT@');sys.path[:0]=[str(ROOT/'work/ldsa_support/deps'),str(ROOT/'work')]
import pandas as pd,numpy as np
from build_site_tables import delta_log10
from analyze_cleaned_ldsa import monthly
W=ROOT/'work/quality_control/workbook';O=ROOT/'outputs/quality_control/workbook';A=O/'analysis';B=Path('@INPUT_ROOT@/data_sources');PREV=ROOT/'outputs/pnsd_cleaning/analysis';sp=pd.read_csv(PREV/'site_spatial_statistics.csv').set_index('Site');oldms=pd.read_csv(PREV/'monthly_series.csv',parse_dates=['date']);flags=json.loads((ROOT/'outputs/pnsd_cleaning/PNSD/QC_protocol.json').read_text())['sentinel_values'];rev={'SCH_RB':'SCH_SUB','WAL_RB':'WAL_SUB','ZUG_RB':'ZUG_SUB'}
from joined_sources import joined_template, joined_parts
sp=sp[~sp.index.str.startswith('PRA_')].copy()
for joined_site in ['PAR_SUB','PRA_SUB']:sp.loc[joined_site,'Period']='2015-2019'
for d in ['pnsd','meta_raw']: (W/d).mkdir(exist_ok=True)
qc=[];sources=[];channelqc=[];configs=[];checks=[];stats=[];sens=[];elig=[];breaks=[];groups=[];mseries=[];fullms=[];weekly=[];modesmanifest=[];comp=[]
def rawmonthly(t,v):
 d=pd.DataFrame({'date':t,'value':v}).replace([np.inf,-np.inf],np.nan).dropna().sort_values('date').drop_duplicates('date',keep='first');return d.groupby(d.date.dt.to_period('M').dt.to_timestamp()).value.mean().rename_axis('date').reset_index()
def hourly_bp(t,v,min_months):
 if min_months<20:return None
 d=pd.DataFrame({'date':t,'value':v}).replace([np.inf,-np.inf],np.nan).dropna().sort_values('date');y=d.value.to_numpy(float)
 if len(y)<4 or np.std(y)==0:return None
 y=y-y.mean();s=np.cumsum(y);k=np.arange(2,len(y)-1);cost=s[k-1]**2/k+(s[-1]-s[k-1])**2/(len(y)-k);i=int(k[np.argmax(cost)]);return pd.Timestamp(d.date.iloc[i])
polbook=pd.ExcelFile(B/'predictive_models/pollutant_data/all_sites/paired_date_value_measurements.xlsx',engine='calamine')
checkpoint=W/'prepare_checkpoint.pkl'
if checkpoint.exists():
 saved=pickle.load(open(checkpoint,'rb'));globals().update(saved)
done={x['Site'] for x in qc}
for site in sorted(sp.index):
 if site in done:continue
 oldsite=rev.get(site,site);physical='PRA_SUB' if re.match(r'PRA_\d',site) else site;lo=31.65 if site=='LND_TR' else 17
 template=joined_template(site) if site in ['PAR_SUB','PRA_SUB'] else pickle.load(open(ROOT/'work/pnsd_cleaning/pnsd'/(site+'.pkl'),'rb'));dp=template['diameters'];cf=template['coefficients'];oldaux=template['aux'];parts=[]
 if site in ['PAR_SUB','PRA_SUB']:
  parts=joined_parts(site)
 elif site=='LAN_UB':
  meta=pickle.load(open(ROOT/'work/pnsd_qc_cache/LAN_UB.pkl','rb'))
  for p in meta:parts.append((Path(p['path']),p['sheet'],p['dcols']))
 else:
  srcsite='PRA_SUB_'+site.split('_')[1] if re.match(r'PRA_\d',site) else oldsite;path=B/'source_data'/srcsite/(srcsite+'.xlsx');parts=[(path,0,None)]
 bs=[];st=[];ts=[];sid=[];sq=dict(Site=site,Physical_site=physical,rows=0,channels=len(dp),negative=0,sentinel=0,nonfinite_or_nonnumeric=0,zero_retained=0,statistical_extremes_removed=0)
 dest=O/'PNSD'/site;dest.mkdir(exist_ok=True)
 for partidx,(path,sheet,dcols) in enumerate(parts):
  r=path.copy() if isinstance(path,pd.DataFrame) else pd.read_excel(path,sheet_name=sheet,engine='calamine');source_path=Path(r.attrs['source_path']) if isinstance(path,pd.DataFrame) else path;dcols=dcols or [(i,float(str(c).removeprefix('D_'))) for i,c in enumerate(r.columns) if i>0 and (isinstance(c,(float,int,np.number)) or re.fullmatch(r'D_[0-9.]+',str(c)))]
  native=np.array([d for i,d in dcols]);ix=np.abs(dp[:,None]-native).argmin(axis=0);assert np.allclose(dp[ix],native,atol=.00051,rtol=0),(site,'native diameters')
  time=pd.to_datetime(r.iloc[:,0],errors='coerce');valid=time.notna()&time.dt.year.between(1990,2030);time=time[valid].reset_index(drop=True);x=r.iloc[valid.to_numpy(),[i for i,d in dcols]].apply(pd.to_numeric,errors='coerce').to_numpy(float)
  if site=='ROC_UB':time=pd.Series([t.to_pydatetime().replace(tzinfo=ZoneInfo('America/New_York')).astimezone(timezone.utc).replace(tzinfo=None) for t in time])
  bad=~np.isfinite(x)|(x<0)|np.isin(x,flags);clean=np.where(bad,np.nan,x);assert np.array_equal(clean[~bad],x[~bad]);assert (clean[x==0]==0).all();bins=np.full((len(x),len(dp)),np.nan);struct=np.ones_like(bins,bool);bins[:,ix]=clean*delta_log10(native);struct[:,ix]=False
  if site=='MUL_SUB':
   before=time.lt('2018-01-29 15:30').to_numpy();struct[before]|=dp<14.075;struct[~before]|=dp>495.807;assert np.isnan(bins[struct]).all()
  ii,jj=np.where(np.isfinite(x)&bad);pd.DataFrame({'Source_excel_row':np.flatnonzero(valid.to_numpy())[ii]+int(r.attrs.get('source_row_start',2)),'Time':time.iloc[ii].to_numpy(),'diameter_nm':native[jj],'original_dN_dlog10Dp':x[ii,jj],'reason':np.where(x[ii,jj]<0,'negative','sentinel')}).to_csv(dest/f'invalid_cells_source_{partidx+1}.csv.gz',index=False)
  for key,count in [('negative',(np.isfinite(x)&(x<0)).sum()),('sentinel',np.isin(x,flags).sum()),('nonfinite_or_nonnumeric',(~np.isfinite(x)).sum()),('zero_retained',(x==0).sum())]:sq[key]+=int(count)
  for j,d in enumerate(native):channelqc.append(dict(Site=site,Config=partidx+1,diameter_nm=d,delta_log10=delta_log10(native)[j],negative=int((np.isfinite(x[:,j])&(x[:,j]<0)).sum()),sentinel=int(np.isin(x[:,j],flags).sum()),zero_retained=int((x[:,j]==0).sum())))
  selected=native[(native>=lo)&(native<=500)];configs.append(dict(Site=site,Physical_site=physical,Config=f'{site}__{partidx+1}',start=time.min(),end=time.max(),nominal_low_nm=lo,nominal_high_nm=500,first_center_nm=selected.min(),last_center_nm=selected.max(),centers=len(selected)))
  sources.append(dict(Site=site,source=str(source_path),sheet=sheet,sha256=hashlib.file_digest(source_path.open('rb'),'sha256').hexdigest(),source_rows=len(r),observation_rows=len(time),nonobservation_rows=int((~valid).sum())))
  bs.append(bins);st.append(struct);ts.append(time);sid.extend([partidx+1]*len(time))
 bins=np.concatenate(bs);struct=np.concatenate(st);tt=pd.concat(ts,ignore_index=True);sq['rows']=len(tt);qc.append(sq)
 masks=pickle.load(open(ROOT/'work/quality_control/old_qc_masks.pkl','rb'));mask=masks.get(site,np.zeros(bins.shape,bool));assert mask.shape==bins.shape
 if mask.any():
  prev=pickle.load(open(ROOT/'work/workbook_processing/pnsd'/(site+'.pkl'),'rb'));assert np.allclose(bins,prev['bins'],equal_nan=True,rtol=1e-12,atol=1e-9);assert tt.equals(prev['aux'].Time.reset_index(drop=True));bins[mask]=np.nan
 sq['official_invalid_cells_removed']=int(mask.sum());sq['official_invalid_rows_affected']=int(mask.any(axis=1).sum())
 # Ancillary observations retain their observed values. Match interval identity for source midpoint/start differences.
 if np.array_equal(tt.to_numpy(),oldaux.Time.to_numpy()):s=oldaux.copy()
 else:
  ax=oldaux.drop(columns='Time').copy();ax.index=pd.DatetimeIndex(oldaux.Time).round('min').floor('h');ax=ax.groupby(level=0).mean();s=ax.reindex(pd.DatetimeIndex(tt).round('min').floor('h')).reset_index(drop=True);s.insert(0,'Time',tt)
 if site=='HEL_TR':
  for c in s.columns:
   if 'SO2' in str(c):s[c]=np.nan
 def integrate(v,a,b,inc=True,strict=False):
  mask=(dp>=a)&((dp<=b) if inc else (dp<b))
  if not mask.any():return np.full(len(tt),np.nan)
  sel=v[:,mask];out=np.sum(np.where(struct[:,mask],0,sel),axis=1) if strict else np.nansum(sel,axis=1);out[np.isnan(sel).all(axis=1)]=np.nan;return out
 pnsd=pd.DataFrame({'Time':tt,'Source_part':sid});pnsd=pd.concat([pnsd,pd.DataFrame(bins,columns=[f'D_{d:g}nm_cm-3' for d in dp])],axis=1);pnsd.to_csv(dest/(site+'_PNSD_bins_step1_only.csv.gz'),index=False);cf.to_csv(dest/(site+'_DF_coefficients.csv'),index=False)
 pickle.dump(dict(site=site,physical=physical,bins=bins,structural=struct,diameters=dp,aux=s.copy(),coefficients=cf),open(W/'pnsd'/(site+'.pkl'),'wb'),protocol=5)
 for region in ['HA','TB','ALV','Total']:
  v=bins*np.pi*(dp/1000)**2*cf[region].to_numpy()
  for name,a,b in [('primary',lo,500),('17_500nm',17,500),('31.65_500nm',31.65,500),('10_500nm',10,500),('10_800nm',10,800)]:s[f'{region}_LDSA_{name}']=integrate(v,a,b)
 s['LDSA']=s.Total_LDSA_primary;s['LDSA_strict_primary']=integrate(bins*np.pi*(dp/1000)**2*cf.Total.to_numpy(),lo,500,strict=True)
 for name,a,b,inc in [('Ntotal',lo,500,True),('Nucleation',lo,25,False),('Aitken',max(lo,25),100,False),('Accumulation',100,500,True)]:
  s[f'PNC_{name}_primary']=integrate(bins,a,b,inc);d=dp[(dp>=a)&((dp<=b) if inc else (dp<b))];modesmanifest.append(dict(Site=site,Metric=name,lower_nm=a,upper_nm=b,actual_min_center_nm=d.min() if len(d) else np.nan,actual_max_center_nm=d.max() if len(d) else np.nan,channel_centers=len(d),available=bool(len(d))))
 for i in np.linspace(0,len(tt)-1,101,dtype=int):
  sel=(dp>=lo)&(dp<=500)&np.isfinite(bins[i]);value=sum(float(bins[i,j])*np.pi*(float(dp[j])/1000)**2*float(cf.Total.iloc[j]) for j in np.flatnonzero(sel)) if sel.any() else np.nan;assert np.isclose(value,s.LDSA.iloc[i],rtol=1e-12,atol=1e-8,equal_nan=True)
 checks.append(dict(Site=site,rows=len(s),valid_primary=int(s.LDSA.notna().sum()),valid_strict_primary=int(s.LDSA_strict_primary.notna().sum()),partial_rows_retained=int((s.LDSA.notna()&s.LDSA_strict_primary.isna()).sum()),independent_formula_rows=101,formula_verified=True))
 s.to_pickle(W/'series'/(site+'_series.pkl'));sd=O/'sites'/site;sd.mkdir(exist_ok=True);s[['Time','LDSA']+[c for c in s if c.startswith(('HA_','TB_','ALV_','Total_','PNC_'))]].to_csv(sd/(site+'_LDSA_PNC_reference.csv.gz'),index=False)
 s[['Time','LDSA','Total_LDSA_17_500nm','Total_LDSA_31.65_500nm','LDSA_strict_primary']].to_csv(W/'meta_raw'/(site+'.csv.gz'),index=False)
 typ=site.rsplit('_',1)[1];yl,yh=map(int,sp.loc[site,'Period'].split('-'));ss=s[s.Time.dt.year.between(yl,yh)];v=ss.LDSA.dropna();stats.append(dict(Site=site,Physical_site=physical,Type=typ,Period=f'{yl}-{yh}',lower_nm=lo,upper_nm=500,n=len(v),Mean=v.mean(),Median=v.median(),P25=v.quantile(.25),P75=v.quantile(.75),P95=v.quantile(.95),**{reg:ss.loc[ss.LDSA.notna(),reg+'_LDSA_primary'].mean() for reg in ['HA','TB','ALV']}))
 mask=ss[['Total_LDSA_'+n for n in ['primary','10_500nm','10_800nm']]].notna().all(axis=1)
 for label,key in [('Primary','primary'),('10-500','10_500nm'),('10-800','10_800nm')]:
  for reg in ['Total','HA','TB','ALV']:
   v=ss.loc[mask,reg+'_LDSA_'+key].dropna();sens.append(dict(Site=site,Type=typ,Period=f'{yl}-{yh}',Range=label,primary_low_nm=lo,Region=reg,n=len(v),Mean=v.mean(),Median=v.median(),P25=v.quantile(.25),P75=v.quantile(.75),P95=v.quantile(.95)))
 raws={'LDSA':pd.DataFrame({'date':s.Time,'value':s.LDSA})}
 for name in ['Ntotal','Nucleation','Aitken','Accumulation']:raws[name]=pd.DataFrame({'date':s.Time,'value':s['PNC_'+name+'_primary']})
 polsite='PRA_SUB' if re.match(r'PRA_\d',site) else oldsite
 if polsite in polbook.sheet_names:
  pol=polbook.parse(polsite)
  for i in range(0,len(pol.columns)-1,2):
   name=str(pol.columns[i+1]).strip()
   if name in ['PM2.5','PM10','PM1','BC','NO2','SO2','O3','CO','NO']:raws[name]=pd.DataFrame({'date':pd.to_datetime(pol.iloc[:,i],errors='coerce'),'value':pd.to_numeric(pol.iloc[:,i+1],errors='coerce')})
 if re.match(r'PRA_\d',site):
  y0,y1=map(int,site.split('_')[1].split('-'))
  for name in list(raws):
   if name not in ['LDSA','Ntotal','Nucleation','Aitken','Accumulation']:raws[name]=raws[name][raws[name].date.dt.year.between(y0,y1)]
 # LAN ancillary monthly records have a separate archived source; preserve their dates and values.
 for name,g in oldms[oldms.Site.eq(site)&~oldms.Metric.isin(raws)].groupby('Metric'):raws[name]=g[['date','value']].copy()
 if site=='HEL_TR':raws.pop('SO2',None)
 from qc_pollutants import apply_current_pollutants
 raws=apply_current_pollutants(site,raws)
 mm={name:monthly(g.date,g.value) for name,g in raws.items()}
 if site=='LAN_UB':
  for name,g in oldms[oldms.Site.eq(site)&~oldms.Metric.isin(['LDSA','Ntotal','Nucleation','Aitken','Accumulation'])].groupby('Metric'):mm[name]=g[['date','value']].copy()
 anchor=hourly_bp(s.Time,s.LDSA,len(mm['LDSA']));anchor_month=anchor.to_period('M').to_timestamp() if anchor is not None else None
 n=len(mm['LDSA']);dates=mm['LDSA'].date;span=(dates.max().year-dates.min().year)*12+dates.max().month-dates.min().month+1 if n else 0;elig.append(dict(Site=site,Physical_site=physical,Type=typ,valid_months=n,paired_months_14days=n,span_months=span,coverage=n/span if span else 0,meta_eligible=False,breakpoint=anchor_month,breakpoint_timestamp=anchor,primary_low_nm=lo,primary_high_nm=500))
 for metric,d in mm.items():
  z=d.copy();z['Site']=site;z['Metric']=metric;mseries.append(z);raw=raws[metric];b=hourly_bp(raw.date,raw.value,len(d));bm=b.to_period('M').to_timestamp() if b is not None else None
  if b is not None:breaks.append(dict(Site=site,Type=typ,Physical_site=physical,Metric=metric,Breakpoint=bm,Breakpoint_timestamp=b,n_months=len(d),offset_months=(bm.year-anchor_month.year)*12+bm.month-anchor_month.month if anchor_month is not None else np.nan))
  for phase in ['Full','Before','After']:
   phase_cut=b if metric in ['LDSA','Ntotal','Nucleation','Aitken','Accumulation'] else anchor_month
   if phase!='Full' and phase_cut is None:continue
   rr=raw if phase=='Full' else raw[raw.date.lt(phase_cut) if phase=='Before' else raw.date.ge(phase_cut)];g=rawmonthly(rr.date,rr.value)
   if phase=='Full':zz=g.copy();zz['Site']=site;zz['Metric']=metric;fullms.append(zz)
   if len(g)<3:continue
   g['Site']=site;g['Type']=typ;g['Metric']=metric;g['Phase']=phase;g['group_id']=site+'|'+metric+'|'+phase;groups.append(g)
  if metric!='LDSA':
   pair=mm['LDSA'].merge(d,on='date',suffixes=('_ldsa','_other'))
   if len(pair)>=3:
    for role,col in [('LDSA','value_ldsa'),(metric,'value_other')]:
     g=pair[['date',col]].rename(columns={col:'value'});g['Site']=site;g['Type']=typ;g['Metric']=role;g['Phase']='Overlap_'+metric;g['group_id']=site+'|'+role+'|Overlap_'+metric;groups.append(g)
 numeric=s.set_index('Time').select_dtypes('number');numeric.index=numeric.index.to_period('W').to_timestamp();wk=numeric.groupby(level=0).mean().rename_axis('date').reset_index();wk['Site']=site;wk['Type']=typ;weekly.append(wk)
 print(site,'rows',len(s),'partial',checks[-1]['partial_rows_retained'],'paired months',n,'hourly breakpoint',anchor,flush=True)
 pd.DataFrame(qc).to_csv(O/'PNSD/QC_summary.csv',index=False);pickle.dump({k:globals()[k] for k in ['qc', 'sources', 'channelqc', 'configs', 'checks', 'stats', 'sens', 'elig', 'breaks', 'groups', 'mseries', 'fullms', 'weekly', 'modesmanifest', 'comp']},open(checkpoint,'wb'),protocol=5);del s,bins,template,pnsd,bs,st;gc.collect()
for name,rows in [('site_spatial_statistics',stats),('size_range_sensitivity',sens),('site_eligibility',elig),('breakpoints',breaks),('size_window_manifest',modesmanifest),('formula_verification',checks),('actual_primary_configuration_ranges',configs)]:pd.DataFrame(rows).to_csv(A/(name+'.csv'),index=False)
pd.DataFrame(sources).to_csv(O/'PNSD/source_provenance.csv',index=False);pd.DataFrame(channelqc).to_csv(O/'PNSD/channel_conversion_and_QC.csv',index=False);pd.DataFrame(configs).to_csv(O/'PNSD/actual_primary_configuration_ranges.csv',index=False);pd.concat(groups).to_csv(W/'trend_input_base.csv',index=False);pd.concat(mseries).to_csv(A/'monthly_series.csv',index=False);pd.concat(fullms).to_csv(A/'full_monthly_series.csv',index=False);pd.concat(weekly).to_pickle(W/'weekly.pkl')
policy={'ancillary_exclusions':{'HEL_TR':{'SO2':'Nearby station measurement; excluded from all analyses and model inputs'}},'source':'Original native PNSD Excel files, including BCN observations and full DRE_TR record; source hashes recorded','primary_nm':[17,500],'site_override_nm':{'LND_TR':[31.65,500]},'QC':'Step1 only: nonnumeric/nonfinite/negative/exact sentinel missing; retain zero and all other finite positive values; no statistical cutoff','sentinel_values':flags,'integration':'Partial sum over available bins; all-missing stays NA; valid zero retained; no interpolation or renormalization','diameter_and_DF':'Use native diameters and verified diameter-to-deposition-fraction mappings; disjoint mode boundaries','breakpoints':'Single linear-kernel change on raw observation sequence, min_size2; retain >=20 valid daily-based months eligibility; calendar-month boundary for phase fits, month shown on figures','full_and_phase_trends':'Direct arithmetic monthly means, no additional monthly completeness criterion; phases split at the start of the breakpoint month','meta':'openair::timeAverage monthly data.thresh=10; require >=60 finite valid months and usable CI','paired_comparison':'Equally weighted daily means aggregated by month >=14 valid days; >=60 exact paired months and uniform >=50% coverage; sensitivity >=36 months','weekly':'Direct weekly mean; no 126-hour filter','PRA':'Native configurations computed with their own diameters, widths and DF, then joined chronologically as one physical-site record; no separate configuration effect sizes','PAR':'Read from the provider workbook with both explicit native headers and matching PAR-D deposition fractions; configurations joined as one physical-site record','limitations':'Partial sums underestimate complete-window concentration when missing channels have nonzero contribution; global sentinel list not source-by-source validated'}
(O/'analysis_policy.json').write_text(json.dumps(policy,ensure_ascii=False,indent=2));(O/'PNSD/QC_protocol.json').write_text(json.dumps({k:policy[k] for k in ['source','QC','sentinel_values','limitations']},ensure_ascii=False,indent=2));print('PREPARED',len(stats),'independent records',flush=True)
