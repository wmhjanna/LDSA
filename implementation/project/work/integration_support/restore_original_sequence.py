from pathlib import Path
import sys,json
R=Path('@PROJECT_ROOT@');sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import numpy as np,pandas as pd
from analyze_cleaned_ldsa import monthly
W=R/'work/integration_support';A=R/'outputs/integration_support/analysis'
native={'LDSA':'LDSA','Ntotal':'PNC_Ntotal_primary','Nucleation':'PNC_Nucleation_primary','Aitken':'PNC_Aitken_primary','Accumulation':'PNC_Accumulation_primary'}
oldinput=pd.read_csv(W/'trend_input.csv',parse_dates=['date']);ms=pd.read_csv(A/'monthly_series.csv',parse_dates=['date']);oldms=pd.read_csv(R/'outputs/pnsd_cleaning/analysis/monthly_series.csv',parse_dates=['date']);e=pd.read_csv(A/'site_eligibility.csv');prevbp=pd.read_csv(A/'breakpoints.csv');book=pd.ExcelFile('@INPUT_ROOT@/data_sources/predictive_models/pollutant_data/all_sites/paired_date_value_measurements.xlsx',engine='calamine')
rev={'SCH_RB':'SCH_SUB','WAL_RB':'WAL_SUB','ZUG_RB':'ZUG_SUB'};parts=[];bounds=[];bps=[];full=[]
def clean(d):return d.replace([np.inf,-np.inf],np.nan).dropna(subset=['date','value']).reset_index(drop=True)
def index_cp(d,n):
 if n<20 or len(d)<4 or d.value.std()==0:return None
 y=d.value.to_numpy(float);y=y-y.mean();cs=np.cumsum(y);k=np.arange(2,len(y)-1);return int(k[np.argmax(cs[k-1]**2/k+(cs[-1]-cs[k-1])**2/(len(y)-k))])
for site in e.Site:
 s=pd.read_pickle(W/'series'/(site+'_series.pkl'));raws={name:clean(pd.DataFrame({'date':s.Time,'value':s[col]})) for name,col in native.items()};polsite='PRA_SUB' if site.startswith('PRA_') else rev.get(site,site)
 if polsite in book.sheet_names:
  pol=book.parse(polsite)
  for i in range(0,len(pol.columns)-1,2):
   name=str(pol.columns[i+1]).strip()
   if name not in ['PM2.5','PM10','PM1','BC','NO2','SO2','O3','CO','NO']:continue
   d=clean(pd.DataFrame({'date':pd.to_datetime(pol.iloc[:,i],errors='coerce'),'value':pd.to_numeric(pol.iloc[:,i+1],errors='coerce')}))
   if site.startswith('PRA_'):
    y0,y1=map(int,site.split('_')[1].split('-'));d=d[d.date.dt.year.between(y0,y1)]
   raws[name]=d.reset_index(drop=True)
 if site=='LAN_UB':
  lan=pd.read_excel('@INPUT_ROOT@/data_sources/LAN/LANbreakpoint/data/LAN.xlsx',engine='calamine')
  lan_months=[]
  for name in lan.columns[1:]:
   raws[name]=clean(pd.DataFrame({'date':pd.to_datetime(lan.date),'value':pd.to_numeric(lan[name],errors='coerce')}));z=monthly(raws[name].date,raws[name].value);z['Site']=site;z['Metric']=name;lan_months.append(z)
  ms=pd.concat([ms[~(ms.Site.eq(site)&~ms.Metric.isin(native))]]+lan_months,ignore_index=True)
 for name,g in oldms[oldms.Site.eq(site)&~oldms.Metric.isin(raws)].groupby('Metric'):raws[name]=clean(g[['date','value']])
 if site=='HEL_TR':raws.pop('SO2',None)
 
 if site=='LAN_UB':
  for name in native:raws[name]=raws[name].sort_values('date',kind='stable').reset_index(drop=True)
 ns=ms[ms.Site.eq(site)].groupby('Metric').size();indices={name:index_cp(d,int(ns.get(name,0))) for name,d in raws.items()};ai=indices['LDSA'];anchor=raws['LDSA'].date.iloc[ai] if ai is not None else None;am=anchor.to_period('M').start_time if anchor is not None else None
 e.loc[e.Site.eq(site),'breakpoint']=str(am) if am is not None else np.nan;e.loc[e.Site.eq(site),'breakpoint_timestamp']=str(anchor) if anchor is not None else np.nan
 for name,d in raws.items():
  k=indices[name];b=d.date.iloc[k] if k is not None else None;bm=b.to_period('M').start_time if b is not None else None
  if b is not None:bps.append(dict(Site=site,Type=site.rsplit('_',1)[1],Physical_site='PRA_SUB' if site.startswith('PRA_') else site,Metric=name,Breakpoint=bm,Breakpoint_timestamp=b,n_months=int(ns.get(name,0)),offset_months=(bm.year-am.year)*12+bm.month-am.month if am is not None else np.nan))
  for phase in ['Full','Before','After']:
   if phase=='Full':rr=d.drop_duplicates('date',keep='first')
   elif name in native:
    if k is None:continue
    rr=d.iloc[:k] if phase=='Before' else d.iloc[k:]
   else:
    if am is None:continue
    cut=anchor if site=='LAN_UB' else am
    rr=d[d.date.lt(cut) if phase=='Before' else d.date.ge(cut)].drop_duplicates(['date','value'])
   g=rr.groupby(rr.date.dt.to_period('M').dt.to_timestamp()).value.mean().rename_axis('date').reset_index()
   if phase=='Full':z=g.copy();z['Site']=site;z['Metric']=name;full.append(z)
   if len(g)<3:continue
   if phase!='Full' and name in native:bounds.append(dict(Site=site,Metric=name,Phase=phase,start=rr.date.min(),end=rr.date.max()))
   g['Site']=site;g['Type']=site.rsplit('_',1)[1];g['Metric']=name;g['Phase']=phase;g['group_id']=site+'|'+name+'|'+phase;parts.append(g)
ms.to_csv(A/'monthly_series.csv',index=False)
ld=ms[ms.Site.eq('LAN_UB')&ms.Metric.eq('LDSA')][['date','value']]
for name,d in ms[ms.Site.eq('LAN_UB')&~ms.Metric.isin(native)].groupby('Metric'):
 pair=ld.merge(d[['date','value']],on='date',suffixes=('_ldsa','_other'))
 if len(pair)<3:continue
 for role,col in [('LDSA','value_ldsa'),(name,'value_other')]:
  z=pair[['date',col]].rename(columns={col:'value'});z['Site']='LAN_UB';z['Type']='UB';z['Metric']=role;z['Phase']='Overlap_'+name;z['group_id']='LAN_UB|'+role+'|Overlap_'+name;parts.append(z)
oldkeep=oldinput[~oldinput.Phase.isin(['Full','Before','After'])&~(oldinput.Site.eq('LAN_UB')&oldinput.Phase.str.startswith('Overlap_'))]
new=pd.concat([oldkeep]+parts,ignore_index=True);new.to_csv(W/'trend_input.csv',index=False);new[new.Phase.ne('MetaFull')].to_csv(W/'trend_input_base.csv',index=False);pd.concat(full).to_csv(A/'full_monthly_series.csv',index=False);e.to_csv(A/'site_eligibility.csv',index=False);pd.DataFrame(bps).to_csv(A/'breakpoints.csv',index=False);pd.DataFrame(bounds).to_csv(W/'native_phase_ranges.csv',index=False)
changed=[]
for key,g in new.groupby('group_id'):
 old=oldinput[oldinput.group_id.eq(key)].sort_values('date');g=g.sort_values('date')
 if len(g)!=len(old) or not np.array_equal(g.date.values,old.date.values) or not np.allclose(g.value,old.value,equal_nan=True,rtol=1e-13,atol=1e-10):changed.append(key)
(W/'changed_fit_ids.json').write_text(json.dumps(changed))
for name in ['trends.csv','invalid_percent_CI_raw.csv']:
 p=A/name
 if p.exists():
  t=pd.read_csv(p);t=t[~(t.Site+'|'+t.Metric+'|'+t.Phase).isin(changed)];t.to_csv(p,index=False)
prevbp.merge(pd.DataFrame(bps),on=['Site','Metric'],suffixes=('_sorted','_original')).to_csv(A/'original_sequence_breakpoint_audit.csv',index=False)
policy=json.loads((A.parent/'analysis_policy.json').read_text());policy['breakpoints']='Original retained observation row sequence (LAN multiple-source records are ordered chronologically to match its original integrated workbook), single linear KernelCPD min_size=2; preserve >=20 valid daily-based months. Native LDSA/PNC split at their own detected row index; synchronized pollutants split at the start of the LDSA breakpoint month.';policy['LAN_ancillary_source']='@INPUT_ROOT@/data_sources/LAN/LANbreakpoint/data/LAN.xlsx';policy['LAN_phase_boundary']='Use exact LDSA timestamp, as specified in the original LAN R script; paired analyses still use common valid calendar months.';policy['duplicate_timestamps']='Full-record trend: first original row per timestamp (original R distinct then arrange). Native phase trends retain original rows. Synchronized pollutant phases remove exact duplicate date/value rows (original sync drop_duplicates). Weekly means retain original rows.';(A.parent/'analysis_policy.json').write_text(json.dumps(policy,ensure_ascii=False,indent=2));print('OBSERVATION SEQUENCE PREPARED',len(changed),'changed fits',new.group_id.nunique(),'total')
