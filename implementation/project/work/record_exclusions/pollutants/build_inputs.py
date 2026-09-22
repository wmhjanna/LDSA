from pathlib import Path
import sys,json,shutil,ast
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent
sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import numpy as np,pandas as pd
from analyze_cleaned_ldsa import monthly
BASE=R/'outputs/record_exclusions';OLD=BASE/'analysis';A=R/'outputs/record_exclusions/pollutants/analysis'
A.mkdir(parents=True,exist_ok=True)
for p in OLD.iterdir():
 if p.is_file():shutil.copy2(p,A/p.name)
for node in ast.parse((R/'work/record_exclusions/workbook/restore_original_sequence.py').read_text()).body:
 if isinstance(node,ast.FunctionDef) and node.name=='index_cp':exec(compile(ast.Module(body=[node],type_ignores=[]),'reference_index_cp','exec'))
raw=pd.read_pickle(W/'selected.pkl');ms=pd.read_csv(OLD/'monthly_series.csv',parse_dates=['date']);fm=pd.read_csv(OLD/'full_monthly_series.csv',parse_dates=['date']);bps=pd.read_csv(OLD/'breakpoints.csv');el=pd.read_csv(OLD/'site_eligibility.csv').set_index('Site')
ms['date']=pd.to_datetime(ms.date,format='mixed');fm['date']=pd.to_datetime(fm.date,format='mixed')
parts=[];mo=[];full=[];bp=[];availability=[];overlap=[]
for (site,metric),d in raw.groupby(['Site','metric']):
 d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=['date','value']).sort_values('date').drop_duplicates('date').reset_index(drop=True)
 if d.empty:continue
 m=monthly(d.date,d.value);m['Site']=site;m['Metric']=metric;mo.append(m)
 anchor=pd.to_datetime(el.loc[site,'breakpoint'],errors='coerce');n=len(m);k=index_cp(d,n)
 if k is not None:
  bt=d.date.iloc[k];bm=bt.to_period('M').start_time
  bp.append(dict(Site=site,Type=site.rsplit('_',1)[1],Physical_site=site,Metric=metric,Breakpoint=bm,Breakpoint_timestamp=bt,n_months=n,offset_months=(bm.year-anchor.year)*12+bm.month-anchor.month if pd.notna(anchor) else np.nan))
 availability.append(dict(Site=site,Metric=metric,n_observations=len(d),first=d.date.min(),last=d.date.max(),valid_daily_months=n))
 for phase in ['Full','Before','After']:
  if phase!='Full' and pd.isna(anchor):continue
  q=d if phase=='Full' else d[d.date.lt(anchor) if phase=='Before' else d.date.ge(anchor)]
  g=q.groupby(q.date.dt.to_period('M').dt.to_timestamp()).value.mean().rename_axis('date').reset_index()
  if phase=='Full':z=g.copy();z['Site']=site;z['Metric']=metric;full.append(z)
  if len(g)<3:continue
  g['Site']=site;g['Type']=site.rsplit('_',1)[1];g['Metric']=metric;g['Phase']=phase;g['group_id']=site+'|'+metric+'|'+phase;parts.append(g)
 ld=ms[ms.Site.eq(site)&ms.Metric.eq('LDSA')][['date','value']];pair=ld.merge(m[['date','value']],on='date',suffixes=('_ldsa','_other')).sort_values('date');pn=len(pair)
 span=(pair.date.max().year-pair.date.min().year)*12+pair.date.max().month-pair.date.min().month+1 if pn else 0
 overlap.append(dict(Site=site,Metric=metric,common_months=pn,span_months=span,coverage=pn/span if span else 0,start=pair.date.min(),end=pair.date.max(),primary_included=pn>=60 and pn/span>=.5 if span else False,sensitivity_included=pn>=36 and pn/span>=.5 if span else False))
 if pn>=3:
  for role,col in [('LDSA','value_ldsa'),(metric,'value_other')]:
   g=pair[['date',col]].rename(columns={col:'value'});g['Site']=site;g['Type']=site.rsplit('_',1)[1];g['Metric']=role;g['Phase']='Overlap_'+metric;g['group_id']=site+'|'+role+'|Overlap_'+metric;parts.append(g)
new=pd.concat(parts,ignore_index=True);prior=pd.read_csv(R/'work/record_exclusions/trend_input.csv',parse_dates=['date'])
prior['date']=pd.to_datetime(prior.date,format='mixed')
assert set(new.group_id).isdisjoint(set(prior.group_id))
combined=pd.concat([prior,new],ignore_index=True);assert not combined.duplicated(['group_id','date']).any()
combined.to_csv(W/'trend_input.csv',index=False);new.to_csv(W/'new_trend_input.csv',index=False)
pd.concat([ms,*mo],ignore_index=True).to_csv(A/'monthly_series.csv',index=False)
pd.concat([fm,*full],ignore_index=True).to_csv(A/'full_monthly_series.csv',index=False)
pd.concat([bps,pd.DataFrame(bp)],ignore_index=True).to_csv(A/'breakpoints.csv',index=False)
pd.DataFrame(availability).to_csv(A/'new_pollutant_availability.csv',index=False);pd.DataFrame(overlap).to_csv(A/'new_pollutant_screening.csv',index=False)
policy=json.loads((BASE/'analysis_policy.json').read_text());policy['ancillary_comparison_scope']='All 65 network sites screened where downloaded matching native pollutant data exist; PAR_SUB/PRA_SUB excluded from conventional comparison as before. HEL_TR SO2 excluded.'
policy['new_pollutant_source']='Archived paper65 EBAS download; netcdf only, 1996-2023. No workbook records replaced.'
policy['new_pollutant_priority']='Exact total component, time mean; highest data level, shortest native resolution, newest revision, provider mass unit copy, ordinary over intensive campaign copy. Overlapping lower-ranked intervals excluded even if preferred observation invalid; no resurrection of flagged observations.'
policy['new_pollutant_dimensions']='Use one provided eBC channel: 880 nm preferred, else closest provided wavelength; minimum inlet height per configuration for tower records. Record every configuration; no cross-wavelength absorption-to-mass conversion.'
policy['new_pollutant_units']='BC ng m-3, CO mg m-3, remaining ug m-3. Elemental N/S converted to molecular mass. Molar-only gases converted at source standard T/P, or stated project 293.15 K/1013.25 hPa if absent; not ambient.'
policy['new_pollutant_QC']='Native observation flags: I/H/M invalid unless 100 override; zero and V warnings retained. Finite negative invalid; no statistical outlier cutoff. Sampling interval retained without upsampling.'
(A.parent/'analysis_policy.json').write_text(json.dumps(policy,ensure_ascii=False,indent=2))
print('NEW FITS',new.group_id.nunique(),'TOTAL',combined.group_id.nunique(),flush=True)
print(pd.DataFrame(overlap).groupby('Metric')[['primary_included','sensitivity_included']].sum().to_string(),flush=True)
