from pathlib import Path
import json,hashlib,shutil,zipfile,sys
sys.path.insert(0,'@PROJECT_ROOT@/work/ldsa_support/deps')
import numpy as np,pandas as pd
from prepare import R,W,O,METRICS,UNITS,no_overlap
BOOK=Path('@INPUT_ROOT@/data_sources/predictive_models/pollutant_data/all_sites/paired_date_value_measurements.xlsx')
SITES=json.load(open(O/'station_mapping.json'));new=pd.read_pickle(O/'EBAS_selected.pkl')
rev={'SCH_RB':'SCH_SUB','WAL_RB':'WAL_SUB','ZUG_RB':'ZUG_SUB'}
book=pd.ExcelFile(BOOK,engine='calamine');rows=[];summ=[];olds=[]
(O/'analysis_inputs').mkdir(exist_ok=True);(O/'site_csv').mkdir(exist_ok=True)
for s in SITES:
 site=s['site'];oldsite=rev.get(site,site);oldparts=[]
 if oldsite in book.sheet_names:
  raw=book.parse(oldsite)
  for i in range(0,len(raw.columns)-1,2):
   m=str(raw.columns[i+1]).strip()
   if m not in METRICS:continue
   q=pd.DataFrame({'date':pd.to_datetime(raw.iloc[:,i],errors='coerce'),'value':pd.to_numeric(raw.iloc[:,i+1],errors='coerce')}).replace([np.inf,-np.inf],np.nan).dropna()
   if q.empty:continue
   q=q.assign(metric=m,unit=UNITS[m],source_id='reference_other_database',source='other_databases',Site=site)
   q['end']=q.date+pd.Timedelta(microseconds=1);oldparts.append(q)
 old=pd.concat(oldparts,ignore_index=True) if oldparts else pd.DataFrame(columns=['date','end','value','metric','unit','source_id','source','Site'])
 olds.append(old);combined=[]
 for m in METRICS:
  a=new[new.code.eq(s['code'])&new.metric.eq(m)].copy()
  lo=pd.Timestamp(s['start']).floor('D');hi=pd.Timestamp(s['end']).ceil('D');a=a[a.date.ge(lo)&a.date.lt(hi)].copy();a['Site']=site;a['source']='EBAS'
  b=old[old.metric.eq(m)].copy();keep=no_overlap(b,a)
  z=pd.concat([a,keep],ignore_index=True).sort_values('date');combined.append(z)
  summ.append(dict(Site=site,EBAS_code=s['code'],Station=s['name'],Metric=m,Unit=UNITS[m],EBAS_n=len(a),EBAS_start=str(a.date.min()) if len(a) else '',EBAS_end=str(a.end.max()) if len(a) else '',Old_n=len(b),Old_retained=len(keep),Old_replaced=len(b)-len(keep),Updated_n=len(z),Status='merge_available_sources' if len(a) and len(keep) else 'EBAS_only' if len(a) else 'complementary_sources_only' if len(b) else 'no_available_observations'))
 out=pd.concat(combined,ignore_index=True);out.to_pickle(O/'analysis_inputs'/f'{site}.pkl');rows.append(out)
 # Paired date/value columns preserve independent sampling intervals; Excel-readable CSV.
 wide=[]
 for m in METRICS:
  q=out[out.metric.eq(m)][['date','value']].reset_index(drop=True);q.columns=[f'date_{m}',f'{m} ({UNITS[m]})'];wide.append(q)
 pd.concat(wide,axis=1).to_csv(O/'site_csv'/f'{site}.csv',index=False,encoding='utf-8-sig')
allout=pd.concat(rows,ignore_index=True);pd.concat(olds,ignore_index=True).to_pickle(O/'complementary_observations.pkl')
summary=pd.DataFrame(summ);summary.to_csv(O/'pollutant_source_summary.csv',index=False,encoding='utf-8-sig')
monthly=[];full=[];weekly=[]
for (site,m),q in allout.groupby(['Site','metric']):
 d=q.groupby(q.date.dt.floor('D')).value.mean();g=d.resample('MS').agg(['mean','count']);g=g[g['count']>=14]['mean'].rename('value').rename_axis('date').reset_index();g['Site']=site;g['Metric']=m;monthly.append(g)
 g=q.groupby(q.date.dt.to_period('M').dt.start_time).value.mean().rename_axis('date').reset_index();g['Site']=site;g['Metric']=m;full.append(g)
 g=q.groupby(q.date.dt.to_period('W').dt.start_time).value.mean().rename_axis('date').reset_index();g['Site']=site;g['Metric']=m;weekly.append(g)
pd.concat(monthly,ignore_index=True).to_csv(O/'analysis_inputs/monthly_series.csv',index=False)
pd.concat(full,ignore_index=True).to_csv(O/'analysis_inputs/full_monthly_series.csv',index=False)
pd.concat(weekly,ignore_index=True).to_csv(O/'analysis_inputs/weekly_pollutants.csv',index=False)
policy=dict(version='20260920',directory=str(O),sites=[s['site'] for s in SITES],metrics=METRICS,units=UNITS,reference_workbook=str(BOOK),reference_sha256=hashlib.file_digest(BOOK.open('rb'),'sha256').hexdigest(),new_priority='Use EBAS valid measurement intervals first; retain legacy observations outside these intervals. No new-data coverage means unchanged legacy series.',raw_time='EBAS UTC sample start/end; legacy dates retained as used in existing project, original timezone unverified',cleaning='EBAS published finite nonnegative arithmetic means; zero retained; QC codes retained; no statistical outlier filtering; no repeated expansion of daily samples',source_priority='Highest data level, shortest native interval, latest revision; lower-priority overlapping intervals excluded rather than double counted',gas_conversion='Use metadata units; ug N/m3 to molecular NO/NO2 and ug S/m3 to SO2; molar-ratio conversion uses source standard conditions when present, otherwise documented reference 293.15 K and 101325 Pa',analysis_status='Exports native, monthly and weekly pollutant inputs for downstream trend analysis.',reference_preservation='Other-database legacy values retained outside new valid intervals, including all unavailable EBAS station/metric series. Original workbook remains unchanged.')
(O/'update_policy.json').write_text(json.dumps(policy,ensure_ascii=False,indent=2))
active=R/'outputs/current_pollutants.json';backup=O/'reference_metadata';backup.mkdir(exist_ok=True)
if active.exists():shutil.copy2(active,backup/active.name)
active.write_text(json.dumps(policy,ensure_ascii=False,indent=2))
changed=summary[summary.EBAS_n.gt(0)]
text=f'''# Station-matched pollutant data

EBAS observations are combined with available complementary records for {changed.Site.nunique()} stations and {len(changed)} station–metric combinations.

Valid EBAS sampling intervals take precedence. Complementary observations outside those intervals are retained. Source timestamps, sampling intervals, quality flags and metadata are preserved in the input records.

BC is expressed in ng m⁻³, CO in mg m⁻³ and other mass concentrations in μg m⁻³. Elemental nitrogen and sulphur concentrations are converted to the corresponding molecular mass. Gas mixing-ratio conversions use source reference conditions where available, otherwise 293.15 K and 101325 Pa.

Monthly paired-comparison inputs require at least 14 valid daily means. Full-record monthly and weekly means are also exported. See update_policy.json for the processing definitions, download_manifest.json for file provenance and sources.json for source metadata.

Data service: https://thredds.nilu.no/thredds/catalog/ebas/catalog.html
'''
(O/'output_notes.md').write_text(text)
print(summary.groupby('Status').size().to_string());print('ACTIVE',active,'new rows',int(summary.EBAS_n.sum()),'legacy retained',int(summary.Old_retained.sum()),flush=True)
