from pathlib import Path
import sys,json,subprocess,os
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;M=W/'model';A=R/'outputs/network_analysis/analysis';sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import pandas as pd,numpy as np
from qc_pollutants import load_site
features=['BC','CO','NO','NO2','O3','PM2.5','SO2','TEMP','RH','PRESS','RAD','WS','WD','RAIN'];polfeatures=features[:7]
prior=pd.read_pickle(M/'weekly_model_reference.pkl');prior=prior[~prior.Type.eq('RB')];wk=pd.read_pickle(W/'weekly.pkl');bps=pd.read_csv(A/'breakpoints.csv');mapping={x['new']:x for x in json.loads((R/'work/station_registry/mapping.json').read_text())};met=pd.read_pickle(W/'meteorology/weekly.pkl');oldbook=Path('@INPUT_ROOT@/data_sources/predictive_models/LG/weekly_observations.xlsx');book=pd.ExcelFile(oldbook,engine='calamine');rev={'SCH_RB':'SCH_SUB','WAL_RB':'WAL_SUB','ZUG_RB':'ZUG_SUB'};rows=[prior];audit=[];ranges=pd.read_csv(M/'native_phase_ranges.csv');newranges=[]
for site in sorted(wk.loc[wk.Type.eq('RB'),'Site'].unique()):
 z=wk[wk.Site.eq(site)][['date','Site','Type','LDSA']].copy()
 for f in features:z[f]=np.nan
 oldsite=rev.get(site,site);oldmetrics=[];newmetrics=[];metmetrics=[]
 if site not in mapping and oldsite in book.sheet_names:
  d=book.parse(oldsite);datecol=None
  for c in d:
   if str(c).lower().strip().startswith(('date','time')):datecol=c;continue
   name=str(c).split('(')[0].split('（')[0].strip()
   if datecol is None or name not in features:continue
   q=pd.DataFrame({'date':pd.to_datetime(d[datecol],errors='coerce'),name:pd.to_numeric(d[c],errors='coerce')}).replace([np.inf,-np.inf],np.nan).dropna();q.loc[q[name].isin([-999,-9999,-99999,-999999,99999,999999,9999999]),name]=np.nan;q=q.dropna()
   if q.empty:continue
   z[name]=z.date.map(q.groupby(q.date.dt.to_period('W').dt.start_time)[name].mean());oldmetrics.append(name)
 raw=load_site(site) if site not in mapping else (pd.read_pickle(R/'work/netcdf_pollutants/native'/f'{site}.pkl') if (R/'work/netcdf_pollutants/native'/f'{site}.pkl').exists() else None)
 if raw is not None:
  for metric,q in raw.groupby('metric'):
   if metric not in polfeatures:continue
   q=q.replace([np.inf,-np.inf],np.nan).dropna(subset=['date','value']);z[metric]=z.date.map(q.groupby(q.date.dt.to_period('W').dt.start_time).value.mean());newmetrics.append(metric)
 for metric,q in met[met.Site.eq(site)].groupby('Metric'):
  z[metric]=z.date.map(q.set_index('date').value);metmetrics.append(metric)
 rows.append(z);audit.append(dict(Site=site,old_workbook_features=';'.join(oldmetrics),current_native_pollutants=';'.join(newmetrics),new_native_meteorology=';'.join(metmetrics),weeks_with_predictors=int(z[features].notna().any(axis=1).sum())))
 if site in mapping:
  b=bps[bps.Site.eq(site)&bps.Metric.eq('LDSA')]
  if len(b):
   bt=pd.Timestamp(b.Breakpoint_timestamp.iloc[0]);d=pd.read_csv(R/'outputs/network_analysis/netcdf/hourly_ldsa'/f"{mapping[site]['old']}_LDSA.csv.gz",parse_dates=['Time_UTC']).dropna(subset=['Total_LDSA_um2_cm-3'])
   for phase,sel in [('Before',d.Time_UTC.lt(bt)),('After',d.Time_UTC.ge(bt))]:
    q=d[sel];newranges.append(dict(Site=site,Metric='LDSA',Phase=phase,start=q.Time_UTC.min(),end=q.Time_UTC.max()))
combined=pd.concat(rows,ignore_index=True);assert not combined.duplicated(['Site','date']).any();combined.to_pickle(M/'weekly_model_reference.pkl');pd.concat([ranges,pd.DataFrame(newranges)],ignore_index=True).drop_duplicates(['Site','Metric','Phase'],keep='last').to_csv(M/'native_phase_ranges.csv',index=False);pd.DataFrame(audit).to_csv(M/'RB_source_screening.csv',index=False)
s=(M/'model.py').read_text().replace("['UB','TR']", "['UB','TR','RB']");(M/'model.py').write_text(s)
print(pd.DataFrame(audit).to_string(index=False),flush=True)
env=os.environ.copy();env['DYLD_LIBRARY_PATH']='/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib'
with (M/'model.log').open('w') as f:subprocess.run([sys.executable,str(M/'model.py')],cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
print('MODELS_DONE',flush=True)
