from pathlib import Path
import sys,json,ast,hashlib
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;M=W/'model';M.mkdir(exist_ok=True);A=R/'outputs/network_analysis/analysis';BASE=R/'outputs/record_exclusions/analysis';sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import numpy as np,pandas as pd
from qc_pollutants import load_site
features=['BC','CO','NO','NO2','O3','PM2.5','SO2','TEMP','RH','PRESS','RAD','WS','WD','RAIN'];polfeatures=features[:7]
oldw=pd.read_pickle(R/'work/quality_control/workbook/weekly_model_reference.pkl');oldwk=pd.read_pickle(R/'work/record_exclusions/weekly.pkl');bps=pd.read_csv(BASE/'breakpoints.csv');mapping=json.loads((R/'work/station_registry/mapping.json').read_text());newsites={x['new'] for x in mapping};oldsites=set(pd.read_csv(R/'outputs/quality_control/workbook/analysis/site_eligibility.csv').Site)
rows=[];audit=[];changes=[];ranges=pd.read_csv(R/'work/quality_control/workbook/native_phase_ranges.csv');newranges=[]
# Retain available meteorology and combine German pollutant observations with the other-database records.
for site in list(oldw.Site.unique()):
 old=oldw[oldw.Site.eq(site)].copy();wk=oldwk[oldwk.Site.eq(site)][['date','Site','Type','LDSA']];z=wk.merge(old[['date']+features],on='date',how='left',validate='one_to_one');raw=load_site(site)
 if raw is not None:
  for metric,d in raw.groupby('metric'):
   if metric not in polfeatures:continue
   d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=['date','value']);g=d.groupby(d.date.dt.to_period('W').dt.start_time).value.mean();before=z[metric].copy();z[metric]=z.date.map(g);n=int((~np.isclose(before,z[metric],equal_nan=True)).sum());changes.append(dict(Site=site,Metric=metric,changed_week_values=n,source='qc_pollutants latest native observations, preserving merged external database records'))
 z.loc[z.Site.eq('HEL_TR'),'SO2']=np.nan;rows.append(z)
for x in mapping:
 site=x['new'];typ=x['type']
 if typ not in ['UB','TR']:continue
 wk=oldwk[oldwk.Site.eq(site)][['date','Site','Type','LDSA']];p=R/'work/netcdf_pollutants/native'/f'{site}.pkl';z=wk.copy()
 for f in features:z[f]=np.nan
 if p.exists():
  raw=pd.read_pickle(p)
  for metric,d in raw.groupby('metric'):
   if metric not in polfeatures:continue
   d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=['date','value']);z[metric]=z.date.map(d.groupby(d.date.dt.to_period('W').dt.start_time).value.mean())
 # Download inventory confirms no meteorology for these four UB stations; leave missing, never substitute RB stations.
 audit.append(dict(Site=site,available_weeks=int(z.LDSA.notna().sum()),weeks_any_predictor=int(z[features].notna().any(axis=1).sum()),available_features=';'.join(f for f in features if z[f].notna().any()),meteorology_downloaded=False))
 b=bps[bps.Site.eq(site)&bps.Metric.eq('LDSA')]
 if len(b):
  bt=pd.Timestamp(b.Breakpoint_timestamp.iloc[0]);d=pd.read_csv(R/'outputs/record_exclusions/netcdf/hourly_ldsa'/f"{x['old']}_LDSA.csv.gz",parse_dates=['Time_UTC']).dropna(subset=['Total_LDSA_um2_cm-3'])
  for phase,sel in [('Before',d.Time_UTC.lt(bt)),('After',d.Time_UTC.ge(bt))]:
   q=d[sel];newranges.append(dict(Site=site,Metric='LDSA',Phase=phase,start=q.Time_UTC.min(),end=q.Time_UTC.max()))
 rows.append(z)
combined=pd.concat(rows,ignore_index=True);assert not combined.duplicated(['Site','date']).any();combined.to_pickle(M/'weekly_model_reference.pkl');pd.concat([ranges,pd.DataFrame(newranges)],ignore_index=True).to_csv(M/'native_phase_ranges.csv',index=False);pd.DataFrame(audit).to_csv(M/'new_station_screening.csv',index=False);pd.DataFrame(changes).to_csv(M/'existing_station_pollutant_updates.csv',index=False)
met=pd.read_csv('@DELIVERY_ROOT@/netcdf_meteorology/netcdf_meteorology_summary.csv');met['Station_class']=met.Site.str.rsplit('_',n=1).str[-1];met['original_UB_TR_model_scope']=met.Station_class.isin(['UB','TR']);met.to_csv(M/'meteorology_scope.csv',index=False)
s=(R/'work/quality_control/workbook/model.py').read_text().replace("STAGE=ROOT/'work/quality_control/workbook'", "STAGE=ROOT/'work/network_analysis/model'").replace("OUT=ROOT/'outputs/quality_control/workbook'", "OUT=ROOT/'outputs/network_analysis'");(M/'model.py').write_text(s)
print('MODEL_INPUT',len(combined));print(pd.DataFrame(audit).to_string(index=False));print(pd.DataFrame(changes).query('changed_week_values>0').to_string(index=False));print(met.groupby('Station_class').ambient_series.sum())
