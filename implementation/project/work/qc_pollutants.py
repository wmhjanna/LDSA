from pathlib import Path
import pandas as pd
from current_pollutants import load_site as german_load
R=Path('@PROJECT_ROOT@')
def load_site(site):
 cleaned=R/'work/quality_control/verification/pollutants'/f'{site}.pkl'
 if cleaned.exists():
  d=pd.read_pickle(cleaned);d['date']=pd.to_datetime(d['date'],format='mixed');return d
 if site!='IPR_RB':return german_load(site)
 rows=[]
 for metric in ['PM2.5','NO2']:
  x=pd.read_csv(R/'outputs/ipr_processing'/f'IPR_RB_{metric}_observations.csv.gz',parse_dates=['date']);x['date']=pd.to_datetime(x['date'],format='mixed');x['metric']=metric;rows.append(x)
 return pd.concat(rows,ignore_index=True)
def apply_current_pollutants(site,raws):
 d=load_site(site)
 if d is not None:
  for name,g in d.groupby('metric'):raws[name]=g[['date','value']].dropna().reset_index(drop=True)
 return raws
