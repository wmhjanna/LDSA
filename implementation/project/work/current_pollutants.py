"""Authoritative pollutant update overlay; native LDSA/PNC and other sites stay unchanged."""
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def policy():
 p=ROOT/'outputs/current_pollutants.json'
 return json.loads(p.read_text()) if p.exists() else {}
def load_site(site):
 p=policy()
 if site not in p.get('sites',[]):return None
 return pd.read_pickle(Path(p['directory'])/'analysis_inputs'/f'{site}.pkl')
def apply_current_pollutants(site,raws):
 d=load_site(site)
 if d is None:return raws
 for m,g in d.groupby('metric'):
  raws[m]=g[['date','value']].sort_values('date').reset_index(drop=True)
 return raws
