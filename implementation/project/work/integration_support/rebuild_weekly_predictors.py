from pathlib import Path
import sys,ast,json,hashlib
ROOT=Path('@PROJECT_ROOT@');sys.path[:0]=[str(ROOT/'work/ldsa_support/deps')]
import numpy as np,pandas as pd
W=ROOT/'work/integration_support';A=ROOT/'outputs/integration_support/analysis'
p=Path('@INPUT_ROOT@/data_sources/predictive_models/LG/weekly_observations.xlsx');book=pd.ExcelFile(p,engine='calamine');wk=pd.read_pickle(W/'weekly.pkl');features=['BC','CO','NO','NO2','O3','PM2.5','SO2','TEMP','RH','PRESS','RAD','WS','WD','RAIN'];rows=[]
for oldsite in book.sheet_names:
 site={'SCH_SUB':'SCH_RB','WAL_SUB':'WAL_RB','ZUG_SUB':'ZUG_RB','LAN':'LAN_UB'}.get(oldsite,oldsite)
 if site not in set(wk.Site) or not site.endswith(('_UB','_TR')):continue
 d=book.parse(oldsite);datecol=None;cols=[]
 for c in d:
  if str(c).lower().strip().startswith(('date','time')):datecol=c;continue
  name=str(c).split('(')[0].split('（')[0].strip()
  if datecol is None or name not in features or (site=='HEL_TR' and name=='SO2'):continue
  q=pd.DataFrame({'date':pd.to_datetime(d[datecol],errors='coerce'),name:pd.to_numeric(d[c],errors='coerce')}).replace([np.inf,-np.inf],np.nan).dropna()
  if q.empty:continue
  q.loc[q[name].isin([-999,-9999,-99999,-999999,99999,999999,9999999]),name]=np.nan;q=q.dropna()
  q.date=q.date.dt.to_period('W').dt.to_timestamp();cols.append(q.groupby('date')[name].mean())
 if not cols:continue
 predictors=pd.concat(cols,axis=1).reset_index();z=wk[wk.Site.eq(site)][['date','Site','Type','LDSA']].merge(predictors,on='date',how='inner')
 for name in features:
  if name not in z:z[name]=np.nan
 rows.append(z)
pd.concat(rows,ignore_index=True).to_pickle(W/'weekly_model_reference.pkl')

print('WEEKLY_PREDICTORS_REBUILT')
