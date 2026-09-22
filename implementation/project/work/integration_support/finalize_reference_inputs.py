from pathlib import Path
import sys,ast,json,hashlib
ROOT=Path('@PROJECT_ROOT@');sys.path[:0]=[str(ROOT/'work/ldsa_support/deps')]
import numpy as np,pandas as pd
W=ROOT/'work/integration_support';A=ROOT/'outputs/integration_support/analysis'
for f in ast.parse((W/'prepare.py').read_text()).body:
 if isinstance(f,ast.FunctionDef) and f.name=='rawmonthly':exec(compile(ast.Module(body=[f],type_ignores=[]),'helper','exec'))
bp=pd.read_csv(A/'breakpoints.csv',parse_dates=['Breakpoint_timestamp']);base=pd.read_csv(W/'trend_input_base.csv',parse_dates=['date']);native=['LDSA','Ntotal','Nucleation','Aitken','Accumulation'];base=base[~(base.Metric.isin(native)&base.Phase.isin(['Before','After']))];add=[];ranges=[]
for b in bp[bp.Metric.isin(native)].itertuples():
 s=pd.read_pickle(W/'series'/(b.Site+'_series.pkl'));col='LDSA' if b.Metric=='LDSA' else 'PNC_'+b.Metric+'_primary'
 for phase in ['Before','After']:
  raw=s[s.Time.lt(b.Breakpoint_timestamp) if phase=='Before' else s.Time.ge(b.Breakpoint_timestamp)][['Time',col]].dropna();g=rawmonthly(raw.Time,raw[col])
  if len(g)<3:continue
  ranges.append(dict(Site=b.Site,Metric=b.Metric,Phase=phase,start=raw.Time.min(),end=raw.Time.max()))
  g['Site']=b.Site;g['Type']=b.Type;g['Metric']=b.Metric;g['Phase']=phase;g['group_id']=b.Site+'|'+b.Metric+'|'+phase;add.append(g)
base=pd.concat([base]+add,ignore_index=True);base.to_csv(W/'trend_input_base.csv',index=False)
meta=pd.read_csv(A/'meta_monthly_series.csv',parse_dates=['date']);pd.concat([base,meta[base.columns]],ignore_index=True).to_csv(W/'trend_input.csv',index=False);pd.DataFrame(ranges).to_csv(W/'native_phase_ranges.csv',index=False)
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
policy=json.loads((A.parent/'analysis_policy.json').read_text());policy.update(breakpoints='Original observation-sequence linear KernelCPD, min_size=2; >=20 valid daily-based months. Native LDSA/PNC phases use their own exact detected timestamp; synchronized co-pollutants use the start of the LDSA breakpoint month.',full_and_phase_trends='Direct monthly means without an additional completeness gate; retain the distinct original native and synchronized pollutant phase boundaries.',weekly='Original independent weekly means of predictors and cleaned LDSA, inner joined by Monday. Phase start/end are the original native-phase month labels converted to Monday, inclusive at both ends; the boundary week may occur in both phases, as in the original script.',weekly_predictor_source=str(p),weekly_predictor_sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest())
(A.parent/'analysis_policy.json').write_text(json.dumps(policy,ensure_ascii=False,indent=2))
print('FINALIZED',base.group_id.nunique(),'trend base groups',len(rows),'weekly source sites',flush=True)
