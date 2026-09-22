from pathlib import Path
import sys,json,ast,hashlib
R=Path('@PROJECT_ROOT@')
sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import numpy as np,pandas as pd
from analyze_cleaned_ldsa import monthly
W=R/'work/record_exclusions/network';A=R/'outputs/record_exclusions/netcdf/analysis'
OLD=R/'outputs/workbook_processing/analysis';B=Path('@DELIVERY_ROOT@')
for p in [A,W/'raw']:p.mkdir(parents=True,exist_ok=True)
for node in ast.parse((R/'work/workbook_processing/restore_original_sequence.py').read_text()).body:
 if isinstance(node,ast.FunctionDef) and node.name=='index_cp':exec(compile(ast.Module(body=[node],type_ignores=[]),'reference_index_cp','exec'))
mapping=json.loads((R/'work/station_registry/mapping.json').read_text())
parts=[];elig=[];bps=[];week=[];provenance=[]
for r in mapping:
 site=r['new'];typ=r['type'];code=r['old'];path=R/'outputs/record_exclusions/netcdf/hourly_ldsa'/f'{code}_LDSA.csv.gz'
 d=pd.read_csv(path,usecols=['Time_UTC','Total_LDSA_um2_cm-3'],parse_dates=['Time_UTC']).rename(columns={'Time_UTC':'date','Total_LDSA_um2_cm-3':'value'})
 assert d.date.is_monotonic_increasing and d.date.is_unique
 d=d.replace([np.inf,-np.inf],np.nan).dropna().reset_index(drop=True);assert len(d)>0
 d.to_csv(W/'raw'/f'{site}.csv.gz',index=False)
 mm=monthly(d.date,d.value);n=len(mm);k=index_cp(d,n)
 bp=d.date.iloc[k] if k is not None else None;bm=bp.to_period('M').start_time if bp is not None else None
 elig.append(dict(Site=site,Physical_site=site,EBAS_code=code,Type=typ,valid_months=np.nan,paired_months_14days=n,span_months=(mm.date.max().year-mm.date.min().year)*12+mm.date.max().month-mm.date.min().month+1 if n else 0,meta_eligible=False,breakpoint=bm,breakpoint_timestamp=bp,primary_low_nm=17,primary_high_nm=500))
 if bp is not None:bps.append(dict(Site=site,Type=typ,Physical_site=site,Metric='LDSA',Breakpoint=bm,Breakpoint_timestamp=bp,n_months=n,offset_months=0))
 for phase in ['Full','Before','After']:
  if phase!='Full' and k is None:continue
  q=d if phase=='Full' else (d.iloc[:k] if phase=='Before' else d.iloc[k:])
  g=q.groupby(q.date.dt.to_period('M').dt.to_timestamp()).value.mean().rename_axis('date').reset_index()
  if len(g)<3:continue
  g['Site']=site;g['Type']=typ;g['Metric']='LDSA';g['Phase']=phase;g['group_id']=site+'|LDSA|'+phase;parts.append(g)
 g=d.groupby(d.date.dt.to_period('W').dt.to_timestamp()).value.mean().rename('LDSA').rename_axis('date').reset_index();g['Site']=site;g['Type']=typ;week.append(g)
 provenance.append(dict(Site=site,EBAS_code=code,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),valid_observations=len(d),start=str(d.date.min()),end=str(d.date.max()),daily_based_months=n,breakpoint_row=k))
 print(site,len(d),'hours',n,'daily-based months','breakpoint',bp,flush=True)
pd.concat(parts,ignore_index=True).to_csv(W/'trend_input_base.csv',index=False)
pd.DataFrame(elig).to_csv(A/'site_eligibility.csv',index=False)
pd.DataFrame(bps).to_csv(A/'breakpoints.csv',index=False)
pd.concat(week,ignore_index=True).to_pickle(W/'weekly.pkl')
(W/'source_manifest.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2))
pd.DataFrame(mapping)[['old','new','type','station_name']].to_csv(A/'new_station_names.csv',index=False)
