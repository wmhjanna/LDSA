from pathlib import Path
import shutil,json
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;O=R/'outputs/quality_control';B=Path('@DELIVERY_ROOT@/provider_quality_control_analysis')
for p in [W/'workbook',W/'workbook/series',W/'workbook/pnsd',W/'workbook/meta_raw',W/'netcdf',W/'network',O/'workbook/analysis',O/'workbook/PNSD',O/'workbook/sites',O/'netcdf',O/'analysis',B/'PNG',B/'tables',B/'validation']:p.mkdir(parents=True,exist_ok=True)
for name in ['prepare.py','joined_sources.py','meta_monthly.R','finalize_reference_inputs.py','restore_original_sequence.py','trends.R','sensitivity.R','final_tables.py','model.py','redraw.py']:
 s=(R/'work/workbook_processing'/name).read_text().replace('work/workbook_processing','work/quality_control/workbook').replace('outputs/workbook_processing','outputs/quality_control/workbook')
 if name=='prepare.py':
  old="bins=np.concatenate(bs);struct=np.concatenate(st);tt=pd.concat(ts,ignore_index=True);sq['rows']=len(tt);qc.append(sq)"
  new=old+"\n masks=pickle.load(open(ROOT/'work/quality_control/old_qc_masks.pkl','rb'));mask=masks.get(site,np.zeros(bins.shape,bool));assert mask.shape==bins.shape\n if mask.any():\n  prev=pickle.load(open(ROOT/'work/workbook_processing/pnsd'/(site+'.pkl'),'rb'));assert np.allclose(bins,prev['bins'],equal_nan=True,rtol=1e-12,atol=1e-9);assert tt.equals(prev['aux'].Time.reset_index(drop=True));bins[mask]=np.nan\n sq['official_invalid_cells_removed']=int(mask.sum());sq['official_invalid_rows_affected']=int(mask.any(axis=1).sum())"
  assert old in s;s=s.replace(old,new)
  # Use a single canonical pollutant input including German observations and the IPR observations.
  s=s.replace('from current_pollutants import apply_current_pollutants','from qc_pollutants import apply_current_pollutants')
 if name=='restore_original_sequence.py':s=s.replace('from current_pollutants import','from qc_pollutants import')
 (W/'workbook'/name).write_text(s)
# module available via the already used ROOT/work path; no existing module modified
(R/'work/qc_pollutants.py').write_text('''from pathlib import Path
import pandas as pd
from current_pollutants import load_site as german_load
R=Path('@PROJECT_ROOT@')
def load_site(site):
 if site!='IPR_RB':return german_load(site)
 rows=[]
 for metric in ['PM2.5','NO2']:
  x=pd.read_csv(R/'outputs/ipr_processing'/f'IPR_RB_{metric}_observations.csv.gz',parse_dates=['date']);x['metric']=metric;rows.append(x)
 return pd.concat(rows,ignore_index=True)
def apply_current_pollutants(site,raws):
 d=load_site(site)
 if d is not None:
  for name,g in d.groupby('metric'):raws[name]=g[['date','value']].dropna().reset_index(drop=True)
 return raws
''')
s=(R/'work/network_trends/recalculate/calculate.py').read_text()
s=s.replace("O=ROOT/'outputs/network_trends/corrected_LDSA'","O=ROOT/'outputs/quality_control/netcdf'")
s=s.replace("(W/'workbook_payload.json').write_text", "(ROOT/'work/quality_control/netcdf/workbook_payload.json').write_text")
s=s.replace("registry=[];", "ft=pd.read_csv(ROOT/'work/provider_flag_audit/ebas_flags.csv',header=None,quotechar=\"'\",names=['flag','description','status']);invalid_codes=ft[ft.status.isin(['I','H','M'])].flag.to_numpy()\nregistry=[];")
needle="raw=np.asarray(np.ma.filled(v[:],np.nan),float).T[mask];hasraw=np.isfinite(raw).any(axis=1)"
s=s.replace(needle,needle+"\n            q=ds[rec['variable']+'_qc'];qa=np.asarray(np.ma.filled(q[:],0));qa=np.transpose(qa,[q.dimensions.index('time'),q.dimensions.index(da)]+[i for i,d in enumerate(q.dimensions) if d not in ['time',da]])[mask];official_bad=np.isin(qa,invalid_codes).any(axis=2)&~(qa==100).any(axis=2)")
s=s.replace("bad=~np.isfinite(rr)|(rr<0)|np.isin(rr,flags)","bad=~np.isfinite(rr)|(rr<0)|np.isin(rr,flags)|official_bad[:,idx][active]")
s=s.replace("mat=weights[:,None]*np.pi", "integration_weights=np.ones_like(weights) if (code=='FI0096G' and len(registry)+1==13) else weights\n                mat=integration_weights[:,None]*np.pi")
s=s.replace("clean[rr==0]==0", "clean[(rr==0)&~bad]==0")
s=s.replace("z*weights[select]", "z*integration_weights[select]")
s=s.replace("negative_cells=int(", "official_invalid_cells=int((official_bad[:,idx][active]&np.isfinite(rr)).sum()),negative_cells=int(")
(W/'netcdf/calculate.py').write_text(s)
(W/'scope.json').write_text(json.dumps({'network':65,'LDSA_methods':'unchanged','netcdf_flags':'per-channel official I/H/M excluded; 100 override; V/unknown kept','workbook_flags':'only spectrum-confirmed native-cell matches','PAL':'source13 interpreted as DeltaN','coanalysis_model':'retain established workbook cohort; PAR/PRA omitted only in conventional pollutant comparison','pollutants':'keep other-database observations; use accepted German EBAS updates and IPR update; no HEL_TR SO2','figures':'PNG, Times New Roman bold, established layout'},ensure_ascii=False,indent=2))
print(B)
