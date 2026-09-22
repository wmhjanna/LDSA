from pathlib import Path
import sys,subprocess,pandas as pd,numpy as np,json,shutil
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;O=R/'outputs/record_exclusions';PY=sys.executable

def run(cmd,name):
 print(name,flush=True)
 with (W/(name+'.log')).open('w') as f:subprocess.run(cmd,cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
run(['@RSCRIPT@',str(W/'network/monthly.R')],'monthly')
run(['@RSCRIPT@',str(W/'workbook/trends.R'),str(W/'network/trend_input.csv'),str(O/'netcdf/analysis')],'new_trends')
run([PY,str(W/'merge_network.py')],'merge')
policy=json.loads((O/'analysis_policy.json').read_text());policy['user_case_exclusions']=json.loads((W/'exclusions.json').read_text());policy['user_case_exclusion_note']='Manual study-specific interval exclusions applied to corresponding native LDSA/PNC records before aggregation, separately from provider quality control.';(O/'analysis_policy.json').write_text(json.dumps(policy,ensure_ascii=False,indent=2))
P=W/'pollutants';P.mkdir(exist_ok=True);oldP=R/'work/netcdf_pollutants';A=O/'pollutants/analysis'
(P/'selected.pkl').symlink_to(oldP/'selected.pkl') if not (P/'selected.pkl').exists() else None
s=(oldP/'build_inputs.py').read_text().replace('quality_control','record_exclusions').replace('outputs/netcdf_pollutants','outputs/record_exclusions/pollutants');(P/'build_inputs.py').write_text(s)
run([PY,str(P/'build_inputs.py')],'poll_inputs')
s=(oldP/'fit.py').read_text().replace('quality_control','record_exclusions').replace('outputs/netcdf_pollutants','outputs/record_exclusions/pollutants');(P/'fit.py').write_text(s)
d=pd.read_csv(P/'new_trend_input.csv');d=d[d.Phase.str.startswith('Overlap_')];prior=pd.read_csv(oldP/'new_trend_input.csv');prior={k:g.sort_values('date').reset_index(drop=True) for k,g in prior.groupby('group_id')};same=[]
for k,g in d.groupby('group_id'):
 a=prior[k];b=g.sort_values('date').reset_index(drop=True)
 if a.date.equals(b.date) and np.array_equal(a.value,b.value,equal_nan=True):same.append(k)
old=pd.read_csv(oldP/'new_overlap_trends.csv');keys=sorted(d.group_id.unique());keepkeys=old.Site+'|'+old.Metric+'|'+old.Phase
for i in range(3):
 p=P/f'fit_{i}';p.mkdir(exist_ok=True);old[keepkeys.isin(set(keys[i::3])&set(same))].to_csv(p/'trends.csv',index=False)
(P/'reuse.json').write_text(json.dumps({'total':len(keys),'same':len(same),'refit':sorted(set(keys)-set(same))},indent=2))
run([PY,str(P/'fit.py')],'poll_fit')
s=(oldP/'render_fig3.py').read_text().replace('outputs/netcdf_pollutants','outputs/record_exclusions/pollutants').replace('netcdf_pollutant_figures','record_exclusion_analysis/PNG');(P/'render_fig3.py').write_text(s)
run([PY,str(P/'render_fig3.py')],'fig3')
run([PY,str(W/'layout/render.py')],'fig2_s5')
run([PY,str(W/'render_s6.py')],'s6')
print('DONE_CORE',flush=True)
