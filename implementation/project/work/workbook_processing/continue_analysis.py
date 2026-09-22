from pathlib import Path
import subprocess,os,sys,time,json
import pandas as pd,numpy as np
R=Path('@PROJECT_ROOT@');W=R/'work/workbook_processing';A=R/'outputs/workbook_processing/analysis';PY='@PYTHON@'
def run(args,name,env=None):
 print('RUN',name,flush=True)
 with (W/(name+'.log')).open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True,env=env,cwd=R)
while not (W/'prepare.log').read_text().rstrip().endswith('PREPARED input records'):
 if 'Traceback (most recent call last)' in (W/'prepare.log').read_text():raise RuntimeError('prepare failed')
 time.sleep(2)
run(['@RSCRIPT@',str(W/'meta_monthly.R')],'meta_monthly')
run([PY,str(W/'finalize_reference_inputs.py')],'finalize_inputs')
run([PY,str(W/'restore_original_sequence.py')],'restore_sequence')
old=pd.read_csv(R/'work/integration_support/trend_input.csv');new=pd.read_csv(W/'trend_input.csv');saved=pd.read_csv(R/'outputs/integration_support/analysis/trends.csv');identical=[]
for k,g in new.groupby('group_id'):
 a=old[old.group_id.eq(k)].sort_values('date').reset_index(drop=True);b=g.sort_values('date').reset_index(drop=True)
 if len(a)==len(b) and a.date.equals(b.date) and np.array_equal(a.value.to_numpy(),b.value.to_numpy(),equal_nan=True):identical.append(k)
keep=(saved.Site+'|'+saved.Metric+'|'+saved.Phase).isin(identical);saved[keep].to_csv(A/'trends.csv',index=False)
invalid=pd.read_csv(R/'outputs/integration_support/analysis/invalid_percent_CI_raw.csv');invalid[(invalid.Site+'|'+invalid.Metric+'|'+invalid.Phase).isin(identical)].to_csv(A/'invalid_percent_CI_raw.csv',index=False)
(A/'fit_reuse_verification.json').write_text(json.dumps({'unchanged_input_groups':len(identical),'total_groups':new.group_id.nunique(),'changed_or_new':sorted(set(new.group_id)-set(identical))},indent=2))
env=os.environ.copy();env['DYLD_LIBRARY_PATH']='/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib'
f=(W/'model.log').open('w');m=subprocess.Popen([PY,str(W/'model.py')],stdout=f,stderr=subprocess.STDOUT,env=env,cwd=R)
run(['@RSCRIPT@',str(W/'trends.R'),str(W/'trend_input.csv'),str(A)],'trends')
assert m.wait()==0;f.close()
run(['@RSCRIPT@',str(W/'sensitivity.R')],'sensitivity')
run([PY,str(W/'final_tables.py')],'final_tables')
run([PY,str(W/'redraw.py')],'redraw')
print('ANALYSIS_AND_FIGURES_COMPLETE',flush=True)
