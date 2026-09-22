from pathlib import Path
import subprocess,sys,os,time,json,shutil
import pandas as pd,numpy as np
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;OW=W/'workbook';A=R/'outputs/quality_control/workbook/analysis';PY=sys.executable

def run(args,log):
 print('RUN',log,flush=True)
 with (W/(log+'.log')).open('w') as f:subprocess.run(args,cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True,env=os.environ)
while True:
 t=(OW/'prepare.log').read_text()
 if 'PREPARED input records' in t:break
 if 'Traceback (most recent call last)' in t:raise RuntimeError(t[-2000:])
 time.sleep(2)
run(['@RSCRIPT@',str(OW/'meta_monthly.R')],'meta_monthly')
run([PY,str(OW/'finalize_reference_inputs.py')],'finalize')
run([PY,str(OW/'restore_original_sequence.py')],'restore')
old=pd.read_csv(R/'work/ipr_processing/trend_input.csv');new=pd.read_csv(OW/'trend_input.csv');saved=pd.read_csv(R/'outputs/ipr_processing/analysis/trends.csv');prior={k:g.sort_values('date').reset_index(drop=True) for k,g in old.groupby('group_id')};same=[]
for k,g in new.groupby('group_id'):
 a=prior.get(k);b=g.sort_values('date').reset_index(drop=True)
 if a is not None and len(a)==len(b) and a.date.equals(b.date) and np.array_equal(a.value.to_numpy(),b.value.to_numpy(),equal_nan=True):same.append(k)
keep=(saved.Site+'|'+saved.Metric+'|'+saved.Phase).isin(same);saved[keep].to_csv(A/'trends.csv',index=False)
bad=pd.read_csv(R/'outputs/workbook_processing/analysis/invalid_percent_CI_raw.csv');bad[(bad.Site+'|'+bad.Metric+'|'+bad.Phase).isin(same)].to_csv(A/'invalid_percent_CI_raw.csv',index=False)
(A/'fit_reuse_verification.json').write_text(json.dumps({'rebuilt_input_groups':new.group_id.nunique(),'exactly_unchanged':int(keep.sum()),'refit':new.group_id.nunique()-int(keep.sum()),'same_keys':same},indent=2))
env=os.environ.copy();env['DYLD_LIBRARY_PATH']='/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib'
f=(W/'model.log').open('w');model=subprocess.Popen([PY,str(OW/'model.py')],cwd=R,stdout=f,stderr=subprocess.STDOUT,env=env)
run(['@RSCRIPT@',str(OW/'trends.R'),str(OW/'trend_input.csv'),str(A)],'workbook_trends')
assert model.wait()==0;f.close()
run(['@RSCRIPT@',str(OW/'sensitivity.R')],'workbook_sensitivity')
run([PY,str(OW/'final_tables.py')],'workbook_tables')
(W/'workbook_done').write_text('Complete')
print('workbook_COMPLETE',flush=True)
