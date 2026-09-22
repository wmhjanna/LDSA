from pathlib import Path
import subprocess,time,os
w=Path('work/workbook_processing');a=Path('outputs/workbook_processing/analysis');py='@PYTHON@'
while not (w/'prepare.log').read_text().rstrip().endswith('PREPARED input records') or not (w/'trend_input.csv').exists():time.sleep(2)
def run(args,name,env=None):
 with (w/(name+'.log')).open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True,env=env)
run([py,str(w/'finalize_reference_inputs.py')],'finalize_inputs')
env=os.environ.copy();env['DYLD_LIBRARY_PATH']='/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib'
f=(w/'model.log').open('w');m=subprocess.Popen([py,str(w/'model.py')],stdout=f,stderr=subprocess.STDOUT,env=env)
run(['@RSCRIPT@',str(w/'trends.R'),str(w/'trend_input.csv'),str(a)],'trends')
assert m.wait()==0;f.close()
run(['@RSCRIPT@',str(w/'sensitivity.R')],'sensitivity')
run([py,str(w/'final_tables.py')],'final_tables')
print('ANALYSES_DONE',flush=True)
