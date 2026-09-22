from pathlib import Path
import subprocess,sys,time
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent
while not (W/'workbook_done').exists():
 if 'Traceback' in (W/'run_analysis.log').read_text():raise RuntimeError((W/'run_analysis.log').read_text()[-1500:])
 time.sleep(2)
def run(name):
 print('RUN',name,flush=True)
 with (W/(name+'.log')).open('w') as f:subprocess.run([sys.executable,str(W/(name+'.py'))],cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
run('merge_network')
# Independent image groups write distinct files in the same delivery directory.
f=(W/'render_remaining.log').open('w');p=subprocess.Popen([sys.executable,str(W/'render_remaining.py')],cwd=R,stdout=f,stderr=subprocess.STDOUT)
run('render_fig2');run('render_s3_s6')
assert p.wait()==0;f.close();(W/'figures_done').write_text('Complete')
print('MAIN_FIGURES_COMPLETE',flush=True)
