from pathlib import Path
import json,subprocess,concurrent.futures,hashlib,shutil,time
W=Path('work/ebas_download');O=Path('@DELIVERY_ROOT@/EBAS_observations')
records={r['file']:r for r in map(json.loads,(O/'download_log.jsonl').read_text().splitlines())}
r=next(r for r in records.values() if r['file'].startswith('NP0001G.') and r['status']=='failed');dest=O/r['relative_path'];part=dest.with_suffix('.nc.part');start=part.stat().st_size;total=1542479876;chunk=128*1024**2
chunks=[];(W/'recovery_chunks').mkdir(exist_ok=True)
for i,st in enumerate(range(start,total,chunk)):chunks.append((st,min(st+chunk,total)-1,W/'recovery_chunks'/f'NP0001G_{i:02d}.part'))
def get(t):
 st,en,p=t
 for attempt in range(25):
  have=p.stat().st_size if p.exists() else 0
  if have==en-st+1:print('CHUNK',p.name,have,flush=True);return p
  assert have<en-st+1
  nextend=min(en,st+have+32*1024**2-1);q=p.with_suffix('.next')
  cp=subprocess.run(['curl','--http1.1','--fail','--location','--silent','--show-error','--connect-timeout','30','--max-time','300','--range',f'{st+have}-{nextend}','--output',str(q),r['url']],capture_output=True)
  if q.exists() and q.stat().st_size:
   assert q.stat().st_size<=nextend-st-have+1
   with p.open('ab') as f,q.open('rb') as src:shutil.copyfileobj(src,f)
   q.unlink()
  if cp.returncode:print('CHUNK_RETRY',p.name,attempt,cp.stderr.decode()[-200:],flush=True);time.sleep(3)
 raise RuntimeError(p)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:got=list(ex.map(get,chunks))
with part.open('ab') as f:
 for p in got:
  with p.open('rb') as src:shutil.copyfileobj(src,f,1024*1024)
assert part.stat().st_size==total
part.replace(dest)
with dest.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
r.pop('error',None);r.update(status='downloaded',bytes=total,sha256=digest,recovered='HTTP byte ranges after transfer interruption')
with (O/'download_log.jsonl').open('a') as f:f.write(json.dumps(r,ensure_ascii=False)+'\n')
print('RECOVERED',dest.name,total,flush=True)
