from pathlib import Path
import subprocess,json,re,time,hashlib,shutil,collections
W=Path('work/ebas_download');O=Path('@DELIVERY_ROOT@/EBAS_observations');attempts=collections.Counter()
def recover(r):
 dest=O/r['relative_path'];part=dest.with_suffix('.nc.part')
 cp=subprocess.run(['curl','--http1.1','--fail','--silent','--show-error','--head','--connect-timeout','25','--max-time','60',r['url']],capture_output=True)
 sizes=re.findall(r'content-length:\s*(\d+)',cp.stdout.decode(),re.I)
 if cp.returncode or not sizes:raise RuntimeError('HEAD failed '+cp.stderr.decode())
 total=int(sizes[-1]);have=part.stat().st_size if part.exists() else 0
 for i in range(20):
  if have==total:break
  en=min(total-1,have+32*1024**2-1);q=W/'recovery_chunks'/(dest.name+'.next');h=q.with_suffix('.headers');q.unlink(missing_ok=True)
  cp=subprocess.run(['curl','--http1.1','--fail','--location','--silent','--show-error','--connect-timeout','25','--max-time','180','--range',f'{have}-{en}','--dump-header',str(h),'--output',str(q),r['url']],capture_output=True)
  ranges=re.findall(r'content-range:\s*bytes (\d+)-(\d+)/(\d+)',h.read_text(),re.I)
  if q.exists() and q.stat().st_size and ranges and int(ranges[-1][0])==have and int(ranges[-1][2])==total:
   assert q.stat().st_size<=en-have+1
   with part.open('ab') as f,q.open('rb') as src:shutil.copyfileobj(src,f)
   have=part.stat().st_size;q.unlink()
  if cp.returncode:print('RANGE_RETRY',r['station_code'],i,cp.stderr.decode()[-150:],flush=True);time.sleep(2)
 assert have==total,(have,total)
 with part.open('rb') as f:assert f.read(8).startswith((b'CDF',b'\x89HDF'))
 part.replace(dest)
 with dest.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
 r.pop('error',None);r.update(status='downloaded',bytes=total,sha256=digest,recovered='HTTP byte ranges')
 with (O/'download_log.jsonl').open('a') as f:f.write(json.dumps(r,ensure_ascii=False)+'\n')
 print('RECOVERED',r['file'],total,flush=True)
while True:
 d={r['file']:r for r in map(json.loads,(O/'download_log.jsonl').read_text().splitlines())}
 pending=[r for r in d.values() if r['status']!='downloaded' and not r['file'].startswith('NP0001G.') and attempts[r['file']]<3]
 for r in pending:
  attempts[r['file']]+=1
  try:recover(r)
  except Exception as e:print('RECOVERY_ERROR',r['file'],str(e),flush=True)
 if (O/'download_manifest.json').exists() and not pending:break
 time.sleep(15)
print('RETRY_WORKER_DONE',dict(attempts),flush=True)
