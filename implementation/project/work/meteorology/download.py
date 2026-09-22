from pathlib import Path
import json,subprocess,concurrent.futures,hashlib,shutil,threading,time
ROOT=Path('@PROJECT_ROOT@')
W=ROOT/'work/meteorology'
BASE=Path('@DELIVERY_ROOT@')
O=BASE/'netcdf_meteorology';O.mkdir(exist_ok=True)
for p in [W/'logs',O/'source_netcdf']:p.mkdir(exist_ok=True)
old=BASE/'EBAS_observations'
cache={r['file']:old/r['relative_path'] for r in json.loads((old/'download_manifest.json').read_text()) if r.get('status')=='downloaded'}
items=json.loads((W/'direct_candidates.json').read_text())
lock=threading.Lock();count=0
def fetch(r):
 global count
 name=r['name'];code=name.split('.')[0];p=O/'source_netcdf'/code/name;p.parent.mkdir(exist_ok=True)
 result={'station_code':code,'file':name,'url':'https://thredds.nilu.no/thredds/fileServer/'+r['urlPath'],'relative_path':str(p.relative_to(O))}
 try:
  if not p.exists():
   if name in cache and cache[name].exists():shutil.copy2(cache[name],p)
   else:
    part=p.with_suffix('.nc.part');log=W/'logs'/(name+'.log')
    with log.open('w') as f:cp=subprocess.run(['curl','--http1.1','--fail','--location','--silent','--show-error','--connect-timeout','20','--max-time','300','--retry','2','--retry-delay','2','--retry-all-errors','-o',str(part),result['url']],stdout=f,stderr=f)
    if cp.returncode:raise RuntimeError(log.read_text()[-500:])
    if not part.read_bytes()[:8].startswith((b'\x89HDF',b'CDF')):raise RuntimeError('Not NetCDF')
    part.replace(p)
  result.update(status='downloaded',bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 except Exception as e:result.update(status='failed',error=str(e))
 with lock:
  count+=1
  with (W/'download.jsonl').open('a') as f:f.write(json.dumps(result)+'\n')
  print(count,len(items),result['status'],name[:100],flush=True)
 return result
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:results=list(pool.map(fetch,items))
(O/'download_manifest.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
print('DONE',len(results),sum(r['status']=='downloaded' for r in results),flush=True)
