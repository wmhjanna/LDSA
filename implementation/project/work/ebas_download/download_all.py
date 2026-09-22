from pathlib import Path
import json,re,datetime,concurrent.futures,subprocess,hashlib,time,threading,shutil
W=Path('work/ebas_download')
O=Path('@DELIVERY_ROOT@/EBAS_observations')
O.mkdir(exist_ok=True);(W/'http_logs').mkdir(exist_ok=True);(W/'dds').mkdir(exist_ok=True)
items=json.loads((W/'catalog_65sites.json').read_text());lock=threading.Lock()
stats={'downloaded':0,'failed':0,'bytes':0,'started':time.time()}
target={'pm25_mass':'PM2.5','pm10_mass':'PM10','pm1_mass':'PM1','nitrogen_dioxide':'NO2','nitrogen_monoxide':'NO','sulphur_dioxide':'SO2','ozone':'O3','carbon_monoxide':'CO','equivalent_black_carbon':'BC','black_carbon':'BC'}
def eligible(x):
 p=x['name'].split('.')
 start=datetime.datetime.strptime(p[1],'%Y%m%d%H%M%S')
 if start.year>2023:return False
 m=re.fullmatch(r'(\d+)(y|mo|w|d|h)',p[6])
 if m:
  days=int(m[1])*{'y':366,'mo':31,'w':7,'d':1,'h':1/24}[m[2]]
  if start+datetime.timedelta(days=days)<datetime.datetime(1996,1,1):return False
 return True
def matches(n):return sorted({label for var,label in target.items() if re.search(r'\b'+var+r'(?:_|\[)',n)})
def fetch(x,category,variables):
 name=x['name'];code=name.split('.')[0];dest=O/'source_netcdf'/code/category/name
 dest.parent.mkdir(parents=True,exist_ok=True);url='https://thredds.nilu.no/thredds/fileServer/'+x['urlPath']
 result={'station_code':code,'file':name,'category':category,'requested_variables':variables,'url':url,'relative_path':str(dest.relative_to(O)),'catalog_size_bytes':x['size_bytes']}
 try:
  if not dest.exists():
   cached=Path('@DELIVERY_ROOT@/EBAS_IPR_PM25_NO2/source_netcdf')/name
   if cached.exists():shutil.copy2(cached,dest)
   else:
    part=dest.with_suffix('.nc.part');log=W/'http_logs'/(name+'.log')
    with log.open('w') as f:
     cp=subprocess.run(['curl','--http1.1','--fail','--location','--silent','--show-error','--connect-timeout','25','--max-time','2400','--speed-limit','1024','--speed-time','90','--retry','2','--retry-all-errors','--retry-delay','3','--continue-at','-','--output',str(part),url],stdout=f,stderr=f)
    if cp.returncode:raise RuntimeError(log.read_text()[-500:])
    with part.open('rb') as f:head=f.read(8)
    if not(head.startswith(b'\x89HDF') or head.startswith(b'CDF')):raise ValueError('Not a NetCDF file')
    part.replace(dest)
  with dest.open('rb') as f:result['sha256']=hashlib.file_digest(f,'sha256').hexdigest()
  result.update(status='downloaded',bytes=dest.stat().st_size)
 except Exception as e:result.update(status='failed',error=str(e))
 with lock:
  with (O/'download_log.jsonl').open('a') as f:f.write(json.dumps(result,ensure_ascii=False)+'\n')
  stats['downloaded' if result['status']=='downloaded' else 'failed']+=1;stats['bytes']+=result.get('bytes',0)
  (W/'progress.json').write_text(json.dumps(stats))
  print('FETCH',stats['downloaded'],stats['failed'],round(stats['bytes']/1e9,3),code,category,name[:60],flush=True)
 return result
def ddscheck(x):
 name=x['name'];p=W/'dds'/(name+'.dds');url='https://thredds.nilu.no/thredds/dodsC/'+x['urlPath']+'.dds'
 if not p.exists():
  cp=subprocess.run(['curl','--fail','--location','--silent','--show-error','--connect-timeout','15','--max-time','65','--retry','1','--output',str(p),url],capture_output=True)
  if cp.returncode:return {'file':name,'status':'metadata_failed','error':cp.stderr.decode()[-300:]}
 text=p.read_text();found=matches(text)
 return {'file':name,'status':'metadata_checked','variables':found}
if __name__=='__main__':
 pnsd=[];direct=[];candidates=[]
 for x in items:
  if not eligible(x):continue
  p=x['name'].split('.');instrument=p[3];comp=p[4]
  if instrument in ['smps','dmps','tdmps','tsmps'] and '.lev1.' not in x['name']:
   pnsd.append(x)
  elif comp in target:direct.append((x,[target[comp]]))
  elif not comp and instrument not in ['wet_only_sampler','bulk_sampler','nephelometer','sun_tracking_filter_radiometer','aps','v-smps','v-tdmps','smps','dmps','cpc','pollen_monitor','HTDMA','optical_particle_size_spectrometer','optical_particle_counter','CCNC','DMPS-CCNC','Hg_mon','PTR-MS'] and p[5] not in ['precip','wet','cloud']:
   candidates.append(x)
 (W/'selection_plan.json').write_text(json.dumps({'PNSD':pnsd,'direct_pollutants':[x for x,v in direct],'metadata_candidates':candidates},indent=2))
 print('PLAN',len(pnsd),len(direct),len(candidates),'PNSD_GB',sum(x['size_bytes'] for x in pnsd)/1e9,flush=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=10) as dl,concurrent.futures.ThreadPoolExecutor(max_workers=4) as md:
  futures=[dl.submit(fetch,x,'PNSD',['PNSD']) for x in sorted(pnsd,key=lambda x:-x['size_bytes'])]
  futures +=[dl.submit(fetch,x,'pollutants',v) for x,v in direct]
  checks=[]
  for x,r in zip(candidates,md.map(ddscheck,candidates)):
   checks.append(r)
   if r.get('variables'):futures.append(dl.submit(fetch,x,'pollutants',r['variables']))
   if len(checks)%25==0:print('METADATA',len(checks),'/',len(candidates),flush=True)
  (W/'metadata_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
  result=[f.result() for f in concurrent.futures.as_completed(futures)]
 (O/'download_manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
 print('DONE',stats,flush=True)
