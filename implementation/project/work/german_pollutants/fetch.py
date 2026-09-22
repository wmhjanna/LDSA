from pathlib import Path
import json,sys,re,subprocess,hashlib,shutil,concurrent.futures,xml.etree.ElementTree as E
import pandas as pd
R=Path('@PROJECT_ROOT@');W=R/'work/german_pollutants'
O=Path('@DELIVERY_ROOT@/german_pollutant_inputs');O.mkdir(exist_ok=True)
C=R/'work/ebas_download';CACHE=O.parent/'EBAS_observations'
target={'pm25_mass':'PM2.5','pm10_mass':'PM10','pm1_mass':'PM1','nitrogen_dioxide':'NO2','nitrogen_monoxide':'NO','sulphur_dioxide':'SO2','ozone':'O3','carbon_monoxide':'CO','equivalent_black_carbon':'BC','black_carbon':'BC'}
sites=[]
for r in json.load(open(C/'local_comparison.json')):
 if not r['code'].startswith('DE') or not r['current_analysis_sites']:continue
 for s in r['current_analysis_sites']:
  d=pd.read_pickle(R/'work/workbook_processing/series'/f'{s}_series.pkl')
  sites.append(dict(code=r['code'],site=s,name=r['name'],start=str(d.Time.min()),end=str(d.Time.max())))
(O/'station_mapping.json').write_text(json.dumps(sites,ensure_ascii=False,indent=2))
codes={r['code'] for r in sites};ns={'t':'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0'}
root=E.parse(W/'catalog.xml');items=[]
for e in root.findall('.//t:dataset',ns):
 n=e.get('name','');p=n.split('.')
 if not e.get('urlPath') or p[0] not in codes or not n.endswith('.nc'):continue
 items.append(dict(name=n,urlPath=e.get('urlPath')))
(W/'catalog_german.json').write_text(json.dumps(items,indent=2))
def check(x):
 n=x['name'];p=n.split('.');code=p[0]
 if int(p[1][:4])>max(int(s['end'][:4]) for s in sites if s['code']==code):return None
 duration=re.fullmatch(r'(\d+)(y|mo|w|d|h)',p[6])
 if duration:
  start=pd.to_datetime(p[1],format='%Y%m%d%H%M%S')
  days=int(duration[1])*{'y':366,'mo':31,'w':7,'d':1,'h':1/24}[duration[2]]
  if start+pd.Timedelta(days=days)<min(pd.Timestamp(s['start']).floor('D') for s in sites if s['code']==code):return None
 if '.lev1.' in n or p[5] in ['precip','wet','cloud']:return None
 if p[3] in ['wet_only_sampler','bulk_sampler','smps','dmps','tdmps','tsmps','v-smps','v-tdmps','nephelometer','aps','HTDMA','cpc','CCNC']:return None
 if p[4] in target:return dict(x,metrics=[target[p[4]]])
 # All other files: inspect variables, including multicomponent files.
 old=C/'dds'/(n+'.dds');local=W/'dds'/(n+'.dds');local.parent.mkdir(exist_ok=True)
 if old.exists():text=old.read_text()
 else:
  cp=subprocess.run(['curl','-fLsS','--connect-timeout','15','--max-time','75','--retry','1','https://thredds.nilu.no/thredds/dodsC/'+x['urlPath']+'.dds','-o',str(local)],capture_output=True)
  if cp.returncode:return dict(x,error=cp.stderr.decode()[-200:])
  text=local.read_text()
 found=[m for k,m in target.items() if re.search(r'\b'+k+r'(?:_|\[)',text)]
 return dict(x,metrics=sorted(set(found))) if found else None
def download(x):
 n=x['name'];dest=O/'source_netcdf'/n.split('.')[0]/n;dest.parent.mkdir(parents=True,exist_ok=True)
 cached=CACHE/'source_netcdf'/n.split('.')[0]/'pollutants'/n
 url='https://thredds.nilu.no/thredds/fileServer/'+x['urlPath']
 try:
  if not dest.exists():
   if cached.exists():shutil.copy2(cached,dest)
   else:
    part=dest.with_suffix('.part');cp=subprocess.run(['curl','-fLsS','--connect-timeout','20','--max-time','600','--retry','2',url,'-o',str(part)],capture_output=True)
    if cp.returncode:raise RuntimeError(cp.stderr.decode()[-300:])
    part.replace(dest)
  with dest.open('rb') as f:
   assert f.read(8).startswith((b'CDF',b'\x89HDF'))
  return dict(x,url=url,path=str(dest),sha256=hashlib.file_digest(dest.open('rb'),'sha256').hexdigest(),status='ok',reused_download=cached.exists())
 except Exception as e:return dict(x,url=url,status='failed',error=str(e))
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
  checked=[]
  for i,x in enumerate(pool.map(check,items)):
   if x:checked.append(x)
   if i%50==0:print('metadata',i,len(items),flush=True)
 (W/'metadata_checks.json').write_text(json.dumps(checked,indent=2))
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:results=list(pool.map(download,[x for x in checked if x.get('metrics')]))
 (O/'download_manifest.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
 print('COMPLETE',len(results),'failures',sum(x['status']!='ok' for x in results),'metadata errors',sum('error' in x for x in checked),flush=True)
