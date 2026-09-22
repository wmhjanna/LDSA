from pathlib import Path
import json,urllib.request,concurrent.futures,hashlib,time
W=Path('work/ipr_download');O=Path('@DELIVERY_ROOT@/EBAS_IPR_PM25_NO2');(O/'source_netcdf').mkdir(parents=True,exist_ok=True)
items=[]
for d in json.loads((W/'ipr_catalog.json').read_text()):
 n=d['name'];yr=int(n.split('.')[1][:4])
 if yr>2019:continue
 no2='chemiluminescence_molybdenum' in n and ('.nitrogen_dioxide.' in n or '..air.' in n) and '.lev1.' not in n
 pm='.pm25.' in n and ('.pm25_mass.' in n or ('..pm25.' in n and ('filter_1pack' in n or 'mass_01' in n)))
 if no2 or pm:items.append(d)
(W/'selected.json').write_text(json.dumps(items,indent=2))
def get(d):
 n=d['name'];url='https://thredds.nilu.no/thredds/fileServer/'+d['urlPath'];p=O/'source_netcdf'/n
 try:
  if not p.exists():
   with urllib.request.urlopen(url,timeout=55) as r:b=r.read()
   p.write_bytes(b)
  b=p.read_bytes();assert b.startswith(b'\x89HDF') or b.startswith(b'CDF')
  return dict(file=n,url=url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
 except Exception as e:return dict(file=n,url=url,error=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
 result=[]
 for r in ex.map(get,items):result.append(r);print(r['file'],r.get('bytes',r.get('error')),flush=True)
(O/'download_manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print('DONE',len(result),flush=True)
