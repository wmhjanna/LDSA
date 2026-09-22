from pathlib import Path
import sys,json,re
sys.path.insert(0,'@PROJECT_ROOT@/work/ipr_download/deps')
import netCDF4
ROOT=Path('@PROJECT_ROOT@');W=ROOT/'work/meteorology'
B=Path('@DELIVERY_ROOT@/EBAS_observations')
cat={r['name']:r for r in json.loads((W/'station_catalog.json').read_text())}
manifest=json.loads((B/'download_manifest.json').read_text())
def get(md,key):
 canon=lambda s:re.sub('[^a-z0-9]','',s.lower()).removeprefix('ebas')
 return next((v for k,v in md.items() if canon(k)==canon(key)),'')
def match(k):return any(x in k.lower() for x in ['temperature','humidity','pressure','wind','radiation','precipitation'])
rows=[];seen=0
for r in manifest:
 if r['file'] not in cat or r.get('status')!='downloaded':continue
 p=B/r['relative_path'];seen+=1
 try:
  with netCDF4.Dataset(p) as ds:
   for key,v in ds.variables.items():
    if not key.endswith('_ebasmetadata') or not match(key):continue
    md=json.loads(str(v[tuple(0 for _ in v.dimensions)]));matrix=get(md,'Matrix');comp=get(md,'Component');inst=get(md,'Instrument type')
    accepted=matrix.lower()=='met' or (comp=='precipitation_amount' and matrix.lower()=='precip')
    rows.append(dict(file=r['file'],relative_path=r['relative_path'],variable=key[:-13],matrix=matrix,component=comp,instrument=inst,accepted=accepted,units=get(md,'Unit')))
 except Exception as e:rows.append(dict(file=r['file'],error=str(e),accepted=False))
 if seen%100==0:print('SCANNED',seen,flush=True)
(W/'cached_met_metadata_scan.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
extra=sorted({r['file'] for r in rows if r['accepted']})
(W/'extra_local_candidates.json').write_text(json.dumps([cat[k] for k in extra],indent=2))
print('DONE',seen,'variables',len(rows),'ambient_files',len(extra));print(extra)
