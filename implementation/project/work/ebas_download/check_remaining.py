from pathlib import Path
import json,concurrent.futures
from download_all import ddscheck,fetch
W=Path('work/ebas_download');items=json.loads((W/'additional_metadata_candidates.json').read_text());items=[x for x in items if x['name'].split('.')[3] not in ['wet_only_sampler','bulk_sampler'] and not x['name'].split('.')[5].startswith('precip')]
print('EXTRA_METADATA',len(items),flush=True)
checks=[];futures=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as md,concurrent.futures.ThreadPoolExecutor(max_workers=3) as dl:
 for x,r in zip(items,md.map(ddscheck,items)):
  checks.append(r)
  if r.get('variables'):futures.append(dl.submit(fetch,x,'pollutants',r['variables']))
  if len(checks)%25==0:print('EXTRA_SCAN',len(checks),flush=True)
 results=[f.result() for f in futures]
(W/'additional_metadata_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
(W/'additional_download_manifest.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
print('EXTRA_DONE',len(checks),'MATCHED',len(results),flush=True)
