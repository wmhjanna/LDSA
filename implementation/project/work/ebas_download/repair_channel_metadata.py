from pathlib import Path
import json
from audit_downloads import audit
O=Path('@DELIVERY_ROOT@/EBAS_observations')
a={r['file']:r for r in map(json.loads,(O/'verified_metadata.jsonl').read_text().splitlines())}
for r in list(a.values()):
 if r['status']=='verified' and 'PNSD' in r.get('available_metrics',[]) and any(v.get('channel_dimensions') and not v.get('native_D_values') for v in r.get('variables',[]) if v.get('metric')=='PNSD'):
  updated=audit(O/r['relative_path'])
  with (O/'verified_metadata.jsonl').open('a') as f:f.write(json.dumps(updated,ensure_ascii=False)+'\n')
  a[r['file']]=updated
  print('CHANNEL_COORDINATES_CHECKED',r['file'],flush=True)
(O/'verified_metadata.json').write_text(json.dumps(list(a.values()),ensure_ascii=False,indent=2))
