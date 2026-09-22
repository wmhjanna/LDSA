from pathlib import Path
import json,zipfile
O=Path('@DELIVERY_ROOT@/EBAS_observations');Z=O.parent/'EBAS_observations_raw.zip';manifest=json.loads((O/'download_manifest.json').read_text())
with zipfile.ZipFile(Z,'a',compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
 existing=set(z.namelist())
 for r in manifest:
  p=O/r['relative_path'];name=r['relative_path']
  if name not in existing:z.write(p,name);existing.add(name)
 for p in sorted(O.iterdir()):
  if p.is_file() and not p.name.endswith('.inspect.ndjson') and p.name not in existing:z.write(p,p.name);existing.add(p.name)
 print('ZIP_ENTRIES',len(z.namelist()),flush=True)
 assert len(existing)==len(z.namelist())
with zipfile.ZipFile(Z) as z:
 for r in manifest:assert z.getinfo(r['relative_path']).file_size==r['bytes']
 assert z.testzip() is None
final=Z.with_name('EBAS_observations.zip');Z.replace(final)
print('ZIP_VERIFIED',str(final),'bytes',final.stat().st_size,flush=True)
