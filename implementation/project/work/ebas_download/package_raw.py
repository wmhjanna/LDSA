from pathlib import Path
import zipfile,time,json
O=Path('@DELIVERY_ROOT@/EBAS_observations')
Z=O.parent/'EBAS_observations_raw.zip'
done=set();total=0
with zipfile.ZipFile(Z,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
 while True:
  for p in sorted((O/'source_netcdf').rglob('*.nc')):
   if p.name in done:continue
   z.write(p,str(p.relative_to(O)));done.add(p.name);total+=p.stat().st_size
   if len(done)%25==0:print('PACK',len(done),round(total/1e9,3),flush=True)
  if (O/'download_manifest.json').exists():break
  time.sleep(10)
print('RAW_ZIP_DONE',len(done),Z.stat().st_size,flush=True)
