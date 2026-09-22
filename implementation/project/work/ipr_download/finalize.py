from pathlib import Path
import sys,json,shutil,datetime,hashlib
sys.path.insert(0,'work/ipr_download/deps')
import netCDF4,numpy as np
from zipfile import ZipFile,ZIP_DEFLATED
W=Path('work/ipr_download');O=Path('@DELIVERY_ROOT@/EBAS_IPR_PM25_NO2')
raw=O/'source_netcdf';other=W/'unselected_sources';other.mkdir(exist_ok=True)
manifest=json.loads((O/'download_manifest.json').read_text());byfile={x['file']:x for x in manifest};summary=[]
for p in sorted(raw.glob('*.nc')):
 with netCDF4.Dataset(p) as d:
  names=[n for n in d.variables if (n=='pm25_mass' or n.startswith('nitrogen_dioxide')) and not any(s in n for s in ['_qc','metadata','_max','_min','_stddev'])]
  start=np.array(netCDF4.num2date(d['time_bnds'][:,0],d['time'].units,only_use_cftime_datetimes=False))
  mask=np.array([2009<=t.year<=2019 for t in start])
  if not names or not mask.any():keep=False
  else:
   keep=True;details=[]
   for n in names:
    v=d[n];a=v[:];valid=~np.ma.getmaskarray(a)&np.isfinite(np.asarray(a))&mask
    years={str(y):int(sum(valid&np.array([t.year==y for t in start]))) for y in range(2009,2020)}
    details.append({'variable':n,'units':getattr(v,'units',''),'ebas_unit':getattr(v,'ebas_unit',''),'valid_observations_2009_2019':int(sum(valid)),'annual_valid_observations':years,'first_valid_start_UTC':str(start[valid][0]) if valid.any() else None,'last_valid_start_UTC':str(start[valid][-1]) if valid.any() else None})
   attrs={k:str(d.getncattr(k)) for k in d.ncattrs()}
   (O/'source_metadata').mkdir(exist_ok=True)
   (O/'source_metadata'/(p.stem+'.json')).write_text(json.dumps({'global':attrs,'variables':{n:{k:str(d[n].getncattr(k)) for k in d[n].ncattrs()} for n in names}},ensure_ascii=False,indent=2))
   row=byfile[p.name]|{'station':'Ispra','station_code':'IT0004R','GAW_ID':'IPR','time_resolution':getattr(d,'ebas_resolution_code','') or p.name.split('.')[7],'full_file_start':getattr(d,'time_coverage_start',''),'full_file_end':getattr(d,'time_coverage_end',''),'variables':details}
   summary.append(row)
 if not keep:shutil.move(p,other/p.name)
(O/'download_manifest.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
lines=['# Ispra PM₂.₅ and NO₂ observations', '',
'Station: IT0004R / Ispra / IPR_RB. The selected input period is 2009–2019.',
'Source NetCDF files retain provider metadata, quality flags and sampling intervals. Availability by year is recorded in the download manifest.',
'Data service: https://thredds.nilu.no/thredds/catalog/ebas/catalog.html',
'Access information: https://ebas.nilu.no/data-access/otheraccess/thredds/', '',
'Time bounds are expressed in UTC. Use the native sampling intervals when aggregating observations.',
'NO₂ expressed as mass of nitrogen requires multiplication by 46.0055/14.0067 to obtain mass of NO₂. Consult both the units and ebas_unit metadata; gas mixing ratios require a separate conversion using reference temperature and pressure.',
'Missing years are not interpolated. NO observations are not substituted for NO₂.', '',
'## Source files', '', '|Metric|File period|Resolution|File|', '|---|---|---|---|']
for r in summary:
 label='PM₂.₅' if any(v['variable']=='pm25_mass' for v in r['variables']) else 'NO₂'
 lines.append(f"|{label}|{r['full_file_start']} to {r['full_file_end']}|{r['time_resolution']}|{r['file']}|")
lines+=['', 'The download manifest records source URLs, SHA-256 checksums and coverage. Dataset citations, identifiers and acknowledgements are retained in the provider metadata. Follow the provider terms when reusing observations.']
(O/'data_notes.md').write_text('\n'.join(lines),encoding='utf-8')
archive=O.with_suffix('.zip')
with ZipFile(archive,'w',ZIP_DEFLATED) as z:
 for p in sorted(O.rglob('*')):
  if p.is_file():z.write(p,p.relative_to(O))
with ZipFile(archive) as z:assert z.testzip() is None
for r in summary:print(r['time_resolution'],r['file'],[(v['variable'],v['valid_observations_2009_2019']) for v in r['variables']])
print('FILES',len(summary),'ZIP',archive,archive.stat().st_size)
