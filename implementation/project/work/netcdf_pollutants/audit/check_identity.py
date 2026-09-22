from pathlib import Path
import sys,json,re,math
import pandas as pd,numpy as np
R=Path('@PROJECT_ROOT@');A=R/'work/netcdf_pollutants/audit';B=Path('@DELIVERY_ROOT@/EBAS_observations');sys.path.insert(0,str(R/'work/ipr_download/deps'));import netCDF4
M={r['old']:r for r in json.load(open(R/'work/station_registry/mapping.json'))};reg=pd.read_csv(R/'outputs/quality_control/netcdf/calculation_parameters/source_qc_registry.csv');reg=reg[reg.selected_rows.gt(0)];refs=[]
norm=lambda k:re.sub('[^a-z0-9]','',k.lower()).removeprefix('ebas')
def loc(d,p):
 try:return float(d[p+'latitude']),float(d[p+'longitude'])
 except:return None,None
for r in reg.itertuples():
 with netCDF4.Dataset(B/'source_netcdf'/r.code/'PNSD'/r.file) as n:
  v=n[r.variable];base={norm(k):str(n.getncattr(k)) for k in n.ncattrs()};base.update({norm(k):str(v.getncattr(k)) for k in v.ncattrs()});md=n[r.variable+'_ebasmetadata'];idx=tuple(slice(None) if 'time' in d else 0 for d in md.dimensions);allmd=[{**base,**{norm(k):str(z) for k,z in json.loads(str(x)).items()}} for x in np.asarray(md[idx]).ravel()]
  for d in allmd:
   sla,slo=loc(d,'station');mla,mlo=loc(d,'measurement');refs.append(dict(code=r.code,Site=M[r.code]['new'],Source_ID=r.Source_ID,file=r.file,station_lat=sla,station_lon=slo,measurement_lat=mla,measurement_lon=mlo,station_name=d.get('stationname',''),platform=d.get('platformcode',''),height=d.get('measurementheight',''),inlet=d.get('inletdescription','')))
pd.DataFrame(refs).to_csv(A/'PNSD_full_reference_locations.csv',index=False)
l=pd.read_csv(A/'pollutant_PNSD_locations.csv');rr=pd.DataFrame(refs)
def distance(la,lo,x,y):
 if not np.isfinite([la,lo,x,y]).all():return np.nan
 p1,p2=np.deg2rad([la,x]);dl=np.deg2rad(lo-y);z=np.sin((p1-p2)/2)**2+np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2;return 6371000*2*np.arcsin(min(1,np.sqrt(z)))
out=[]
for r in l.itertuples():
 ds=rr[rr.code.eq(r.code)];dist=[]
 for ref in ds.itertuples():
  for pr in ['station','measurement']:
   x=getattr(ref,pr+'_lat');y=getattr(ref,pr+'_lon');v=distance(r.latitude,r.longitude,x,y)
   if np.isfinite(v):dist.append((v,ref.Source_ID,pr,x,y,ref.station_name,ref.platform))
 best=min(dist) if dist else None;d=r._asdict();d.pop('Index',None);d['min_distance_full_PNSD_m']=best[0] if best else None;d['closest_PNSD_detail']=json.dumps(best,ensure_ascii=False);out.append(d)
o=pd.DataFrame(out);o.to_csv(A/'pollutant_PNSD_full_location_comparison.csv',index=False)
# A clear internal station/measurement mismatch that cannot safely be resolved from code alone.
filename='IE0031R.20140101000000.20250902084524.TEOM...6y.1h.IE03L_TEOM_1405_MH..lev2.nc';exclude=[dict(file=filename,reason='Mace Head name/code but measurement coordinate 58.38,8.25 lies at Birkenes; inconsistent with local PNSD Mace Head coordinates by about 1,259 km. Excluded because station identity is unresolved; the metadata do not establish a station relocation.')];(A/'identity_exclusions.json').write_text(json.dumps(exclude,indent=2))
print(o.groupby('Site').min_distance_full_PNSD_m.agg(['min','max']).to_string());print('above500',o[o.min_distance_full_PNSD_m.gt(500)][['Site','latitude','longitude','min_distance_full_PNSD_m']].drop_duplicates().to_string(index=False))
