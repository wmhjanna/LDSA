from pathlib import Path
import json,hashlib,shutil,zipfile,sys
import numpy as np,pandas as pd
from prepare import R,W,O,no_overlap
sys.path.insert(0,str(R/'work'));from current_pollutants import load_site,apply_current_pollutants
manifest=json.load(open(O/'download_manifest.json'));assert len(manifest)==74 and all(r['status']=='ok' for r in manifest)
keep={Path(r['path']).name for r in manifest};unused=W/'out_of_period_downloads';unused.mkdir(exist_ok=True)
for f in (O/'source_netcdf').rglob('*'):
 if f.is_file() and f.name not in keep:shutil.move(f,unused/f.name)
for r in manifest:assert hashlib.file_digest(Path(r['path']).open('rb'),'sha256').hexdigest()==r['sha256']
old=pd.read_pickle(O/'complementary_observations.pkl');summary=pd.read_csv(O/'pollutant_source_summary.csv');sites=json.load(open(O/'station_mapping.json'))
checks=[]
for s in sites:
 d=load_site(s['site']);assert d is not None
 for r in summary[summary.Site.eq(s['site'])].itertuples():
  q=d[d.metric.eq(r.Metric)];a=q[q.source.eq('EBAS')];b=q[q.source.ne('EBAS')];p=old[old.Site.eq(s['site'])&old.metric.eq(r.Metric)]
  assert len(q)==r.Updated_n and len(a)==r.EBAS_n and len(b)==r.Old_retained
  assert len(no_overlap(b,a))==len(b)
  expected=no_overlap(p,a).sort_values(['date','value']).reset_index(drop=True);actual=b.sort_values(['date','value']).reset_index(drop=True)
  assert np.array_equal(expected.date.to_numpy(dtype='datetime64[ns]'),actual.date.to_numpy(dtype='datetime64[ns]'))
  assert np.array_equal(expected.value.to_numpy(float),actual.value.to_numpy(float))
  assert len(p)==len(b)+r.Old_replaced
  if len(a):
   assert np.allclose(a.value,a.raw_value*a.conversion_factor,rtol=1e-13)
   aa=a.sort_values('date');assert not aa.date.duplicated().any()
   assert np.all(aa.date.to_numpy(dtype='datetime64[ns]')[1:]>=aa.end.to_numpy(dtype='datetime64[ns]')[:-1])
  checks.append(dict(Site=s['site'],Metric=r.Metric,reference_values_preserved=True,overlap_checked=True,EBAS_n=len(a),reference_n=len(b)))
 sentinel=pd.DataFrame({'date':[pd.Timestamp('2000-01-01')],'value':[123.]})
 assert apply_current_pollutants(s['site'],{'LDSA':sentinel})['LDSA'].equals(sentinel)
 # Human-readable provenance alongside compact runtime data.
 d.to_csv(O/'site_csv'/f'{s["site"]}_native_observations.csv.gz',index=False)
assert apply_current_pollutants('HEL_TR',{'SO2':sentinel})['SO2'].equals(sentinel) # loader does not change other-site policies
sources=json.load(open(O/'sources.json'));assert not any(x['metric']=='BC' for x in sources)
assert set(summary[summary.EBAS_n.gt(0)].Site)=={'HOH_RB','MEL_RB','NEU_SUB','SCH_RB','WAL_RB','ZUG_RB'}
result=dict(all_checks_passed=True,download_files=len(manifest),mapped_site_entries=len(sites),updated_sites=6,updated_metric_pairs=int(summary.EBAS_n.gt(0).sum()),new_observations=int(summary.EBAS_n.sum()),reference_rows_before=int(summary.Old_n.sum()),reference_rows_retained=int(summary.Old_retained.sum()),reference_rows_replaced=int(summary.Old_replaced.sum()),unavailable_EBAS_metric_pairs_with_reference=int(((summary.EBAS_n==0)&(summary.Old_n>0)).sum()),checks=checks)
(O/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
archive=O/'station_pollutant_inputs.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in sorted(O.rglob('*')):
  if not f.is_file() or f==archive or any(p in f.parts for p in ['source_netcdf','reference_metadata']) or f.suffix=='.pkl':continue
  z.write(f,f.relative_to(O))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print(json.dumps({k:v for k,v in result.items() if k!='checks'},ensure_ascii=False,indent=2));print('ZIP',archive,archive.stat().st_size)
