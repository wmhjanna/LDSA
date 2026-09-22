from pathlib import Path
import json,re
import pandas as pd,numpy as np
W=Path('work/ebas_download');rows=[]
for p in sorted(Path('work/workbook_processing/series').glob('*_series.pkl')):
 site=p.name.removesuffix('_series.pkl');d=pd.read_pickle(p);t=pd.to_datetime(d.Time,errors='coerce')
 for col in d.columns:
  key=str(col).split('（')[0]
  if key not in ['PM2.5','PM10','PM1','NO2','NO','SO2','O3','CO','BC']:continue
  v=pd.to_numeric(d[col],errors='coerce');ok=np.isfinite(v)&t.notna()
  rows.append({'site':site,'metric':key,'column':col,'finite_records':int(ok.sum()),'first':str(t[ok].min()) if ok.any() else '', 'last':str(t[ok].max()) if ok.any() else ''})
(W/'local_current_pollutants.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
print('local pollutant columns checked',len(rows))
