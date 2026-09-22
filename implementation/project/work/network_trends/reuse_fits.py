from pathlib import Path
import sys,json
R=Path('@PROJECT_ROOT@');sys.path.insert(0,str(R/'work/ldsa_support/deps'))
import pandas as pd,numpy as np
W=R/'work/network_trends';A=R/'outputs/network_trends/analysis';P=W/'before_qc_correction'
old=pd.read_csv(P/'trend_input.csv');new=pd.read_csv(W/'trend_input.csv');changed=[]
for key,g in new.groupby('group_id'):
 a=old[old.group_id.eq(key)].sort_values('date');g=g.sort_values('date')
 if len(a)!=len(g) or not np.array_equal(a.date.values,g.date.values) or not np.allclose(a.value,g.value,rtol=1e-13,atol=1e-10):changed.append(key)
for name in ['trends.csv','invalid_percent_CI_raw.csv']:
 t=pd.read_csv(P/name);t=t[~(t.Site+'|'+t.Metric+'|'+t.Phase).isin(changed)];t.to_csv(A/name,index=False)
(W/'changed_fits_after_qc.json').write_text(json.dumps(changed,indent=2));print(changed)
