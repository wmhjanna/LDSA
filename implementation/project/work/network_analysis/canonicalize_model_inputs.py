from pathlib import Path
import pandas as pd,numpy as np,json
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;M=W/'model'
old=pd.read_pickle(R/'work/quality_control/workbook/weekly_model_reference.pkl').set_index(['Site','date']);new=pd.read_pickle(M/'weekly_model_reference.pkl').set_index(['Site','date']);idx=old.index.intersection(new.index);audit=[]
for c in ['BC','CO','NO','NO2','O3','PM2.5','SO2']:
 a=old.loc[idx,c];b=new.loc[idx,c];same=np.isfinite(a)&np.isfinite(b)&np.isclose(a,b,rtol=0,atol=1e-10);changed=same&(a!=b);new.loc[idx[same],c]=a[same];audit.append(dict(feature=c,roundoff_only_values_canonicalized=int(changed.sum()),absolute_tolerance=1e-10))
new.reset_index().to_pickle(M/'weekly_model_reference.pkl');(M/'roundoff_canonicalization.json').write_text(json.dumps(audit,indent=2))
