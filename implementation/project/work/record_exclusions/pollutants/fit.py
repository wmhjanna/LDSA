"""Fit newly integrated groups only; old 65-site results remain unchanged."""
from pathlib import Path
import subprocess,concurrent.futures,json
import pandas as pd
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;A=R/'outputs/record_exclusions/pollutants/analysis';OLD=R/'outputs/record_exclusions/analysis'
d=pd.read_csv(W/'new_trend_input.csv')
# Fit the phase groups selected in the input table.
d=d[d.Phase.str.startswith('Overlap_')]
keys=sorted(d.group_id.unique());workers=3
src=(R/'work/record_exclusions/workbook/trends.R').read_text().split("\ne<-read.csv(file.path(args[2],'site_eligibility.csv'))",1)[0]
(W/'new_trends.R').write_text(src)
def fit(i):
 p=W/f'fit_{i}';p.mkdir(exist_ok=True);inp=p/'input.csv';d[d.group_id.isin(keys[i::workers])].to_csv(inp,index=False)
 with (p/'run.log').open('w') as f:subprocess.run(['@RSCRIPT@',str(W/'new_trends.R'),str(inp),str(p)],cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
 return pd.read_csv(p/'trends.csv')
with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:results=list(pool.map(fit,range(workers)))
new=pd.concat(results,ignore_index=True);old=pd.read_csv(OLD/'trends.csv')
assert set(new.Site+'|'+new.Metric+'|'+new.Phase).isdisjoint(set(old.Site+'|'+old.Metric+'|'+old.Phase))
pd.concat([old,new],ignore_index=True).to_csv(A/'trends.csv',index=False);new.to_csv(W/'new_overlap_trends.csv',index=False)
bad=[pd.read_csv(OLD/'invalid_percent_CI_raw.csv')]+[pd.read_csv(W/f'fit_{i}/invalid_percent_CI_raw.csv') for i in range(workers)]
pd.concat(bad,ignore_index=True).to_csv(A/'invalid_percent_CI_raw.csv',index=False)
(W/'fit_complete.json').write_text(json.dumps(dict(requested=len(keys),fitted=len(new),scope='Exact-overlap pollutant trends; full-record and phase ancillary analyses are outside this stage'),indent=2))
print('FITTED',len(new),'of',len(keys),flush=True)
