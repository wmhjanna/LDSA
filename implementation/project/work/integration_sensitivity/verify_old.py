from pathlib import Path
import pandas as pd,numpy as np,pickle,json,gc
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;a=pd.read_csv(W/'workbook_sources.csv');report=[]
for r in a.itertuples():
 s=pd.read_pickle(r.source);d=pickle.load(open(R/'work/quality_control/workbook/pnsd'/f'{r.Site}.pkl','rb'));dp=d['diameters'];b=d['bins'];cf=d['coefficients'];n=0
 for i in np.linspace(0,len(s)-1,17,dtype=int):
  for name,low,high in [('primary',31.65 if r.Site=='LND_TR' else 17,500),('10_500nm',10,500),('10_800nm',10,800)]:
   jj=np.flatnonzero((dp>=low)&(dp<=high)&np.isfinite(b[i]))
   for region in ['HA','TB','ALV','Total']:
    val=sum(float(b[i,j])*np.pi*(float(dp[j])/1000)**2*float(cf[region].iloc[j]) for j in jj) if len(jj) else np.nan
    assert np.isclose(val,s[region+'_LDSA_'+name].iloc[i],equal_nan=True,rtol=1e-12,atol=1e-8),(r.Site,i,region,name);n+=1
 report.append(dict(Site=r.Site,independent_scalar_checks=n,passed=True));del s,d,b;gc.collect()
pd.DataFrame(report).to_csv(W/'workbook_formula_verification.csv',index=False);print('CHECKED',sum(x['independent_scalar_checks'] for x in report))
