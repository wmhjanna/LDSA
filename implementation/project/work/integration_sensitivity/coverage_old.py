from pathlib import Path
import pandas as pd,numpy as np,pickle,gc
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;stats=pd.read_csv(W/'workbook_sources.csv');rows=[]
for r in stats.itertuples():
 s=pd.read_pickle(r.source);d=pickle.load(open(R/'work/quality_control/workbook/pnsd'/f'{r.Site}.pkl','rb'));t=pd.to_datetime(s.Time);ys=list(map(int,str(r.Period).split('-')));mask=t.dt.year.between(ys[0],ys[-1])&s.Total_LDSA_primary.notna();dp=d['diameters'];struct=d['structural'][mask];configs,counts=np.unique(struct,axis=0,return_counts=True)
 for i,(st,n) in enumerate(zip(configs,counts)):
  native=dp[~st];rows.append(dict(Site=r.Site,Source_ID=i+1,native_low_nm=native.min(),native_high_nm=native.max(),n=n,Period=r.Period,covers_10_500=native.min()<=10+1e-8 and native.max()>=500-1e-8,covers_10_800=native.min()<=10+1e-8 and native.max()>=800-1e-8,source=r.source,sha256=r.sha256))
 del s,d,struct;gc.collect()
pd.DataFrame(rows).to_csv(W/'workbook_coverage.csv',index=False);print(pd.DataFrame(rows).groupby('Site')[['covers_10_500','covers_10_800']].all().sum())
