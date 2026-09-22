from pathlib import Path
import sys,pickle,json,hashlib
ROOT=Path('@PROJECT_ROOT@');sys.path[:0]=[str(ROOT/'work/ldsa_support/deps'),str(ROOT/'work')]
import pandas as pd,numpy as np
from joined_sources import joined_parts
from build_site_tables import delta_log10
W=ROOT/'work/workbook_processing';O=ROOT/'outputs/workbook_processing';A=O/'analysis';flags=json.loads((O/'analysis_policy.json').read_text())['sentinel_values'];rows=[];cfrows=[]
for site in ['PAR_SUB','PRA_SUB']:
 obj=pickle.load(open(W/'pnsd'/(site+'.pkl'),'rb'));s=pd.read_pickle(W/'series'/(site+'_series.pkl'));offset=0;old=ROOT/'work/integration_support/series'
 for k,(path,sh,_) in enumerate(joined_parts(site)):
  r=path if isinstance(path,pd.DataFrame) else pd.read_excel(path,sheet_name=sh,engine='calamine');t=pd.to_datetime(r.iloc[:,0],errors='coerce');valid=t.notna()&t.dt.year.between(1990,2030);r=r[valid];t=t[valid];dp=np.array(r.columns[1:],float);x=r.iloc[:,1:].apply(pd.to_numeric,errors='coerce').to_numpy(float);x[~np.isfinite(x)|(x<0)|np.isin(x,flags)]=np.nan;bins=x*delta_log10(dp);cf=obj['coefficients'].set_index('diameter_nm').loc[dp];mask=(dp>=17)&(dp<=500);v=np.nansum(bins[:,mask]*np.pi*(dp[mask]/1000)**2*cf.Total.to_numpy()[mask],axis=1);v[np.isnan(bins[:,mask]).all(axis=1)]=np.nan
  selected=s.iloc[offset:offset+len(r)];assert np.array_equal(selected.Time.to_numpy(),t.to_numpy());assert np.allclose(selected.LDSA.to_numpy(),v,rtol=1e-12,atol=1e-9,equal_nan=True)
  if site=='PRA_SUB':
   oldsite=['PRA_2012-2016_SUB','PRA_2017-2019_SUB'][k];oldvals=pd.read_pickle(old/(oldsite+'_series.pkl'));assert np.array_equal(oldvals.Time.to_numpy(),t.to_numpy());assert np.allclose(oldvals.LDSA.to_numpy(),v,rtol=1e-12,atol=1e-9,equal_nan=True)
  z=pd.DataFrame({'date':t.to_numpy(),'value':v});z=z.dropna();mm=z.groupby(z.date.dt.to_period('M')).value.mean()
  rows.append(dict(Site=site,Config=k+1,source_rows=len(r),native_channels=len(dp),valid_LDSA=int(np.isfinite(v).sum()),actual_low_nm=dp[mask].min(),actual_high_nm=dp[mask].max(),all_rows_formula_verified=True,monthly_mean=mm.mean(),first_valid=z.date.min(),last_valid=z.date.max()))
  offset+=len(r)
 assert offset==len(s)
pd.DataFrame(rows).to_csv(A/'joined_config_independent_verification.csv',index=False)
DF=Path('@INPUT_ROOT@/data_sources/Mortality-LDSA/PNSD_inputs/wetransfer_pnsd-data-xiansheng_2026-03-11_1136/PNSD data Xiansheng/LDSA/DF1.xlsx')
d=pd.read_excel(DF,engine='calamine').iloc[:,110:115].apply(pd.to_numeric,errors='coerce').dropna();d.to_csv(O/'PNSD/PAR_SUB/PAR_second_configuration_native_DF.csv',index=False)
(O/'PNSD/PAR_SUB/DF_source.json').write_text(json.dumps({'source':str(DF),'sheet':'Sheet1','range':'DG1:DK122','sha256':hashlib.file_digest(DF.open('rb'),'sha256').hexdigest(),'channel_count':121,'exact_diameter_match':True,'interpolation':False},ensure_ascii=False,indent=2))
print(pd.DataFrame(rows).to_string(index=False))
