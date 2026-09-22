from pathlib import Path
import json,pickle
import pandas as pd,numpy as np
ROOT=Path('@PROJECT_ROOT@');B=Path('@INPUT_ROOT@/data_sources')
SRC=B/'Mortality-LDSA/PNSD_inputs/PNSD data Xiansheng/Paris_PNSD_2012_2019.xlsx'
DF=B/'Mortality-LDSA/PNSD_inputs/wetransfer_pnsd-data-xiansheng_2026-03-11_1136/PNSD data Xiansheng/LDSA/DF1.xlsx'
cache={}
def joined_parts(site):
 if site in cache:return cache[site]['parts']
 dfs=[];parts=[];aux=[]
 if site=='PAR_SUB':
  raw=pd.read_excel(SRC,sheet_name='Data',header=None,engine='calamine')
  headers=[i for i in range(len(raw)) if i==0 or raw.iloc[i,20]=='CHANGE IN SIZE BINS']
  assert headers==[0,35065]
  tables=[pd.read_excel(B/'source_data/PAR_SUB/DF10.27-807.76.xlsx',engine='calamine'),pd.read_excel(DF,engine='calamine').iloc[:,110:115]]
  for k,h in enumerate(headers):
   rr=raw.iloc[h];dp=[]
   for v in rr.iloc[21:142]:
    if not isinstance(v,(int,float,np.number)) or not np.isfinite(v):break
    dp.append(float(v))
   end=headers[k+1] if k+1<len(headers) else len(raw)
   d=raw.iloc[h+1:end,[2]+list(range(21,21+len(dp)))].copy();d.columns=['Time']+dp;d.attrs={'source_path':str(SRC),'source_row_start':h+2}
   table=tables[k].apply(pd.to_numeric,errors='coerce').dropna();assert len(table)==len(dp);assert np.allclose(table.iloc[:,0]*1000,dp,rtol=0,atol=1e-10)
   cf=pd.DataFrame(table.iloc[:,1:5].to_numpy(),columns=['HA','TB','ALV','Total']);cf.insert(0,'diameter_nm',dp);dfs.append(cf);parts.append((d,'Data',None));aux.append(pd.DataFrame({'Time':pd.to_datetime(d.Time,errors='coerce')}))
 elif site=='PRA_SUB':
  for period in ['2012-2016','2017-2019']:
   old='PRA_'+period+'_SUB';p=pickle.load(open(ROOT/'work/pnsd_cleaning/pnsd'/(old+'.pkl'),'rb'));cf=p['coefficients'];dfs.append(cf);aux.append(p['aux'])
   path=B/'source_data'/('PRA_SUB_'+period)/('PRA_SUB_'+period+'.xlsx');parts.append((path,0,None))
 else:raise ValueError(site)
 coeff=pd.concat(dfs,ignore_index=True)
 for d,g in coeff.groupby('diameter_nm'):
  assert np.allclose(g.iloc[:,1:].to_numpy(),g.iloc[0,1:].to_numpy(),rtol=0,atol=1e-12),(site,'conflicting shared DF',d)
 coeff=coeff.drop_duplicates('diameter_nm').sort_values('diameter_nm').reset_index(drop=True)
 oldaux=pd.concat(aux,ignore_index=True);oldaux=oldaux[oldaux.Time.notna()].reset_index(drop=True)
 if site=='PRA_SUB':assert not oldaux.Time.duplicated().any(),site
 cache[site]=dict(parts=parts,diameters=coeff.diameter_nm.to_numpy(),coefficients=coeff,aux=oldaux)
 return parts
def joined_template(site):
 joined_parts(site);return cache[site]
