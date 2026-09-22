from pathlib import Path
import sys
import numpy as np
import pandas as pd
ROOT=Path('@PROJECT_ROOT@');STAGE=ROOT/'work/ldsa_support';sys.path.insert(0,str(STAGE/'deps'))
for p in sorted(STAGE.glob('*_series.pkl')):
    s=pd.read_pickle(p)
    if 'PNC_Ntotal_17_500nm' in s:continue
    site=p.name.removesuffix('_series.pkl');d=ROOT/'outputs/ldsa_cleaned_reanalysis/sites'/site
    source=ROOT/f'outputs/site_data_pnsd_cm3_cleaned/{site}/{site}_standardized.xlsx'
    raw=pd.read_excel(source,engine='calamine');bins=raw[[c for c in raw if str(c).startswith('D_')]].to_numpy(float)
    dps=pd.read_csv(d/f'{site}_deposition_coefficients.csv').diameter_nm.to_numpy()
    for label,lo,hi,inclusive in [('Nucleation',17,25,False),('Aitken',25,100,False),('Accumulation',100,500,True),('Ntotal',17,500,True)]:
        mask=(dps>=lo)&((dps<=hi) if inclusive else (dps<hi));v=np.nansum(bins[:,mask],axis=1)
        c=f'Total_{"LDSA" if label=="Ntotal" else label}_{lo}_{hi}nm';v[s[c].isna()]=np.nan;s[f'PNC_{label}_{lo}_{hi}nm']=v
    assert np.allclose(s.PNC_Ntotal_17_500nm,s.PNC_Nucleation_17_25nm+s.PNC_Aitken_25_100nm+s.PNC_Accumulation_100_500nm,equal_nan=True)
    s.to_pickle(p);print(site,flush=True)
