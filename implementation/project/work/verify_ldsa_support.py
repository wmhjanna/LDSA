from pathlib import Path
import sys,json,gzip,csv
import numpy as np
import pandas as pd
ROOT=Path('@PROJECT_ROOT@');STAGE=ROOT/'work/ldsa_support';sys.path.insert(0,str(STAGE/'deps'))
OUT=ROOT/'outputs/ldsa_cleaned_reanalysis';rows=[]
for d in sorted((OUT/'sites').iterdir()):
    if not (d/'complete.json').exists():continue
    site=d.name;s=pd.read_pickle(STAGE/f'{site}_series.pkl');info=json.loads((d/'complete.json').read_text())
    assert len(s)==info['rows']
    for reg in ['HA','TB','ALV','Total']:
        total=s[f'{reg}_LDSA_17_500nm'];parts=s[f'{reg}_Nucleation_17_25nm']+s[f'{reg}_Aitken_25_100nm']+s[f'{reg}_Accumulation_100_500nm']
        assert np.allclose(total,parts,equal_nan=True,rtol=1e-12,atol=1e-10),(site,reg)
    p=ROOT/f'outputs/site_data_pnsd_cm3_cleaned/{site}/{site}_standardized.xlsx'
    inp=pd.read_excel(p,engine='calamine',nrows=30);coeff=pd.read_csv(d/f'{site}_deposition_coefficients.csv');bins=inp[[c for c in inp if str(c).startswith('D_')]].to_numpy(float)
    mask=coeff.diameter_nm.between(17,500).to_numpy();manual=np.sum(bins[:,mask]*np.pi*(coeff.diameter_nm.to_numpy()[mask]/1000)**2*coeff.Total.to_numpy()[mask],axis=1)
    # First-period structural channels can be blank for LAN/MUL; compare only complete sampled rows here.
    finite=np.isfinite(manual);assert np.allclose(manual[finite],s.LDSA.iloc[:len(inp)].to_numpy()[finite],rtol=1e-12,atol=1e-10)
    xlsx=d/f'{site}_LDSA_standardized.xlsx';verified=False
    if xlsx.exists():
        book=pd.read_excel(xlsx,engine='calamine');ref=pd.read_csv(d/f'{site}_LDSA_standardized.csv.gz')
        assert book.shape==ref.shape and list(book.columns)==list(ref.columns),(site,book.shape,ref.shape)
        assert np.allclose(book.iloc[:,4:].apply(pd.to_numeric,errors='coerce').to_numpy(float),ref.iloc[:,4:].apply(pd.to_numeric,errors='coerce').to_numpy(float),equal_nan=True,rtol=1e-12,atol=1e-9),site
        assert book.Site.eq(site).all();assert book.Type.eq('Regional background' if site.endswith('_RB') else {'UB':'Urban background','TR':'Traffic','SUB':'Suburban background'}[site.rsplit('_',1)[-1]]).all()
        verified=True
    rows.append(dict(Site=site,rows=len(s),valid_primary=int(s.LDSA.notna().sum()),mode_sum_verified=True,sampled_formula_verified=True,xlsx_all_values_verified=verified))
    print(site,verified,flush=True)
pd.DataFrame(rows).to_csv(OUT/'analysis/verification.csv',index=False)
