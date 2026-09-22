from pathlib import Path
import sys,pickle,json
import pandas as pd,numpy as np
R=Path('@PROJECT_ROOT@');W=R/'work/quality_control';V=W/'verification';masks=pickle.load(open(W/'old_qc_masks.pkl','rb'));rows=[]
assert len(list((W/'workbook/series').glob('*_series.pkl')))==34
for fn in sorted((W/'workbook/pnsd').glob('*.pkl')):
 site=fn.stem;o=pickle.load(open(R/'work/workbook_processing/pnsd'/fn.name,'rb'));n=pickle.load(open(fn,'rb'));mask=masks.get(site,np.zeros(o['bins'].shape,bool));assert np.array_equal(o['diameters'],n['diameters']);assert np.array_equal(o['structural'],n['structural']);pd.testing.assert_frame_equal(o['coefficients'],n['coefficients']);assert o['aux'].Time.equals(n['aux'].Time)
 expected=o['bins'].copy();expected[mask]=np.nan;assert np.array_equal(n['bins'],expected,equal_nan=True),site
 os=pd.read_pickle(R/'work/workbook_processing/series'/f'{site}_series.pkl');ns=pd.read_pickle(W/'workbook/series'/f'{site}_series.pkl');cols=[c for c in os if 'LDSA' in c or c.startswith('PNC_')];changed=np.zeros(len(os),bool);outside=np.zeros(len(os),bool)
 for c in cols:
  delta=~np.isclose(os[c],ns[c],equal_nan=True,rtol=1e-12,atol=1e-9);changed|=delta;outside|=delta&~mask.any(1)
 assert not outside.any(),site
 prim_removed=os.LDSA.notna()&ns.LDSA.isna();rows.append(dict(Site=site,rows=len(os),bins_removed=int(mask.sum()),masked_rows=int(mask.any(1).sum()),primary_removed=int(prim_removed.sum()),any_LDSA_PNC_changed=int(changed.sum()),all_unmasked_unchanged=True))
q=pd.DataFrame(rows);q.to_csv(V/'workbook_recalculation_comparison.csv',index=False);assert len(q)==34;assert q.primary_removed.sum()==68;result=dict(workbook_verified=True,unaffected_sites=int(q.masked_rows.eq(0).sum()),masked_rows_total=int(q.masked_rows.sum()),primary_finite_to_missing_total=int(q.primary_removed.sum()),affected=q[q.masked_rows.gt(0)].to_dict('records'));(V/'workbook_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
