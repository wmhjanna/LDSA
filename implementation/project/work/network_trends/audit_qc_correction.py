from pathlib import Path
import sys,json
R=Path('@PROJECT_ROOT@');sys.path.insert(0,str(R/'work/ldsa_support/deps'))
import numpy as np,pandas as pd
old=Path('@DELIVERY_ROOT@/netcdf_ldsa_inputs')
new=R/'outputs/network_trends/corrected_LDSA';W=R/'work/network_trends'
records=[]
for f in sorted((old/'hourly_ldsa').glob('*.csv.gz')):
 a=pd.read_csv(f);b=pd.read_csv(new/'hourly_ldsa'/f.name)
 assert a.Time_UTC.equals(b.Time_UTC) and a.Source_ID.equals(b.Source_ID)
 cols=['Head_LDSA_um2_cm-3','P_LDSA_um2_cm-3','TB_LDSA_um2_cm-3','Total_LDSA_um2_cm-3','PNC_cm-3']
 change=~np.isclose(a[cols].to_numpy(),b[cols].to_numpy(),equal_nan=True,rtol=1e-11,atol=1e-9).all(axis=1)
 records.append(dict(code=f.name[:7],changed_rows=int(change.sum()),new_missing_rows=int((a['Total_LDSA_um2_cm-3'].notna()&b['Total_LDSA_um2_cm-3'].isna()).sum()),first_changed=a.loc[change,'Time_UTC'].min() if change.any() else None,last_changed=a.loc[change,'Time_UTC'].max() if change.any() else None,source_ids=list(map(int,a.loc[change,'Source_ID'].unique()))))
 if change.any():
  detail=a.loc[change,['Time_UTC','Source_ID','Total_LDSA_um2_cm-3','Valid_bins']].copy();detail['corrected_Total_LDSA_um2_cm-3']=b.loc[change,'Total_LDSA_um2_cm-3'];detail['corrected_Valid_bins']=b.loc[change,'Valid_bins'];detail.to_csv(new/f'{f.name[:7]}_quality_control_changes.csv.gz',index=False)
pd.testing.assert_frame_equal(pd.read_csv(old/'site_spatial_statistics.csv'),pd.read_csv(new/'site_spatial_statistics.csv'),check_exact=False,rtol=1e-12,atol=1e-10)
report=dict(added_exact_missing_code=999999.99999,all31_period_statistics_unchanged=True,stations=records)
(W/'qc_correction_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));print(json.dumps([r for r in records if r['changed_rows']],ensure_ascii=False,indent=2))
