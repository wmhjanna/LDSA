from pathlib import Path
import shutil,json,sys
import pandas as pd,numpy as np
R=Path('@PROJECT_ROOT@');W=R/'work/record_exclusions';O=R/'outputs/record_exclusions';D=Path('@DELIVERY_ROOT@/record_exclusion_analysis')
oldW=R/'work/quality_control';oldO=R/'outputs/quality_control'
for p in [O,D,D/'PNG',W/'network',W/'layout',W/'s3_s6']:p.mkdir(parents=True,exist_ok=True)
# Independent outputs; immutable large workbook data are referenced, not rewritten.
for name in ['workbook']:
 (W/name).symlink_to(oldW/name,target_is_directory=True) if not (W/name).exists() else None
 (O/name).symlink_to(oldO/name,target_is_directory=True) if not (O/name).exists() else None
(W/'workbook_done').write_text('Workbook input reused unchanged')
for name in ['analysis','netcdf/analysis','netcdf/calculation_parameters']:shutil.copytree(oldO/name,O/name,dirs_exist_ok=True)
hour=O/'netcdf/hourly_ldsa';hour.mkdir(parents=True,exist_ok=True)
cases={'MT0001R':('GIO_RB',['2013-02-18','2013-03-25','2013-04-01']),'CZ0003R':('KOS_RB',['2019-11-18']),'FI0096G':('PAL_RB',['2014-01-13'])};audit=[];removed=[]
for p in (oldO/'netcdf/hourly_ldsa').glob('*.csv.gz'):
 code=p.name.split('_LDSA')[0];dest=hour/p.name
 if code not in cases:
  if not dest.exists():dest.symlink_to(p)
  continue
 d=pd.read_csv(p,parse_dates=['Time_UTC']);site,weeks=cases[code];mask=pd.Series(False,index=d.index)
 for week in weeks:
  start=pd.Timestamp(week);end=start+pd.Timedelta(days=7);sel=d.Time_UTC.ge(start)&d.Time_UTC.lt(end);finite=sel&d['Total_LDSA_um2_cm-3'].notna();mask|=sel
  audit.append(dict(Site=site,EBAS=code,start_UTC=str(start),end_UTC_exclusive=str(end),rows=int(sel.sum()),additional_finite_records_excluded=int(finite.sum()),weekly_mean_before=d.loc[finite,'Total_LDSA_um2_cm-3'].mean(),reason='Manual study-specific record exclusion; applied separately from provider quality control'))
 z=d[mask].copy();z.insert(0,'Site',site);removed.append(z)
 cols=[c for c in d if 'LDSA_' in c or c=='PNC_cm-3'];d.loc[mask,cols]=np.nan;d.loc[mask,'Valid_bins']=0;d['User_case_excluded']=mask
 d.to_csv(dest,index=False)
pd.DataFrame(audit).to_csv(D/'record_exclusion_summary.csv',index=False,encoding='utf-8-sig');pd.concat(removed).to_csv(O/'excluded_records_before_mask.csv.gz',index=False)
(W/'exclusions.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2));print(pd.DataFrame(audit).to_string(index=False))
for name in ['network/prepare.py','network/monthly.R','merge_network.py','network_sensitivity.R','render_fig2.py','render_s3_s6.py']:
 s=(oldW/name).read_text().replace('quality_control','record_exclusions').replace('provider_quality_control_analysis','record_exclusion_analysis');(W/name).write_text(s)
# Cache only sites outside the explicit exclusion set.
for f in ['trends.csv','invalid_percent_CI_raw.csv']:
 d=pd.read_csv(O/'netcdf/analysis'/f);d[~d.Site.isin([v[0] for v in cases.values()])].to_csv(O/'netcdf/analysis'/f,index=False)
for name in ['render.py','align_layout.py']:
 s=(R/'work/station_diagnostic_layout'/name).read_text().replace('quality_control','record_exclusions').replace('station_layout_figures','record_exclusion_analysis/PNG');(W/'layout'/name).write_text(s)
s=(R/'work/station_diagnostic_export/run.py').read_text().split("pointer_path=",1)[0].replace('quality_control','record_exclusions').replace('station_diagnostic_figures','record_exclusion_analysis');(W/'render_s6.py').write_text(s+"\n(TASK/'s6_verification.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))\n")
