from pathlib import Path
import json,hashlib,copy
import pandas as pd,numpy as np,openpyxl
ROOT=Path('@PROJECT_ROOT@');W=ROOT/'work/quality_control/spatial'
BASE=Path('@DELIVERY_ROOT@');O=BASE/'provider_quality_control_analysis'
p=json.loads((W/'payload.json').read_text());checks=[]
for n in [65,49]:
 source=BASE/f'Figure1_bc_{n}_station_tables/location_LDSA-{n}sites_site_abbreviationexport.xlsx'
 path=O/'tables'/f'location_LDSA-{n}sites_quality_control.xlsx'
 before=openpyxl.load_workbook(source,data_only=False);after=openpyxl.load_workbook(path,data_only=False)
 assert before.sheetnames==after.sheetnames==['location']
 b=before['location'];a=after['location'];assert a.max_row==n+1 and a.max_column==27
 for r in range(1,n+2):
  for c in range(1,28):
   aa=a.cell(r,c);bb=b.cell(r,c)
   if r==1 or c<=11:assert aa.value==bb.value,(n,r,c,aa.value,bb.value)
   else:
    exp=p['statistics'][p['sites'].index(a.cell(r,2).value)][c-12]
    if isinstance(exp,str):assert aa.value==exp
    else:assert isinstance(aa.value,(int,float)) and np.isclose(aa.value,exp,rtol=1e-13,atol=1e-12),(n,r,c)
   assert aa.font.name==bb.font.name and aa.font.sz==bb.font.sz and aa.font.b==bb.font.b,(n,r,c,'font')
   assert copy.copy(aa.alignment)==copy.copy(bb.alignment),(n,r,c,'alignment')
   assert aa.number_format==bb.number_format,(n,r,c,'format',aa.number_format,bb.number_format)
   assert aa.fill.patternType==bb.fill.patternType,(n,r,c,'fill pattern')
   if aa.fill.patternType=='solid':assert aa.fill.fgColor==bb.fill.fgColor,(n,r,c,'visible solid fill')
 checks.append(dict(workbook=str(path),rows=n,columns=27,metadata_unchanged=True,numeric_statistics_verified=True,font_alignment_fill_number_formats_preserved=True,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
periods=pd.read_csv(O/'tables/spatial_analysis_periods.csv');attrib=pd.read_csv(O/'tables/spatial_statistics_comparison.csv');means=attrib[attrib.field.eq('Total_Mean')]
report='# Spatial-statistics verification\n\n'
report+=f"Stations using 2015–2019: {int((~periods.fallback).sum())}; stations using the nearest available year: {int(periods.fallback.sum())}.\n\n| Station | Analysis period | Valid records |\n|---|---:|---:|\n"
for r in periods[periods.fallback].itertuples():report+=f'| {r.Site} | {r.Period} | {r.valid_records} |\n'
report+='\n| Station | Reference mean | Harmonised mean | QC mean | Harmonisation difference | QC difference |\n|---|---:|---:|---:|---:|---:|\n'
for r in means.itertuples():
 if not(np.isclose(r.source_harmonization_delta,0,atol=1e-7) and np.isclose(r.official_QC_delta,0,atol=1e-7)):
  report+=f'| {r.Site} | {r.reference_table:.4f} | {r.latest_source_before_QC:.4f} | {r.after_QC:.4f} | {r.source_harmonization_delta:+.4f} | {r.official_QC_delta:+.4f} |\n'
report+='\nAll workbook checks completed.\n'
(O/'tables/spatial_validation.md').write_text(report)
(W/'workbook_verification.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
print('WORKBOOKS_VERIFIED');print(means[['Site','reference_table','latest_source_before_QC','after_QC','source_harmonization_delta','official_QC_delta']].to_string(index=False))
