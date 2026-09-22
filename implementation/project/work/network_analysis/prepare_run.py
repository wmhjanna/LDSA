from pathlib import Path
import sys,json,shutil,subprocess
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;O=R/'outputs/network_analysis';A=O/'base_analysis';P=W/'pollutants';PA=O/'analysis';B=Path('@DELIVERY_ROOT@');D=B/'ancillary_analysis_figures'
sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import numpy as np,pandas as pd
from analyze_cleaned_ldsa import monthly
for p in [O,A,P,PA,D,D/'PNG',D/'analysis_results',W/'model']:p.mkdir(parents=True,exist_ok=True)
OLD=R/'outputs/record_exclusions';OW=R/'work/record_exclusions';shutil.copytree(OLD/'analysis',A,dirs_exist_ok=True)
shutil.copy2(OW/'trend_input.csv',W/'trend_input.csv');shutil.copy2(OW/'weekly.pkl',W/'weekly.pkl')
# Keep all established masks and additionally use the accepted GIO four-hour exclusion.
hour=O/'netcdf/hourly_ldsa';hour.mkdir(parents=True,exist_ok=True)
for p in (OLD/'netcdf/hourly_ldsa').glob('*.csv.gz'):
 origin=R/'work/gio_record_filter/MT0001R_LDSA.csv.gz' if p.name=='MT0001R_LDSA.csv.gz' else p
 if not (hour/p.name).exists():(hour/p.name).symlink_to(origin)
gio=pd.read_csv(hour/'MT0001R_LDSA.csv.gz',parse_dates=['Time_UTC']);raw=gio[['Time_UTC','Total_LDSA_um2_cm-3']].rename(columns={'Time_UTC':'date','Total_LDSA_um2_cm-3':'value'}).dropna();raw.to_csv(W/'GIO_raw.csv',index=False)
for filename,g in [('monthly_series.csv',monthly(raw.date,raw.value)),('full_monthly_series.csv',raw.groupby(raw.date.dt.to_period('M').dt.start_time).value.mean().rename_axis('date').reset_index())]:
 old=pd.read_csv(A/filename);old=old[~(old.Site.eq('GIO_RB')&old.Metric.eq('LDSA'))];g=g.assign(Site='GIO_RB',Metric='LDSA');pd.concat([old,g],ignore_index=True).to_csv(A/filename,index=False)
w=pd.read_pickle(W/'weekly.pkl');ng=pd.read_pickle(R/'work/gio_record_filter/weekly.pkl');pd.concat([w[~w.Site.eq('GIO_RB')],ng],ignore_index=True).to_pickle(W/'weekly.pkl')
s=f'''suppressPackageStartupMessages(library(openair))
d<-read.csv('{W}/GIO_raw.csv');d$date<-as.POSIXct(d$date,tz='UTC')
z<-as.data.frame(timeAverage(d,avg.time='month',data.thresh=10,progress=FALSE));z<-z[is.finite(z$value),c('date','value')];z$date<-format(z$date,'%Y-%m-%d');z$Site<-'GIO_RB';z$Type<-'RB';z$Metric<-'LDSA';z$Phase<-'MetaFull';z$group_id<-'GIO_RB|LDSA|MetaFull';write.csv(z,'{W}/gio_meta.csv',row.names=FALSE)
''';(W/'gio_monthly.R').write_text(s);subprocess.run(['@RSCRIPT@',str(W/'gio_monthly.R')],cwd=R,check=True)
new=pd.concat([pd.read_csv(R/'work/gio_record_filter/trend_input.csv'),pd.read_csv(W/'gio_meta.csv')]);inp=pd.read_csv(W/'trend_input.csv');inp=pd.concat([inp[~(inp.Site.eq('GIO_RB')&inp.Metric.eq('LDSA'))],new],ignore_index=True);inp.to_csv(W/'trend_input.csv',index=False)
el=pd.read_csv(A/'site_eligibility.csv');nm=len(pd.read_csv(W/'gio_meta.csv'));el.loc[el.Site.eq('GIO_RB'),'valid_months']=nm;el.loc[el.Site.eq('GIO_RB'),'meta_eligible']=nm>=60;el.to_csv(A/'site_eligibility.csv',index=False)
meta=pd.read_csv(A/'meta_monthly_series.csv');pd.concat([meta[~meta.Site.eq('GIO_RB')],pd.read_csv(W/'gio_meta.csv')]).to_csv(A/'meta_monthly_series.csv',index=False)
tr=pd.read_csv(A/'trends.csv');tr=tr[~(tr.Site.eq('GIO_RB')&tr.Metric.eq('LDSA'))];tr=pd.concat([tr,pd.read_csv(R/'work/gio_record_filter/trends.csv')]);tr.to_csv(A/'trends.csv',index=False)
with (W/'base_fit.log').open('w') as f:subprocess.run(['@RSCRIPT@',str(OW/'workbook/trends.R'),str(W/'trend_input.csv'),str(A)],cwd=R,check=True,stdout=f,stderr=subprocess.STDOUT)
policy=json.loads((OLD/'analysis_policy.json').read_text());policy['additional_GIO_2011_exclusion']=json.loads((R/'work/gio_record_filter/verification.json').read_text())['excluded_times_UTC'];(O/'analysis_policy.json').write_text(json.dumps(policy,ensure_ascii=False,indent=2))
oldP=R/'work/netcdf_pollutants'
if not (P/'selected.pkl').exists():(P/'selected.pkl').symlink_to(oldP/'selected.pkl')
s=(oldP/'build_inputs.py').read_text().replace("BASE=R/'outputs/quality_control';OLD=BASE/'analysis';A=R/'outputs/netcdf_pollutants/analysis'", "BASE=R/'outputs/network_analysis';OLD=BASE/'base_analysis';A=BASE/'analysis'").replace("R/'work/quality_control/trend_input.csv'", "R/'work/network_analysis/trend_input.csv'");(P/'build_inputs.py').write_text(s)
with (W/'build_pollutants.log').open('w') as f:subprocess.run([sys.executable,str(P/'build_inputs.py')],cwd=R,check=True,stdout=f,stderr=subprocess.STDOUT)
d=pd.read_csv(P/'new_trend_input.csv');prior=pd.read_csv(OW/'pollutants/new_trend_input.csv');prior={k:g.sort_values('date').reset_index(drop=True) for k,g in prior.groupby('group_id')};same=[]
for k,g in d.groupby('group_id'):
 a=prior.get(k);b=g.sort_values('date').reset_index(drop=True)
 if a is not None and a.date.equals(b.date) and np.array_equal(a.value,b.value,equal_nan=True):same.append(k)
old=pd.read_csv(OW/'pollutants/new_overlap_trends.csv');keys=sorted(d.group_id.unique());oldkeys=old.Site+'|'+old.Metric+'|'+old.Phase
for i in range(3):
 p=P/f'fit_{i}';p.mkdir(exist_ok=True);old[oldkeys.isin(set(keys[i::3])&set(same))].to_csv(p/'trends.csv',index=False)
s=(oldP/'fit.py').read_text().replace("A=R/'outputs/netcdf_pollutants/analysis';OLD=R/'outputs/quality_control/analysis'", "A=R/'outputs/network_analysis/analysis';OLD=R/'outputs/network_analysis/base_analysis'").replace("d=d[d.Phase.str.startswith('Overlap_')]", "# Fit all Full, Before, After and Overlap groups.").replace('full-record and phase ancillary analyses are outside this stage','full-record and phase ancillary analyses included');(P/'fit.py').write_text(s)
(W/'fit_reuse.json').write_text(json.dumps(dict(total=len(keys),reused=int(oldkeys.isin(same).sum()),refit=len(keys)-int(oldkeys.isin(same).sum())),indent=2))
with (W/'all_pollutant_fit.log').open('w') as f:subprocess.run([sys.executable,str(P/'fit.py')],cwd=R,check=True,stdout=f,stderr=subprocess.STDOUT)
(W/'trends_done').write_text('Complete')
print('ALL_TRENDS_DONE',flush=True)
