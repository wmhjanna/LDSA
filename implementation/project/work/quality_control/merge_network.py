from pathlib import Path
import sys,shutil,json,subprocess
import pandas as pd,numpy as np
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;O=R/'outputs/quality_control';A=O/'analysis';old=O/'workbook/analysis';new=O/'netcdf/analysis'
assert (W/'workbook_done').exists()
for p in old.iterdir():
 if p.is_file() and p.name!='site_spatial_statistics.csv':shutil.copy2(p,A/p.name)
for f in ['trends.csv','invalid_percent_CI_raw.csv','site_eligibility.csv','breakpoints.csv','meta_monthly_series.csv']:
 a=pd.read_csv(old/f);b=pd.read_csv(new/f);pd.concat([a,b],ignore_index=True).to_csv(A/f,index=False)
shutil.copy2(new/'new_station_names.csv',A/'new_station_names.csv')
inp=pd.concat([pd.read_csv(W/'workbook/trend_input.csv'),pd.read_csv(W/'network/trend_input.csv')],ignore_index=True);assert not inp.duplicated(['group_id','date']).any();inp.to_csv(W/'trend_input.csv',index=False)
week=pd.concat([pd.read_pickle(W/'workbook/weekly.pkl'),pd.read_pickle(W/'network/weekly.pkl')],ignore_index=True);week.to_pickle(W/'weekly.pkl')
sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
from analyze_cleaned_ldsa import monthly
mo=[];full=[]
for p in sorted((W/'network/raw').glob('*.csv.gz')):
 d=pd.read_csv(p,parse_dates=['date']);site=p.name.removesuffix('.csv.gz');g=monthly(d.date,d.value);g['Site']=site;g['Metric']='LDSA';mo.append(g)
 g=d.groupby(d.date.dt.to_period('M').dt.to_timestamp()).value.mean().rename_axis('date').reset_index();g['Site']=site;g['Metric']='LDSA';full.append(g)
for name,rows in [('monthly_series',mo),('full_monthly_series',full)]:pd.concat([pd.read_csv(old/(name+'.csv')),*rows],ignore_index=True).to_csv(A/(name+'.csv'),index=False)
with (W/'combined_meta.log').open('w') as f:subprocess.run(['@RSCRIPT@',str(W/'workbook/trends.R'),str(W/'trend_input.csv'),str(A)],cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
policy=json.loads((O/'workbook/analysis_policy.json').read_text());policy['QC']='Nonfinite, negative and exact sentinels missing; official invalid flags excluded only with verified channel/source correspondence; zero and valid warning observations retained; no statistical cutoff';policy['official_QC']='Apply official I/H/M per native cell; 100 overrides invalid; valid warnings and unknown 360 retained. Workbook values masked only after native spectrum match.';policy['network']=65;policy['ancillary_model_scope']='Model and pollutant eligibility are evaluated separately from the spatial network; consult the network_analysis stage for final eligible stations.';policy['netcdf_DF']='MPPD deposition fractions; native channel centres and log10 bin-width weights; PAL source 13 interpreted as DeltaN';policy['pollutants']='Latest German EBAS additions plus retained other-database records; IPR PM2.5 and NO2 update; official invalid values missing, HEL_TR SO2 excluded';(O/'analysis_policy.json').write_text(json.dumps(policy,ensure_ascii=False,indent=2))
print('MERGED',len(pd.read_csv(A/'site_eligibility.csv')),inp.group_id.nunique(),len(pd.read_csv(A/'meta_sites.csv')),flush=True)
