from pathlib import Path
import sys,os,ast,json,subprocess,hashlib
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;D=Path('@DELIVERY_ROOT@/GIO_station_diagnostics');D.mkdir(exist_ok=True)
sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')];os.environ['MPLCONFIGDIR']=str(W/'mplconfig')
import pandas as pd,numpy as np
from analyze_cleaned_ldsa import monthly
source=R/'outputs/record_exclusions/netcdf/hourly_ldsa/MT0001R_LDSA.csv.gz';d=pd.read_csv(source,parse_dates=['Time_UTC']);timestamps=pd.date_range('2011-07-22 10:00','2011-07-22 13:00',freq='h');mask=d.Time_UTC.isin(timestamps);assert mask.sum()==4 and d.loc[mask,'Total_LDSA_um2_cm-3'].notna().all();removed=d[mask].copy();removed.to_csv(W/'excluded_four_hours.csv',index=False)
cols=[c for c in d if 'LDSA_' in c or c=='PNC_cm-3'];d.loc[mask,cols]=np.nan;d.loc[mask,'Valid_bins']=0;d.loc[mask,'User_case_excluded']=True;d.to_csv(W/'MT0001R_LDSA.csv.gz',index=False)
z=d[['Time_UTC','Total_LDSA_um2_cm-3']].rename(columns={'Time_UTC':'date','Total_LDSA_um2_cm-3':'value'}).dropna().reset_index(drop=True)
for node in ast.parse((R/'work/workbook_processing/restore_original_sequence.py').read_text()).body:
 if isinstance(node,ast.FunctionDef) and node.name=='index_cp':exec(compile(ast.Module(body=[node],type_ignores=[]),'existing_cp','exec'))
mm=monthly(z.date,z.value);k=index_cp(z,len(mm));bt=z.date.iloc[k];bm=bt.to_period('M').start_time
parts=[]
for phase,q in [('Full',z),('Before',z.iloc[:k]),('After',z.iloc[k:])]:
 g=q.groupby(q.date.dt.to_period('M').dt.start_time).value.mean().rename_axis('date').reset_index();g=g.assign(Site='GIO_RB',Type='RB',Metric='LDSA',Phase=phase,group_id='GIO_RB|LDSA|'+phase);parts.append(g)
pd.concat(parts).to_csv(W/'trend_input.csv',index=False)
s=(R/'work/quality_control/workbook/trends.R').read_text().split("\ne<-read.csv(file.path(args[2],'site_eligibility.csv'))",1)[0];(W/'fit.R').write_text(s)
with (W/'fit.log').open('w') as f:subprocess.run(['@RSCRIPT@',str(W/'fit.R'),str(W/'trend_input.csv'),str(W)],cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
# Preserve the exact current S6 template, only substitute this site's recalculated inputs.
s=(R/'work/record_exclusions/render_s6.py').read_text().split('exec(compile(src,\'s6_only\'',1)[0]
# Everything before the final draw creates template functions only.
ns={'__file__':str(R/'work/record_exclusions/render_s6.py')};exec(compile(s,'latest_S6_style','exec'),ns);m=ns['m']
m.DEST=D;m.sp=m.sp[m.sp.Site.eq('GIO_RB')].copy();m.tr=pd.read_csv(W/'trends.csv');m.bp=pd.DataFrame([dict(Site='GIO_RB',Type='RB',Metric='LDSA',Breakpoint=bm,Breakpoint_timestamp=bt)])
w=z.groupby(z.date.dt.to_period('W').dt.start_time).value.mean().rename('LDSA').rename_axis('date').reset_index().assign(Site='GIO_RB',Type='RB');m.wk=w;w.to_pickle(W/'weekly.pkl');m.bp.to_csv(W/'breakpoints.csv',index=False)
exec(compile(ns['src'],'single_S6','exec'),m.__dict__);m.s456()
p=D/'Figure_S6/Figure_S6_GIO_RB.png';p.rename(D/p.name);(D/'Figure_S6').rmdir()
assert np.isclose(w.loc[w.date.eq('2011-07-18'),'LDSA'].iloc[0],67.42043701851463)
verify=dict(excluded_times_UTC=[str(x) for x in timestamps],source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),new_week_mean=float(w.loc[w.date.eq('2011-07-18'),'LDSA'].iloc[0]),breakpoint=str(bt),figure=str(D/p.name),scope='GIO_RB S6 only; other figures and full-network results not recomputed',annotation=ns['audit'])
(W/'verification.json').write_text(json.dumps(verify,ensure_ascii=False,indent=2));print(json.dumps(verify,ensure_ascii=False,indent=2))
from PIL import Image
im=Image.open(D/p.name);im.thumbnail((1400,1100));im.save(W/'preview.png')
