from pathlib import Path
import sys,json,ast,shutil,datetime,hashlib
R=Path('@PROJECT_ROOT@');sys.path[:0]=[str(R/'work/ipr_download/deps'),str(R/'work/ldsa_support/deps'),str(R/'work')]
import numpy as np,pandas as pd,netCDF4
from analyze_cleaned_ldsa import monthly
W=Path(__file__).parent;A=R/'outputs/pal_channel_conversion/analysis';A.mkdir(parents=True,exist_ok=True)
OLD=R/'outputs/network_trends/analysis';OW=R/'work/network_trends';B=Path('@DELIVERY_ROOT@')
for p in OLD.glob('*'):shutil.copy2(p,A/p.name)
oldpath=R/'outputs/network_trends/corrected_LDSA/hourly_ldsa/FI0096G_LDSA.csv.gz'
d=pd.read_csv(oldpath,parse_dates=['Time_UTC']);original=d.copy()
c=pd.read_csv(R/'outputs/network_trends/corrected_LDSA/calculation_parameters/diameter_deposition_weights.csv');c=c[c.Source_ID.eq(13)].sort_values('channel')
reg=pd.read_csv(R/'outputs/network_trends/corrected_LDSA/calculation_parameters/source_qc_registry.csv');rr=reg[reg.Source_ID.eq(13)].iloc[0]
audits=json.loads((OW/'recalculate/source_audit.json').read_text());a=next(a for a in audits if a['record']['file']==rr.file)
p=B/'EBAS_observations'/a['audit']['relative_path']
flags=json.loads((R/'outputs/pnsd_cleaning/PNSD/QC_protocol.json').read_text())['sentinel_values']+[999999.99999]
with netCDF4.Dataset(p) as ds:
 raw=np.asarray(np.ma.filled(ds[rr.variable][:],np.nan),float).T
 tv=ds['time'];ep=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units,calendar=tv.calendar)
 time=pd.to_datetime(np.rint((np.asarray(ds['time_bnds'][:])[:,0]-ep)*86400).astype('int64'),unit='s').round('h')
 raw=np.where(~np.isfinite(raw)|(raw<0)|np.isin(raw,flags),np.nan,raw)
dp=c.diameter_nm.to_numpy();select=(dp>=17)&(dp<=500);z=raw[:,select];valid=np.isfinite(z).sum(axis=1)
surface=np.pi*(dp[select]/1000)**2;df=c[['DF_head','DF_ALV','DF_TB','DF_total']].to_numpy()[select]
mat=surface[:,None]*df
values=np.nan_to_num(z,nan=0)@mat;values[valid==0]=np.nan
legacy=np.nan_to_num(z,nan=0)@(mat*c.delta_log10.to_numpy()[select,None]);legacy[valid==0]=np.nan
cols=['Head_LDSA_um2_cm-3','P_LDSA_um2_cm-3','TB_LDSA_um2_cm-3','Total_LDSA_um2_cm-3']
mask=d.Source_ID.eq(13);idx=pd.Index(time).get_indexer(d.loc[mask,'Time_UTC']);assert (idx>=0).all()
assert np.allclose(d.loc[mask,cols],legacy[idx],equal_nan=True,rtol=1e-10,atol=1e-9)
assert np.array_equal(d.loc[mask,'Valid_bins'],valid[idx])
d.loc[mask,cols]=values[idx]
pnc=np.nansum(z,axis=1);pnc[valid==0]=np.nan;d.loc[mask,'PNC_cm-3']=pnc[idx]
pd.testing.assert_frame_equal(d.loc[~mask],original.loc[~mask]);assert np.allclose(values[:,:3].sum(axis=1),values[:,3],equal_nan=True)
assert np.array_equal(d[cols].isna(),original[cols].isna())
out=A.parent/'FI0096G_LDSA_early_deltaN.csv.gz';d.to_csv(out,index=False,float_format='%.12g')
g=d[['Time_UTC','Total_LDSA_um2_cm-3']].rename(columns={'Time_UTC':'date','Total_LDSA_um2_cm-3':'value'}).dropna().reset_index(drop=True)
g.to_csv(W/'PAL_RB.csv.gz',index=False)
for node in ast.parse((R/'work/workbook_processing/restore_original_sequence.py').read_text()).body:
 if isinstance(node,ast.FunctionDef) and node.name=='index_cp':exec(compile(ast.Module(body=[node],type_ignores=[]),'reference_cp','exec'))
mm=monthly(g.date,g.value);k=index_cp(g,len(mm));bp=g.date.iloc[k];bm=bp.to_period('M').start_time
e=pd.read_csv(A/'site_eligibility.csv');e.loc[e.Site.eq('PAL_RB'),['breakpoint','breakpoint_timestamp']]=[str(bm),str(bp)];e.to_csv(A/'site_eligibility.csv',index=False)
b=pd.read_csv(A/'breakpoints.csv');b.loc[b.Site.eq('PAL_RB'),['Breakpoint','Breakpoint_timestamp']]=[str(bm),str(bp)];b.to_csv(A/'breakpoints.csv',index=False)
parts=[]
for phase,q in [('Full',g),('Before',g.iloc[:k]),('After',g.iloc[k:])]:
 z=q.groupby(q.date.dt.to_period('M').dt.to_timestamp()).value.mean().rename_axis('date').reset_index();z['Site']='PAL_RB';z['Type']='RB';z['Metric']='LDSA';z['Phase']=phase;z['group_id']='PAL_RB|LDSA|'+phase;parts.append(z)
t=pd.read_csv(OW/'trend_input.csv');pd.concat([t[t.Site.ne('PAL_RB')],*parts],ignore_index=True).to_csv(W/'trend_input_base.csv',index=False)
wk=pd.read_pickle(OW/'weekly.pkl');nw=g.groupby(g.date.dt.to_period('W').dt.to_timestamp()).value.mean().rename('LDSA').rename_axis('date').reset_index();nw['Site']='PAL_RB';nw['Type']='RB';pd.concat([wk[wk.Site.ne('PAL_RB')],nw],ignore_index=True).to_pickle(W/'weekly.pkl')
for name in ['trends.csv','invalid_percent_CI_raw.csv']:
 t=pd.read_csv(A/name);t[t.Site.ne('PAL_RB')].to_csv(A/name,index=False)
summary={'assumption':'Source 13 (2000–2006) uses the study assumption of per-bin ΔN; the provider unit definition is unconfirmed.','formula':'LDSA = sum(DeltaN_i * pi * (Dp_i/1000)^2 * DF_i); no delta_log10 factor.','source_file':p.name,'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'updated_rows':int(mask.sum()),'valid_updated_rows':int(d.loc[mask,cols[-1]].notna().sum()),'before_mean':float(original.loc[mask,cols[-1]].mean()),'after_mean':float(d.loc[mask,cols[-1]].mean()),'new_breakpoint':str(bp),'later_rows_identical':True,'missing_pattern_identical':True,'range_nm':[17,500],'source_id':13}
(W/'conversion_audit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));print(json.dumps(summary,ensure_ascii=False,indent=2))
