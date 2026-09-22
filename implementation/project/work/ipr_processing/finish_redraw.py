from pathlib import Path
import importlib.util,json,shutil,hashlib
import pandas as pd,numpy as np
R=Path('@PROJECT_ROOT@');W=R/'work/ipr_processing';A=R/'outputs/ipr_processing/analysis';OLD=R/'outputs/workbook_processing/analysis';U=Path('@DELIVERY_ROOT@');D=U/'IPR_pollutant_figures';D.mkdir(exist_ok=True)
old=pd.read_csv(OLD/'trends.csv');new=pd.read_csv(W/'fits/trends.csv')
remove=old.Site.eq('IPR_RB')&(old.Metric.isin(['PM2.5','NO2'])|old.Phase.isin(['Overlap_PM2.5','Overlap_NO2']))
assert set(new.Site)=={'IPR_RB'}
combined=pd.concat([old[~remove],new],ignore_index=True);assert not combined.duplicated(['Site','Metric','Phase']).any();combined.to_csv(A/'trends.csv',index=False)
keep=['trends.csv','monthly_series.csv','full_monthly_series.csv','breakpoints.csv','site_eligibility.csv','site_spatial_statistics.csv','size_range_sensitivity.csv','meta_summary.csv','meta_sites.csv','meta_subgroup_test.csv','meta_monthly_series.csv','size_window_manifest.csv']
for f in A.glob('*.csv'):
 if f.name not in keep:f.unlink()
spec=importlib.util.spec_from_file_location('current_redraw',R/'work/workbook_processing/redraw.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.ANA=A;m.WORK=W;m.DEST=D;m.tr=combined;m.ms=pd.read_csv(A/'monthly_series.csv',parse_dates=['date']);m.bp=pd.read_csv(A/'breakpoints.csv');m.bp['Breakpoint']=pd.to_datetime(m.bp.Breakpoint,format='mixed')
pd.testing.assert_frame_equal(combined[combined.Metric.eq('LDSA')&~combined.Phase.str.startswith('Overlap_')].reset_index(drop=True),old[old.Metric.eq('LDSA')&~old.Phase.str.startswith('Overlap_')].reset_index(drop=True))
for name in ['meta_summary.csv','meta_sites.csv','meta_subgroup_test.csv']:
 assert (A/name).read_bytes()==(OLD/name).read_bytes()
save=m.save
def onlyrequested(fig,name,*args,**kwargs):
 if name in ['Figure_2','Figure_3']:save(fig,name,*args,**kwargs)
 else:m.plt.close(fig)
m.save=onlyrequested;m.main12()
m.tr=m.tr[~m.tr.Site.isin(['PRA_SUB','PAR_SUB'])].copy();m.ms=m.ms[~m.ms.Site.isin(['PRA_SUB','PAR_SUB'])].copy();m.main3()
for name in ['exact_overlap_pairs','classification','exact_overlap_breakpoints_min12','correlation_temporal_sensitivity']:
 shutil.copy2(W/(name+'.csv'),A/(name+'.csv'));shutil.copy2(W/(name+'.csv'),D/(name+'.csv'))
newpairs=pd.read_csv(W/'exact_overlap_pairs.csv');previous=pd.read_csv(U/'pollutant_reference_comparisons/exact_overlap_pairs.csv')
pd.testing.assert_frame_equal(newpairs[newpairs.Site.ne('IPR_RB')].reset_index(drop=True),previous[previous.Site.ne('IPR_RB')].reset_index(drop=True),check_exact=False,rtol=1e-12,atol=1e-12)
co=pd.read_csv(W/'correlation_temporal_sensitivity.csv');before=pd.read_csv(U/'pollutant_reference_comparisons/correlation_temporal_sensitivity.csv');compare=before.merge(co,on=['Metric','Strategy'],suffixes=('_before','_after'));compare.to_csv(D/'correlation_before_after.csv',index=False)
ipr=newpairs[newpairs.Site.eq('IPR_RB')&newpairs.Pollutant.isin(['PM2.5','NO2'])];ipr.to_csv(D/'IPR_updated_pairs.csv',index=False)
assert ipr.primary_included.all();assert not newpairs.Site.isin(['PRA_SUB','PAR_SUB']).any()
(W/'render_audit.json').write_text(json.dumps(m.manifest,ensure_ascii=False,indent=2))
oldfig=U/'workbook_analysis_figures/Figure_2.png'
same=oldfig.read_bytes()==(D/'Figure_2.png').read_bytes()
(W/'verification.json').write_text(json.dumps({'Fig2_LDSA_inputs_unchanged':True,'Fig2_identical_PNG':same,'all_non_IPR_pairs_unchanged':True,'excluded_sites':['PRA_SUB','PAR_SUB'],'updated_fit_groups':len(new),'IPR_primary_pairs':ipr[['Pollutant','n_months','coverage']].to_dict('records')},indent=2))
print(compare[compare.Strategy.eq('Exact overlap ≥60 months')].to_string(index=False));print(ipr.to_string(index=False));print('Figure 2 identical PNG:',same)
