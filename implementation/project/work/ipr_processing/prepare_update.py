from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd
R=Path('@PROJECT_ROOT@');W=R/'work/ipr_processing';OLD=R/'outputs/workbook_processing/analysis';A=R/'outputs/ipr_processing/analysis';A.mkdir(parents=True,exist_ok=True)
SRC=Path('@INPUT_ROOT@/documents/IPR_EBAS_PM25_NO2_2009-2019.xlsx')
COPY_FILES=['trends.csv','monthly_series.csv','full_monthly_series.csv','breakpoints.csv','site_eligibility.csv','site_spatial_statistics.csv','size_range_sensitivity.csv','meta_summary.csv','meta_sites.csv','meta_subgroup_test.csv','meta_monthly_series.csv','size_window_manifest.csv']
for name in COPY_FILES:shutil.copy2(OLD/name,A/name)
for f in ['analysis_policy.json']:shutil.copy2(OLD.parent/f,A.parent/f)
ms=pd.read_csv(OLD/'monthly_series.csv',parse_dates=['date']);full=pd.read_csv(OLD/'full_monthly_series.csv',parse_dates=['date']);bp=pd.read_csv(OLD/'breakpoints.csv');oldinput=pd.read_csv(R/'work/workbook_processing/trend_input.csv',parse_dates=['date'])
def selected(d):return d.Site.eq('IPR_RB')&d.Metric.isin(['PM2.5','NO2'])
ms_keep=ms[~selected(ms)];full_keep=full[~selected(full)];bp_keep=bp[~selected(bp)]
ld=ms[ms.Site.eq('IPR_RB')&ms.Metric.eq('LDSA')][['date','value']]
anchor=pd.Timestamp(bp[bp.Site.eq('IPR_RB')&bp.Metric.eq('LDSA')].Breakpoint.iloc[0])
def monthly(d):
 day=d.groupby(d.date.dt.floor('D')).value.mean();g=day.resample('MS').agg(['mean','count']);return g[g['count']>=14]['mean'].rename('value').rename_axis('date').reset_index()
parts=[];mo=[];fu=[];bprows=[];audit=[]
for metric in ['PM2.5','NO2']:
 d=pd.read_pickle(W/(metric+'.pkl'));t=pd.to_datetime(d.iloc[:,0]);v=pd.to_numeric(d.iloc[:,2],errors='coerce');bad=~np.isfinite(v)|(v<0);v=v.mask(bad);factor=46/14 if metric=='NO2' else 1
 raw=pd.DataFrame({'date':t,'end':pd.to_datetime(d.iloc[:,1]),'original_value':pd.to_numeric(d.iloc[:,2],errors='coerce'),'value':v*factor,'original_unit':'µg N m-3' if metric=='NO2' else 'µg m-3','unit':'µg m-3','quality_flag':d.iloc[:,4 if metric=='NO2' else 3],'source_id':d.iloc[:,-1]})
 raw.to_csv(A.parent/(f'IPR_RB_{metric}_observations.csv.gz'),index=False)
 x=raw[['date','value']].dropna();assert not x.date.duplicated().any()
 m=monthly(x);m['Site']='IPR_RB';m['Metric']=metric;mo.append(m)
 for phase in ['Full','Before','After']:
  z=x if phase=='Full' else x[x.date.lt(anchor) if phase=='Before' else x.date.ge(anchor)]
  z=z.groupby(z.date.dt.to_period('M').dt.to_timestamp()).value.mean().rename_axis('date').reset_index()
  if phase=='Full':q=z.copy();q['Site']='IPR_RB';q['Metric']=metric;fu.append(q)
  if len(z)<3:continue
  z=z.assign(Site='IPR_RB',Type='RB',Metric=metric,Phase=phase,group_id='IPR_RB|'+metric+'|'+phase);parts.append(z)
 pair=ld.merge(m[['date','value']],on='date',suffixes=('_ldsa','_other'))
 for role,col in [('LDSA','value_ldsa'),(metric,'value_other')]:
  z=pair[['date',col]].rename(columns={col:'value'}).assign(Site='IPR_RB',Type='RB',Metric=role,Phase='Overlap_'+metric,group_id='IPR_RB|'+role+'|Overlap_'+metric);parts.append(z)
 if len(m)>=20:
  y=x.value.to_numpy(float);y-=y.mean();cs=np.cumsum(y);k=np.arange(2,len(y)-1);idx=int(k[np.argmax(cs[k-1]**2/k+(cs[-1]-cs[k-1])**2/(len(y)-k))]);stamp=x.date.iloc[idx];bm=stamp.to_period('M').to_timestamp();bprows.append(dict(Site='IPR_RB',Type='RB',Physical_site='IPR_RB',Metric=metric,Breakpoint=bm,Breakpoint_timestamp=stamp,n_months=len(m),offset_months=(bm.year-anchor.year)*12+bm.month-anchor.month))
 audit.append(dict(Metric=metric,source_rows=len(raw),valid_rows=len(x),negative_removed=int((raw.original_value<0).sum()),unit_conversion_factor=factor,paired_months=len(pair),paired_start=str(pair.date.min()),paired_end=str(pair.date.max()),valid_months=len(m),quality_flag_policy='Published values retained except missing, nonfinite and negative concentrations; nonzero flags not automatically rejected'))
new=pd.concat(parts,ignore_index=True);new.to_csv(W/'changed_trend_input.csv',index=False)
mask=oldinput.Site.eq('IPR_RB')&((oldinput.Metric.isin(['PM2.5','NO2']))|(oldinput.Phase.isin(['Overlap_PM2.5','Overlap_NO2'])))
pd.concat([oldinput[~mask],new],ignore_index=True).to_csv(W/'trend_input.csv',index=False)
pd.concat([ms_keep]+mo,ignore_index=True).to_csv(A/'monthly_series.csv',index=False)
pd.concat([full_keep]+fu,ignore_index=True).to_csv(A/'full_monthly_series.csv',index=False)
pd.concat([bp_keep,pd.DataFrame(bprows)],ignore_index=True).to_csv(A/'breakpoints.csv',index=False)
source={'source_workbook':str(SRC),'sha256':hashlib.file_digest(SRC.open('rb'),'sha256').hexdigest(),'site':'IPR_RB','source_station':'IT0004R','aggregation':'Use UTC sample start; retain native resolution; equal-weight valid daily means for paired monthly data; minimum 14 valid days per month; no interpolation or expansion of daily data','primary_pair_rule':'At least 60 common valid months and 50% coverage','NO2_gaps':'2009, 2012 and 2013 remain missing; no previous-source backfill','Fig3_exclusions':['PAR_SUB','PRA_SUB'],'fig2':'All LDSA inputs unchanged','metrics':audit}
(A.parent/'source_update.json').write_text(json.dumps(source,ensure_ascii=False,indent=2));print(json.dumps(audit,ensure_ascii=False,indent=2))
policy=json.loads((A.parent/'analysis_policy.json').read_text());policy['IPR_PM25_NO2_update']=source;(A.parent/'analysis_policy.json').write_text(json.dumps(policy,ensure_ascii=False,indent=2))
# Fit changed groups using the specified openair options and per-group seeds.
r=(R/'work/workbook_processing/trends.R').read_text().split("e<-read.csv(file.path(args[2],'site_eligibility.csv'))")[0]
(W/'fit_changed.R').write_text(r);(W/'fits').mkdir(exist_ok=True)
