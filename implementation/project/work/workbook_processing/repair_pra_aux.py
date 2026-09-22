from pathlib import Path
import sys,ast,json
ROOT=Path('@PROJECT_ROOT@');sys.path[:0]=[str(ROOT/'work/ldsa_support/deps'),str(ROOT/'work')]
import pandas as pd,numpy as np
from analyze_cleaned_ldsa import monthly
W=ROOT/'work/workbook_processing';A=ROOT/'outputs/workbook_processing/analysis';p=W/'prepare.py';src=p.read_text()
for f in ast.parse(src).body:
 if isinstance(f,ast.FunctionDef) and f.name in ['rawmonthly','hourly_bp']:exec(compile(ast.Module(body=[f],type_ignores=[]),'prepare_helpers','exec'),globals())
src=src.replace(" if oldsite in polbook.sheet_names:\n  pol=polbook.parse(oldsite)"," polsite='PRA_SUB' if site.startswith('PRA_') else oldsite\n if polsite in polbook.sheet_names:\n  pol=polbook.parse(polsite)")
needle=" # LAN ancillary monthly records have a separate archived source; preserve their dates and values."
src=src.replace(needle," if site.startswith('PRA_'):\n  y0,y1=map(int,site.split('_')[1].split('-'))\n  for name in list(raws):\n   if name not in ['LDSA','Ntotal','Nucleation','Aitken','Accumulation']:raws[name]=raws[name][raws[name].date.dt.year.between(y0,y1)]\n"+needle)
p.write_text(src)
pol=pd.read_excel('@INPUT_ROOT@/data_sources/predictive_models/pollutant_data/all_sites/paired_date_value_measurements.xlsx',sheet_name='PRA_SUB',engine='calamine');native=['LDSA','Ntotal','Nucleation','Aitken','Accumulation'];base=pd.read_csv(W/'trend_input_base.csv',parse_dates=['date']);ms=pd.read_csv(A/'monthly_series.csv',parse_dates=['date']);full=pd.read_csv(A/'full_monthly_series.csv',parse_dates=['date']);bps=pd.read_csv(A/'breakpoints.csv',parse_dates=['Breakpoint','Breakpoint_timestamp']);adds=[];madds=[];fadds=[];badds=[]
mask=base.Site.str.startswith('PRA_')&(~base.Metric.isin(native)|base.Phase.str.startswith('Overlap_'));removed=set(base.loc[mask,'group_id']);base=base[~mask]
for site in ['PRA_2012-2016_SUB','PRA_2017-2019_SUB']:
 y0,y1=map(int,site.split('_')[1].split('-'));anchor=bps.loc[bps.Site.eq(site)&bps.Metric.eq('LDSA'),'Breakpoint_timestamp'].iloc[0];am=anchor.to_period('M').to_timestamp();ld=ms[ms.Site.eq(site)&ms.Metric.eq('LDSA')][['date','value']]
 for i in range(0,len(pol.columns)-1,2):
  metric=pol.columns[i+1];raw=pd.DataFrame({'date':pd.to_datetime(pol.iloc[:,i],errors='coerce'),'value':pd.to_numeric(pol.iloc[:,i+1],errors='coerce')});raw=raw[raw.date.dt.year.between(y0,y1)];d=monthly(raw.date,raw.value);z=d.copy();z['Site']=site;z['Metric']=metric;madds.append(z);b=hourly_bp(raw.date,raw.value,len(d));bm=b.to_period('M').to_timestamp() if b is not None else None
  if b is not None:badds.append(dict(Site=site,Type='SUB',Physical_site='PRA_SUB',Metric=metric,Breakpoint=bm,Breakpoint_timestamp=b,n_months=len(d),offset_months=(bm.year-am.year)*12+bm.month-am.month))
  for phase in ['Full','Before','After']:
   rr=raw if phase=='Full' else raw[raw.date.lt(anchor) if phase=='Before' else raw.date.ge(anchor)];g=rawmonthly(rr.date,rr.value)
   if phase=='Full':z=g.copy();z['Site']=site;z['Metric']=metric;fadds.append(z)
   if len(g)<3:continue
   g['Site']=site;g['Type']='SUB';g['Metric']=metric;g['Phase']=phase;g['group_id']=site+'|'+metric+'|'+phase;adds.append(g)
  pair=ld.merge(d,on='date',suffixes=('_ldsa','_other'))
  if len(pair)>=3:
   for role,col in [('LDSA','value_ldsa'),(metric,'value_other')]:
    g=pair[['date',col]].rename(columns={col:'value'});g['Site']=site;g['Type']='SUB';g['Metric']=role;g['Phase']='Overlap_'+metric;g['group_id']=site+'|'+role+'|Overlap_'+metric;adds.append(g)
for df,added,path in [(ms,madds,A/'monthly_series.csv'),(full,fadds,A/'full_monthly_series.csv'),(bps,[pd.DataFrame(badds)],A/'breakpoints.csv')]:
 df=df[~(df.Site.str.startswith('PRA_')&~df.Metric.isin(native))];pd.concat([df]+added,ignore_index=True).to_csv(path,index=False)
base=pd.concat([base]+adds,ignore_index=True);base.to_csv(W/'trend_input_base.csv',index=False);meta=pd.read_csv(A/'meta_monthly_series.csv',parse_dates=['date']);pd.concat([base,meta[base.columns]],ignore_index=True).to_csv(W/'trend_input.csv',index=False)
for name in ['trends.csv','invalid_percent_CI_raw.csv']:
 t=pd.read_csv(A/name);k=t.Site+'|'+t.Metric+'|'+t.Phase;t[~k.isin(removed)].to_csv(A/name,index=False)
print('PRA repaired',len(removed),'previous fits removed',len(pd.concat(adds).group_id.unique()),'current fit groups')
p=W/'trends.R';s=p.read_text();needle="write.csv(result[bad,],file.path(args[2],'invalid_percent_CI_raw.csv'),row.names=FALSE)";replacement="""rawbad<-result[bad,]
rawpath<-file.path(args[2],'invalid_percent_CI_raw.csv')
if(file.exists(rawpath)){
 previousbad<-read.csv(rawpath);previousbad$percent_CI_valid<-NULL
 keys<-paste(result$Site,result$Metric,result$Phase,sep='|');oldkeys<-paste(previousbad$Site,previousbad$Metric,previousbad$Phase,sep='|')
 keep<-oldkeys %in% keys[!is.finite(result$CI_Lower_Percent)&!is.finite(result$CI_Upper_Percent)]
 rawbad<-rbind(rawbad,previousbad[keep,names(rawbad)])
}
write.csv(rawbad,rawpath,row.names=FALSE)""";s=s.replace(needle,replacement);p.write_text(s)
