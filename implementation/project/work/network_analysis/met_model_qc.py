from pathlib import Path
import pandas as pd,numpy as np,json
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;M=W/'model';d=pd.read_pickle(M/'weekly_model_reference.pkl');audit=[]
# Ispra's source workbook has no units, and mixes temperature around 280 with values around 10.
# No defensible Kelvin conversion can be asserted from this workbook alone: quarantine its TEMP feature.
mask=d.Site.eq('IPR_RB')&d.TEMP.notna();audit.append(dict(Site='IPR_RB',Metric='TEMP',excluded_weeks=int(mask.sum()),reason='Legacy workbook TEMP has no unit metadata and mixes two incompatible scales; whole IPR TEMP predictor quarantined, no inferred Kelvin conversion.'));d.loc[mask,'TEMP']=np.nan
# Do source-level screening before weekly means. Preserve raw native records and provider flags separately.
s=pd.read_pickle(W/'meteorology/selected.pkl');s['additional_physical_invalid']=False
rules={'RAD':lambda x:x.lt(0),'WS':lambda x:x.lt(0),'RAIN':lambda x:x.lt(0),'WD':lambda x:x.lt(0)|x.gt(360),'RH':lambda x:x.lt(0)|x.gt(100)}
for metric,rule in rules.items():
 bad=s.metric.eq(metric)&s.value.notna()&rule(s.value)
 for site,g in s[bad].groupby('Site'):
  audit.append(dict(Site=site,Metric=metric,excluded_native_records=len(g),reason='Additional model input domain screen: RAD/WS/RAIN >=0, WD 0..360, RH 0..100; separate from provider flags. Relative humidity >100 is excluded as outside nominal RH range, not asserted to be a provider-invalid record.'))
 s.loc[bad,'additional_physical_invalid']=True
s[s.additional_physical_invalid].to_csv(W/'meteorology/additional_model_exclusions.csv.gz',index=False)
s.loc[s.additional_physical_invalid,'value']=np.nan
for (site,metric),q in s.groupby(['Site','metric']):
 v=q.groupby(q.date.dt.to_period('W').dt.start_time).value.mean();idx=d.Site.eq(site);d.loc[idx,metric]=d.loc[idx,'date'].map(v)
d.to_pickle(M/'weekly_model_reference.pkl');pd.DataFrame(audit).to_csv(M/'additional_meteorology_QC.csv',index=False);print(pd.DataFrame(audit).to_string(index=False))
