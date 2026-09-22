from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
ROOT=Path('@PROJECT_ROOT@')
STAGE=ROOT/'work/ldsa_support';sys.path.insert(0,str(STAGE/'deps'))
import ruptures as rpt
OUT=ROOT/'outputs/ldsa_cleaned_reanalysis';ANA=OUT/'analysis';ANA.mkdir(parents=True,exist_ok=True)
BASE=Path('@INPUT_ROOT@/data_sources')
REV={'SCH_RB':'SCH_SUB','WAL_RB':'WAL_SUB','ZUG_RB':'ZUG_SUB'}
METRICS={'LDSA':'LDSA','PNC_Nucleation_17_25nm':'Nucleation','PNC_Aitken_25_100nm':'Aitken','PNC_Accumulation_100_500nm':'Accumulation','PNC_Ntotal_17_500nm':'Ntotal'}
POLL=['PM2.5','PM10','PM1','BC','NO2','SO2','O3','CO','NO']

def monthly(t,v):
    d=pd.DataFrame({'date':pd.to_datetime(t,errors='coerce'),'value':pd.to_numeric(v,errors='coerce')}).replace([np.inf,-np.inf],np.nan).dropna()
    day=d.groupby(d.date.dt.floor('D')).value.mean()
    m=day.resample('MS').agg(['mean','count']);m.loc[m['count']<14,'mean']=np.nan
    return m['mean'].dropna().rename('value').rename_axis('date').reset_index()

def bp(d):
    if len(d)<20:return None
    k=rpt.KernelCPD(kernel='linear',min_size=2,jump=1).fit(d.value.to_numpy()).predict(n_bkps=1)[0]
    return pd.Timestamp(d.date.iloc[k])

def main():
    files=sorted(STAGE.glob('*_series.pkl'))
    spatial=[];sensitivity=[];breaks=[];elig=[];groups=[];weekly_all=[];monthly_all=[]
    periods=pd.read_csv(ROOT/'outputs/german_RB_revision/site_spatial_statistics.csv').set_index('Site').Period
    polbook=pd.ExcelFile(BASE/'predictive_models/pollutant_data/all_sites/paired_date_value_measurements.xlsx',engine='calamine')
    for file in files:
        site=file.name.removesuffix('_series.pkl');old=REV.get(site,site);typ=site.rsplit('_',1)[-1]
        s=pd.read_pickle(file).sort_values('Time');s=s.loc[s.Time.notna()]
        lo,hi=[int(x) for x in str(periods[site]).replace('–','-').split('-')]
        spatial_s=s[s.Time.dt.year.between(lo,hi)]
        primary=spatial_s['LDSA'].dropna()
        spatial.append(dict(Site=site,Type=typ,Period=f'{lo}-{hi}',n=len(primary),Mean=primary.mean(),Median=primary.median(),P25=primary.quantile(.25),P75=primary.quantile(.75),P95=primary.quantile(.95),**{r:spatial_s.loc[spatial_s.LDSA.notna(),f'{r}_LDSA_17_500nm'].mean() for r in ['HA','TB','ALV']}))
        common=spatial_s[[f'Total_LDSA_{a}_{b}nm' for a,b in [(17,500),(10,500),(10,800)]]].notna().all(axis=1)
        sensitivity_s=spatial_s.loc[common]
        for low,high in [(17,500),(10,500),(10,800)]:
            for reg in ['Total','HA','TB','ALV']:
                c=f'{reg}_LDSA_{low}_{high}nm';x=sensitivity_s[c].dropna()
                row=dict(Site=site,Type=typ,Period=f'{lo}-{hi}',Range=f'{low}-{high}',Region=reg,n=len(x),Mean=x.mean(),Median=x.median(),P25=x.quantile(.25),P75=x.quantile(.75),P95=x.quantile(.95))
                sensitivity.append(row)
        ms={name:monthly(s.Time,s[c]) for c,name in METRICS.items()}
        if old in polbook.sheet_names:
            p=polbook.parse(old)
            for i in range(0,len(p.columns)-1,2):
                metric=str(p.columns[i+1]).strip()
                if metric in POLL:ms[metric]=monthly(p.iloc[:,i],p.iloc[:,i+1])
        anchor=bp(ms['LDSA'])
        n=len(ms['LDSA']);span=0 if n==0 else (ms['LDSA'].date.max().year-ms['LDSA'].date.min().year)*12+ms['LDSA'].date.max().month-ms['LDSA'].date.min().month+1
        elig.append(dict(Site=site,Type=typ,valid_months=n,span_months=span,coverage=n/span if span else 0,meta_eligible=n>=60,breakpoint=anchor))
        for metric,d in ms.items():
            d=d.copy();d['Site']=site;d['Metric']=metric;monthly_all.append(d)
            b=bp(d)
            if b is not None:breaks.append(dict(Site=site,Type=typ,Metric=metric,Breakpoint=b,n_months=len(d),offset_months=((b.year-anchor.year)*12+b.month-anchor.month) if anchor is not None else np.nan))
            for phase in ['Full','Before','After']:
                if phase!='Full' and anchor is None:continue
                g=d if phase=='Full' else d[d.date.lt(anchor) if phase=='Before' else d.date.ge(anchor)]
                if len(g)<3:continue
                g=g[['date','value']].copy();g['Site']=site;g['Type']=typ;g['Metric']=metric;g['Phase']=phase;g['group_id']=f'{site}|{metric}|{phase}';groups.append(g)
            if metric!='LDSA':
                pair=ms['LDSA'].merge(d[['date','value']],on='date',suffixes=('_ldsa','_other'))
                if len(pair)>=3:
                    for role,c in [('LDSA','value_ldsa'),(metric,'value_other')]:
                        g=pair[['date',c]].rename(columns={c:'value'});g['Site']=site;g['Type']=typ;g['Metric']=role;g['Phase']='Overlap_'+metric;g['group_id']=f'{site}|{role}|Overlap_{metric}';groups.append(g)
        # Count unique hourly means, avoiding duplicate timestamps inflating completeness.
        hourly=s.set_index('Time').select_dtypes('number').resample('h').mean()
        weekly=hourly.resample('W-MON',label='left',closed='left').mean()
        counts=hourly.resample('W-MON',label='left',closed='left').count()
        weekly=weekly.where(counts>=126)
        weekly=weekly.reset_index().rename(columns={'Time':'date'});weekly['Site']=site;weekly['Type']=typ
        weekly_all.append(weekly)
        print(site,'months',n,'breakpoint',anchor,'weeks',weekly.LDSA.notna().sum(),flush=True)
    pd.DataFrame(spatial).to_csv(ANA/'site_spatial_statistics.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(ANA/'size_range_sensitivity.csv',index=False)
    pd.DataFrame(breaks).to_csv(ANA/'breakpoints.csv',index=False)
    pd.DataFrame(elig).to_csv(ANA/'site_eligibility.csv',index=False)
    pd.concat(groups,ignore_index=True).to_csv(STAGE/'trend_input.csv',index=False)
    pd.concat(weekly_all,ignore_index=True).to_pickle(STAGE/'weekly.pkl')
    pd.concat(monthly_all,ignore_index=True).to_csv(ANA/'monthly_series.csv',index=False)

if __name__=='__main__':main()
