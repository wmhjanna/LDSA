from pathlib import Path
import sys,json,pickle,hashlib
import pandas as pd,numpy as np
ROOT=Path('@PROJECT_ROOT@')
BASE=Path('@DELIVERY_ROOT@')
W=ROOT/'work/quality_control/spatial';W.mkdir(parents=True,exist_ok=True)
O=BASE/'provider_quality_control_analysis';T=O/'tables';T.mkdir(parents=True,exist_ok=True)
S=O/'diameter_coverage';S.mkdir(exist_ok=True)
A=ROOT/'outputs/quality_control/analysis';A.mkdir(parents=True,exist_ok=True)
N=ROOT/'outputs/quality_control/netcdf'
OLD=ROOT/'work/quality_control/workbook'
src=BASE/'Figure1_65_station_tables/location_LDSA-65sites_site_abbreviationexport.xlsx'
tab=pd.read_excel(src,sheet_name='location');assert len(tab)==65
names=json.loads((ROOT/'work/station_registry/mapping.json').read_text());codes={x['new']:x['old'] for x in names}
reg=pd.read_csv(N/'calculation_parameters/source_qc_registry.csv').set_index('Source_ID')
assert len(list((OLD/'series').glob('*_series.pkl')))==34,'Wait for complete workbook preparation'
statistics=[];periods=[];detail=[];changes=[];sources=[];coverage=[];configs=[];spatial=[];attribution=[]
for i,row in tab.iterrows():
    site=row.Acronym;code=codes.get(site,site);typ={'U':'UB','F':'RB','M':'RB','P':'RB'}.get(row.type,row.type)
    if i<34:
        file=OLD/'series'/f'{site}_series.pkl';s=pd.read_pickle(file);t=pd.to_datetime(s.Time)
        cols=[f'{r}_LDSA_primary' for r in ['HA','ALV','TB','Total']]
        df_description='workbook-format deposition fractions'
    else:
        file=N/'hourly_ldsa'/f'{code}_LDSA.csv.gz'
        cols=[f'{r}_LDSA_um2_cm-3' for r in ['Head','P','TB','Total']]
        s=pd.read_csv(file,usecols=['Time_UTC','Source_ID']+cols);t=pd.to_datetime(s.Time_UTC)
        df_description='netCDF-format inputvalidationMPPD DF'
    values=s[cols].replace([np.inf,-np.inf],np.nan);valid=t.notna()&values.iloc[:,3].notna()
    mask=valid&t.dt.year.between(2015,2019);fallback=not mask.any()
    if fallback:
        years=sorted(t[valid].dt.year.unique());assert years
        year=min(years,key=lambda y:(max(2015-y,y-2019,0),-y));mask=valid&t.dt.year.eq(year);period=str(year)
    else:period='2015-2019'
    x=values.loc[mask];assert x.notna().all().all(),site;out=[];wide={}
    for col,region in zip(cols,['HA','ALV','TB','Total']):
        v=x[col].to_numpy();q=np.quantile(v,[.25,.5,.75,.95]);stats=[float(v.mean()),float(q[1]),f'{q[0]:.2f}-{q[2]:.2f}',float(q[3])];out.extend(stats)
        wide[region]=stats[0]
        detail.append(dict(Site=site,region=region,Period=period,n=len(v),Mean=stats[0],Median=stats[1],P25=float(q[0]),P75=float(q[2]),P95=stats[3]))
        if region=='Total':wide.update(Mean=stats[0],Median=stats[1],P25=float(q[0]),P75=float(q[2]),P95=stats[3])
    dates=t[mask];periods.append(dict(Site=site,code=code,Period=period,fallback=bool(fallback),actual_start=str(dates.min()),actual_end=str(dates.max()),available_years=','.join(map(str,sorted(dates.dt.year.unique()))),valid_records=int(mask.sum()),valid_months=int(dates.dt.to_period('M').nunique()),DF=df_description,source=str(file)))
    statistics.append(out)
    if i<34:
        previous=ROOT/'work/workbook_processing/series'/f'{site}_series.pkl';before=pd.read_pickle(previous);bt=pd.to_datetime(before.Time)
    else:
        previous=ROOT/'outputs/pal_channel_conversion/FI0096G_LDSA_early_deltaN.csv.gz' if code=='FI0096G' else ROOT/'outputs/network_trends/corrected_LDSA/hourly_ldsa'/f'{code}_LDSA.csv.gz'
        before=pd.read_csv(previous,usecols=['Time_UTC']+cols);bt=pd.to_datetime(before.Time_UTC)
    bv=before[cols].replace([np.inf,-np.inf],np.nan);bm=bt.notna()&bv.iloc[:,3].notna()&bt.dt.year.between(2015,2019)
    if not bm.any():
        years=sorted(bt[bv.iloc[:,3].notna()].dt.year.unique());by=min(years,key=lambda y:(max(2015-y,y-2019,0),-y));bm=bv.iloc[:,3].notna()&bt.dt.year.eq(by)
    bstats=[]
    for col in cols:
        v=bv.loc[bm,col].to_numpy();q=np.quantile(v,[.25,.5,.75,.95]);bstats.extend([float(v.mean()),float(q[1]),f'{q[0]:.2f}-{q[2]:.2f}',float(q[3])])
    for k in [0,1,3,4,5,7,8,9,11,12,13,15]:
        old=float(row.iloc[11+k]);matched=float(bstats[k]);new=float(out[k]);qc_delta=new-matched;harm_delta=matched-old
        reason=[]
        if not np.isclose(old,matched,rtol=1e-9,atol=1e-7):reason.append('source_harmonisation')
        if not np.isclose(matched,new,rtol=1e-9,atol=1e-7):reason.append('provider_quality_control')
        attribution.append(dict(Site=site,field=tab.columns[11+k],reference_table=old,latest_source_before_QC=matched,after_QC=new,source_harmonization_delta=harm_delta,official_QC_delta=qc_delta,reason='; '.join(reason) if reason else 'Unchanged within numerical tolerance',comparison_source=str(previous)))
    del before,bv
    for k,new in enumerate(out):
        old=row.iloc[11+k];same=old==new if isinstance(new,str) else np.isclose(float(old),new,rtol=1e-9,atol=1e-7)
        changes.append(dict(Site=site,field=tab.columns[11+k],before=old,after=new,changed=not bool(same),Period=period))
    grids=[]
    if i<34:
        with open(OLD/'pnsd'/f'{site}.pkl','rb') as f:z=pickle.load(f)
        masks,counts=np.unique(z['structural'][mask],axis=0,return_counts=True)
        for m,n in zip(masks,counts):
            dp=np.asarray(z['diameters'])[~m];assert len(dp)>0
            grids.append((float(dp.min()),float(dp.max()),int(n),'native structural grid'))
        del z
    else:
        for sid,n in s.loc[mask].Source_ID.value_counts().items():
            z=reg.loc[sid];assert z.code==code
            grids.append((float(z.native_low_nm),float(z.native_high_nm),int(n),f'Source_ID={sid}'))
    keep=all(a<=17+1e-9 and b>=500-1e-9 for a,b,n,source in grids) and site!='LND_TR'
    reasons=[]
    if site=='LND_TR':reasons.append('LND_TR integration lower bound is 31.65 nm')
    if any(a>17+1e-9 for a,b,n,source in grids):reasons.append('Native lower diameter exceeds 17 nm')
    if any(b<500-1e-9 for a,b,n,source in grids):reasons.append('Native upper diameter is below 500 nm')
    ranges='; '.join(f'{a:g}–{b:g}' for a,b in sorted(set((a,b) for a,b,n,source in grids)))
    coverage.append(dict(Site=site,code=code,Period=period,native_ranges_nm=ranges,keep=bool(keep),reason='; '.join(reasons) if reasons else 'Native configurations cover 17–500 nm'))
    for a,b,n,source in grids:configs.append(dict(Site=site,code=code,Period=period,source=source,native_low_nm=a,native_high_nm=b,valid_records=n,configuration_covers_17_500=bool(a<=17+1e-9 and b>=500-1e-9)))
    spatial.append(dict(Site=site,Type=typ,Period=period,n=len(x),Physical_site=site,lower_nm=31.65 if site=='LND_TR' else 17,upper_nm=500,Country=row.Country,Record_Start=str(t[valid].min()),Record_End=str(t[valid].max()),native_ranges_nm=ranges,**wide))
    sources.append(dict(path=str(file),sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
    print(site,period,len(x),round(out[12],6),'strict',keep,flush=True)
assert len(statistics)==65
pd.DataFrame(spatial).to_csv(A/'site_spatial_statistics.csv',index=False)
pd.DataFrame(periods).to_csv(T/'spatial_analysis_periods.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(detail).to_csv(T/'spatial_statistics_details.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(changes).to_csv(T/'spatial_statistics_changes.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(attribution).to_csv(T/'spatial_statistics_comparison.csv',index=False,encoding='utf-8-sig')
cov=pd.DataFrame(coverage);cfg=pd.DataFrame(configs)
cov.to_csv(S/'coverage_selection.csv',index=False,encoding='utf-8-sig');cfg.to_csv(S/'configuration_coverage.csv',index=False,encoding='utf-8-sig')
keep=cov.loc[cov.keep,'Site'].tolist();assert len(keep)==49
upper_excluded=set(cfg.loc[cfg.native_high_nm<500-1e-9,'Site']);assert len(upper_excluded)==12
payload=dict(statistics=statistics,sites=tab.Acronym.tolist(),keep=keep,source=str(src))
(W/'payload.json').write_text(json.dumps(payload,ensure_ascii=False,allow_nan=False))
(W/'coverage_selection.json').write_text(json.dumps(dict(keep=keep,excluded=cov.loc[~cov.keep,'Site'].tolist(),upper_excluded=sorted(upper_excluded)),ensure_ascii=False,indent=2))
(W/'provenance.json').write_text(json.dumps(dict(source_workbook=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),sites=65,fallback_sites=[r for r in periods if r['fallback']],sources=sources,changed_sites=sorted(set(c['Site'] for c in changes if c['changed']))),ensure_ascii=False,indent=2))
print('SPATIAL_READY',len(keep),flush=True)
