from pathlib import Path
import sys,json,hashlib,collections,datetime,shutil
ROOT=Path('@PROJECT_ROOT@')
sys.path[:0]=[str(ROOT/'work/ipr_download/deps'),str(ROOT/'work/ldsa_support/deps')]
import numpy as np,pandas as pd,netCDF4
W=ROOT/'work/network_trends/recalculate'
D=Path('@DELIVERY_ROOT@/EBAS_observations')
O=ROOT/'outputs/quality_control/netcdf';O.mkdir(exist_ok=True)
R=ROOT/'work/deposition/mppd'
for sub in ['hourly_ldsa','calculation_parameters']: (O/sub).mkdir(exist_ok=True)
sources=json.loads((W/'source_audit.json').read_text())
details={(r['file'],r['var']):r for r in json.loads((W/'metadata_details.json').read_text())}
subgrids={(r['file'],r['var']):r['groups'] for r in json.loads((W/'subgrid_groups.json').read_text()) if r['file'].startswith(('FR0027U','FR0030R'))}
stations=[r for r in json.loads((R/'inputs/stations.json').read_text()) if r['scope']=='europe31']
df=pd.read_csv(R/'results/direct_europe31/DF_by_station_grid.csv')
flags=json.loads((ROOT/'outputs/pnsd_cleaning/PNSD/QC_protocol.json').read_text())['sentinel_values']+[999999.99999]
ft=pd.read_csv(ROOT/'work/provider_flag_audit/ebas_flags.csv',header=None,quotechar="'",names=['flag','description','status']);invalid_codes=ft[ft.status.isin(['I','H','M'])].flag.to_numpy()
registry=[];excluded=[];coefficients=[];checks=[];statrows=[];overview=[];location=[]
manifest={r['file']:r for r in json.loads((D/'download_manifest.json').read_text())}
regions=['Head','P','TB','Total'];dfnames=['DF_head','DF_ALV','DF_TB','DF_total']
types={'F':'Forest','P':'Pristine','M':'Mountain','U':'Urban','RB':'Rural background','RB, Coast':'Rural background, coastal'}
def mdget(md,key):
    # Older EBAS metadata uses e.g. 'Volume std. pressure'; newer files use ebas_volume_std_pressure.
    # Ignore punctuation for the lookup without changing the underlying concentration data.
    import re
    canon=lambda s:re.sub('[^a-z0-9]','',s.lower()).removeprefix('ebas')
    target=canon(key)
    return next((v for k,v in md.items() if canon(k)==target),'')
def delta(dp):
    assert np.all(np.diff(dp)>0)
    d=np.diff(np.log10(dp));return np.r_[d,d[-1]]
for st in stations:
    code=st['code'];parts=[]
    for x in [r for r in sources if r['record']['code']==code]:
        rec=x['record'];au=x['audit'];md=details[(rec['file'],rec['variable'])]['metadata'];lev=str(mdget(md,'Data level'));matrix=mdget(md,'Matrix')
        if lev!='2' or 'non_volatile' in matrix or 'nonvol' in mdget(md,'Method ref'):
            excluded.append(dict(code=code,file=rec['file'],variable=rec['variable'],reason='Non-volatile spectrum' if ('non' in matrix or 'nonvol' in mdget(md,'Method ref')) else 'NRT / not confirmed Level 2'));continue
        path=D/au['relative_path']
        with netCDF4.Dataset(path) as ds:
            tv=ds['time'];cal=getattr(tv,'calendar','standard');tb=ds['time_bnds'][:].filled(np.nan) if np.ma.isMaskedArray(ds['time_bnds'][:]) else np.asarray(ds['time_bnds'][:])
            fac={'days':86400,'hours':3600,'minutes':60,'seconds':1}[tv.units.split()[0]];ep=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units,calendar=cal)
            seconds=np.rint((tb-ep)*fac).astype('int64');rawtime=pd.to_datetime(seconds[:,0],unit='s');time=rawtime.round('h');mask=(time>='1996-01-01')&(time<'2024-01-01');time=time[mask];rawtime=rawtime[mask]
            v=ds[rec['variable']];dp=np.array(rec['diameters_nm']);da=v.dimensions[0];assert np.array_equal(ds[da][:],dp)
            raw=np.asarray(np.ma.filled(v[:],np.nan),float).T[mask];hasraw=np.isfinite(raw).any(axis=1)
            q=ds[rec['variable']+'_qc'];qa=np.asarray(np.ma.filled(q[:],0));qa=np.transpose(qa,[q.dimensions.index('time'),q.dimensions.index(da)]+[i for i,d in enumerate(q.dimensions) if d not in ['time',da]])[mask];official_bad=np.isin(qa,invalid_codes).any(axis=2)&~(qa==100).any(axis=2)
            # Published 1h records have decimal-day rounding; Preila 2013 additionally drifts up to 947 s.
            assert np.max(np.abs(np.asarray((rawtime-time).total_seconds())[hasraw]),initial=0)<=1000,(rec['file'],'hour alignment exceeds inspected bound')
            assert not time[hasraw].duplicated().any(),(rec['file'],'finite source rows collapse on hourly alignment')
            groups=subgrids.get((rec['file'],rec['variable']),[dict(indices=list(range(len(dp))),metadata=md)])
            if code=='FR0030R' and len(groups)>1:
                best=max(groups,key=lambda g:int(mdget(g['metadata'],'Version') or 1))
                for g in groups:
                    if g is not best:excluded.append(dict(code=code,file=rec['file'],variable=rec['variable'],reason=f"Superseded within-file channel group, version {mdget(g['metadata'],'Version')}, {len(g['indices'])} bins"))
                groups=[best]
            for gi,g in enumerate(groups):
                idx=np.array(g['indices']);native=dp[idx];r=raw[:,idx];dmeta=g['metadata'];weights=delta(native)
                # VDA 2022 contains two sequential grids: only rows in that grid's observation span belong to it.
                active=np.ones(len(time),bool)
                if code=='FR0027U' and len(groups)>1:
                    anyraw=np.isfinite(r).any(axis=1)
                    assert anyraw.any()
                    first,last=np.flatnonzero(anyraw)[[0,-1]];active[:first]=False;active[last+1:]=False
                rr=r[active];tt=time[active];bad=~np.isfinite(rr)|(rr<0)|np.isin(rr,flags)|official_bad[:,idx][active];clean=np.where(bad,np.nan,rr)
                cc=df[(df.site_code==code)&(df.grid_id==rec['grid'])].sort_values('channel_index').iloc[idx]
                assert np.allclose(cc.diameter_nm,native,rtol=0,atol=1e-10)
                select=(native>=17)&(native<=500);assert select.any()
                dcoef=cc[dfnames].to_numpy();assert np.allclose(dcoef[:,:3].sum(axis=1),dcoef[:,3],atol=1e-12)
                integration_weights=np.ones_like(weights) if (code=='FI0096G' and len(registry)+1==13) else weights
                mat=integration_weights[:,None]*np.pi*(native[:,None]/1000)**2*dcoef
                z=clean[:,select];value=np.nan_to_num(z,nan=0)@mat[select];validbins=np.isfinite(z).sum(axis=1);value[validbins==0]=np.nan
                assert np.array_equal(clean[~bad],rr[~bad]);assert (clean[(rr==0)&~bad]==0).all()
                assert np.allclose(value[:,:3].sum(axis=1),value[:,3],equal_nan=True,atol=1e-9)
                sid=len(registry)+1
                rev=str(mdget(dmeta,'Revision date') or rec['file'].split('.')[2]);version=int(mdget(dmeta,'Version') or 1)
                priority=(int(rev),version,len(time),-sid)
                s=pd.DataFrame(value,columns=regions);s.insert(0,'Time_UTC',tt.to_numpy());s['Raw_start_UTC']=rawtime[active].to_numpy();s['Raw_duration_seconds']=(seconds[mask,1]-seconds[mask,0])[active];s['Source_ID']=sid;s['Valid_bins']=validbins;s['Window_bins']=int(select.sum());s['First_center_nm']=native[select].min();s['Last_center_nm']=native[select].max();s['PNC_cm-3']=np.nansum(z*integration_weights[select],axis=1);s.loc[validbins==0,'PNC_cm-3']=np.nan
                s['_priority']=[priority]*len(s);parts.append(s)
                for n in np.linspace(0,len(s)-1,min(25,len(s)),dtype=int):
                    vals=[sum(float(clean[n,j])*float(mat[j,k]) for j in np.flatnonzero(select&np.isfinite(clean[n]))) if validbins[n] else np.nan for k in range(4)]
                    assert np.allclose(vals,value[n],rtol=1e-12,atol=1e-8,equal_nan=True)
                registry.append(dict(Source_ID=sid,code=code,file=rec['file'],variable=rec['variable'],parent_grid=rec['grid'],subgrid=gi+1,version=version,revision=rev,rows=len(s),native_low_nm=native.min(),native_high_nm=native.max(),native_channels=len(native),first_center_nm=native[select].min(),last_center_nm=native[select].max(),window_channels=int(select.sum()),official_invalid_cells=int((official_bad[:,idx][active]&np.isfinite(rr)).sum()),negative_cells=int((np.isfinite(rr)&(rr<0)).sum()),sentinel_cells=int(np.isin(rr,flags).sum()),nonfinite_cells=int((~np.isfinite(rr)).sum()),zero_retained=int((rr==0).sum()),volume_temperature=mdget(dmeta,'Volume std temperature'),volume_pressure=mdget(dmeta,'Volume std pressure'),instrument=mdget(dmeta,'Instrument name'),method=mdget(dmeta,'Method ref'),url=rec['url'],sha256=manifest[rec['file']]['sha256']))
                for j,d in enumerate(native):coefficients.append(dict(Source_ID=sid,code=code,channel=j+1,diameter_nm=d,delta_log10=weights[j],selected_17_500=bool(select[j]),**dict(zip(dfnames,dcoef[j]))))
    allrows=pd.concat(parts,ignore_index=True);allrows.sort_values(['Time_UTC','_priority'],ascending=[True,False],inplace=True,kind='stable');duplicates=int(allrows.duplicated('Time_UTC').sum());s=allrows.drop_duplicates('Time_UTC',keep='first').drop(columns='_priority').sort_values('Time_UTC').reset_index(drop=True)
    assert not s.Time_UTC.duplicated().any();assert s.Time_UTC.dt.minute.eq(0).all() and s.Time_UTC.dt.second.eq(0).all(),(code,'nonaligned')
    for row in registry:
        if row['code']==code:row['selected_rows']=int(s.Source_ID.eq(row['Source_ID']).sum());row['selected_valid_rows']=int((s.Source_ID.eq(row['Source_ID'])&s.Total.notna()).sum())
    for reg in regions:s.rename(columns={reg:reg+'_LDSA_um2_cm-3'},inplace=True)
    s.to_csv(O/'hourly_ldsa'/f'{code}_LDSA.csv.gz',index=False,float_format='%.12g')
    v=s.dropna(subset=['Total_LDSA_um2_cm-3']);assert len(v)
    core=v[v.Time_UTC.dt.year.between(2015,2019)];fallback=core.empty
    if fallback:
        years=sorted(v.Time_UTC.dt.year.unique());yr=min(years,key=lambda y:(max(2015-y,y-2019,0),-y));core=v[v.Time_UTC.dt.year.eq(yr)];period=str(yr)
    else:period='2015-2019'
    calc={}
    for reg in regions:
        vv=core[reg+'_LDSA_um2_cm-3'];q=vv.quantile([.25,.75,.95]);calc[reg]=dict(Mean=float(vv.mean()),Median=float(vv.median()),P25=float(q.loc[.25]),P75=float(q.loc[.75]),P95=float(q.loc[.95]))
        statrows.append(dict(code=code,name=st['name'],region=reg,period=period,n=len(core),start=str(core.Time_UTC.min()),end=str(core.Time_UTC.max()),**calc[reg]))
    used=sorted(core.Source_ID.unique());regused=[r for r in registry if r['Source_ID'] in used];native_ranges=sorted({(r['native_low_nm'],r['native_high_nm']) for r in regused});actual_ranges=sorted({(r['first_center_nm'],r['last_center_nm']) for r in regused})
    fmt=lambda ranges:'; '.join(f'{a:g}-{b:g}' for a,b in ranges)
    fields=[]
    for reg in regions:
        k=calc[reg];fields += [k['Mean'],k['Median'],f"{k['P25']:.2f}-{k['P75']:.2f}",k['P95']]
    loc=[34+len(location),code,code,st['name'],code[:2],types[st['paper_type']],code,st['paper_type'].split(',')[0],fmt(native_ranges),st['lat'],st['lon']]+fields;assert len(loc)==27;location.append(loc)
    record=dict(code=code,name=st['name'],paper_type=st['paper_type'],period=period,fallback=bool(fallback),statistics_start=str(core.Time_UTC.min()),statistics_end=str(core.Time_UTC.max()),statistics_valid_hours=len(core),statistics_months=int(core.Time_UTC.dt.to_period('M').nunique()),total_mean=calc['Total']['Mean'],full_start=str(v.Time_UTC.min()),full_end=str(v.Time_UTC.max()),full_valid_hours=len(v),full_months=int(v.Time_UTC.dt.to_period('M').nunique()),overlap_rows_removed=duplicates,native_ranges=fmt(native_ranges),actual_selected_ranges=fmt(actual_ranges),partial_hours=int((core.Valid_bins<core.Window_bins).sum()),source_ids=list(map(int,used)),short_native_upper=any(b<450 for a,b in native_ranges))
    overview.append(record);checks.append(dict(code=code,valid_hours=len(v),duplicate_timestamps_after=0,region_sum_verified=True,independent_formula_verified=True));print(code,'period',period,'valid',len(core),'mean',round(calc['Total']['Mean'],3),'overlaps',duplicates,'actual',fmt(actual_ranges),flush=True)
pd.DataFrame(registry).to_csv(O/'calculation_parameters/source_qc_registry.csv',index=False)
pd.DataFrame(excluded).to_csv(O/'calculation_parameters/excluded_input_records.csv',index=False)
pd.DataFrame(coefficients).to_csv(O/'calculation_parameters/diameter_deposition_weights.csv',index=False)
pd.DataFrame(statrows).to_csv(O/'site_spatial_statistics.csv',index=False)
pd.DataFrame(overview).to_csv(O/'netcdf_station_summary.csv',index=False)
(O/'calculation_parameters/integration_validation.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
(ROOT/'work/quality_control/netcdf/workbook_payload.json').write_text(json.dumps(dict(rows=location,overview=overview),ensure_ascii=False,allow_nan=False,indent=2))
shutil.copy2(R/'model_profile.json',O/'calculation_parameters/mppd_model_profile.json')
shutil.copy2(__file__,O/'calculation_parameters/calculate.py')
print('COMPLETE',len(overview),sum(x['full_valid_hours'] for x in overview),flush=True)
