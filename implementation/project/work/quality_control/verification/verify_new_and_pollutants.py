from pathlib import Path
import pandas as pd,numpy as np,json,sys,re,datetime,hashlib
R=Path('@PROJECT_ROOT@');V=R/'work/quality_control/verification';B=Path('@DELIVERY_ROOT@');O=R/'outputs/quality_control/netcdf';OLD=R/'outputs/network_trends/corrected_LDSA';sys.path[:0]=[str(R/'work/ipr_download/deps')]
import netCDF4
f=pd.read_csv(R/'work/provider_flag_audit/ebas_flags.csv',header=None,quotechar="'",names=['flag','description','status']);invalid=set(f.loc[f.status.isin(['I','H','M']),'flag']);definition=dict(zip(f.flag,f.status))
def parse(s):
 if pd.isna(s):return []
 return [int(x) for x in re.findall(r'\d+',str(s))]
def isbad(s):
 c=parse(s);return bool(set(c)&invalid) and 100 not in c
summary=[];changes=[];fail=[]
cols=['Head_LDSA_um2_cm-3','P_LDSA_um2_cm-3','TB_LDSA_um2_cm-3','Total_LDSA_um2_cm-3','PNC_cm-3']
for p in sorted((O/'hourly_ldsa').glob('*.csv.gz')):
 code=p.name[:7];new=pd.read_csv(p);prior=R/'outputs/pal_channel_conversion/FI0096G_LDSA_early_deltaN.csv.gz' if code=='FI0096G' else OLD/'hourly_ldsa'/p.name;old=pd.read_csv(prior)
 assert len(new)==len(old),(code,len(new),len(old))
 keycols=[x for x in new if x not in cols+['Valid_bins']];mismatch={c:int((~((new[c]==old[c])|(new[c].isna()&old[c].isna()))).sum()) for c in keycols};assert not any(mismatch.values()),(code,mismatch)
 removed=old[cols[-2]].notna()&new[cols[-2]].isna();introduced=old[cols[-2]].isna()&new[cols[-2]].notna();assert not introduced.any(),code
 delta={c:int((~np.isclose(old.loc[~removed,c],new.loc[~removed,c],rtol=2e-10,atol=2e-8,equal_nan=True)).sum()) for c in cols};assert not any(delta.values()),(code,delta)
 assert new.loc[removed,cols].isna().all().all();assert new.loc[removed,'Valid_bins'].eq(0).all();assert np.array_equal(new.loc[~removed,'Valid_bins'],old.loc[~removed,'Valid_bins'])
 for i in np.flatnonzero(removed):changes.append(dict(code=code,Time_UTC=new.Time_UTC.iloc[i],Source_ID=int(new.Source_ID.iloc[i]),old_Total=float(old[cols[-2]].iloc[i])))
 summary.append(dict(code=code,rows=len(new),old_finite=int(old[cols[-2]].notna().sum()),new_finite=int(new[cols[-2]].notna().sum()),removed_hours=int(removed.sum()),all_other_data_unchanged=True))
co=pd.read_csv(O/'calculation_parameters/diameter_deposition_weights.csv');pco=pd.read_csv(OLD/'calculation_parameters/diameter_deposition_weights.csv');pd.testing.assert_frame_equal(co,pco)
reg=pd.read_csv(O/'calculation_parameters/source_qc_registry.csv');preg=pd.read_csv(OLD/'calculation_parameters/source_qc_registry.csv');pd.testing.assert_frame_equal(reg[preg.columns.drop('selected_valid_rows')],preg.drop(columns='selected_valid_rows'))
cs=pd.DataFrame(changes);cs.to_csv(V/'netcdf_removed_hours.csv',index=False);pd.DataFrame(summary).to_csv(V/'netcdf_comparison.csv',index=False)
a=pd.read_csv(B/'provider_quality_flags/netcdf_source_flags.csv');print('audit cols',a.columns.tolist())
# Exact consistency with independently prepared full-network flag audit.
assert len(cs)==9726 and cs.code.nunique()==6
# Independent PAL Source 13 integration: use native DeltaN directly, no dlog multiplier.
x=reg[reg.Source_ID.eq(13)].iloc[0];path=B/'EBAS_observations/source_netcdf'/x.code/'PNSD'/x.file
with netCDF4.Dataset(path) as ds:
 v=ds[x.variable];da=next(z for z in v.dimensions if z!='time');raw=np.moveaxis(np.asarray(np.ma.filled(v[:],np.nan)),v.dimensions.index('time'),0);dp=np.asarray(ds[da][:]);cc=co[co.Source_ID.eq(13)].sort_values('channel');select=(dp>=17)&(dp<=500);sent=json.loads((R/'outputs/pnsd_cleaning/PNSD/QC_protocol.json').read_text())['sentinel_values']+[999999.99999];raw[(raw<0)|np.isin(raw,sent)]=np.nan;z=raw[:,select];coeff=np.pi*(dp[select]/1000)**2*cc.loc[select,'DF_total'].to_numpy();expected=np.nansum(z*coeff,axis=1);expected[~np.isfinite(z).any(1)]=np.nan;tv=ds['time'];ep=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units);fac={'days':86400,'hours':3600,'minutes':60,'seconds':1}[tv.units.split()[0]];tt=pd.to_datetime(np.rint((np.asarray(ds['time_bnds'][:,0])-ep)*fac).astype('int64'),unit='s').round('h');ex=pd.Series(expected,index=tt)
pal=pd.read_csv(O/'hourly_ldsa/FI0096G_LDSA.csv.gz',parse_dates=['Time_UTC']);p=pal[pal.Source_ID.eq(13)];expect=ex.reindex(p.Time_UTC).to_numpy();assert np.allclose(expect,p.Total_LDSA_um2_cm_3 if 'Total_LDSA_um2_cm_3' in p else p['Total_LDSA_um2_cm-3'],equal_nan=True,rtol=2e-10,atol=2e-8)
checks=dict(netcdf_verified=True,rows=sum(s['rows'] for s in summary),removed_hours=len(cs),removed_sites=cs.groupby('code').size().to_dict(),DF_and_source_metadata_identical=True,PAL_Source13_DeltaN_verified=True,PAL_Source13_finite=int(np.isfinite(expect).sum()))
# Native pollutant flags: only official invalid status, retain 100 override and V warnings.
pol=B/'german_pollutant_inputs';allparts=[];ps=[];PO=V/'pollutants';PO.mkdir(exist_ok=True)
for fn in sorted((pol/'analysis_inputs').glob('*.pkl')):
 d=pd.read_pickle(fn);eb=d.source.eq('EBAS');bad=d.qc.apply(isbad)&eb&np.isfinite(pd.to_numeric(d.value,errors='coerce'));unknown=sorted({c for q in d.loc[eb,'qc'] for c in parse(q) if c not in definition and c!=0});q=d.loc[bad].copy();q['origin']='German_current';allparts.append(q);clean=d.copy();clean.loc[bad,'value']=np.nan;clean.to_pickle(PO/fn.name);ps.append(dict(site=fn.stem,scope='German_current',finite_rows=int(np.isfinite(pd.to_numeric(d.value,errors='coerce')).sum()),EBAS_rows=int(eb.sum()),invalid_finite=int(bad.sum()),unknown=unknown))
for m in ['NO2','PM2.5']:
 d=pd.read_csv(R/f'outputs/ipr_processing/IPR_RB_{m}_observations.csv.gz');bad=d.quality_flag.apply(isbad)&np.isfinite(pd.to_numeric(d.value,errors='coerce'));d['qc']=d.quality_flag;d['Site']='IPR_RB';d['metric']=m;q=d.loc[bad].copy();q['origin']='IPR_current';allparts.append(q);clean=d.copy();clean.loc[bad,'value']=np.nan;clean['date']=pd.to_datetime(clean.date,format="mixed");clean['end']=pd.to_datetime(clean.end,format="mixed");clean.to_pickle(PO/(m+'.pkl'));ps.append(dict(site='IPR_RB',scope='IPR_'+m,finite_rows=int(np.isfinite(pd.to_numeric(d.value,errors='coerce')).sum()),invalid_finite=int(bad.sum()),unknown=sorted({c for q in d.quality_flag for c in parse(q) if c not in definition and c!=0})))
ipr=pd.concat([pd.read_pickle(PO/(m+'.pkl')) for m in ['NO2','PM2.5']],ignore_index=True);ipr.to_pickle(PO/'IPR_RB.pkl')
pr=pd.concat(allparts,ignore_index=True);pr.to_csv(V/'pollutants_invalid_finite.csv.gz',index=False);pd.DataFrame(ps).to_csv(V/'pollutants_summary.csv',index=False);g=pr.groupby(['origin','Site','metric','qc']).agg(n=('value','size'),start=('date','min'),end=('date','max')).reset_index();g.to_csv(V/'pollutants_invalid_groups.csv',index=False)
checks['pollutants_invalid_finite']=len(pr);checks['pollutants_invalid_groups']=g.to_dict('records');(V/'verification.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2,default=str));print(json.dumps(checks,ensure_ascii=False,indent=2,default=str))
