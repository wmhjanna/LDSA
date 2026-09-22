exec(open('work/quality_control/verification/verify_old_masks.py').read().split('for site,g in a.groupby')[0])
site='LND_UB';s=pickle.load(open(R/'work/workbook_processing/pnsd'/f'{site}.pkl','rb'));dp=s['diameters'];bins=s['bins'];ts=pd.DatetimeIndex(s['aux'].Time);width=np.r_[np.diff(np.log10(dp)),np.diff(np.log10(dp))[-1]];orig=bins/width;r=a[a.Site.eq(site)].iloc[0]
with netCDF4.Dataset(D/r.code/'PNSD'/r.file) as ds:
 v=ds[r.variable];da=next(x for x in v.dimensions if x!='time');ndp=np.asarray(ds[da][:]);q=ds[r.variable+'_qc'];qa=np.asarray(np.ma.filled(q[:],0));qa=np.transpose(qa,[q.dimensions.index('time'),q.dimensions.index(da)]+[i for i,x in enumerate(q.dimensions) if x not in ['time',da]]);bad=np.isin(qa,invalid).any(2)&~(qa==100).any(2);raw=np.moveaxis(np.asarray(np.ma.filled(v[:],np.nan)),v.dimensions.index('time'),0);tv=ds['time'];ep=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units);fac={'days':86400,'hours':3600,'seconds':1,'minutes':60}[tv.units.split()[0]];tt=pd.to_datetime(np.rint((np.asarray(ds['time_bnds'][:,0])-ep)*fac).astype('int64'),unit='s').round('h')
ix=np.abs(dp[:,None]-ndp).argmin(1);matchdp=np.isclose(dp,ndp[ix],rtol=0,atol=.001);mi=pd.Index(tt).get_indexer(ts.round('min').floor('h'));idx=np.flatnonzero(mi>=0);rr=raw[mi[idx]][:,ix];xx=orig[idx];finite=np.isfinite(xx)&np.isfinite(rr)&matchdp;flagged=bad[mi[idx]][:,ix]&matchdp&np.isfinite(rr);bc=flagged&finite
comp=[]
for fact in [1,32]:
 for atol in [.001,.016,.16,1.61]:
  ag=finite&np.isclose(xx,rr*fact,rtol=0,atol=atol);n=finite.sum(1);rat=ag.sum(1)/np.maximum(n,1);exact=(n>=30)&(rat==1);selected=bc&exact[:,None];comp.append(dict(factor=fact,atol=atol,matching_spectra=int(exact.sum()),new_bad_rows=int(selected.any(1).sum()),new_bad_cells=int(selected.sum())))
  if fact==32 and atol==1.61:
   out=pd.DataFrame(dict(row=idx,Time=ts[idx],finite_comparable=n,ratio=rat,all_agree=exact,flagged_finite_cells=bc.sum(1),max_abs_diff=np.nanmax(np.where(finite,np.abs(xx-rr*fact),np.nan),axis=1)));out.to_csv(V/'LND_factor32_diagnostic.csv',index=False)
print(pd.DataFrame(comp).to_string(index=False));(V/'LND_factor32_diagnostic_summary.json').write_text(json.dumps(comp,indent=2))
