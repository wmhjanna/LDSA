exec(open('@PROJECT_ROOT@/work/provider_flag_audit/audit.py').read().split('results=[];')[0])
table=pd.read_csv(W/'ebas_flags.csv',header=None,quotechar="'",names=['flag','description','status']);badflags=table[table.status.isin(['I','H','M'])].flag.to_numpy()
data=pd.read_csv(P/'hourly_ldsa/MT0001R_LDSA.csv.gz',parse_dates=['Time_UTC']);data['corrected']=data['Total_LDSA_um2_cm-3'];audit=[]
for sid,g in data.groupby('Source_ID'):
 rr=reg[reg.Source_ID.eq(sid)].iloc[0];c=cf[cf.Source_ID.eq(sid)].sort_values('channel');dp=c.diameter_nm.to_numpy();sel=(dp>=17)&(dp<=500);weight=c.delta_log10.to_numpy()*np.pi*(dp/1000)**2*c.DF_total.to_numpy()
 with netCDF4.Dataset(D/'MT0001R/PNSD'/rr.file) as n:
  tb=np.asarray(n['time_bnds'][:,0]);tv=n['time'];ep=netCDF4.date2num(datetime.datetime(1970,1,1),tv.units);times=pd.to_datetime(np.rint((tb-ep)*86400).astype('int64'),unit='s').round('h');idx=pd.Index(times).get_indexer(g.Time_UTC);assert (idx>=0).all()
  raw=np.asarray(np.ma.filled(n[rr.variable][:,idx],np.nan),float).T;z=np.where(~np.isfinite(raw)|(raw<0)|np.isin(raw,flags),np.nan,raw)
  qc=np.asarray(n[rr.variable+'_qc'][:,:,idx].filled(0));bad=np.isin(qc,badflags).any(axis=1).T;override=(qc==100).any(axis=1).T;bad&=~override
  finite=np.isfinite(z)&sel;affected=(bad&finite).any(axis=1);audit.append(dict(Source_ID=int(sid),finite_cells_invalidated=int((bad&finite).sum()),affected_hours=int(affected.sum()),flags=sorted(set(map(int,qc[np.isin(qc,badflags)])))))
  z[bad]=np.nan;calc=np.nansum(z[:,sel]*weight[sel],axis=1);calc[~np.isfinite(z[:,sel]).any(axis=1)]=np.nan;data.loc[g.index,'corrected']=calc
data['week']=data.Time_UTC.dt.to_period('W').dt.start_time;out=[]
for week in cases['GIO_RB'][1]:
 g=data[data.week.eq(pd.Timestamp(week))];out.append(dict(week=week,before=g['Total_LDSA_um2_cm-3'].mean(),after=g.corrected.mean(),before_hours=g['Total_LDSA_um2_cm-3'].count(),after_hours=g.corrected.count()))
d=data[['Time_UTC','corrected']].rename(columns={'Time_UTC':'date','corrected':'value'}).dropna().reset_index(drop=True);m=monthly(d.date,d.value);k=index_cp(d,len(m));result=dict(source_audit=audit,flag_table_source='https://raw.githubusercontent.com/metno/pyaerocom/master/pyaerocom/data/ebas_flags.csv',new_breakpoint=str(d.date.iloc[k]),new_valid_hours=len(d),original_valid_hours=int(data['Total_LDSA_um2_cm-3'].count()),weekly_comparison=out,scope='Diagnostic copy only; no active data or figures changed')
(O/'quality_control_sensitivity.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,default=str));print(json.dumps(result,ensure_ascii=False,indent=2,default=str))
data.to_csv(O/'quality_control_sensitivity.csv.gz',index=False)
