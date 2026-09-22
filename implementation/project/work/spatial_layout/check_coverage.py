"""Filter whole sites by native diameter coverage in the actual plotted period."""
from pathlib import Path
import sys,json,pickle
sys.path.insert(0,'@PROJECT_ROOT@/work/ldsa_support/deps')
import numpy as np
import pandas as pd
ROOT=Path('@PROJECT_ROOT@')
BASE=Path('@DELIVERY_ROOT@')
W=ROOT/'work/spatial_layout'
O=BASE/'Figure1_49_station_coverage';O.mkdir(exist_ok=True)
N=BASE/'netcdf_ldsa_inputs'
periods=pd.read_csv(BASE/'spatial_reference_tables/spatial_analysis_periods.csv')
reg=pd.read_csv(N/'calculation_parameters/source_qc_registry.csv').set_index('Source_ID')
rows=[];configs=[]
for i,r in periods.iterrows():
    site=r.site
    lo,hi=(2015,2019) if r.period=='2015-2019' else (int(r.period),int(r.period))
    grids=[]
    if i<34:
        with open(ROOT/'work/workbook_processing/pnsd'/f'{site}.pkl','rb') as f:x=pickle.load(f)
        s=pd.read_pickle(ROOT/'work/workbook_processing/series'/f'{site}_series.pkl')
        mask=s.Time.dt.year.between(lo,hi)&s.LDSA.notna()
        assert mask.sum()==r.valid_records,(site,mask.sum(),r.valid_records)
    # Structural masks identify each native grid independently of missing observations.
        masks,counts=np.unique(x['structural'][mask],axis=0,return_counts=True)
        for m,n in zip(masks,counts):
            dp=np.asarray(x['diameters'])[~m];assert len(dp)>0
            grids.append((float(dp.min()),float(dp.max()),int(n),'native structural grid'))
    else:
        s=pd.read_csv(N/'hourly_ldsa'/f'{site}_LDSA.csv.gz',usecols=['Time_UTC','Total_LDSA_um2_cm-3','Source_ID'])
        t=pd.to_datetime(s.Time_UTC)
        mask=t.dt.year.between(lo,hi)&s['Total_LDSA_um2_cm-3'].notna()
        assert mask.sum()==r.valid_records,(site,mask.sum(),r.valid_records)
        for sid,n in s.loc[mask].Source_ID.value_counts().items():
            z=reg.loc[sid];assert z.code==site
            grids.append((float(z.native_low_nm),float(z.native_high_nm),int(n),f'Source_ID={sid}'))
    keep=all(a<=17+1e-9 and b>=500-1e-9 for a,b,n,source in grids) and site!='LND_TR'
    reasons=[]
    if site=='LND_TR':reasons.append('Integration lower bound is 31.65 nm; full 17–500 nm coverage is unavailable.')
    if any(a>17+1e-9 for a,b,n,source in grids):reasons.append('Native lower diameter exceeds 17 nm')
    if any(b<500-1e-9 for a,b,n,source in grids):reasons.append('Native upper diameter is below 500 nm')
    ranges='; '.join(f'{a:g}–{b:g}' for a,b in sorted(set((a,b) for a,b,n,source in grids)))
    rows.append(dict(site=site,period=r.period,native_ranges_nm=ranges,keep=keep,reason='; '.join(reasons) if reasons else 'Native configurations cover 17–500 nm'))
    for a,b,n,source in grids:configs.append(dict(site=site,period=r.period,source=source,native_low_nm=a,native_high_nm=b,valid_records=n,configuration_covers_17_500=bool(a<=17+1e-9 and b>=500-1e-9)))
    print(site,keep,ranges,flush=True)
out=pd.DataFrame(rows);out.to_csv(O/'coverage_selection.csv',index=False,encoding='utf-8-sig')
out[~out.keep].to_csv(O/'excluded_stations.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(configs).to_csv(O/'configuration_coverage.csv',index=False,encoding='utf-8-sig')
selection={'keep':out.loc[out.keep,'site'].tolist(),'excluded':out.loc[~out.keep,'site'].tolist(),'rows':rows}
(W/'coverage_selection.json').write_text(json.dumps(selection,ensure_ascii=False,indent=2))
text='For the study protocol and reproduction instructions, see docs/METHODS.md and docs/WORKFLOW.md.\n'
text+=f"Retained {len(selection['keep'])} stations; excluded {len(selection['excluded'])} stations.\n\n| Excluded station | Native diameter range (nm) |\n|---|---|\n"
for r in rows:
    if not r['keep']:text+=f"| {r['site']} | {r['native_ranges_nm']} |\n"
(O/'selectionnotes.md').write_text(text)
print('KEPT',len(selection['keep']),'EXCLUDED',len(selection['excluded']))
