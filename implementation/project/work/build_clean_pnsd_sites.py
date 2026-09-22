from pathlib import Path
import sys, json, pickle, gzip, csv, gc, shutil
from datetime import datetime, timedelta, timezone
from zoneinfo import ZoneInfo
sys.path.insert(0,str(Path(__file__).parent/'pnsd_qc_dependencies'))
import numpy as np
import pandas as pd
import build_site_tables as old

CACHE=Path('work/pnsd_qc_cache')
OUT=Path('outputs/site_data_pnsd_cm3_cleaned')
RENAME={'SCH_SUB':'SCH_RB','WAL_SUB':'WAL_RB','ZUG_SUB':'ZUG_RB'}

def structural_missing(site,df,dps):
    mask=np.zeros((len(df),len(dps)),dtype=bool)
    if site=='MUL_SUB':
        times=pd.to_datetime(df.iloc[:,0],errors='coerce')
        before=(times<pd.Timestamp('2018-01-29 15:30:00')).to_numpy()
        mask[before,:]=dps<14.075
        mask[~before,:]=dps>495.807
    return mask

def main():
    OUT.mkdir(exist_ok=True)
    protocol=json.loads((CACHE/'protocol.json').read_text())
    flags=protocol['sentinel_values']; summary=[]
    sites=sys.argv[1:] or old.SITES
    if sys.argv[1:]:summary=pd.read_csv(OUT/'QC_summary.csv').to_dict('records')
    for site in sites:
        name=RENAME.get(site,site); dest=OUT/name;dest.mkdir(exist_ok=True)
        with (CACHE/f'{site}.pkl').open('rb') as f:parts=pickle.load(f)
        union=sorted({round(dp,6) for p in parts for _,dp in p['dcols']})
        positions={dp:j for j,dp in enumerate(union)}
        headers=old.BASE_HEADERS+[old.format_dp(dp) for dp in union]+old.MODE_HEADERS
        wb,ws,_,_=old.create_workbook(dest/f'{name}_standardized.xlsx',headers)
        nrow=1; interval=timedelta(hours=1)
        counts=dict(site=name,source_site=site,rows=0,bins=len(union),negative=0,sentinel=0,extreme=0,nonfinite_or_missing=0,affected_rows=0)
        with gzip.open(dest/'excluded_PNSD_cells.csv.gz','wt',encoding='utf-8-sig',newline='') as logfile:
            log=csv.writer(logfile);log.writerow(['site','source_file','sheet','excel_row','UTC_start','diameter_nm','original_dN_dlog10Dp','reason','upper_threshold'])
            for part,p in enumerate(parts):
                df=p['df'];dps=np.array([dp for _,dp in p['dcols']]); idx=[i for i,_ in p['dcols']]
                raw=df.iloc[:,idx].apply(pd.to_numeric,errors='coerce').to_numpy(float)
                threshold=np.asarray(p['thresholds'])
                reason=np.zeros(raw.shape,dtype=np.uint8)
                reason[~np.isfinite(raw)]=1;reason[np.isfinite(raw)&(raw<0)]=2
                reason[np.isin(raw,flags)]=3
                reason[(reason==0)&(raw>threshold)]=4
                clean=np.where(reason==0,raw,np.nan)
                bins=clean*old.delta_log10(dps)
                structural=structural_missing(site,df,dps)
                assert not np.isfinite(raw[structural]).any()
                counts['structural_empty_bins']=counts.get('structural_empty_bins',0)+int(structural.sum())
                integration_bins=np.where(structural,0,bins)
                assert np.all(np.isnan(bins)|(bins>=0))
                modes=[]
                for lo,hi in [(10,25),(25,100),(100,800.000001)]:
                    select=(dps>=lo)&(dps<hi)
                    # Ordinary sum deliberately propagates missing/removed measured bins.
                    result=np.sum(integration_bins[:,select],axis=1) if select.any() else np.full(len(df),np.nan)
                    if select.any():result[np.all(structural[:,select],axis=1)]=np.nan
                    modes.append(result)
                modes += [modes[1]+modes[2],modes[0]+modes[1]+modes[2]]
                mapping={}
                for j,c in enumerate(df.columns):
                    k=old.output_header_key(c)
                    if k and k not in mapping:mapping[k]=j
                if p['kind']=='enriched':
                    times=pd.to_datetime(df.iloc[:,mapping['UTC start']],errors='coerce')
                else:
                    tc=next((c for c in df.columns if old.norm_header(c) in {'TIME','DATETIME','DATE'}),df.columns[0])
                    times=pd.to_datetime(df[tc],errors='coerce')
                    diffs=times.dropna().sort_values().diff().dropna();diffs=diffs[(diffs>pd.Timedelta(0))&(diffs<=pd.Timedelta(days=1))]
                    if len(diffs):interval=diffs.median().to_pytimedelta()
                valid_time=times.notna().to_numpy()
                r=reason[valid_time]
                for key,code in [('nonfinite_or_missing',1),('negative',2),('sentinel',3),('extreme',4)]:counts[key]+=int(np.sum(r==code))
                counts['affected_rows']+=int(np.any(r>=2,axis=1).sum())
                for j,arr in enumerate(modes):counts[f'valid_{j+1}']=counts.get(f'valid_{j+1}',0)+int(np.isfinite(arr[valid_time]).sum())
                pos=[positions[round(dp,6)] for dp in dps]
                for i,row in enumerate(df.itertuples(index=False,name=None)):
                    if not valid_time[i]:continue
                    start=times.iloc[i]
                    if p['kind']=='enriched':
                        base=[row[mapping[h]] if h in mapping else None for h in old.BASE_HEADERS]
                        base[2]=start
                        if site=='ROC_UB':
                            for t in (2,3):
                                if isinstance(base[t],datetime):base[t]=base[t].replace(tzinfo=ZoneInfo('America/New_York')).astimezone(timezone.utc).replace(tzinfo=None)
                    else:
                        base=[None]*20;base[2]=start;base[3]=start+interval-timedelta(seconds=1)
                        for h,j in [('TEMP (°C)',13),('PRESS (hPa)',15)]:
                            if h in mapping:
                                val=pd.to_numeric(row[mapping[h]],errors='coerce')
                                base[j]=val-273.15 if j==13 and pd.notna(val) and val>150 else val
                    base[0]=name;base[1]='Regional background' if name.endswith('_RB') else old.TYPE_MAP[name.rsplit('_',1)[1]]
                    b=[None]*len(union)
                    for j,target in enumerate(pos):b[target]=bins[i,j]
                    old.write_row(ws,nrow,base+b+[x[i] for x in modes]);nrow+=1
                    for j in np.flatnonzero(reason[i]>=2):
                        log.writerow([name,p['path'],p['sheet'],i+2,base[2],dps[j],raw[i,j],{2:'negative',3:'sentinel',4:'statistical_extreme'}[int(reason[i,j])],threshold[j]])
                del raw,clean,bins,reason,r;gc.collect()
        ws.autofilter(0,0,nrow-1,len(headers)-1);wb.close()
        counts['rows']=nrow-1;summary=[x for x in summary if x['site']!=name];summary.append(counts)
        pd.DataFrame(summary).to_csv(OUT/'QC_summary.csv',index=False,encoding='utf-8-sig')
        print(counts,flush=True)
        del parts,df;gc.collect()
    shutil.copy2(CACHE/'channel_thresholds.csv',OUT/'channel_thresholds_source_IDs.csv')
    shutil.copy2(CACHE/'protocol.json',OUT/'QC_protocol.json')
    (OUT/'README.txt').write_text('For the study protocol and reproduction instructions, see docs/METHODS.md and docs/WORKFLOW.md.\n',encoding='utf-8')

if __name__=='__main__':main()
