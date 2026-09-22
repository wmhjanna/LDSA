"""Recompute LDSA exclusively from exported cleaned bin number concentrations."""
from pathlib import Path
import sys, re, json, pickle, hashlib, gc
import numpy as np
import pandas as pd
ROOT=Path('@PROJECT_ROOT@')
sys.path.insert(0,str(ROOT/'work/ldsa_support/deps'))
BASE=Path('@INPUT_ROOT@/data_sources')
OUT=ROOT/'outputs/ldsa_cleaned_reanalysis'
STAGE=ROOT/'work/ldsa_support'
REV={'SCH_RB':'SCH_SUB','WAL_RB':'WAL_SUB','ZUG_RB':'ZUG_SUB'}
REGIONS=['HA','TB','ALV','Total']
RANGES=[(17,500),(10,500),(10,800)]

def coefficients(old,dps):
    if old=='LAN_UB':paths=[BASE/'LAN/deposition_fractions.xlsx']
    elif old=='PRA_SUB':paths=sorted((BASE/'source_data').glob('PRA_SUB_*/DF*.xlsx'))
    else:paths=list((BASE/'source_data'/old).glob('DF*.xlsx'))
    pool=[]
    for p in paths:
        for sh,d in pd.read_excel(p,sheet_name=None).items():
            d=d.iloc[:,:5].apply(pd.to_numeric,errors='coerce').dropna()
            if old=='LEI_TR':
                ref=pd.read_excel(BASE/'source_data/ANB_UB/DF10-800.xlsx')
                assert np.array_equal(d.iloc[:,1:].values,ref.iloc[:,1:].values)
                d.iloc[:,0]=ref.iloc[:,0].values
            if old=='HEL_TR':
                # Retain historical channel correspondence, correct surface-area diameters.
                raw=pd.read_excel(BASE/'source_data/HEL_TR/HEL_TR.xlsx',nrows=0)
                assert len(d)==len(raw.columns)-1
                d.iloc[:,0]=np.array(raw.columns[1:],float)/1000
            for row in d.to_numpy(float):pool.append((row[0]*1000,row[1:],str(p),str(sh)))
    result=[]; audit=[]
    for dp in dps:
        candidates=sorted(pool,key=lambda x:abs(dp-x[0]));r=candidates[0]
        assert abs(dp-r[0])<=max(.06,dp*.005),(old,dp,r[0])
        result.append(r[1]);audit.append(dict(diameter_nm=dp,DF_diameter_nm=r[0],HA=r[1][0],TB=r[1][1],ALV=r[1][2],Total=r[1][3],source=r[2],sheet=r[3]))
    return np.array(result),audit

def main():
    for name in ['sites','analysis','figures']: (OUT/name).mkdir(parents=True,exist_ok=True)
    chosen=sys.argv[1:]
    for path in sorted((ROOT/'outputs/site_data_pnsd_cm3_cleaned').glob('*/*.xlsx')):
        site=path.parent.name
        if chosen and site not in chosen:continue
        dest=OUT/'sites'/site;dest.mkdir(exist_ok=True)
        if (dest/'complete.json').exists():continue
        old=REV.get(site,site)
        df=pd.read_excel(path,engine='calamine')
        times=pd.to_datetime(df['UTC start'],errors='coerce')
        cols=[c for c in df if str(c).startswith('D_')]
        dps=np.array([float(re.search(r'D_(.*?) nm',c).group(1)) for c in cols])
        bins=df[cols].to_numpy(float)
        assert not np.isinf(bins).any() and np.all(np.isnan(bins)|(bins>=0))
        structural=np.ones_like(bins,dtype=bool)
        with (ROOT/f'work/pnsd_qc_cache/{old}.pkl').open('rb') as f:parts=pickle.load(f)
        offset=0
        precise={}
        for p in parts:
            # PNSD export preserves valid-time row order within each input part.
            tc=next((c for c in p['df'] if 'start' in str(c).lower()),p['df'].columns[0])
            pt=pd.to_datetime(p['df'][tc],errors='coerce').dropna()
            n=len(pt);pdps=np.array([x[1] for x in p['dcols']])
            for dp in pdps:
                j=int(np.argmin(abs(dps-dp)));assert abs(dps[j]-dp)<.00051
                precise[j]=dp;structural[offset:offset+n,j]=False
            offset+=n
        assert offset==len(df),(site,offset,len(df))
        dps=np.array([precise[j] for j in range(len(cols))])
        if old=='MUL_SUB':
            before=(times<pd.Timestamp('2018-01-29 15:30')).to_numpy()
            structural[before,:]|=dps<14.075
            structural[~before,:]|=dps>495.807
        assert np.isnan(bins[structural]).all()
        cf,audit=coefficients(old,dps)
        pd.DataFrame(audit).to_csv(dest/f'{site}_deposition_coefficients.csv',index=False)
        meta=df.iloc[:,:20].copy()
        series=pd.DataFrame({'Time':times})
        for reg,k in zip(REGIONS,range(4)):
            values=bins*(np.pi*(dps/1000)**2*cf[:,k])
            detail=pd.DataFrame(values,columns=[f'LDSA_{dp:.6g}nm_um2_cm3' for dp in dps])
            detail.insert(0,'Time',times)
            detail.to_csv(dest/f'{site}_{reg}_bins.csv.gz',index=False)
            integration=np.where(structural,0,values)
            for lo,hi in RANGES:
                for label,a,b,inclusive in [('LDSA',lo,hi,True),('Nucleation',lo,25,False),('Aitken',25,100,False),('Accumulation',100,hi,True)]:
                    mask=(dps>=a)&((dps<=b) if inclusive else (dps<b))
                    v=integration[:,mask].sum(axis=1) if mask.any() else np.full(len(df),np.nan)
                    if mask.any():v[structural[:,mask].all(axis=1)]=np.nan
                    key=f'{reg}_{label}_{a}_{b}nm'
                    series[key]=v
            del detail,values,integration
        # Aliases for the primary response and number metrics are explicit.
        series['LDSA']=series['Total_LDSA_17_500nm']
        for label,lo,hi,inclusive in [('Nucleation',17,25,False),('Aitken',25,100,False),('Accumulation',100,500,True),('Ntotal',17,500,True)]:
            select=(dps>=lo)&((dps<=hi) if inclusive else (dps<hi))
            value=np.nansum(bins[:,select],axis=1)
            reference=f'Total_{"LDSA" if label=="Ntotal" else label}_{lo}_{hi}nm'
            value[series[reference].isna()]=np.nan
            series[f'PNC_{label}_{lo}_{hi}nm']=value
        for c in df.columns[4:20]:series[c]=pd.to_numeric(df[c],errors='coerce')
        for c in df.columns[-5:]:series[c]=pd.to_numeric(df[c],errors='coerce')
        series.to_pickle(STAGE/f'{site}_series.pkl')
        full=pd.concat([meta,series.iloc[:,1:37]],axis=1)
        # Select every distinct derived region/range/mode column (duplicate modes appear only once).
        derived=[c for c in series if any(c.startswith(r+'_') for r in REGIONS)]
        full=pd.concat([meta,series[derived]],axis=1)
        full.to_csv(dest/f'{site}_LDSA_standardized.csv.gz',index=False)
        summary={'site':site,'rows':len(df),'bins':len(dps),'valid_primary':int(series.LDSA.notna().sum()),'source':str(path),'source_sha256':hashlib.file_digest(path.open('rb'),'sha256').hexdigest(),'columns':list(full.columns)}
        (dest/'complete.json').write_text(json.dumps(summary,indent=2))
        print(summary['site'],summary['rows'],summary['valid_primary'],flush=True)
        del df,bins,series,full,meta,parts;gc.collect()

if __name__=='__main__':main()
