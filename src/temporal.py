"""Temporal preprocessing used for the study's station diagnostics."""
import numpy as np
import pandas as pd

def paired_monthly(time,values):
    """Average daily means equally; require at least 14 observed days per month."""
    d=pd.DataFrame({'date':pd.to_datetime(time,errors='coerce'),'value':pd.to_numeric(values,errors='coerce')}).replace([np.inf,-np.inf],np.nan).dropna()
    day=d.groupby(d.date.dt.floor('D')).value.mean()
    m=day.resample('MS').agg(['mean','count']);m.loc[m['count']<14,'mean']=np.nan
    return m['mean'].dropna().rename('value').rename_axis('date').reset_index()

def observation_breakpoint(time,values,valid_months):
    """One linear-kernel change on ordered observations; total eligibility >=20 months."""
    if valid_months<20:return None
    d=pd.DataFrame({'date':time,'value':values}).replace([np.inf,-np.inf],np.nan).dropna().sort_values('date')
    y=d.value.to_numpy(float)
    if len(y)<4 or np.std(y)==0:return None
    y=y-y.mean();s=np.cumsum(y);k=np.arange(2,len(y)-1)
    cost=s[k-1]**2/k+(s[-1]-s[k-1])**2/(len(y)-k)
    return pd.Timestamp(d.date.iloc[int(k[np.argmax(cost)])])
