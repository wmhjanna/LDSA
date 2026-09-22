from pathlib import Path
import os,sys,json
import numpy as np
import pandas as pd
ROOT=Path('@PROJECT_ROOT@');STAGE=ROOT/'work/quality_control/workbook';sys.path.insert(0,str(ROOT/'work/ldsa_support/deps'));sys.path.insert(0,str(ROOT/'work'))
import lightgbm as lgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score,mean_squared_error,mean_absolute_error
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'outputs/quality_control/workbook';ANA=OUT/'analysis';FIG=OUT/'figures'
FEATURES={'BC (ng/m³)':'BC','CO (mg/m³)':'CO','NO (µg/m³)':'NO','NO2 (µg/m³)':'NO2','O3 (µg/m³)':'O3','PM2.5 (µg/m³)':'PM2.5','SO2 (µg/m³)':'SO2','TEMP (°C)':'TEMP','RH (%)':'RH','PRESS (hPa)':'PRESS','RAD (W/m²)':'RAD','WS (m/s)':'WS','WD (°)':'WD','RAIN (mm)':'RAIN'}
LABELS={v:k for k,v in FEATURES.items()}
plt.rcParams.update({'font.family':['Times New Roman','DejaVu Serif'],'font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
from ldsa_reference_style import configure,apply_style
configure()
def save(fig,name):
    apply_style(fig,name)
    for ext in ['png','pdf','svg']:fig.savefig(FIG/f'{name}.{ext}',dpi=300,bbox_inches='tight')
    plt.close(fig)
def main():
    w=pd.read_pickle(STAGE/'weekly_model_reference.pkl').rename(columns=FEATURES);bps=pd.read_csv(ANA/'breakpoints.csv',parse_dates=['Breakpoint']);audit=[];eligible=[]
    w.loc[w.Site.eq('HEL_TR'),'SO2']=np.nan
    for site,g in w[w.Type.isin(['UB','TR'])].groupby('Site',sort=False):
        b=bps[bps.Site.eq(site)&bps.Metric.eq('LDSA')]
        if b.empty:continue
        g=g.dropna(subset=['LDSA']);g=g[g[list(FEATURES.values())].notna().any(axis=1)].copy()
        ranges=pd.read_csv(STAGE/'native_phase_ranges.csv',parse_dates=['start','end'])
        parts=[]
        for phase in ['Before','After']:
            rr=ranges[ranges.Site.eq(site)&ranges.Metric.eq('LDSA')&ranges.Phase.eq(phase)]
            if rr.empty:continue
            start=pd.Timestamp(rr.start.iloc[0]).to_period('M').start_time.to_period('W').start_time
            end=pd.Timestamp(rr.end.iloc[0]).to_period('M').start_time.to_period('W').start_time
            z=g[g.date.between(start,end)].copy();z['Phase']=phase;parts.append(z)
        if len(parts)!=2:continue
        g=pd.concat(parts);counts=g.Phase.value_counts();keep=counts.get('Before',0)>=52 and counts.get('After',0)>=52
        audit.append(dict(Site=site,before_weeks=counts.get('Before',0),after_weeks=counts.get('After',0),included=keep))
        if keep:eligible.append(g)
    pd.DataFrame(audit).to_csv(ANA/'model_site_eligibility.csv',index=False)
    if not eligible:raise RuntimeError('No sites meet both 52-week phase criteria')
    data=pd.concat(eligible);data[['Site','date','Type','Phase','LDSA']+list(FEATURES.values())].to_csv(ANA/'weekly_model_input.csv',index=False)
    metrics=[];importance=[];models={}
    for typ in ['UB','TR']:
        for phase in ['Before','After']:
            d=data[data.Type.eq(typ)&data.Phase.eq(phase)];X=d[list(FEATURES.values())];y=d.LDSA
            if len(d)<20:continue
            train,test=train_test_split(np.arange(len(d)),test_size=.2,random_state=42);train,valid=train_test_split(train,test_size=.25,random_state=42)
            model=lgb.train({'objective':'regression','metric':'rmse','learning_rate':.05,'num_leaves':31,'max_depth':6,'feature_fraction':.8,'seed':42,'verbosity':-1,'num_threads':4,'deterministic':True,'force_col_wise':True},lgb.Dataset(X.iloc[train],label=y.iloc[train]),num_boost_round=1000,valid_sets=[lgb.Dataset(X.iloc[valid],label=y.iloc[valid])],callbacks=[lgb.early_stopping(50,verbose=False)])
            pred=model.predict(X.iloc[test]);contrib=model.predict(X.iloc[test],pred_contrib=True);shap=contrib[:,:-1]
            assert np.allclose(contrib.sum(axis=1),pred,rtol=1e-6,atol=1e-6)
            model.save_model(str(ANA/f'LightGBM_{typ}_{phase}.txt'))
            m=dict(Type=typ,Phase=phase,sites=d.Site.nunique(),weeks=len(d),train=len(train),validation=len(valid),test=len(test),R2=r2_score(y.iloc[test],pred),RMSE=np.sqrt(mean_squared_error(y.iloc[test],pred)),MAE=mean_absolute_error(y.iloc[test],pred),trees=model.best_iteration);metrics.append(m);print(m,flush=True)
            for f,val in zip(X.columns,np.abs(shap).mean(axis=0)):importance.append(dict(Type=typ,Phase=phase,Feature=f,mean_abs_SHAP=val))
            held=d.iloc[test][['Site','date','LDSA']].copy();held['Prediction']=pred
            for j,f in enumerate(X.columns):held[f]=X.iloc[test][f].values;held[f+'_SHAP']=shap[:,j]
            held.to_csv(ANA/f'heldout_SHAP_{typ}_{phase}.csv',index=False);models[typ,phase]=held
    pd.DataFrame(metrics).to_csv(ANA/'model_performance.csv',index=False);imp=pd.DataFrame(importance);imp.to_csv(ANA/'SHAP_importance.csv',index=False)

if __name__=='__main__':main()
