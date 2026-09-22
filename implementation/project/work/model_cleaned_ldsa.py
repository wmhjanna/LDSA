from pathlib import Path
import os,sys,json
import numpy as np
import pandas as pd
ROOT=Path('@PROJECT_ROOT@');STAGE=ROOT/'work/ldsa_support';sys.path.insert(0,str(STAGE/'deps'))
import lightgbm as lgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score,mean_squared_error,mean_absolute_error
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'outputs/ldsa_cleaned_reanalysis';ANA=OUT/'analysis';FIG=OUT/'figures'
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
    w=pd.read_pickle(STAGE/'weekly.pkl').rename(columns=FEATURES);bps=pd.read_csv(ANA/'breakpoints.csv',parse_dates=['Breakpoint']);audit=[];eligible=[]
    for site,g in w[w.Type.isin(['UB','TR'])].groupby('Site'):
        b=bps[bps.Site.eq(site)&bps.Metric.eq('LDSA')]
        if b.empty:continue
        g=g.dropna(subset=['LDSA']);g=g[g[list(FEATURES.values())].notna().any(axis=1)].copy()
        g['Phase']=np.where(g.date<b.Breakpoint.iloc[0],'Before','After');counts=g.Phase.value_counts();keep=counts.get('Before',0)>=52 and counts.get('After',0)>=52
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
    fig,axs=plt.subplots(2,2,figsize=(12,10),layout='constrained')
    for i,typ in enumerate(['UB','TR']):
        for j,phase in enumerate(['Before','After']):
            d=imp[imp.Type.eq(typ)&imp.Phase.eq(phase)].sort_values('mean_abs_SHAP');axs[i,j].barh(d.Feature,d.mean_abs_SHAP,color='#1f77b4' if phase=='Before' else '#ff7f0e');axs[i,j].set_title(f'{typ} {phase}',loc='left');axs[i,j].set_xlabel('Mean |SHAP| (μm² cm⁻³)')
    save(fig,'Figure_S10')
    # Predictor-response displays use the held-out observations and polynomial confidence bands.
    def dependence(ax,feature,typ):
        for phase,color in [('Before','#1f77b4' if typ=='UB' else '#2ca02c'),('After','#ff7f0e' if typ=='UB' else '#d62728')]:
            if (typ,phase) not in models:continue
            d=models[typ,phase][[feature,feature+'_SHAP']].dropna();x=d[feature].to_numpy();y=d[feature+'_SHAP'].to_numpy();ax.scatter(x,y,c=color,s=50,alpha=.3,edgecolors='white',label=phase)
            if len(x)>10 and np.unique(x).size>4:
                center=x.mean();scale=x.std();xx=np.linspace(np.quantile(x,.01),np.quantile(x,.99),100);z=(x-center)/scale;zz=(xx-center)/scale;A=np.vander(z,4);B=np.vander(zz,4);coef=np.linalg.lstsq(A,y,rcond=None)[0];v=np.sum((y-A@coef)**2)/(len(x)-4);cov=v*np.linalg.pinv(A.T@A);se=np.sqrt(np.maximum(0,np.sum((B@cov)*B,axis=1)));fit=B@coef;ax.plot(xx,fit,c=color,lw=3.0);ax.fill_between(xx,fit-1.96*se,fit+1.96*se,color=color,alpha=.12)
        ax.set_title(typ+' · '+feature,loc='left');ax.set_xlabel(LABELS[feature]);ax.set_ylabel('SHAP (μm² cm⁻³)');ax.grid(alpha=.12);ax.legend(frameon=False,fontsize=7,loc='best')
    for typ in ['UB','TR']:
        fig,axs=plt.subplots(4,4,figsize=(22,20),layout='constrained')
        for ax,f in zip(axs.flat,FEATURES.values()):dependence(ax,f,typ)
        for ax in list(axs.flat)[14:]:ax.axis('off')
        save(fig,f'Figure_S11_{typ}')
    fig=plt.figure(figsize=(22,15));gs=fig.add_gridspec(3,4,height_ratios=[1.35,1,1],hspace=.55,wspace=.65)
    for i,typ in enumerate(['UB','TR']):
        ax=fig.add_subplot(gs[0,i*2:(i+1)*2]);d=imp[imp.Type.eq(typ)].pivot(index='Feature',columns='Phase',values='mean_abs_SHAP');v=(d.After-d.Before).sort_values();ax.barh(v.index,v,color=np.where(v>=0,'#ff7f0e','#1f77b4'),edgecolor='black',linewidth=.5);ax.axvline(0,c='#777',lw=.7);ax.set_title(typ+' change in predictor contribution',loc='left');ax.set_xlabel('After − before mean |SHAP| (μm² cm⁻³)')
    for i,typ in enumerate(['UB','TR']):
        for j,f in enumerate(['BC','NO2','O3','WS']):dependence(fig.add_subplot(gs[i+1,j]),f,typ)
    save(fig,'Figure_4')

if __name__=='__main__':main()
