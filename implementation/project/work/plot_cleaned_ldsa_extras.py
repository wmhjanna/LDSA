from pathlib import Path
import sys,json
ROOT=Path('@PROJECT_ROOT@');sys.path.insert(0,str(ROOT/'work/ldsa_support/deps'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_cleaned_ldsa import scatter,head,save,manifest
OUT=ROOT/'outputs/ldsa_cleaned_reanalysis';ANA=OUT/'analysis';FIG=OUT/'figures';STAGE=ROOT/'work/ldsa_support'
tr=pd.read_csv(ANA/'trends.csv');ms=pd.read_csv(ANA/'monthly_series.csv',parse_dates=['date']);bp=pd.read_csv(ANA/'breakpoints.csv',parse_dates=['Breakpoint'])
# Regulatory and fleet measurements remain the archived contextual inputs; only the LDSA anchor changes.
date=bp.loc[bp.Site.eq('LND_TR')&bp.Metric.eq('LDSA'),'Breakpoint'].iloc[0]
source=(ROOT/'work/redraw_regulatory_traffic_context.py').read_text()
source=source.replace('/outputs/regulatory_traffic_context','/outputs/ldsa_cleaned_reanalysis/figures').replace('breakpoint = 2013 + 2 / 12',f'breakpoint = {date.year+(date.month-1)/12}').replace('Mar 2013',date.strftime('%b %Y')).replace('ax_b.axvline(2013,','ax_b.axvline(breakpoint,').replace('Regulatory_timing_traffic_fleet_context','Figure_S3')
(STAGE/'redraw_regulatory_updated.py').write_text(source)
exec(compile(source,str(STAGE/'redraw_regulatory_updated.py'),'exec'),{'__name__':'regulatory_replot'})
plt.rcParams['font.family']=['Times New Roman','DejaVu Serif']
import ast
import matplotlib.dates as mdates
from scipy.interpolate import UnivariateSpline
from ldsa_reference_style import configure
configure()
BSPLINE_SMOOTHING_FRACTION=.045;BSPLINE_POINTS_PER_MONTH=8;BSPLINE_MAX_GAP_DAYS=100
oldtree=ast.parse((ROOT/'work/redraw_lnd_sto_breakpoint_context_smooth.py').read_text());funcs=[n for n in oldtree.body if isinstance(n,ast.FunctionDef) and n.name=='bspline_segments'];exec(compile(ast.Module(body=funcs,type_ignores=[]),'original_bspline_display','exec'))
line_colors={'LDSA':'#202020','Ntotal':'#009E73','BC':'#7B3294','PM2.5':'#D55E00','PM10':'#E69F00','NO2':'#0072B2','Aitken':'#E69F00','Accumulation':'#7B3294'}
summaries=[]
fig,axs=plt.subplots(2,2,figsize=(14,9),layout='constrained')
for i,site in enumerate(['LND_TR','STO_TR']):
    d=ms[ms.Site.eq(site)].pivot(index='date',columns='Metric',values='value');b=bp.loc[bp.Site.eq(site)&bp.Metric.eq('LDSA'),'Breakpoint'].iloc[0]
    norm=d.copy()
    for c in d:
        baseline=d.loc[(d.index>=b-pd.DateOffset(months=24))&(d.index<b),c].mean();norm[c]=100*d[c]/baseline if baseline!=0 else np.nan
        pre=d.loc[(d.index>=b-pd.DateOffset(months=36))&(d.index<b),c];post=d.loc[(d.index>=b)&(d.index<b+pd.DateOffset(months=36)),c]
        summaries.append(dict(Site=site,Metric=c,Breakpoint=b,pre_months=pre.notna().sum(),post_months=post.notna().sum(),pre_mean=pre.mean(),post_mean=post.mean(),change_percent=100*(post.mean()/pre.mean()-1)))
    norm.to_csv(ANA/f'{site}_normalized_monthly.csv')
    for j,metrics in enumerate([['LDSA','PM2.5','BC','NO2'],['LDSA','Ntotal','Aitken','Accumulation']]):
        ax=axs[i,j]
        for c in metrics:
            if c in norm:
                for k,(xx,yy) in enumerate(bspline_segments(pd.Series(norm.index),norm[c].reset_index(drop=True))):ax.plot(xx,yy,label=c if k==0 else '_nolegend_',lw=1.9,color=line_colors[c])
        ax.axvline(b,c='#333333',ls='--',lw=1.25);ax.axhline(100,c='#9A9A9A',ls=(0,(2,2)),lw=.9);ax.xaxis.set_major_locator(mdates.YearLocator(2));ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'));ax.set_ylabel('Index (pre-breakpoint 24 months = 100)');ax.legend(frameon=False,ncol=2,fontsize=8);head(ax,chr(97+i*2+j),site)
save(fig,'LND_STO_synchronized_trajectories','monthly_series.csv; breakpoints.csv')
summary=pd.DataFrame(summaries);summary.to_csv(ANA/'city_36month_changes.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(13,5.7),layout='constrained')
order=['LDSA','Ntotal','Nucleation','Aitken','Accumulation'];labels=['LDSA','Total PNC','Nucleation PNC','Aitken PNC','Accumulation PNC']
for ax,site,city in zip(axs,['LND_TR','STO_TR'],['London','Stockholm']):
    d=summary[summary.Site.eq(site)].set_index('Metric');v=d.loc[order,'change_percent'];y=np.arange(5);ax.scatter(v,y,s=155,color='#2C7FB8',edgecolors='#2C7FB8',lw=.8,zorder=3)
    for x,yy in zip(v,y):ax.annotate(f'{x:.1f}%',(x,yy),xytext=(12,0),textcoords='offset points',va='center',fontsize=13,fontweight='bold')
    ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlim(min(-52,v.min()-8),max(3,v.max()+12));ax.axvline(0,c='#6E6E6E',lw=1.4);ax.grid(axis='x',color='#D9D9D9',lw=.9,alpha=.75);ax.set_xlabel('Post- vs pre-transition mean change (%)');b=bp.loc[bp.Site.eq(site)&bp.Metric.eq('LDSA'),'Breakpoint'].iloc[0];head(ax,'a' if site=='LND_TR' else 'b',f'{site} ({city}; {b:%b %Y} LDSA transition)');ax.yaxis.grid(False)
save(fig,'Size_mode_pre_post_changes','city_36month_changes.csv')
paired=[];coverage_rows=[]
fig,axs=plt.subplots(2,3,figsize=(14,8),layout='constrained')
for ax,metric in zip(axs.flat,['PM2.5','PM10','BC','NO2','SO2']):
    z=tr[tr.Phase.eq('Overlap_'+metric)].pivot(index=['Site','Type'],columns='Metric',values='Slope_Percent').reset_index();z['Metric']=metric
    included=[]
    for site in z.Site:
        q=ms[ms.Site.eq(site)&ms.Metric.isin(['LDSA',metric])].pivot(index='date',columns='Metric',values='value').dropna()
        n=len(q);span=(q.index.max().year-q.index.min().year)*12+q.index.max().month-q.index.min().month+1 if n else 0
        keep=n>=60 and n/span>=.6;included.append(keep);coverage_rows.append(dict(Site=site,Metric=metric,n_months=n,span_months=span,coverage=n/span if span else 0,primary_included=keep))
    z['primary_included']=included;paired.append(z.rename(columns={metric:'Other_slope','LDSA':'LDSA_slope'}));scatter(ax,z[z.primary_included],'LDSA',metric);ax.set(xlabel='LDSA slope (% yr⁻¹)',ylabel=metric+' slope (% yr⁻¹)');head(ax,chr(97+list(axs.flat).index(ax)),metric)
axs.flat[-1].axis('off');axs.flat[-1].legend(handles=[__import__('matplotlib').lines.Line2D([],[],color=c,marker=m,linestyle='none',label=t) for t,c,m in [('UB','#1683C7','o'),('TR','#F05A14','s'),('SUB','#6755B5','^'),('RB','#3C9A35','D')]],loc='center',frameon=False,ncol=2,fontsize=11);save(fig,'Exact_overlap_trend_comparison','trends.csv; monthly_series.csv');pd.concat(paired,ignore_index=True).to_csv(ANA/'exact_overlap_pairs.csv',index=False)
pd.DataFrame(coverage_rows).to_csv(ANA/'exact_overlap_coverage.csv',index=False)
order=['BC','PM10','PM2.5','NO2','NO','SO2','CO','O3']
phase=tr[tr.Phase.isin(['Before','After'])];median=phase[phase.Metric.isin(order)].groupby(['Metric','Phase']).Slope_Percent.median().unstack()
paired_phase=phase[phase.Metric.isin(order)].merge(phase[phase.Metric.eq('LDSA')][['Site','Phase','Slope_Percent']],on=['Site','Phase'],suffixes=('','_LDSA')).dropna(subset=['Slope_Percent','Slope_Percent_LDSA']);paired_phase['Same_sign']=paired_phase.Slope_Percent*paired_phase.Slope_Percent_LDSA>0
fraction=paired_phase.groupby(['Metric','Phase']).Same_sign.mean().unstack();fraction.to_csv(ANA/'phase_direction_fraction.csv')
fig,axs=plt.subplots(1,2,figsize=(12.6,5.4),layout='constrained')
for j,(ax,table,xlabel,reference) in enumerate(zip(axs,[median,fraction],['Median co-pollutant slope (% yr⁻¹)','Fraction of sites with same-sign LDSA slope'],[0,.5])):
    for yi,metric in enumerate(order):
        before=table.loc[metric,'Before'];after=table.loc[metric,'After'];ax.plot([before,after],[yi-.10,yi+.10],color='#8f8f8f',lw=1.25,zorder=1);ax.scatter(before,yi-.10,s=48,color='#2c7fb8',label='Before' if yi==0 else None,zorder=3);ax.scatter(after,yi+.10,s=48,color='#f47a2a',label='After' if yi==0 else None,zorder=3)
    ax.set_yticks(range(len(order)),order);ax.invert_yaxis();ax.set_xlabel(xlabel);ax.axvline(reference,color='#a6a6a6',ls='--',lw=1);ax.grid(axis='x',color='#dedede',lw=.8,alpha=.7);ax.text(-.08,1.015,chr(97+j),transform=ax.transAxes,fontsize=15,fontweight='bold')
axs[0].legend(loc='lower left',frameon=False,fontsize=11);axs[1].set_xlim(-.02,1.02)
save(fig,'Phase_copollutant_trends','trends.csv; phase_direction_fraction.csv')
if (ANA/'correlation_sensitivity.csv').exists():
    c=pd.read_csv(ANA/'correlation_sensitivity.csv');c=c[c.Subset.eq('All sites')].sort_values('Pearson_r');fig,ax=plt.subplots(figsize=(9,7));y=np.arange(len(c));ax.scatter(c.Pearson_r,y,label='Pearson r',c='#347C9E');ax.scatter(c.Spearman_rho,y+.16,label='Spearman ρ',c='#D59149');ax.hlines(y,c.leave_one_out_r_min,c.leave_one_out_r_max,color='#347C9E',alpha=.4);ax.set_yticks(y,c.Metric);ax.set(xlim=(-1.05,1.05),xlabel='Correlation coefficient');ax.axvline(0,c='#aaa',ls='--');ax.legend(frameon=False);head(ax,'','Full-record correlations and leave-one-out range');save(fig,'Correlation_robustness','correlation_sensitivity.csv')
previous=json.loads((ANA/'figure_manifest.json').read_text());previous+=manifest+[dict(figure='Figure_S3',source='archived regulatory and traffic values; updated breakpoints.csv')]
(ANA/'figure_manifest.json').write_text(json.dumps(list({x['figure']:x for x in previous}.values()),indent=2))
