from pathlib import Path
import json,sys
sys.path.insert(0,'@PROJECT_ROOT@/work/ldsa_support/deps')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from scipy.stats import pearsonr
ROOT=Path('@PROJECT_ROOT@');OUT=ROOT/'outputs/ldsa_cleaned_reanalysis';ANA=OUT/'analysis';FIG=OUT/'figures';STAGE=ROOT/'work/ldsa_support'
from ldsa_reference_style import C,MARKERS,configure,apply_style
ORDER=list(C)
plt.rcParams.update({'font.family':['Times New Roman','DejaVu Serif'],'font.size':10,'axes.labelweight':'bold','axes.titleweight':'bold','pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'stix'})
configure()
manifest=[]
def save(fig,name,source):
    apply_style(fig,name)
    fig.savefig(FIG/f'{name}.png',dpi=300,bbox_inches='tight',facecolor='white')
    for ext in ['pdf','svg']:fig.savefig(FIG/f'{name}.{ext}',bbox_inches='tight',facecolor='white')
    plt.close(fig);manifest.append(dict(figure=name,source=source));print(name,flush=True)
def head(ax,letter,title):
    ax.set_title(title,loc='left',pad=12);ax.text(-.09,1.03,letter,transform=ax.transAxes,fontweight='bold',fontsize=15);ax.grid(alpha=.16)
def boxes(ax,groups,labels,colors=None):
    vals=[np.asarray(v)[np.isfinite(v)] for v in groups]
    ax.boxplot(vals,positions=range(len(vals)),showfliers=False,widths=.5)
    for i,v in enumerate(vals):ax.scatter(i+np.linspace(-.10,.10,len(v)),v,s=18,c=(colors or ['#477DA1']*len(vals))[i],zorder=3)
    ax.set_xticks(range(len(vals)),[f'{l}\n(n={len(v)})' for l,v in zip(labels,vals)])
def scatter(ax,d,x,y):
    d=d.dropna(subset=[x,y]);a=d[x].to_numpy();b=d[y].to_numpy()
    for t,g in d.groupby('Type'):ax.scatter(g[x],g[y],label=t,c=C[t],marker=MARKERS[t],s=35,alpha=.85)
    if len(d)>2 and np.std(a)>0 and np.std(b)>0:
        r,p=pearsonr(a,b);q=np.polyfit(a,b,1);xx=np.array([min(a),max(a)]);ax.plot(xx,np.polyval(q,xx),'--',color='#777');text=f'r = {r:.3f}; P = {p:.3g}; n = {len(d)}'
    else:text=f'n = {len(d)}; correlation unavailable'
    ax.margins(y=.18);ax.text(.03,.96,text,transform=ax.transAxes,va='top',fontsize=9);ax.grid(alpha=.18)

def main():
    sp=pd.read_csv(ANA/'site_spatial_statistics.csv');sens=pd.read_csv(ANA/'size_range_sensitivity.csv');tr=pd.read_csv(ANA/'trends.csv');bp=pd.read_csv(ANA/'breakpoints.csv',parse_dates=['Breakpoint']);wk=pd.read_pickle(STAGE/'weekly.pkl');meta=pd.read_csv(ANA/'meta_summary.csv');sites=pd.read_csv(ANA/'meta_sites.csv')
    fig,ax=plt.subplots(figsize=(6.5,4.7),layout='constrained')
    vals=[sp.loc[sp.Type.eq(t),'Mean'].dropna().to_numpy() for t in ORDER]
    bp1=ax.boxplot(vals,positions=np.arange(4),patch_artist=True,showfliers=False,widths=.5)
    for patch,t in zip(bp1['boxes'],ORDER):patch.set(facecolor=C[t],alpha=.25)
    for i,(t,v) in enumerate(zip(ORDER,vals)):ax.scatter(i+np.linspace(-.10,.10,len(v)),v,color=C[t],marker=MARKERS[t],s=35,zorder=3)
    ax.set_xticks(range(4),[f'{t}\n(n={len(v)})' for t,v in zip(ORDER,vals)]);ax.grid(axis='y',alpha=.2);ax.set_xlabel('Station type');ax.set_ylabel(r'Site-mean LDSA ($\mu$m$^2$ cm$^{-3}$)');head(ax,'b','Distribution by station type');ax.xaxis.grid(False);save(fig,'Figure_1_b','site_spatial_statistics.csv')
    fig,ax=plt.subplots(figsize=(7,10),layout='constrained');d=sp.sort_values('Mean')
    for i,r in enumerate(d.itertuples()):ax.hlines(i,r.P25,r.P75,color=C[r.Type],lw=2);ax.scatter(r.Mean,i,c=C[r.Type],s=38,zorder=3);ax.scatter(r.P95,i,facecolors='white',edgecolors=C[r.Type],s=30,zorder=3)
    ax.set_yticks(range(len(d)),d.Site,fontsize=9);ax.set_xlabel(r'LDSA ($\mu$m$^2$ cm$^{-3}$)');ax.grid(axis='x',alpha=.2);head(ax,'c','Ranked site statistics');ax.yaxis.grid(False)
    ax.legend(handles=[Line2D([],[],color='#555',lw=2,label='P25–P75'),Line2D([],[],color='#555',marker='o',lw=0,label='Mean'),Line2D([],[],color='#555',marker='o',mfc='white',lw=0,label='P95')],loc='lower right',fontsize=9,frameon=False);save(fig,'Figure_1_c','site_spatial_statistics.csv')

    def timepanel(ax,site):
        d=wk[wk.Site.eq(site)].dropna(subset=['LDSA']);ax.scatter(d.date,d.LDSA,c='#95a5a6',s=16,alpha=.5)
        b=bp[bp.Site.eq(site)&bp.Metric.eq('LDSA')]
        if len(b):
            date=b.Breakpoint.iloc[0];ax.axvline(date,c='#333',ls='--');ax.text(.03,.94,date.strftime('Breakpoint: %Y-%m'),transform=ax.transAxes,va='top',fontsize=8)
            for phase,col in [('Before','#E74C3C'),('After','#2980B9')]:
                g=d[d.date.lt(date) if phase=='Before' else d.date.ge(date)];z=tr[tr.Site.eq(site)&tr.Metric.eq('LDSA')&tr.Phase.eq(phase)]
                if len(g)>1 and len(z):
                    x=(g.date-g.date.min()).dt.total_seconds()/86400/365.25;y=x*z.Slope_Unit.iloc[0];y+=g.LDSA.mean()-y.mean();ax.plot(g.date,y,c=col,lw=3)
        ax.set_title(site,loc='left');ax.set_ylabel('LDSA (μm² cm⁻³)');ax.grid(alpha=.15)
        locator=mdates.AutoDateLocator(minticks=3,maxticks=5);ax.xaxis.set_major_locator(locator);ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator));ax.tick_params(axis='x',labelsize=8)
    def forest(fig,spec):
        grid=spec.subgridspec(1,4,width_ratios=[.25,.23,.10,.42],wspace=.025);axs=[fig.add_subplot(grid[0,i]) for i in range(4)];rows=[]
        names={'UB':'Urban-Background','SUB':'Suburban-Background','TR':'Urban-Traffic','RB':'Regional-Background'}
        for t in ['UB','SUB','TR','RB']:
            rows.append(dict(kind='heading',label=names[t]));group=sites[sites.Type.eq(t)].sort_values('Site')
            for r in group.itertuples():rows.append(dict(kind='study',label=r.Site,effect=r.Slope_Percent,low=r.CI_Lower_Percent,high=r.CI_Upper_Percent,weight=r.weight_percent))
            r=meta[meta.Type.eq(t)].iloc[0];rows.append(dict(kind='summary',label='Random effects model',effect=r.effect,low=r.lower,high=r.upper,weight=group.weight_percent.sum()));rows.append(dict(kind='space'))
        r=meta[meta.Type.eq('Overall')].iloc[0];rows.append(dict(kind='summary',label='Random effects model',effect=r.effect,low=r.lower,high=r.upper,weight=100))
        study,ci,weight,ax=axs
        for a in axs:a.set_ylim(-.5,len(rows)+2)
        for a in axs[:3]:a.set_xlim(0,1);a.axis('off')
        study.text(0,len(rows)+1,'Study',fontsize=14,va='center');ci.text(.5,len(rows)+1,'[95%-CI]',fontsize=14,ha='center',va='center');weight.text(0,len(rows)+1,'Weight',fontsize=14,va='center')
        for i,r in enumerate(rows):
            y=len(rows)-i-1
            if r['kind']=='space':continue
            study.text(0,y,r['label'],fontsize=13 if r['kind']=='heading' else 12,va='center')
            if r['kind']=='heading':continue
            ci.text(.5,y,f"{r['effect']:.2f} [{r['low']:.2f}; {r['high']:.2f}]",fontsize=11.5,ha='center',va='center');weight.text(0,y,f"{r['weight']:.1f}%",fontsize=12,va='center')
            if r['kind']=='study':ax.hlines(y,r['low'],r['high'],color='#777777',lw=1.05);ax.scatter(r['effect'],y,s=28+22*r['weight'],marker='s',facecolor='#f2b01e',edgecolor='#9a6b00',lw=.9,zorder=3)
            else:ax.add_patch(Polygon([(r['low'],y),(r['effect'],y+.25),(r['high'],y),(r['effect'],y-.25)],closed=True,facecolor='#b7c7c9',edgecolor='#4f5b5d',lw=1))
        ax.set_yticks([]);ax.axvline(0,c='#777',lw=1.1);ax.axvline(meta.loc[meta.Type.eq('Overall'),'effect'].iloc[0],c='#999',ls=':',lw=1.1);ax.spines[['left','right','top']].set_visible(False);ax.set_xlabel(r'% yr$^{-1}$',fontsize=13,labelpad=10)
        test=pd.read_csv(ANA/'meta_subgroup_test.csv').iloc[0];study.text(0,-.055,f"Subgroup difference among UB, SUB, TR and RB: Q = {test.Q:.2f}, P = {test.p:.3f}",transform=study.transAxes,fontsize=12,va='top')
        study.text(-.035,1.025,'c',transform=study.transAxes,fontsize=16,va='bottom');study.text(.085,1.025,'Random-effects synthesis of full-record LDSA slopes by station class',transform=study.transAxes,fontsize=14,va='bottom')
    fig=plt.figure(figsize=(22,12.5));gs=fig.add_gridspec(2,2,left=.055,right=.985,bottom=.11,top=.92,width_ratios=[.27,.65],hspace=.30,wspace=.16)
    ax=fig.add_subplot(gs[0,0]);timepanel(ax,'LND_TR');head(ax,'a','');ax.set_title('LND_TR (LDSA (17–500 nm))',fontsize=16,fontweight='bold');ax.xaxis.set_major_locator(mdates.YearLocator());ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'));ax.grid(linestyle=':',alpha=.4)
    line_data=wk[wk.Site.eq('LND_TR')].dropna(subset=['LDSA']);ax.plot(line_data.date,line_data.LDSA,color='#95a5a6',alpha=.2,lw=1,zorder=0)
    annotation=ax.texts[0];annotation.set_position((.45,.66));phase_slopes=tr[tr.Site.eq('LND_TR')&tr.Metric.eq('LDSA')].set_index('Phase').Slope_Unit
    annotation.set_text(annotation.get_text()+'\n'+r'Pre-slope: '+f'{phase_slopes.Before:.2f}'+r' $\mu$m$^2$ cm$^{-3}$ yr$^{-1}$'+'\n'+r'Post-slope: '+f'{phase_slopes.After:.2f}'+r' $\mu$m$^2$ cm$^{-3}$ yr$^{-1}$')
    ax.legend(handles=[Line2D([],[],color=c,lw=3,label=p) for c,p in zip(['#E74C3C','#2980B9'],['Pre-Breakpoint Trend','Post-Breakpoint Trend'])],loc='upper right',frameon=False,fontsize=9)
    ax=fig.add_subplot(gs[1,0]);groups=[];rng=np.random.default_rng(6)
    for j,per in enumerate(['Before','After']):
        for i,t in enumerate(ORDER):
            v=tr.loc[tr.Type.eq(t)&tr.Phase.eq(per)&tr.Metric.eq('LDSA'),'Slope_Percent'].dropna().to_numpy();groups.append(v);x=i+(j-.5)*.32
            ax.boxplot([v],positions=[x],widths=.28,showfliers=False,medianprops={'color':'black'},manage_ticks=False);ax.scatter(x+rng.uniform(-.06,.06,len(v)),v,c=['#E74C3C','#2980B9'][j],s=20,zorder=3)
    ax.set_xticks(range(4),ORDER);ax.set_xlabel('Station type');ax.axhline(0,c='#888',ls='--');ax.set_ylabel(r'LDSA slope (% yr$^{-1}$)');head(ax,'b','Pre- and post-breakpoint slopes');ax.set_yscale('symlog',linthresh=10);ax.text(.98,.97,'Symmetric log scale',transform=ax.transAxes,ha='right',va='top',fontsize=8)
    ax.legend(handles=[Line2D([],[],color=c,marker='o',lw=0,label=p) for c,p in zip(['#E74C3C','#2980B9'],['Before','After'])],loc='upper left',frameon=False,fontsize=9)
    values=np.concatenate(groups);ax.set_ylim(min(-20,float(values.min())*1.4),max(100,float(values.max())*1.5));forest(fig,gs[:,1]);save(fig,'Figure_2','trends.csv; breakpoints.csv; meta_sites.csv; weekly.pkl')

    full=tr[tr.Phase.eq('Full')].pivot(index=['Site','Type'],columns='Metric',values='Slope_Percent').reset_index()
    correlations=[]
    for metric in [c for c in full if c not in ['Site','Type','LDSA']]:
        d=full[['Site','Type','LDSA',metric]].dropna();r,p=pearsonr(d.LDSA,d[metric]) if len(d)>2 and d.LDSA.std()>0 and d[metric].std()>0 else (np.nan,np.nan)
        correlations.append(dict(Metric=metric,n=len(d),r=r,p=p,same_sign=int((d.LDSA*d[metric]>0).sum()),opposite_sign=int((d.LDSA*d[metric]<0).sum())))
    cor=pd.DataFrame(correlations);cor.to_csv(ANA/'slope_correlations.csv',index=False);full.to_csv(ANA/'full_trend_pairs.csv',index=False)
    fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained');d=cor.sort_values('r');axs[0,0].barh(d.Metric+' (n='+d.n.astype(str)+')',d.r,color='#4F89A7');axs[0,0].set_xlabel('Pearson r');head(axs[0,0],'a','Whole-record slope associations');axs[0,1].barh(d.Metric,d.same_sign,label='Same sign',color='#5A987F');axs[0,1].barh(d.Metric,d.opposite_sign,left=d.same_sign,label='Opposite sign',color='#D7A15F');axs[0,1].set_xlabel('Number of paired sites');axs[0,1].legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.13),ncol=2);head(axs[0,1],'b','Trend-direction concordance')
    for ax,per,l in zip(axs[1],['Before','After'],['c','d']):
        z=tr[tr.Phase.eq(per)].pivot(index=['Site','Type'],columns='Metric',values='Slope_Percent').reset_index();scatter(ax,z,'LDSA','BC');ax.set(xlabel='LDSA slope (% yr⁻¹)',ylabel='BC slope (% yr⁻¹)');head(ax,l,per+' breakpoint')
    save(fig,'Figure_3','trends.csv; slope_correlations.csv')

    import importlib.util
    def old_module(name):
        spec=importlib.util.spec_from_file_location(name,ROOT/f'outputs/{name}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
    s1=old_module('redraw_figure_s1');s1.configure_style();plt.rcParams['font.family']=['Times New Roman','DejaVu Serif'];plt.rcParams['font.weight']='normal'
    fig,axs=plt.subplots(2,3,figsize=(11.6,6.65))
    for row,stat in enumerate(['Mean','Median']):
        wide=sens[sens.Region.eq('Total')].pivot(index='Site',columns='Range',values=stat)
        for j,(a,b) in enumerate(s1.COMPARISONS):
            pair=wide[[a,b]].dropna();pair.columns=['x','y'];lo=min(pair.min())*.9;hi=max(pair.max())*1.10
            s1.draw_panel(axs[row,j],pair,a,b,stat.lower(),chr(97+row*3+j),(lo,hi),False)
    fig.subplots_adjust(left=.075,right=.992,bottom=.095,top=.985,wspace=.24,hspace=.22);save(fig,'Figure_S1','size_range_sensitivity.csv')
    s2=old_module('redraw_figure_s2');s2.configure_style();plt.rcParams['font.family']=['Times New Roman','DejaVu Serif'];fig,axs=plt.subplots(2,2,figsize=(9.3,6))
    for j,stat in enumerate(['Mean','Median']):
        w=sens.pivot(index='Site',columns=['Region','Range'],values=stat)
        vals=[(100*(w[r,'10-800']-w[r,'10-500'])/w[r,'10-800']).dropna().to_numpy() for r in ['HA','TB','ALV','Total']]
        extent=np.concatenate(vals);limits=(min(-2.5,extent.min()-2),max(54.5,extent.max()+2))
        s2.draw_boxplot(axs[0,j],vals,['HA','TB','ALV','Total'],['#1f77b4','#ff7f0e','#2ca02c','#d62728'],chr(97+j),'Site-level '+stat.lower()+'s','Contribution of 500–800 nm (%)',limits,10,101+j)
        vals=[(100*(w[r,'10-800']/w['Total','10-800']-w[r,'10-500']/w['Total','10-500'])).dropna().to_numpy() for r in ['HA','TB','ALV']];extent=np.concatenate(vals);limits=(min(-3.4,extent.min()-.5),max(6.8,extent.max()+.5))
        s2.draw_boxplot(axs[1,j],vals,['HA','TB','ALV'],['#1f77b4','#ff7f0e','#2ca02c'],chr(99+j),'10–800 nm minus 10–500 nm ('+stat.lower()+'s)','Change in fraction (percentage points)',limits,2,103+j)
    fig.subplots_adjust(left=.095,right=.992,bottom=.085,top=.975,wspace=.23,hspace=.24);save(fig,'Figure_S2','size_range_sensitivity.csv');configure()
    # Re-render the archived regulatory data with the newly estimated LND breakpoint.
    source=(ROOT/'work/redraw_regulatory_traffic_context.py').read_text()
    date=bp.loc[bp.Site.eq('LND_TR')&bp.Metric.eq('LDSA'),'Breakpoint'].iloc[0]
    (STAGE/'new_lnd_breakpoint.txt').write_text(str(date))
    chosen=['LND_TR','STO_TR','DRE_TR','LEI_UB'];fig,axs=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for ax,site in zip(axs.flat,chosen):timepanel(ax,site)
    save(fig,'Figure_S4','weekly.pkl; breakpoints.csv; trends.csv')
    b=bp[bp.Metric.eq('LDSA')].copy();b['year']=b.Breakpoint.dt.year+(b.Breakpoint.dt.dayofyear-1)/365.25;b=b.sort_values('year')
    fig,axs=plt.subplots(1,2,figsize=(14,5),layout='constrained');axs[0].scatter(range(len(b)),b.year,c=b.Type.map(C));axs[0].set_xticks(range(len(b)),b.Site,rotation=60,ha='right',fontsize=7);axs[0].set_ylabel('Breakpoint year');boxes(axs[1],[b[b.Type.eq(t)].year for t in ORDER],ORDER,list(C.values()));head(axs[0],'a','Site-specific breakpoints');head(axs[1],'b','Station classes');save(fig,'Figure_S5','breakpoints.csv')
    for page,start in enumerate(range(0,len(b),6),1):
        fig,axs=plt.subplots(3,2,figsize=(14,14),layout='constrained')
        for ax,site in zip(axs.flat,b.Site.tolist()[start:start+6]):timepanel(ax,site)
        for ax in list(axs.flat)[len(b.Site.tolist()[start:start+6]):]:ax.axis('off')
        save(fig,f'Figure_S6_page_{page:02}','weekly.pkl; breakpoints.csv; trends.csv')
    bw=bp.pivot(index=['Site','Type'],columns='Metric',values='Breakpoint').reset_index()
    for c in bw.columns[2:]:bw[c]=bw[c].dt.year+(bw[c].dt.dayofyear-1)/365.25
    fig,axs=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    for ax,metric in zip(axs,['BC','PM2.5','NO2']):
        scatter(ax,bw,'LDSA',metric);lims=[2000,2025];ax.plot(lims,lims,c='#555',lw=1);ax.set(xlabel='LDSA breakpoint year',ylabel=metric+' breakpoint year');head(ax,chr(97+list(axs).index(ax)),metric)
    save(fig,'Figure_S7','breakpoints.csv')
    fig,axs=plt.subplots(2,2,figsize=(13,10),layout='constrained')
    for ax,t in zip(axs.flat,ORDER):
        g=bp[bp.Type.eq(t)&bp.Metric.ne('LDSA')].groupby('Metric').offset_months.agg(['mean','min','max','count']).dropna().sort_values('mean');ax.barh(g.index,g['mean'],color=C[t]);ax.errorbar(g['mean'],g.index,xerr=[g['mean']-g['min'],g['max']-g['mean']],fmt='none',ecolor='#555',capsize=3);ax.axvline(0,c='#444',ls='--');ax.set_xlabel('Offset from LDSA breakpoint (months)');head(ax,chr(97+ORDER.index(t)),t)
    save(fig,'Figure_S8','breakpoints.csv')
    fig,axs=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for ax,metric in zip(axs.flat,['Accumulation','Aitken','BC','PM10','PM2.5','O3']):scatter(ax,full,'LDSA',metric);ax.set(xlabel='LDSA slope (% yr⁻¹)',ylabel=metric+' slope (% yr⁻¹)');head(ax,chr(97+list(axs.flat).index(ax)),metric)
    save(fig,'Figure_S9','trends.csv')
    (ANA/'figure_manifest.json').write_text(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
