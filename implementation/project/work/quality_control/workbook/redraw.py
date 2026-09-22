"""PNG-only redraw from cleaned analysis, with manuscript-embedded figure templates.
Run with figure groups as arguments: main12 main3 main4 s12 s3 s456 s7 s8 s9 s10 s11.
Source data are read-only. New pairwise calculations and rendering audit stay in work/.
"""
from pathlib import Path
import sys, ast, textwrap, math, json, shutil, re, importlib.util
ROOT=Path('@PROJECT_ROOT@')
sys.path[:0]=[str(ROOT/'work/ldsa_support/deps'),str(ROOT/'work')]
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Patch
from scipy.stats import pearsonr
from scipy.interpolate import UnivariateSpline
import seaborn as sns
import ruptures as rpt
ANA=ROOT/'outputs/quality_control/workbook/analysis'
WORK=ROOT/'work/quality_control/workbook'
DEST=Path('@DELIVERY_ROOT@/workbook_analysis_figures')
DEST.mkdir(parents=True,exist_ok=True)
STAGE=ROOT/'work/ldsa_support'
sp=pd.read_csv(ANA/'site_spatial_statistics.csv'); sens=pd.read_csv(ANA/'size_range_sensitivity.csv')
tr=pd.read_csv(ANA/'trends.csv');ms=pd.read_csv(ANA/'monthly_series.csv',parse_dates=['date'])
bp=pd.read_csv(ANA/'breakpoints.csv',parse_dates=['Breakpoint']);wk=pd.read_pickle(WORK/'weekly.pkl')
meta=pd.read_csv(ANA/'meta_summary.csv');sites=pd.read_csv(ANA/'meta_sites.csv')
assert len(sp)==34 and sp.Physical_site.nunique()==34
assert {'PAR_SUB','PRA_SUB'}.issubset(set(sp.Site)) and sp.Site.is_unique
C={'UB':'#1683C7','TR':'#F05A14','SUB':'#6755B5','RB':'#3C9A35'};MARKERS={'UB':'o','TR':'s','SUB':'^','RB':'D'};ORDER=list(C)
manifest=[]
def configure():
    mpl.rcdefaults()
    plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'font.weight':'bold','axes.labelweight':'bold','axes.titleweight':'bold','mathtext.fontset':'stix','axes.unicode_minus':True})
def save(fig,name,source='',dpi=300):
    if name=='Figure_2':
        time_ax=fig.axes[0]
        time_ax.tick_params(axis='x',labelsize=9)
        annotation=time_ax.texts[0]
        annotation.set_fontsize(9)
        fig.canvas.draw()
        legend_box=time_ax.get_legend().get_window_extent(fig.canvas.get_renderer())
        legend_left=min(handle.get_window_extent(fig.canvas.get_renderer()).x0 for handle in time_ax.get_legend().legend_handles)
        legend_bottom=time_ax.transAxes.inverted().transform((legend_left,legend_box.y0))
        annotation.set_position((legend_bottom[0],legend_bottom[1]-.025))
        annotation.set_ha('left');annotation.set_va('top');annotation.set_multialignment('left')
        # Bring the forest axis closer to the already compact text columns.
        forest_ax=fig.axes[-1]
        pos=forest_ax.get_position()
        forest_ax.set_position([pos.x0-.045,pos.y0,pos.width,pos.height])
    if name=='Figure_S1':
        ax=fig.axes[3]
        limits=ax.get_ylim()
        ax.set_yticks([tick for tick in ax.get_yticks() if not np.isclose(tick,80)])
        ax.set_ylim(limits)
    if name=='Figure_S2':
        for ax in fig.axes:
            ax.title.set_fontweight('bold');ax._left_title.set_fontweight('bold');ax._right_title.set_fontweight('bold')
    if name=='Figure_3':
        legs=[ax.get_legend() for ax in fig.axes if ax.get_legend() is not None]
        le,lf=legs[-2:];fig.canvas.draw();ren=fig.canvas.get_renderer();be=le.get_window_extent(ren);bf=lf.get_window_extent(ren);y=(be.y0+bf.y0)/2
        for leg,box in [(le,be),(lf,bf)]:
            anchor=leg.get_bbox_to_anchor().transformed(fig.transFigure.inverted());leg.set_bbox_to_anchor((anchor.x0,anchor.y0+(y-box.y0)/fig.bbox.height),transform=fig.transFigure)
        fig.canvas.draw();assert abs(le.get_window_extent(fig.canvas.get_renderer()).y0-lf.get_window_extent(fig.canvas.get_renderer()).y0)<1
    target=DEST/(name+'.png');target.parent.mkdir(parents=True,exist_ok=True)
    for ax in fig.axes:
        for artist in list(ax.texts):
            if artist.get_text().startswith('2013–2016: Total'):artist.remove()
    for artist in fig.findobj(mpl.text.Text):
        artist.set_text(artist.get_text().replace('PRA_2012-2016_SUB','PRA_SUB').replace('PRA_2017-2019_SUB','PRA_SUB').replace('μm²/cm³','μm² cm⁻³'))
    for ax in fig.axes:
        for axis in ['x','y']:
            labels=getattr(ax,'get_'+axis+'ticklabels')();texts=[t.get_text() for t in labels]
            if any('PRA_' in t for t in texts):
                labels=[t.replace('PRA_2012-2016_SUB','PRA_SUB').replace('PRA_2017-2019_SUB','PRA_SUB') for t in texts]
                getattr(ax,'set_'+axis+'ticklabels')(labels)
    for ax in fig.axes:
        for axis in [ax.xaxis,ax.yaxis]:
            formatter=axis.get_major_formatter()
            if isinstance(formatter,mpl.ticker.FixedFormatter):
                formatter.seq=[str(t).replace('PRA_2012-2016_SUB','PRA_SUB').replace('PRA_2017-2019_SUB','PRA_SUB') for t in formatter.seq]
    fig.savefig(target,dpi=dpi,bbox_inches='tight',facecolor='white')
    manifest.append({'figure':name,'source':source,'dpi':dpi,'axes':len(fig.axes),'text':[{'text':t.get_text(),'pt':t.get_fontsize(),'family':t.get_fontfamily()} for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]})
    plt.close(fig);print('SAVED',name,flush=True)
def head(ax,letter,title):
    ax.set_title(title,loc='left',pad=12);ax.text(-.09,1.03,letter,transform=ax.transAxes,fontweight='bold',fontsize=15);ax.grid(alpha=.16)
def load_module(path):
    spec=importlib.util.spec_from_file_location('template_'+path.stem,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main12(only_fig1=False):
    configure()
    fig,ax=plt.subplots(figsize=(6.5,4.7));plt.rcParams['font.weight']='normal'
    vals=[sp.loc[sp.Type.eq(t),'Mean'].dropna().to_numpy() for t in ORDER]
    b=ax.boxplot(vals,positions=np.arange(4),patch_artist=True,showfliers=False,widths=.5,medianprops={'color':'black','linewidth':1.5},whiskerprops={'color':'#666','linewidth':1.2},capprops={'color':'#666','linewidth':1.2})
    rng=np.random.default_rng(20260909)
    for i,(p,t,v) in enumerate(zip(b['boxes'],ORDER,vals)):
        p.set(facecolor=mpl.colors.to_rgba(C[t],.22),edgecolor=C[t],linewidth=1.5)
        ax.scatter(i+rng.uniform(-.07,.07,len(v)),v,s=35,color=C[t],edgecolors='white',linewidths=.35,zorder=3)
    ax.set_xticks(range(4),ORDER,fontsize=12);ax.set_xlabel('Station type',fontsize=12,fontweight='normal');ax.set_ylabel(r'Site-mean total LDSA ($\mu$m$^2$ cm$^{-3}$)',fontsize=12,fontweight='normal')
    ax.tick_params(labelsize=11);ax.grid(axis='y',ls='--',alpha=.55);ax.set_axisbelow(True);ax.text(-.13,1.02,'b',transform=ax.transAxes,fontsize=18,fontweight='bold')
    for t in ax.get_xticklabels()+ax.get_yticklabels():t.set_fontweight('normal')
    save(fig,'Figure_1_b','site_spatial_statistics.csv',400)
    configure();fig,ax=plt.subplots(figsize=(6.1,10));d=sp.sort_values('Mean',ascending=False)
    for i,r in enumerate(d.itertuples()):
        ax.hlines(i,r.P25,r.P75,color=C[r.Type],lw=1.5);ax.scatter(r.Mean,i,c=C[r.Type],s=33,edgecolors='white',linewidths=.6,zorder=3);ax.scatter(r.P95,i,facecolors='white',edgecolors=C[r.Type],s=31,linewidths=1.2,zorder=3)
    ax.set_yticks(range(len(d)),d.Site,fontsize=9);ax.invert_yaxis();ax.set_xlabel(r'Total LDSA ($\mu$m$^2$ cm$^{-3}$)',fontsize=11);ax.set_title('Ranked site statistics for total LDSA',fontsize=13,pad=12);ax.grid(axis='x',ls='--',alpha=.3);ax.set_axisbelow(True)
    leg=ax.legend(handles=[Line2D([],[],color='#555',lw=1.5,label='P25–P75'),Line2D([],[],color='black',marker='o',lw=0,label='Mean'),Line2D([],[],color='black',marker='o',mfc='white',lw=0,label='P95')],title='Statistic',loc='lower right',bbox_to_anchor=(1,.16),fontsize=9,title_fontsize=10,frameon=False)
    ax.add_artist(leg);ax.legend(handles=[Line2D([],[],color=C[t],marker='o',lw=0,label=t) for t in ORDER],title='Site type',loc='lower right',fontsize=9,title_fontsize=10,frameon=False)
    ax.text(-.15,1.02,'c',transform=ax.transAxes,fontsize=18,fontweight='bold');save(fig,'Figure_1_c','site_spatial_statistics.csv',400)
    configure()
    if only_fig1:return
    code=(ROOT/'work/plot_cleaned_ldsa.py').read_text().split('    def timepanel',1)[1].split('\n    full=',1)[0]
    code=textwrap.dedent('    def timepanel'+code)
    code=code.replace('figsize=(22,12.5)','figsize=(15.8,12.5)').replace("wspace=.16)","wspace=.05)")
    code=code.replace("head(ax,'b','Pre- and post-breakpoint slopes')","head(ax,'b','Pre- and post-breakpoint LDSA slopes')")
    code=code.replace("for i,t in enumerate(ORDER):","for i,t in enumerate(['TR','UB','SUB','RB']):").replace("ax.set_xticks(range(4),ORDER)","ax.set_xticks(range(4),['TR','UB','SUB','RB'])")
    code=code.replace("ax.set_title('LND_TR (LDSA (17–500 nm))',fontsize=16", "ax.set_title('LND_TR (LDSA (17–500 nm))',fontsize=14")
    code=code.replace("ax.set_ylabel('LDSA (μm² cm⁻³)')", "ax.set_ylabel(r'LDSA ($\\mu$m$^2$ cm$^{-3}$)')")
    code=code.replace("label='Random effects model',effect=r.effect,low=r.lower,high=r.upper,weight=group.weight_percent.sum()", "label=('Single site' if len(group)==1 else 'Random effects model'),effect=r.effect,low=r.lower,high=r.upper,weight=group.weight_percent.sum()")
    code=code.replace('LND_TR (LDSA (17–500 nm))','LND_TR (LDSA)')
    code=code.replace('width_ratios=[.25,.23,.10,.42]','width_ratios=[.18,.23,.10,.49]')
    code=code.replace('ci.text(.5,','ci.text(.25,')
    code=code.replace('weight.text(0,','weight.text(-.5,')
    exec(compile(code,'main2_template','exec'),globals().copy())

def fig1():
    main12(only_fig1=True)

def pair_tables():
    allpairs=[];offset=[];corr=[]; sheets={}
    metric_order=['PM2.5','PM10','BC','NO2','SO2','NO','CO','O3']
    for metric in metric_order:
        rows=[]
        for site,g in tr[tr.Phase.eq('Overlap_'+metric)].groupby('Site'):
            if not {'LDSA',metric}.issubset(set(g.Metric)):continue
            g=g.set_index('Metric');q=ms[ms.Site.eq(site)&ms.Metric.isin(['LDSA',metric])].pivot(index='date',columns='Metric',values='value').dropna().sort_index()
            n=len(q);span=(q.index[-1].year-q.index[0].year)*12+q.index[-1].month-q.index[0].month+1 if n else 0
            l=g.loc['LDSA'];p=g.loc[metric]
            def state(r):
                if not np.isfinite([r.CI_Lower_Percent,r.CI_Upper_Percent]).all():return 'Stable/uncertain'
                return 'Declining' if r.CI_Upper_Percent<0 else 'Increasing' if r.CI_Lower_Percent>0 else 'Stable/uncertain'
            row={'Site':site,'Station type':site.rsplit('_',1)[1],'Pollutant':metric,'Pollutant slope (% yr^-1) [X]':p.Slope_Percent,'LDSA slope (% yr^-1) [Y]':l.Slope_Percent,'n_months':n,'coverage':n/span if span else 0,'LDSA_class':state(l),'Pollutant_class':state(p)}
            row['primary_included']=n>=60 and row['coverage']>=.5
            row['sensitivity_included']=n>=36 and row['coverage']>=.5
            rows.append(row);allpairs.append(row)
            if row['primary_included']:
                bdates={}
                for col in ['LDSA',metric]:
                    y=q[col].to_numpy();z=(y-y.mean())/y.std()
                    k=rpt.KernelCPD(kernel='linear',min_size=12,jump=1).fit(z.reshape(-1,1)).predict(n_bkps=1)[0];bdates[col]=q.index[k]
                off=(bdates[metric].year-bdates['LDSA'].year)*12+bdates[metric].month-bdates['LDSA'].month
                offset.append({'Site':site,'Station type':row['Station type'],'Pollutant':metric,'Offset Δt (months) [Y]':off,'LDSA_breakpoint':bdates['LDSA'],'Metric_breakpoint':bdates[metric],'n_months':n,'min_segment_months':12})
        primary=[r for r in rows if r['primary_included']]
        letter={'PM2.5':'a','PM10':'b','BC':'c','NO2':'d'}.get(metric)
        if letter:sheets[f'Fig3{letter}_{metric}']=primary
        full=tr[tr.Phase.eq('Full')&tr.Metric.isin(['LDSA',metric])].pivot(index='Site',columns='Metric',values='Slope_Percent').dropna()
        for strategy,d in [('Full available record',full),('Exact overlap ≥60 months',pd.DataFrame(primary)),('Exact overlap ≥36 months',pd.DataFrame([r for r in rows if r['sensitivity_included']]))]:
            if strategy=='Full available record':a=d[metric];b=d.LDSA
            else:a=d.get('Pollutant slope (% yr^-1) [X]',pd.Series(dtype=float));b=d.get('LDSA slope (% yr^-1) [Y]',pd.Series(dtype=float))
            valid=np.isfinite(a)&np.isfinite(b);a=a[valid];b=b[valid]
            corr.append({'Metric':metric,'Strategy':strategy,'n':len(a),'Pearson_r':pearsonr(a,b).statistic if len(a)>=3 and np.std(a)>0 and np.std(b)>0 else np.nan})
    e=[]
    for m in metric_order[:4]:
        d=[r for r in allpairs if r['Pollutant']==m and r['primary_included']];n=len(d)
        vals={'Same classification (%)':0,'LDSA decline only (%)':0,'Pollutant decline only (%)':0,'Other mismatch (%)':0}
        for r in d:
            l,p=r['LDSA_class'],r['Pollutant_class'];key='Same classification (%)' if l==p else 'LDSA decline only (%)' if l=='Declining' else 'Pollutant decline only (%)' if p=='Declining' else 'Other mismatch (%)';vals[key]+=100/n
        e.append({'Pollutant':m,'n':n,**vals})
    sheets['Fig3e_plot_data']=e;sheets['Fig3f_offsets_raw']=offset
    sheets['Fig3f_summary']=[{'Pollutant':m,'n':len(v),'Median |Δt| (months)':np.median(np.abs(v)) if v else np.nan} for m in metric_order[:5] for v in [[r['Offset Δt (months) [Y]'] for r in offset if r['Pollutant']==m]]]
    for name,rows in [('exact_overlap_pairs',allpairs),('classification',e),('exact_overlap_breakpoints_min12',offset),('correlation_temporal_sensitivity',corr)]:pd.DataFrame(rows).to_csv(WORK/(name+'.csv'),index=False)
    return sheets,pd.DataFrame(corr)

def main3():
    configure();sheets,_=pair_tables();src=(ROOT/'work/redraw_fig3_all_panels.py').read_text();prefix,body=src.split('wb = load_workbook(SOURCE, data_only=True, read_only=True)',1)
    scope=globals().copy();exec(compile(prefix,'fig3_functions','exec'),scope)
    scope.update(wb=None,read_sheet=lambda wb,key:sheets[key]);body=body.split('\nbase = OUTPUT_DIR',1)[0]
    offsets=[r['Offset Δt (months) [Y]'] for r in sheets['Fig3f_offsets_raw'] if r['Pollutant'] in ['PM2.5','PM10','BC','NO2','SO2']]
    bottom=min(-85,min(offsets)-10);top=max(100,max(offsets)+28)
    body=body.replace('ax.set_ylim(-85, 100)',f'ax.set_ylim({bottom}, {top})').replace('        94,',f'        {top-6},').replace('        85,',f'        {top-15},')
    # Label collision placement must use the final enlarged point size, not enlarge afterwards.
    body=body.replace('    text_artist.set_fontsize(text_artist.get_fontsize() + 4)','    if text_artist not in site_label_artists:\n        text_artist.set_fontsize(text_artist.get_fontsize() + 4)')
    body=body.replace('# Panel e: trend-classification agreement.',"site_label_artists=[a for ax in fig.axes for a in ax.texts if isinstance(a, mpl.text.Annotation)]\n# Panel e: trend-classification agreement.")
    # Keep source's enlarged labels and re-run collision detection at the final size.
    scope['place_site_labels']=place_final_site_labels
    exec(compile(body,'fig3_current_data','exec'),scope);save(scope['fig'],'Figure_3','exact_overlap_pairs.csv; exact_overlap_breakpoints_min12.csv',500)

def place_final_site_labels(fig,ax,records,stats):
    # lines, which made crowded groups fall back to overlapping label positions).
    stats.set_fontsize(15.5)
    fig.canvas.draw();renderer=fig.canvas.get_renderer();bounds=ax.get_window_extent(renderer)
    occupied=[stats.get_window_extent(renderer).expanded(1.08,1.1)]
    pointboxes=[mpl.transforms.Bbox.from_bounds(*(ax.transData.transform((r['x'],r['y']))-7),14,14) for r in records]
    candidates=[(7,6),(7,-8),(-7,6),(-7,-8)]
    for radius in [15,23,32,43,56,72,90]:
        candidates.extend([(radius,0),(-radius,0),(0,radius),(0,-radius),(radius,radius),(-radius,radius),(radius,-radius),(-radius,-radius)])
    # Prefer labels in the densest cluster first; distant outliers are easy.
    for rec in sorted(records,key=lambda r:(r['y'],r['x'])):
        best=None
        for dx,dy in candidates:
            ann=ax.annotate(rec['label'],(rec['x'],rec['y']),xytext=(dx,dy),textcoords='offset points',ha='left' if dx>=0 else 'right',va='bottom' if dy>=0 else 'top',fontsize=11.2,fontweight='bold',color='#222222',bbox=dict(facecolor='white',edgecolor='none',alpha=.82,pad=.08),arrowprops=dict(arrowstyle='-',color='#777777',linewidth=.45,shrinkA=1.5,shrinkB=4),zorder=6,clip_on=True)
            ann.update_positions(renderer);ann.update_bbox_position_size(renderer);box=ann.get_bbox_patch().get_window_extent(renderer)
            inside=box.x0>=bounds.x0+2 and box.x1<=bounds.x1-2 and box.y0>=bounds.y0+2 and box.y1<=bounds.y1-2
            collision=any(box.overlaps(b) for b in occupied+pointboxes)
            if inside and not collision:
                occupied.append(box.expanded(1.08,1.15));best=ann;break
            ann.remove()
        if best is None:raise RuntimeError('No collision-free label position: '+rec['label'])
    stats.set_fontsize(11.5)

def main4():
    configure();path=ROOT/'work/redraw_lnd_sto_breakpoint_context_smooth.py';src=path.read_text().replace('import nature_analysis_123 as source','')
    scope=globals().copy();scope['__name__']='figure_template';scope['__file__']=str(path);exec(compile(src,str(path),'exec'),scope)
    scope['BREAKPOINTS']={s:bp.loc[bp.Site.eq(s)&bp.Metric.eq('LDSA'),'Breakpoint'].iloc[0] for s in ['LND_TR','STO_TR']}
    scope['SITE_CONFIG']['LND_TR']['note']='2013–2016: Road VKM −2.4%; diesel-car +6.8%; diesel-LGV +4.3%\nDiesel passenger cars Euro 5/6: 62.4% (2016)'
    scope['set_style']();fig,axes=plt.subplots(2,2,figsize=(12.8,8.2),gridspec_kw={'width_ratios':[1.25,1]})
    for i,site in enumerate(['LND_TR','STO_TR']):
        data=ms[ms.Site.eq(site)].pivot(index='date',columns='Metric',values='value').reset_index()
        scope['draw_time_series'](axes[i,0],site,data);scope['draw_change_panel'](axes[i,1],site,data)
        for artist in list(axes[i,1].texts):
            if artist.get_text()==scope['SITE_CONFIG'][site]['note']:artist.remove()
    fig.subplots_adjust(left=.085,right=.985,bottom=.095,top=.94,wspace=.22,hspace=.24);save(fig,'Figure_4','monthly_series.csv; breakpoints.csv; archived traffic context',600)

def s12():
    configure();src=(ROOT/'work/plot_cleaned_ldsa.py').read_text();code=textwrap.dedent(src.split('    import importlib.util',1)[1].split('    # Re-render the archived regulatory',1)[0]);code=code.replace("s1.configure_style();","s1.COMPARISONS=(('Primary','10-500'),('Primary','10-800'),('10-500','10-800'));s1.configure_style();");code=code.replace('min(-2.5,extent.min()-2),max(54.5,extent.max()+2)','min(-1,extent.min()-1),np.ceil((extent.max()+1)/5)*5').replace('limits,10,101+j','limits,5,101+j').replace('min(-3.4,extent.min()-.5),max(6.8,extent.max()+.5)','np.floor((extent.min()-.2)*2)/2,np.ceil((extent.max()+.2)*2)/2').replace('limits,2,103+j','limits,.5,103+j');exec(compile(code,'s1_s2_original_templates','exec'),globals().copy())

def s3():
    configure();date=bp.loc[bp.Site.eq('LND_TR')&bp.Metric.eq('LDSA'),'Breakpoint'].iloc[0];src=(ROOT/'work/redraw_regulatory_traffic_context.py').read_text()
    src=src.replace('breakpoint = 2013 + 2 / 12',f'breakpoint = {date.year+(date.month-1)/12}').replace('Mar 2013',date.strftime('%b %Y')).replace('ax_b.axvline(2013,','ax_b.axvline(breakpoint,')
    src=src.split('base = OUTDIR',1)[0]
    src=src.replace('ax_c.legend(\n','ax_c.legend(\n').replace('    ncol=2,','    ncol=4,').replace('    columnspacing=1.7,','    columnspacing=1.1,')
    scope=globals().copy();exec(compile(src,'s3_template','exec'),scope);save(scope['fig'],'Figure_S3','archived regulatory/traffic data; updated breakpoint',600)

def s456():
    configure();path=Path('@INPUT_ROOT@/data_sources/plotting/breakpoint_plots/site_diagnostic.py');m=load_module(path);m.configure_style()
    node=next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='draw_site')
    code=ast.unparse(node).replace('fig.savefig(output_file, dpi=dpi, facecolor=\'white\')', '_clear_breakpoint_label(fig, ax)\n    fig.savefig(output_file, dpi=dpi, facecolor=\'white\')')
    m.__dict__['_clear_breakpoint_label']=clear_breakpoint_label;exec(code,m.__dict__)
    for site in sorted(sp.Site):
        w=wk[wk.Site.eq(site)].dropna(subset=['LDSA'])[['date','LDSA']].rename(columns={'date':'Date','LDSA':'Value'})
        dates=bp.loc[bp.Site.eq(site)&bp.Metric.eq('LDSA'),'Breakpoint'];date=dates.iloc[0] if len(dates) else pd.NaT
        rows=tr[tr.Site.eq(site)&tr.Metric.eq('LDSA')&tr.Phase.isin(['Before','After'])].rename(columns={'Phase':'Period'}).copy();rows['End_Time_Range']=date;rows['Start_Time_Range']=date
        target=DEST/'Figure_S6'/f"Figure_S6_{ {'PRA_2012-2016_SUB':'PRA_SUB_1','PRA_2017-2019_SUB':'PRA_SUB_2'}.get(site,site) }.png";target.parent.mkdir(exist_ok=True);m.draw_site(site,w,rows,target,300)
        print('SAVED',target.name,flush=True)
        if site in ['LEI_TR','ZUR_UB','MEL_RB','PAR_SUB','PRA_SUB']:shutil.copy2(target,DEST/f"Figure_S4_{site}.png")
    s5()

def clear_breakpoint_label(fig,ax):
    # Keep the source's frameless annotation, but end a coincident dashed line
    # immediately below the text instead of drawing through the date or slope.
    for artist in fig.findobj(mpl.text.Text):
        artist.set_text(artist.get_text().replace('PRA_2012-2016_SUB','PRA_SUB').replace('PRA_2017-2019_SUB','PRA_SUB').replace('μm²/cm³',r'μm² cm$^{\mathbf{-3}}$'))
    if 'Breakpoint: NA' in ax.texts[0].get_text():
        ax.texts[0].set_text('Breakpoint: NA\nFewer than 20 valid months')
    manifest.append({'figure':'standalone_'+fig.texts[0].get_text(),'text':[{'text':t.get_text(),'pt':t.get_fontsize(),'family':t.get_fontfamily()} for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]})
    fig.canvas.draw();renderer=fig.canvas.get_renderer();box=ax.texts[0].get_window_extent(renderer);ab=ax.get_window_extent(renderer)
    for line in ax.lines:
        if line.get_linestyle()=='--' and line.get_color()=='black':
            x=ax.transData.transform((mdates.date2num(line.get_xdata()[0]),0))[0]
            if box.x0-3<x<box.x1+3:line.set_ydata([0,max(.05,(box.y0-ab.y0-6)/ab.height)])

def s5():
    configure();b=bp[bp.Metric.eq('LDSA')].copy();b['Breakpoint_Year_Decimal']=b.Breakpoint.dt.year+(b.Breakpoint.dt.month-1)/12
    src=Path('@INPUT_ROOT@/data_sources/plotting/breakpoint_plots/breakpoint_distribution.py').read_text();body='type_palette ='+src.split('type_palette =',1)[1].split("fig.savefig(",1)[0]
    scope=globals().copy();scope['plot_data']=b.sort_values('Breakpoint_Year_Decimal');exec(body,scope)
    # Match the embedded supplement's visible right-hand year labels.
    scope['ax_b'].tick_params(axis='y',labelleft=True,left=True,length=5,width=1.2)
    for artist in scope['fig'].findobj(mpl.text.Text):artist.set_fontweight('bold')
    save(scope['fig'],'Figure_S5','breakpoints.csv',600)

def s7():
    configure();_,corr=pair_tables();order=['PM2.5','PM10','BC','NO2','SO2','NO','CO','O3'];labels=['PM$_{2.5}$','PM$_{10}$','BC','NO$_2$','SO$_2$','NO','CO','O$_3$']
    fig,ax=plt.subplots(figsize=(12,4.6));styles=[('Full available record','#1f77b4','o'),('Exact overlap ≥60 months','#ff7f0e','s'),('Exact overlap ≥36 months','#2ca02c','^')]
    for j,(strategy,c,m) in enumerate(styles):
        d=corr[corr.Strategy.eq(strategy)].set_index('Metric').reindex(order);ax.scatter(np.arange(8)+(j-1)*.21,d.Pearson_r,s=75,c=c,marker=m,label=strategy,zorder=3)
    ax.set_xticks(range(8),labels,fontsize=14);ax.set_ylim(-1,1.04);ax.set_yticks(np.arange(-1,1.01,.25));ax.tick_params(labelsize=13);ax.set_ylabel('Pearson r between LDSA and metric trend slopes',fontsize=14);ax.axhline(0,color='#91bde0',ls='--',lw=1);ax.grid(axis='y',color='#ddd',alpha=.8);ax.spines[['top','right']].set_visible(True);ax.spines[:].set_linewidth(1.5);ax.legend(loc='lower center',bbox_to_anchor=(.5,1),ncol=3,frameon=False,fontsize=12.5);save(fig,'Figure_S7','full-record and paired temporal sensitivity',300)

def s8():
    configure();src=(ROOT/'work/redraw_phase_copollutant_dumbbell.py').read_text();scope=globals().copy()
    exec(src.split('# Read the phase-specific',1)[0],scope)
    phase=tr[tr.Phase.isin(['Before','After'])];phase=phase[phase.Metric.isin(['LDSA','BC','PM10','PM2.5','NO2','NO','SO2','CO','O3'])]
    median=phase[phase.Metric.ne('LDSA')].groupby(['Metric','Phase']).Slope_Percent.median().unstack();paired=phase[phase.Metric.ne('LDSA')].merge(phase[phase.Metric.eq('LDSA')][['Site','Phase','Slope_Percent']],on=['Site','Phase'],suffixes=('','_LDSA')).dropna(subset=['Slope_Percent','Slope_Percent_LDSA']);paired['Same_Sign']=paired.Slope_Percent*paired.Slope_Percent_LDSA>0
    same=paired.groupby(['Metric','Phase']).Same_Sign.mean().unstack();scope.update(median_slope=median,same_sign=same)
    body='pollutant_order ='+src.split('pollutant_order =',1)[1].split('\nstem =',1)[0];exec(body,scope);save(scope['fig'],'Figure_S8','trends.csv',300)

def s9():
    configure();src=(ROOT/'work/redraw_fig_s9_size_mode_changes.py').read_text();prefix,body=src.split('metric_order =',1);scope=globals().copy();exec(prefix,scope)
    changes=pd.read_csv(ANA/'city_36month_changes.csv');mapping={'LDSA':'LDSA','Total PNC':'Ntotal','Nucleation PNC':'Nucleation','Aitken PNC':'Aitken','Accumulation PNC':'Accumulation'}
    scope['calculate_site_changes']=lambda site,breakpoint:{k:float(changes.loc[changes.Site.eq(site)&changes.Metric.eq(v),'change_percent'].iloc[0]) for k,v in mapping.items()}
    body='metric_order ='+body.split('\nbase = OUTPUT_DIR',1)[0]
    for site,old in [('LND_TR','2013-03-01'),('STO_TR','2013-10-01')]:
        date=bp.loc[bp.Site.eq(site)&bp.Metric.eq('LDSA'),'Breakpoint'].iloc[0];body=body.replace(old,date.strftime('%Y-%m-%d')).replace('Mar 2013' if site=='LND_TR' else 'Oct 2013',date.strftime('%b %Y'))
    body=body.replace('ax.set_xlim(-52, 3)',"ax.set_xlim(min(-52,np.nanmin(values)-8),max(3,np.nanmax(values)+12))")
    body=body.replace('    for xval, ypos in zip(values, y):', "    for xval, ypos in zip(values, y):\n        if not np.isfinite(xval):\n            ax.text(.05,ypos,'Not covered',transform=ax.get_yaxis_transform(),fontsize=12)\n            continue")
    exec(body,scope);save(scope['fig'],'Figure_S9','city_36month_changes.csv',600)

def s10():
    configure();path=Path('@INPUT_ROOT@/data_sources/predictive_models/LG/shap_importance.py');src=path.read_text();scope=globals().copy();exec(src.split('# Load SHAP comparison tables.',1)[0],scope)
    imp=pd.read_csv(ANA/'SHAP_importance.csv');dfs={}
    for typ in ['UB','TR']:
        d=imp[imp.Type.eq(typ)].pivot(index='Feature',columns='Phase',values='mean_abs_SHAP').reset_index().rename(columns={'Before':'SHAP_Before','After':'SHAP_After'});d.Feature=d.Feature.map(scope['format_feature_name']);dfs[typ]=d
    scope['dfs']=dfs;body='# Configure shared axes.'+src.split('# Configure shared axes.',1)[1].split('save_filename =',1)[0];exec(body,scope);save(scope['fig'],'Figure_S10','SHAP_importance.csv',300)

def s11():
    configure();path=Path('@INPUT_ROOT@/data_sources/predictive_models/LG/shap_dependence.py');src=path.read_text();scope=globals().copy();scope['mlines']=mpl.lines
    prefix=src.split('# Configure fonts and predictor labels.',1)[1].split('# Input file paths.',1)[0];prefix=prefix[prefix.index("plt.rcParams"):];exec(prefix,scope)
    meteorology={'TEMP':('Temperature',r'Temperature ($^\circ$C)'),
                 'RH':('Relative humidity','Relative humidity (%)'),
                 'WS':('Wind speed','Wind speed (m/s)'),
                 'WD':('Wind direction',r'Wind direction ($^\circ$)'),
                 'PRESS':('Atmospheric pressure','Atmospheric pressure (mbar)'),
                 'RAD':('Radiation',r'Radiation (W/m$^2$)'),
                 'RAIN':('Rainfall','Rainfall (mm)')}
    original_label=scope['format_predictor_name'];original_heading=scope['format_predictor_name_no_unit']
    scope['format_predictor_name']=lambda name:meteorology[name][1] if name in meteorology else original_label(name)
    scope['format_predictor_name_no_unit']=lambda name:meteorology[name][0] if name in meteorology else original_heading(name)
    plt.rcParams.update({'xtick.labelsize':14,'ytick.labelsize':14})
    data=[];order=['O3','NO','NO2','SO2','CO','BC','PM2.5','TEMP','RH','WS','WD','PRESS','RAD','RAIN']
    for typ in ['UB','TR']:
        for phase in ['Before','After']:
            d=pd.read_csv(ANA/f'heldout_SHAP_{typ}_{phase}.csv')
            for pred in order:
                q=d[[pred,pred+'_SHAP']].dropna().rename(columns={pred:'Predictor_Value',pred+'_SHAP':'SHAP_Value'});q['Predictor']=pred;q['Category']=typ;q['Period']=phase;data.append(q)
    scope['final_pdp_df']=pd.concat(data,ignore_index=True);scope['pred_order']=order
    body=src.split("    unique_preds =",1)[1].split('    # Predictor importance comparison.',1)[0]
    body=textwrap.dedent('    unique_preds ='+body);body=body.replace("unique_preds = final_pdp_df['Predictor'].unique()","unique_preds = pred_order")
    body=body.split("plt.savefig('Unified_Full_PDP_Combined.png'",1)[0]
    body=body.replace('scatter=False, order=3, color=color,','scatter=False, order=3, seed=20260909, color=color,')
    body=body.replace('fontsize=14,', 'fontsize=18,')
    body=body.replace('y=1.02','y=.992').replace('plt.tight_layout(pad=3.0)',"plt.subplots_adjust(left=.055,right=.992,bottom=.045,top=.95,wspace=.20,hspace=.22)")
    exec(body,scope)
    axes=scope['axes'];legend_pos=axes[14].get_position();last_pos=axes[15].get_position()
    axes[14].set_position([legend_pos.x0,legend_pos.y0,last_pos.x1-legend_pos.x0,legend_pos.height])
    save(scope['fig'],'Figure_S11','heldout_SHAP_{UB,TR}_{Before,After}.csv',300)

if __name__=='__main__':
    groups=sys.argv[1:] or ['main12','main3','main4','s12','s3','s456','s7','s8','s9','s10','s11']
    for name in groups:
        print('START',name,flush=True);globals()[name]()
    (WORK/('audit_'+'_'.join(groups)+'.json')).write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
