"""Recompose Fig2 from its exact current template; export S5a independently."""
from pathlib import Path
import json,hashlib,os,sys
from PIL import Image
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent
D=Path('@DELIVERY_ROOT@/record_exclusion_analysis/PNG');D.mkdir(exist_ok=True)
A=R/'outputs/record_exclusions/analysis'
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in A.iterdir() if p.is_file()}
source=(R/'work/record_exclusions/render_fig2.py').read_text()
source=source.replace("O=Path('@DELIVERY_ROOT@/provider_quality_control_analysis/PNG')",f'O=Path({str(D)!r})')
source=source.replace("im.save(W/'preview.png')",f"im.save(Path({str(W/'Figure_2_preview.png')!r}))")
source=source.replace("(W/'render_fig2_audit.json')",f"Path({str(W/'render_fig2_audit.json')!r})")
injection='''
def draw_distribution(fig, spec):
 ax=fig.add_subplot(spec)
 b=m.bp[m.bp.Metric.eq('LDSA')].copy()
 b['Year']=b.Breakpoint.dt.year+(b.Breakpoint.dt.month-1)/12
 b=b.sort_values('Year')
 palette={'UB':'#2ca02c','TR':'#d62728','SUB':'#1f77b4','RB':'#ff7f0e'}
 order=list(palette)
 m.sns.boxplot(ax=ax,data=b,x='Type',y='Year',order=order,hue='Type',hue_order=order,palette=palette,width=.5,boxprops={'alpha':.45},showfliers=False,legend=False)
 m.sns.swarmplot(ax=ax,data=b,x='Type',y='Year',order=order,hue='Type',hue_order=order,palette=palette,size=5,edgecolor='none',linewidth=0,legend=False)
 ax.set_xlabel('Station class');ax.set_ylabel('Breakpoint year')
 ax.yaxis.set_major_locator(m.mpl.ticker.MultipleLocator(4));ax.yaxis.set_major_formatter(m.mpl.ticker.FormatStrFormatter('%d'))
 m.head(ax,'b','Breakpoint year by station class')
 ax.grid(axis='y',linestyle='-',linewidth=.8,alpha=.3);ax.grid(axis='x',visible=False);ax.set_axisbelow(True)
 for t in ax.get_xticklabels()+ax.get_yticklabels():t.set_fontweight('bold')
 assert sum(len(c.get_offsets()) for c in ax.collections if c.__class__.__name__=='PathCollection')==len(b)==61
 return ax
m.draw_distribution=draw_distribution
source=source.replace("    exec(compile(code,'main2_template','exec'),globals().copy())", "    code=code.replace('fig.add_gridspec(2,2,', 'fig.add_gridspec(3,2,').replace('hspace=.30', 'hspace=.30').replace('ax=fig.add_subplot(gs[1,0]);groups=[]', 'ax=fig.add_subplot(gs[2,0]);groups=[]').replace(\\\"head(ax,'b','Pre- and post-breakpoint LDSA slopes')\\\", \\\"head(ax,'c','Pre- and post-breakpoint LDSA slopes')\\\").replace(\\\"study.text(-.035,1.025,'c',\\\", \\\"study.text(-.035,1.025,'d',\\\").replace(\\\"ax.set_xlabel('Station type')\\\", \\\"ax.set_xlabel('Station class')\\\").replace('forest(fig,gs[:,1]);save', 'draw_distribution(fig,gs[1,0]);forest(fig,gs[:,1]);save')\\n    exec(compile(code,'main2_template','exec'),globals().copy())")
'''
source=source.replace("exec(compile(source,'updated_main12','exec'),m.__dict__)",injection+"\nexec(compile(source,'updated_main12','exec'),m.__dict__)")
exec((Path(__file__).parent/'align_layout.py').read_text(),globals())
source=source.replace(' original_save(fig,name,*args,**kwargs)', ' align_fig2_layout(fig)\n original_save(fig,name,*args,**kwargs)')
if '--s5-only' in sys.argv:source=source.split('m.save=save;m.main12()',1)[0]
exec(compile(source,'current_fig2_relayout','exec'),globals())
if '--fig2-only' in sys.argv:
 assert hashes=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in A.iterdir() if p.is_file()}
 sys.exit(0)

# Independent S5a uses the same dates, colours, ordering and tight station bounds.
m.configure();b=m.bp[m.bp.Metric.eq('LDSA')].copy();b['Year']=b.Breakpoint.dt.year+(b.Breakpoint.dt.month-1)/12;b=b.sort_values('Year')
palette={'UB':'#2ca02c','TR':'#d62728','SUB':'#1f77b4','RB':'#ff7f0e'}
fig,ax=m.plt.subplots(figsize=(14,5.8))
ax.scatter(range(len(b)),b.Year,c=b.Type.map(palette),s=85,edgecolors='none',linewidths=0,zorder=3)
ax.set_xticks(range(len(b)),b.Site,rotation=90,ha='center',fontsize=10);ax.set_xlim(-.6,len(b)-.4)
ax.set_title('Station-specific breakpoint year',fontsize=15,fontweight='bold',pad=12)
ax.set_xlabel('Station',fontsize=13);ax.set_ylabel('Breakpoint year',fontsize=13)
ax.legend(handles=[m.Line2D([0],[0],marker='o',linestyle='',markerfacecolor=c,markeredgecolor='none',markersize=9,label=t) for t,c in palette.items()],title='Station class',loc='upper left',frameon=False)
ax.yaxis.set_major_locator(m.mpl.ticker.MultipleLocator(2));ax.yaxis.set_major_formatter(m.mpl.ticker.FormatStrFormatter('%d'))
ax.tick_params(axis='x',length=5,width=1);ax.tick_params(axis='y',length=6,width=1.2,labelsize=11)
ax.grid(axis='y',linestyle='-',linewidth=.8,alpha=.35);ax.grid(axis='x',visible=False);ax.set_axisbelow(True)
for spine in ax.spines.values():spine.set_color('#555555');spine.set_linewidth(.75)
fig.subplots_adjust(left=.065,right=.99,bottom=.24,top=.91)
for t in fig.findobj(m.mpl.text.Text):t.set_fontproperties(FontProperties(fname=FONT,size=t.get_fontsize(),weight='bold'))
fig.savefig(D/'Figure_S5_a.png',dpi=600,bbox_inches='tight',facecolor='white');m.plt.close(fig)
im=Image.open(D/'Figure_S5_a.png');im.thumbnail((2200,1000));im.save(Path(__file__).parent/'Figure_S5_a_preview.png')
assert hashes=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in A.iterdir() if p.is_file()}
(Path(__file__).parent/'verification.json').write_text(json.dumps(dict(analysis_hashes_unchanged=True,breakpoint_points=len(b),forest_stations=len(m.sites),Figure2_panels={'a':'LND_TR example','b':'Breakpoint year distribution by station class','c':'Before/After slopes','d':'Random-effects forest'},S5a_sites=b.Site.tolist()),indent=2))
print('COMPLETE',D,flush=True)
