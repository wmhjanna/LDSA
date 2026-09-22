from pathlib import Path
import json,sys
sys.path.insert(0,'@PROJECT_ROOT@/work/ldsa_support/deps')
import numpy as np,pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
ROOT=Path('@PROJECT_ROOT@')
W=ROOT/'work/quality_control/spatial'
B=Path('@DELIVERY_ROOT@/provider_quality_control_analysis')
all_data=pd.read_csv(ROOT/'outputs/quality_control/analysis/site_spatial_statistics.csv')
keep=json.loads((W/'coverage_selection.json').read_text())['keep']
C={'UB':'#1683C7','TR':'#F05A14','SUB':'#6755B5','RB':'#3C9A35'};ORDER=list(C)
manifest=[]
def configure():
 mpl.rcdefaults();plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'font.weight':'bold','axes.labelweight':'bold','axes.titleweight':'bold','mathtext.fontset':'stix','axes.unicode_minus':True})
def save(fig,name,d,output):
 for text in fig.findobj(mpl.text.Text):text.set_fontproperties(FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf',size=text.get_fontsize(),weight='bold'))
 fig.canvas.draw();labels=[t.get_text() for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
 assert not any(t in ['b','c'] or 'Primary LDSA' in t or 'physical monitoring sites' in t or '31.65' in t for t in labels)
 if name=='Figure_1_c':
  assert [t.get_text() for t in fig.axes[0].get_yticklabels()]==d.sort_values('Mean',ascending=False,kind='stable').Site.tolist()
  boxes=[t.get_window_extent(fig.canvas.get_renderer()) for t in fig.axes[0].get_yticklabels()]
  assert all(a.y0>b.y1 for a,b in zip(boxes,boxes[1:]))
 fig.savefig(output/(name+'.png'),dpi=400,bbox_inches='tight',facecolor='white')
 fig.savefig(W/(name+f'_{len(d)}_preview.png'),dpi=110,bbox_inches='tight',facecolor='white')
 manifest.append(dict(figure=name,sites=len(d),type_counts=d.Type.value_counts().to_dict(),dpi=400,labels=labels));plt.close(fig)
for count in [65,49]:
 d=all_data.copy() if count==65 else all_data[all_data.Site.isin(keep)].copy();assert len(d)==count
 output=B/'PNG' if count==65 else B/'diameter_coverage'/'Fig1_49_stations';output.mkdir(parents=True,exist_ok=True)
 configure();fig,ax=plt.subplots(figsize=(6.5,4.7));vals=[d.loc[d.Type.eq(t),'Mean'].to_numpy() for t in ORDER]
 bp=ax.boxplot(vals,positions=np.arange(len(ORDER)),patch_artist=True,showfliers=False,widths=.5,medianprops={'color':'black','linewidth':1.5},whiskerprops={'color':'#666','linewidth':1.2},capprops={'color':'#666','linewidth':1.2})
 rng=np.random.default_rng(20260909)
 for i,(patch,typ,v) in enumerate(zip(bp['boxes'],ORDER,vals)):
  patch.set(facecolor=mpl.colors.to_rgba(C[typ],.22),edgecolor=C[typ],linewidth=1.5)
  ax.scatter(i+rng.uniform(-.07,.07,len(v)),v,s=35,color=C[typ],edgecolors='white',linewidths=.35,zorder=3)
 assert sum(len(c.get_offsets()) for c in ax.collections)==len(d)
 ax.set_xticks(range(len(ORDER)),ORDER,fontsize=12);ax.set_xlabel('Station type',fontsize=12)
 ax.set_ylabel(r'Site-mean total LDSA ($\mathbf{\mu m^{2}\,cm^{-3}}$)',fontsize=12)
 ax.tick_params(labelsize=11);ax.grid(axis='y',ls='--',alpha=.55);ax.set_axisbelow(True);save(fig,'Figure_1_b',d,output)
 configure();fig,ax=plt.subplots(figsize=(6.1,18.5*len(d)/65));ranked=d.sort_values('Mean',ascending=False,kind='stable')
 for i,r in enumerate(ranked.itertuples()):
  ax.hlines(i,r.P25,r.P75,color=C[r.Type],lw=1.5)
  ax.scatter(r.Mean,i,c=C[r.Type],s=33,edgecolors='white',linewidths=.6,zorder=3)
  ax.scatter(r.P95,i,facecolors='white',edgecolors=C[r.Type],s=31,linewidths=1.2,zorder=3)
 ax.set_yticks(range(len(d)),ranked.Site,fontsize=9);ax.invert_yaxis()
 ax.set_xlabel(r'Total LDSA ($\mathbf{\mu m^{2}\,cm^{-3}}$)',fontsize=11);ax.set_title('Ranked site statistics for total LDSA',fontsize=13,pad=12)
 ax.grid(axis='x',ls='--',alpha=.3);ax.set_axisbelow(True)
 statleg=ax.legend(handles=[Line2D([],[],color='#555',lw=1.5,label='P25–P75'),Line2D([],[],color='black',marker='o',lw=0,label='Mean'),Line2D([],[],color='black',marker='o',mfc='white',lw=0,label='P95')],title='Statistic',loc='lower right',bbox_to_anchor=(1,.17),fontsize=9,title_fontsize=10,frameon=False);ax.add_artist(statleg)
 typeleg=ax.legend(handles=[Line2D([],[],color=C[t],marker='o',lw=0,label=t) for t in ORDER],title='Site type',loc='lower right',fontsize=9,title_fontsize=10,frameon=False)
 fig.canvas.draw();renderer=fig.canvas.get_renderer();assert not statleg.get_window_extent(renderer).overlaps(typeleg.get_window_extent(renderer))
 for leg in [statleg,typeleg]:
  box=leg.get_window_extent(renderer).expanded(1.06,1.03)
  for i,r in enumerate(ranked.itertuples()):
   for value in [r.P25,r.P75,r.Mean,r.P95]:assert not box.contains(*ax.transData.transform((value,i))),('Legend occlusion',r.Site)
 save(fig,'Figure_1_c',d,output)
(W/'figure_verification.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
print('FIG1_READY')
