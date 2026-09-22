from pathlib import Path
import sys,json
sys.path.insert(0,'@PROJECT_ROOT@/work/ldsa_support/deps')
import numpy as np,pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.font_manager import fontManager, FontProperties

BASE=Path('@DELIVERY_ROOT@')
O=BASE/'record_exclusion_analysis/diameter_coverage'
W=Path('@PROJECT_ROOT@/work/record_exclusions/spatial')
s=pd.read_csv(O/'class_summary.csv');overall=pd.read_csv(O/'network_summary.csv')
for f in ['Times New Roman.ttf','Times New Roman Bold.ttf']:fontManager.addfont('@FONT_ROOT@/'+f)
mpl.rcParams.update({'font.family':'Times New Roman','font.size':7,'font.weight':'bold','axes.labelsize':7,'axes.labelweight':'bold',
 'axes.titlesize':7,'axes.titleweight':'bold','axes.linewidth':.65,'xtick.labelsize':7,'ytick.labelsize':7,
 'xtick.major.width':.6,'ytick.major.width':.6,'xtick.major.size':2.5,'ytick.major.size':0,
 'mathtext.fontset':'custom','mathtext.rm':FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf').get_fontconfig_pattern(),'mathtext.bf':FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf').get_fontconfig_pattern(),'mathtext.it':FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf').get_fontconfig_pattern(),'mathtext.default':'regular','savefig.facecolor':'white','figure.facecolor':'white'})
C={'UB':'#1683C7','TR':'#F05A14','SUB':'#6755B5','RB':'#3C9A35','All sites':'#363636'}
types=['UB','TR','SUB','RB']
schemes=[('Main_65','All sites (65)','o'),('Upper_complete_53','Upper coverage (53)','s'),('Strict_49','Strict coverage (49)','^')]
fig=plt.figure(figsize=(183/25.4,96/25.4))
gs=fig.add_gridspec(1,2,left=.10,right=.985,bottom=.17,top=.79,wspace=.23,width_ratios=[1.12,1])
a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1])
fig.legend(handles=[Line2D([],[],marker=mk,markersize=4.5,lw=0,color='#333333',markerfacecolor='#333333' if j==0 else 'white',markeredgewidth=.9,label=label) for j,(_,label,mk) in enumerate(schemes)],
 loc='upper center',bbox_to_anchor=(.53,.98),ncol=3,frameon=False,columnspacing=2.2,handletextpad=.6,fontsize=7)
checks=[]
for i,t in enumerate(types):
    for j,(name,label,marker) in enumerate(schemes):
        r=s[(s.Type==t)&(s.scenario==name)].iloc[0];y=i+(j-1)*.23
        a.plot([r.P25,r.P75],[y,y],color=C[t],lw=1.3,solid_capstyle='round',alpha=.72,zorder=2)
        a.plot([r.P25,r.P25],[y-.045,y+.045],color=C[t],lw=.65)
        a.plot([r.P75,r.P75],[y-.045,y+.045],color=C[t],lw=.65)
        a.plot(r['mean'],y,marker=marker,ms=4.4,mec=C[t],mfc=C[t] if j==0 else 'white',mew=.95,lw=0,zorder=3)
        checks.append(dict(panel='a',Type=t,scenario=name,n=int(r['n']),mean=float(r['mean']),P25=float(r.P25),P75=float(r.P75)))
a.set_yticks(range(4),types);a.set_ylim(3.62,-.62);a.set_xlim(0,80);a.set_xticks([0,20,40,60,80])
a.set_xlabel(r'Site-mean LDSA ($\mathrm{\mu m^{2}\,cm^{-3}}$)',labelpad=7)
a.set_title('LDSA by site type',loc='left',pad=14)
a.set_ylabel('Station type',labelpad=8)

labels=['All sites']+types;positions=[0,1.3,2.3,3.3,4.3]
for t,y in zip(labels,positions):
    values=[]
    for j,(name,label,marker) in enumerate(schemes[1:]):
        r=overall[overall.scenario==name].iloc[0] if t=='All sites' else s[(s.Type==t)&(s.scenario==name)].iloc[0]
        change=float(r.mean_change_percent);yy=y+(j-.5)*.34;values.append((change,yy))
        b.plot(change,yy,marker=marker,ms=4.5,mec=C[t],mfc='white',mew=1,lw=0,zorder=3)
        # Labels on separate rows avoid implying uncertainty or significance.
        b.annotate(f'{change:+.2f}%',(change,yy),xytext=(5,0),textcoords='offset points',va='center',ha='left',fontsize=6.5,color='#222222',bbox=dict(facecolor='white',edgecolor='none',pad=.5))
        checks.append(dict(panel='b',Type=t,scenario=name,change_percent=change))
    b.plot([v[0] for v in values],[v[1] for v in values],color=C[t],alpha=.30,lw=.9,zorder=1)
b.axvline(0,color='#777777',lw=.8,ls=(0,(3,3)),zorder=0)
b.axhline(.64,color='#D7D7D7',lw=.65)
b.set_yticks(positions,labels);b.set_ylim(4.95,-.65);b.set_xlim(-4.5,12.5);b.set_xticks([-4,0,4,8,12])
b.set_xlabel('Change in mean LDSA (%)',labelpad=7)
b.set_title('Change relative to all sites',loc='left',pad=14)
for ax,letter in [(a,'a'),(b,'b')]:
    for side in ['top','right']:ax.spines[side].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_color('#333333')
    ax.tick_params(axis='y',pad=5)
    ax.set_axisbelow(True)
    ax.text(-.18,1.105,letter,transform=ax.transAxes,fontsize=9,fontweight='bold',ha='left',va='bottom')
a.grid(axis='x',color='#E8E8E8',lw=.5)

fig.canvas.draw()
for t in fig.findobj(mpl.text.Text):
    t.set_fontproperties(FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf',size=t.get_fontsize(),weight='bold'))
fig.canvas.draw();ren=fig.canvas.get_renderer()
# Ensure all text uses the explicit bold font and lies inside the canvas.
assert all(t.get_fontproperties().get_file()=='@FONT_ROOT@/Times New Roman Bold.ttf' for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text())
for t in fig.findobj(mpl.text.Text):
    if not t.get_visible() or not t.get_text():continue
    box=t.get_window_extent(ren)
    assert box.x0>=-1 and box.x1<=fig.bbox.width+1 and box.y0>=-1 and box.y1<=fig.bbox.height+1,(t.get_text(),box.bounds)
fig.savefig(O/'Figure_size_coverage_sensitivity.png',dpi=600)
fig.savefig(W/'preview.png',dpi=180)
(W/'plotted_values.json').write_text(json.dumps(checks,indent=2))
pd.DataFrame(checks).to_csv(O/'plot_data.csv',index=False)
(W/'figure_specifications.json').write_text(json.dumps({'width_mm':183,'height_mm':96,'dpi':600,'font':'Times New Roman','font_weight':'bold','panel_wspace':.23,'font_size_points':[6.5,7,9],'output':'PNG only','panel_a_intervals':'between-site interquartile range, not confidence intervals','panel_a_markers':'arithmetic mean of site means','data_points':len(checks),'guideline':'https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/'},indent=2))
print('SENSITIVITY_FIGURE_READY')
