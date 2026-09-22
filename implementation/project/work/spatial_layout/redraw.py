"""Render the spatial comparison panels from the 65-station location workbook."""
from pathlib import Path
import json, hashlib, sys
sys.path.insert(0,"@PROJECT_ROOT@/work/ldsa_support/deps")
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties

ROOT=Path('@PROJECT_ROOT@')
W=ROOT/'work/spatial_layout'
BASE=Path('@DELIVERY_ROOT@')
SRC=BASE/'spatial_reference_tables/station_spatial_statistics.xlsx'
filtered='--covered' in sys.argv
OUT=BASE/('Figure1_49_station_coverage' if filtered else 'Figure1_65_station_coverage');OUT.mkdir(exist_ok=True)
if filtered:
    selection=json.loads((W/'coverage_selection.json').read_text())
    W=W/'covered';W.mkdir(exist_ok=True)
C={'UB':'#1683C7','TR':'#F05A14','SUB':'#6755B5','RB':'#3C9A35'}
ORDER=list(C)
tab=pd.read_excel(SRC,sheet_name='location')
assert len(tab)==65 and tab.Acronym.is_unique
original_types=tab.type.copy()
tab['type']=tab.type.replace({'U':'UB','F':'RB','M':'RB','P':'RB'})
assert set(tab.type)==set(ORDER)
assert tab.type.value_counts().to_dict()=={'RB':33,'UB':21,'TR':6,'SUB':5}
pd.DataFrame({'Site':tab.Acronym,'Source_type':original_types,'Figure_type':tab.type}).to_csv(W/'type_mapping.csv',index=False)
d=pd.DataFrame({'Site':tab.Acronym,'Type':tab.type,'Mean':tab.Total_Mean,'Median':tab.Total_Medi,'P95':tab.Total_P95})
d[['P25','P75']]=tab.Total_P25_.str.split('-',expand=True).astype(float)
assert d[['Mean','Median','P95','P25','P75']].notna().all().all()
assert (d.P25<=d.Median).all() and (d.Median<=d.P75).all() and (d.P75<=d.P95).all()
assert np.isclose(d.loc[d.Site.eq('DRE_TR'),'Mean'].iloc[0],51.339476792448544)
if filtered:
    d=d[d.Site.isin(selection['keep'])].copy()
    assert len(d)==len(selection['keep']) and set(d.Site)==set(selection['keep'])
if '--renamed' in sys.argv:
    rename_work='station_labels_49' if filtered else 'station_registry'
    mapping=json.loads((ROOT/'work'/rename_work/'mapping.json').read_text())
    names={r['old']:r['new'] for r in mapping}
    d['Site']=d.Site.replace(names)
    assert d.Site.is_unique
    W=ROOT/'work'/rename_work
    OUT=BASE/('Figure1_49_station_tables' if filtered else 'Figure1_65_station_tables');OUT.mkdir(exist_ok=True)
d.to_csv(W/'source_site_statistics.csv',index=False)
manifest=[]

def configure():
    mpl.rcdefaults()
    plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'font.weight':'bold',
        'axes.labelweight':'bold','axes.titleweight':'bold','mathtext.fontset':'stix','axes.unicode_minus':True})

def save(fig,name):
    for text in fig.findobj(mpl.text.Text):
        text.set_fontproperties(FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf',size=text.get_fontsize(),weight='bold'))
    fig.canvas.draw()
    labels=[t.get_text() for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
    assert not any(t in ['b','c'] or 'Primary LDSA' in t or 'physical monitoring sites' in t or '31.65' in t for t in labels)
    if name=='Figure_1_c':
        assert [t.get_text() for t in fig.axes[0].get_yticklabels()]==d.sort_values('Mean',ascending=False,kind='stable').Site.tolist()
        boxes=[t.get_window_extent(fig.canvas.get_renderer()) for t in fig.axes[0].get_yticklabels()]
        assert all(a.y0>b.y1 for a,b in zip(boxes,boxes[1:]))
    fig.savefig(OUT/(name+'.png'),dpi=400,bbox_inches='tight',facecolor='white')
    fig.savefig(W/(name+'_preview.png'),dpi=110,bbox_inches='tight',facecolor='white')
    manifest.append({'figure':name,'sites':len(d),'type_counts':d.Type.value_counts().to_dict(),'dpi':400,'labels':labels})
    plt.close(fig)

configure()
# Retain the prior four-category width, fonts, colours and marks.
fig,ax=plt.subplots(figsize=(6.5,4.7))
vals=[d.loc[d.Type.eq(t),'Mean'].to_numpy() for t in ORDER]
bp=ax.boxplot(vals,positions=np.arange(len(ORDER)),patch_artist=True,showfliers=False,widths=.5,
    medianprops={'color':'black','linewidth':1.5},whiskerprops={'color':'#666','linewidth':1.2},capprops={'color':'#666','linewidth':1.2})
rng=np.random.default_rng(20260909)
for i,(patch,typ,v) in enumerate(zip(bp['boxes'],ORDER,vals)):
    patch.set(facecolor=mpl.colors.to_rgba(C[typ],.22),edgecolor=C[typ],linewidth=1.5)
    ax.scatter(i+rng.uniform(-.07,.07,len(v)),v,s=35,color=C[typ],edgecolors='white',linewidths=.35,zorder=3)
assert sum(len(c.get_offsets()) for c in ax.collections)==len(d)
ax.set_xticks(range(len(ORDER)),ORDER,fontsize=12)
ax.set_xlabel('Station type',fontsize=12)
ax.set_ylabel(r'Site-mean total LDSA ($\mathbf{\mu m^{2}\,cm^{-3}}$)',fontsize=12)
ax.tick_params(labelsize=11);ax.grid(axis='y',ls='--',alpha=.55);ax.set_axisbelow(True)
save(fig,'Figure_1_b')

configure()
# Preserve the former physical row spacing and font sizes with 65 site rows.
fig,ax=plt.subplots(figsize=(6.1,18.5*len(d)/65));ranked=d.sort_values('Mean',ascending=False,kind='stable')
for i,r in enumerate(ranked.itertuples()):
    ax.hlines(i,r.P25,r.P75,color=C[r.Type],lw=1.5)
    ax.scatter(r.Mean,i,c=C[r.Type],s=33,edgecolors='white',linewidths=.6,zorder=3)
    ax.scatter(r.P95,i,facecolors='white',edgecolors=C[r.Type],s=31,linewidths=1.2,zorder=3)
ax.set_yticks(range(len(d)),ranked.Site,fontsize=9);ax.invert_yaxis()
ax.set_xlabel(r'Total LDSA ($\mathbf{\mu m^{2}\,cm^{-3}}$)',fontsize=11)
ax.set_title('Ranked site statistics for total LDSA',fontsize=13,pad=12)
ax.grid(axis='x',ls='--',alpha=.3);ax.set_axisbelow(True)
statleg=ax.legend(handles=[Line2D([],[],color='#555',lw=1.5,label='P25–P75'),
    Line2D([],[],color='black',marker='o',lw=0,label='Mean'),
    Line2D([],[],color='black',marker='o',mfc='white',lw=0,label='P95')],
    title='Statistic',loc='lower right',bbox_to_anchor=(1,.17),fontsize=9,title_fontsize=10,frameon=False)
ax.add_artist(statleg)
typeleg=ax.legend(handles=[Line2D([],[],color=C[t],marker='o',lw=0,label=t) for t in ORDER],
    title='Site type',loc='lower right',fontsize=9,title_fontsize=10,frameon=False)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
assert not statleg.get_window_extent(renderer).overlaps(typeleg.get_window_extent(renderer))
# Confirm legends occupy empty data space, including all P95 markers.
for leg in [statleg,typeleg]:
    box=leg.get_window_extent(renderer).expanded(1.06,1.03)
    for i,r in enumerate(ranked.itertuples()):
        for value in [r.P25,r.P75,r.Mean,r.P95]:
            assert not box.contains(*ax.transData.transform((value,i))),('Legend occlusion',r.Site)
save(fig,'Figure_1_c')
(W/'verification.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
(W/'provenance.json').write_text(json.dumps({'source':str(SRC),'sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),
    'scope':'Fig1 b/c only, statistics copied from updated 65-site workbook. No other analyses changed.',
    'classification':'Station-class mapping: U to UB and F/M/P to RB; retain all other source classes.',
    'period':'2015–2019 priority: 56 sites; nearest available year for 9 sites, as specified in the station table.',
    'native_coverage_filtered':filtered,'retained_sites':d.Site.tolist(),'colors':C},indent=2,ensure_ascii=False))
print('Saved',OUT)
