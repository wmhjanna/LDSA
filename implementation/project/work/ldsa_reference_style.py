"""Typography settings for the study figures."""
from pathlib import Path
import json
import matplotlib.pyplot as plt
from matplotlib.text import Text
ROOT=Path('@PROJECT_ROOT@')
FAMILY=['Times New Roman','DejaVu Serif']
C={'UB':'#1683C7','TR':'#F05A14','SUB':'#6755B5','RB':'#3C9A35'}
MARKERS={'UB':'o','TR':'s','SUB':'^','RB':'D'}
def configure():
    plt.rcParams.update({'font.family':FAMILY,'font.size':11,'font.weight':'bold','axes.labelweight':'bold','axes.titleweight':'bold','mathtext.fontset':'stix','pdf.fonttype':42,'svg.fonttype':'none'})
def apply_style(fig,name):
    if name in ['Figure_S1','Figure_S2']:return
    model=name=='Figure_4' or name.startswith(('Figure_S10','Figure_S11'))
    exact=name=='Exact_overlap_trend_comparison'
    weekly=name=='Figure_S4' or name.startswith('Figure_S6')
    context=name=='LND_STO_synchronized_trajectories'
    size=name=='Size_mode_pre_post_changes'
    phase=name=='Phase_copollutant_trends'
    tick,label,title,note=(12,14,16,11) if model else (10.5,13.5,17,11.5) if exact else (9.6,11.5,11,9) if weekly else (10.8,11.5,13.5,10.8) if context else (13,14,15,13) if size else (11,13,13.2,9) if phase else (11,11,13.2,9)
    weight='normal' if weekly or phase else 'bold'
    for text in fig.findobj(Text):text.set_fontfamily(FAMILY);text.set_fontweight(weight)
    for ax in fig.axes:
        if not ax.axison:continue  # Forest table typography is set explicitly.
        ax.xaxis.label.set_fontsize(label);ax.yaxis.label.set_fontsize(label)
        ax.xaxis.label.set_fontweight('bold');ax.yaxis.label.set_fontweight('bold')
        ax.tick_params(labelsize=tick,width=1.15 if exact else 1.1 if model else .8,length=4.5 if exact else 3.5)
        for tt in [ax.title,ax._left_title,ax._right_title]:tt.set_fontsize(title);tt.set_fontweight('bold')
        for tt in ax.texts:tt.set_fontsize((18 if exact else 16) if len(tt.get_text())==1 else note)
        for spine in ax.spines.values():spine.set_linewidth(2 if model and (name.startswith('Figure_S11') or ax.get_ylabel().startswith('SHAP')) else 1.25 if exact else 1.2 if context else 1.35 if size else .8)
        if model or exact:
            ax.grid(True,linestyle=':',color='#D9D9D9',alpha=.6,linewidth=1.5 if model else .8)
        leg=ax.get_legend()
        if leg:
            for t in leg.get_texts():t.set_fontsize(10 if model else 10.5 if exact else 9);t.set_fontweight(weight)
        if phase:
            ax.tick_params(axis='y',labelsize=12)
            for t in ax.get_yticklabels()+ax.get_xticklabels():t.set_fontweight('bold')
            if leg:
                for t in leg.get_texts():t.set_fontsize(11);t.set_fontweight('bold')
        if name=='Figure_1_c':ax.tick_params(axis='y',labelsize=9)
        if name=='Figure_S5' and ax==fig.axes[0]:ax.tick_params(axis='x',labelsize=8)
    if name=='Figure_2':
        fig.axes[0].title.set_fontsize(16);fig.axes[0].yaxis.label.set_fontsize(12)
        for ax in fig.axes[2:]:
            if ax.axison:ax.xaxis.label.set_fontsize(13)
    fig.canvas.draw()
    audit=[]
    for ax in fig.axes:
        audit.append({'xlabel_pt':ax.xaxis.label.get_fontsize(),'ylabel_pt':ax.yaxis.label.get_fontsize(),'title_pt':ax._left_title.get_fontsize(),'tick_pt':sorted({x.get_fontsize() for x in ax.get_xticklabels()+ax.get_yticklabels()}),'text_pt':sorted({t.get_fontsize() for t in ax.texts})})
    p=ROOT/'work/ldsa_restyle_reference'/f'{name}_style.json';p.write_text(json.dumps(audit,indent=2))
configure()
