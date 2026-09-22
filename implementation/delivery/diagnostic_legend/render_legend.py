from pathlib import Path
import sys,json
sys.path.insert(0,'@PROJECT_ROOT@/work/ldsa_support/deps')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
D=Path(__file__).parent
font=FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf',size=18)
handles=[Line2D([],[],linestyle='none',marker='o',markerfacecolor='#95a5a6',markeredgecolor='#95a5a6',alpha=.5,markersize=7,label='Weekly Mean'),Line2D([],[],color='#e74c3c',lw=3,label='Pre-Breakpoint Trend'),Line2D([],[],color='#2980b9',lw=3,label='Post-Breakpoint Trend'),Line2D([],[],color='black',ls='--',lw=1.8,label='Breakpoint')]
fig=plt.figure(figsize=(8,1.5));legend=fig.legend(handles=[handles[0],handles[3],handles[1],handles[2]],loc='center',ncol=2,frameon=False,prop=font,handlelength=2.4,handletextpad=.65,columnspacing=1.6,borderpad=.15)
fig.canvas.draw();ren=fig.canvas.get_renderer();boxes=[t.get_window_extent(ren) for t in legend.get_texts()];assert not any(a.overlaps(b) for i,a in enumerate(boxes) for b in boxes[i+1:])
bbox=legend.get_window_extent(ren).transformed(fig.dpi_scale_trans.inverted()).expanded(1.02,1.25)
fig.savefig(D/'Figure_S6_legend.png',dpi=600,bbox_inches=bbox,transparent=True)
fig.savefig(D/'preview.png',dpi=120,bbox_inches=bbox,facecolor='white')
plt.close(fig)
(D/'legendnotes.json').write_text(json.dumps({'font':'Times New Roman Bold','dpi':600,'layout':'two columns, two rows','background':'transparent','source':'@INPUT_ROOT@/data_sources/plotting/breakpoint_plots/site_diagnostic.py','items':{'Weekly Mean':'Grey circles #95a5a6, alpha 0.5; weekly observations','Pre-Breakpoint Trend':'Red line #e74c3c','Post-Breakpoint Trend':'Blue line #2980b9','Breakpoint':'Black dashed line'}},ensure_ascii=False,indent=2))
print(D/'Figure_S6_legend.png')
