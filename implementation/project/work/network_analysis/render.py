from pathlib import Path
import sys,os,importlib.util,inspect,json,shutil
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;A=R/'outputs/network_analysis/analysis';D=Path('@DELIVERY_ROOT@/ancillary_analysis_figures');os.environ['MPLCONFIGDIR']=str(W/'mplconfig');sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import pandas as pd,numpy as np
from PIL import Image
from matplotlib.font_manager import fontManager,FontProperties
FONT='@FONT_ROOT@/Times New Roman Bold.ttf'
for n in ['Times New Roman.ttf','Times New Roman Bold.ttf','Times New Roman Italic.ttf','Times New Roman Bold Italic.ttf']:fontManager.addfont('@FONT_ROOT@/'+n)
spec=importlib.util.spec_from_file_location('s7_s11_template',R/'work/workbook_processing/redraw.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.ANA=A;m.WORK=W;m.DEST=D/'PNG';m.tr=pd.read_csv(A/'trends.csv');m.ms=pd.read_csv(A/'monthly_series.csv');m.ms['date']=pd.to_datetime(m.ms.date,format='mixed');m.bp=pd.read_csv(A/'breakpoints.csv');m.bp['Breakpoint']=pd.to_datetime(m.bp.Breakpoint,format='mixed');m.tr=m.tr[~m.tr.Site.isin(['PAR_SUB','PRA_SUB'])].copy();m.ms=m.ms[~m.ms.Site.isin(['PAR_SUB','PRA_SUB'])].copy()
original_save=m.save
checks=[]
def save(fig,name,*args,**kwargs):
 for t in fig.findobj(m.mpl.text.Text):t.set_fontproperties(FontProperties(fname=FONT,size=t.get_fontsize(),weight='bold'))
 fig.canvas.draw()
 if name=='Figure_S10':
  assert len(fig.axes)==6;assert sum('RB stations' in a.get_title() for a in fig.axes)==2
  limits=[]
  for ax in fig.axes:
   bars=[p for p in ax.patches if isinstance(p,m.mpl.patches.Rectangle)]
   low=min(p.get_y() for p in bars);high=max(p.get_y()+p.get_height() for p in bars)
   old=ax.get_ylim();new=(high+.06,low-.06) if ax.yaxis_inverted() else (low-.06,high+.06)
   ax.set_ylim(new);limits.append(dict(title=ax.get_title(),before=list(old),after=list(new),bar_extent=[low,high]))
  fig.canvas.draw()
  (W/'s10_compact/layout_check.json').write_text(json.dumps(dict(change='Within-panel y-axis margins only; original subplot positions retained',axes=limits),indent=2))
 if name=='Figure_S11':
  legend=fig.axes[14].get_legend();assert len(legend.get_texts())==12
  title=fig._suptitle.get_text();assert 'RB' in title
  fig.canvas.draw();box=legend.get_window_extent(fig.canvas.get_renderer());assert box.x0>=0 and box.x1<=fig.bbox.width
 checks.append(dict(figure=name,axes=len(fig.axes),texts=[t.get_text() for t in fig.findobj(m.mpl.text.Text) if t.get_visible() and t.get_text()]))
 original_save(fig,name,*args,**kwargs)
 im=Image.open(m.DEST/(name+'.png'));im.thumbnail((2100,2100));im.save(W/(name+'_preview.png'))
m.save=save
s=inspect.getsource(m.s10).replace("['UB','TR']", "['UB','TR','RB']")
s=s.replace("scope['dfs']=dfs;body=", "scope['dfs']=dfs;scope['categories']=['UB','TR','RB'];src=src.replace('nrows=2, ncols=2, figsize=(16, 12)', 'nrows=3, ncols=2, figsize=(16, 18)');body=")
exec(compile(s,'s10_three_classes','exec'),m.__dict__)
s=inspect.getsource(m.s11).replace("['UB','TR']", "['UB','TR','RB']")
inject='''    body=body.replace("('TR', 'After'): ('#d62728', 'TR - After')", "('TR', 'After'): ('#d62728', 'TR - After'), ('RB', 'Before'): ('#9467bd', 'RB - Before'), ('RB', 'After'): ('#8c564b', 'RB - After')")
    body=body.replace('custom_labels = [h.get_label() for h in custom_handles]', "custom_handles += [mlines.Line2D([], [], color='#9467bd', lw=4, label='Trend: RB - Before'), mlines.Line2D([], [], marker='o', color='w', markerfacecolor='#9467bd', markersize=14, alpha=.5, label='Scatter: RB - Before'), mlines.Line2D([], [], color='#8c564b', lw=4, label='Trend: RB - After'), mlines.Line2D([], [], marker='o', color='w', markerfacecolor='#8c564b', markersize=14, alpha=.5, label='Scatter: RB - After')]\\ncustom_labels = [h.get_label() for h in custom_handles]")
    body=body.replace("loc='center', ncol=2", "loc='center', ncol=3").replace('UB and TR stations', 'UB, TR and RB stations')
'''
s=s.replace('    exec(body,scope)',inject+'    exec(body,scope)');exec(compile(s,'s11_three_classes','exec'),m.__dict__)
for group in (sys.argv[1:] or ['s7','s8','s9','s10','s11']):
 print('START',group,flush=True);getattr(m,group)()
for name in ['exact_overlap_pairs','classification','exact_overlap_breakpoints_min12','correlation_temporal_sensitivity']:
 if (W/(name+'.csv')).exists():shutil.copy2(W/(name+'.csv'),A/(name+'.csv'))
(W/'render_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2));(W/'render_manifest.json').write_text(json.dumps(m.manifest,ensure_ascii=False,indent=2))
print('RENDER_DONE',flush=True)
