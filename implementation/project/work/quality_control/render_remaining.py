import os
os.environ['MPLCONFIGDIR']='@PROJECT_ROOT@/work/quality_control/mplconfig'
from pathlib import Path
import sys,importlib.util,json,shutil
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;A=R/'outputs/quality_control/workbook/analysis'
sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
from matplotlib.font_manager import fontManager,FontProperties
for n in ['Times New Roman.ttf','Times New Roman Bold.ttf','Times New Roman Italic.ttf','Times New Roman Bold Italic.ttf']:fontManager.addfont('@FONT_ROOT@/'+n)
spec=importlib.util.spec_from_file_location('updated',W/'workbook/redraw.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.DEST=Path('@DELIVERY_ROOT@/provider_quality_control_analysis/PNG');m.WORK=W/'workbook';m.DEST.mkdir(exist_ok=True)
orig=m.save
def save(fig,name,*args,**kwargs):
 for t in fig.findobj(m.mpl.text.Text):t.set_fontproperties(FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf',size=t.get_fontsize(),weight='bold'))
 orig(fig,name,*args,**kwargs)
m.save=save
for group in ['s12','main4','s9','s10','s11']:
 print('START',group,flush=True);getattr(m,group)()
m.tr=m.tr[~m.tr.Site.isin(['PAR_SUB','PRA_SUB'])].copy();m.ms=m.ms[~m.ms.Site.isin(['PAR_SUB','PRA_SUB'])].copy()
for group in ['main3','s7','s8']:
 print('START',group,flush=True);getattr(m,group)()
for name in ['exact_overlap_pairs','classification','exact_overlap_breakpoints_min12','correlation_temporal_sensitivity']:
 shutil.copy2(m.WORK/(name+'.csv'),A/(name+'.csv'));shutil.copy2(m.WORK/(name+'.csv'),R/'outputs/quality_control/analysis'/(name+'.csv'))
(W/'render_workbook_audit.json').write_text(json.dumps(m.manifest,ensure_ascii=False,indent=2))
