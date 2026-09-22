import os
from pathlib import Path
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;A=R/'outputs/netcdf_pollutants/analysis'
os.environ['MPLCONFIGDIR']=str(W/'mplconfig')
import sys,importlib.util,json,shutil
sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import pandas as pd
from matplotlib.font_manager import fontManager,FontProperties
for n in ['Times New Roman.ttf','Times New Roman Bold.ttf','Times New Roman Italic.ttf','Times New Roman Bold Italic.ttf']:fontManager.addfont('@FONT_ROOT@/'+n)
spec=importlib.util.spec_from_file_location('new_pollutants_fig3',R/'work/workbook_processing/redraw.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.ANA=A;m.WORK=W;m.tr=pd.read_csv(A/'trends.csv');m.ms=pd.read_csv(A/'monthly_series.csv',parse_dates=['date']);m.bp=pd.read_csv(A/'breakpoints.csv',parse_dates=['Breakpoint'])
m.ms['date']=pd.to_datetime(m.ms.date,format='mixed');m.bp['Breakpoint']=pd.to_datetime(m.bp.Breakpoint,format='mixed')
m.tr=m.tr[~m.tr.Site.isin(['PAR_SUB','PRA_SUB'])].copy();m.ms=m.ms[~m.ms.Site.isin(['PAR_SUB','PRA_SUB'])].copy()
m.DEST=Path('@DELIVERY_ROOT@/netcdf_pollutant_figures');m.DEST.mkdir(exist_ok=True)
orig=m.save
def save(fig,name,*args,**kwargs):
 import numpy as np
 ax=fig.axes[5]
 offsets=pd.read_csv(W/'exact_overlap_breakpoints_min12.csv')
 vals=offsets.loc[offsets.Pollutant.isin(['PM2.5','PM10','BC','NO2','SO2']),'Offset Δt (months) [Y]']
 bottom=min(-85,float(vals.min())-10);top=max(100,float(vals.max())+50)
 ax.set_ylim(bottom,top)
 ax.set_yticks(np.arange(np.ceil(bottom/50)*50,np.floor(top/50)*50+1,50))
 for t in ax.texts:
  if t.get_text().startswith(('n =','med |Δt|')):
   t.set_transform(ax.get_xaxis_transform());t.set_position((t.get_position()[0],.98 if t.get_text().startswith('n =') else .925))
 for t in fig.findobj(m.mpl.text.Text):t.set_fontproperties(FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf',size=t.get_fontsize(),weight='bold'))
 orig(fig,name,*args,**kwargs)
 # Export small image for visual verification, not a replacement for the 500-dpi PNG.
 from PIL import Image
 p=m.DEST/(name+'.png');im=Image.open(p);im.thumbnail((2200,2200));im.save(W/(name+'_preview.png'))
m.save=save
m.main3()
for name in ['exact_overlap_pairs','classification','exact_overlap_breakpoints_min12']:
 shutil.copy2(W/(name+'.csv'),A/(name+'.csv'))
(W/'render_manifest.json').write_text(json.dumps(m.manifest,ensure_ascii=False,indent=2))
print('FIG3_READY',flush=True)
