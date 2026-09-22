from pathlib import Path
import sys,importlib.util,inspect,json,hashlib
R=Path('@PROJECT_ROOT@')
sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import pandas as pd,numpy as np
from matplotlib.font_manager import fontManager,FontProperties
FONT='@FONT_ROOT@/Times New Roman Bold.ttf'
for f in ['Times New Roman.ttf','Times New Roman Bold.ttf','Times New Roman Italic.ttf','Times New Roman Bold Italic.ttf']:
 fontManager.addfont('@FONT_ROOT@/'+f)
W=R/'work/network_trends';A=R/'outputs/network_trends/analysis'
O=Path('@DELIVERY_ROOT@/network_trend_figures');O.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('figure_template',R/'work/workbook_processing/redraw.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
configure=m.configure
def configure_bold():
 configure()
 pattern=FontProperties(fname=FONT).get_fontconfig_pattern()
 m.mpl.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':pattern,'mathtext.bf':pattern,'mathtext.it':pattern,'mathtext.default':'regular'})
m.configure=configure_bold
m.ANA=A;m.WORK=W;m.DEST=O
m.tr=pd.read_csv(A/'trends.csv');m.bp=pd.read_csv(A/'breakpoints.csv');m.bp['Breakpoint']=pd.to_datetime(m.bp.Breakpoint,format='mixed')
m.wk=pd.read_pickle(W/'weekly.pkl');m.meta=pd.read_csv(A/'meta_summary.csv');m.sites=pd.read_csv(A/'meta_sites.csv')
height=max(12.5,(len(m.sites)+15)*.285)
source=inspect.getsource(m.main12)
source=source.replace("    exec(compile(code,'main2_template','exec'),globals().copy())",f"    code=code.replace('figsize=(15.8,12.5)', 'figsize=(15.8,{height})')\n    exec(compile(code,'main2_template','exec'),globals().copy())")
exec(compile(source,'updated_main12','exec'),m.__dict__)
original_save=m.save
def save(fig,name,*args,**kwargs):
 if name!='Figure_2':m.plt.close(fig);return
 # PAL has a very wide published-method percentage CI. Show its estimate, keep the
 # full interval numerically, and use a right-pointing continuation arrow for the CI.
 forest=fig.axes[-1];forest.set_xlim(-12,32);forest.set_xticks([-10,0,10,20,30])
 for i,line in enumerate(forest.collections):
  if hasattr(line,'get_segments'):
   for seg in line.get_segments():
    if len(seg)==2 and seg[:,0].max()>32:
     y=seg[0,1];forest.annotate('',xy=(31.95,y),xytext=(30.5,y),arrowprops={'arrowstyle':'->','color':'#777777','lw':1.05})
 m.mpl.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':FontProperties(fname=FONT).get_fontconfig_pattern(),'mathtext.bf':FontProperties(fname=FONT).get_fontconfig_pattern(),'mathtext.it':FontProperties(fname=FONT).get_fontconfig_pattern(),'mathtext.default':'regular'})
 for a in fig.axes:
  for t in a.texts:
   t.set_text(t.get_text().replace('P = 0.000','P < 0.001'))
   if t.get_text()=='Symmetric log scale':
    t.set_position((.035,.83));t.set_ha('left')
 fig.canvas.draw()
 for t in fig.findobj(m.mpl.text.Text):t.set_fontproperties(FontProperties(fname=FONT,size=t.get_fontsize(),weight='bold'))
 assert all(t.get_fontproperties().get_file()==FONT for t in fig.findobj(m.mpl.text.Text) if t.get_visible() and t.get_text())
 forest_markers=[c for c in forest.collections if hasattr(c,'get_offsets') and c.__class__.__name__=='PathCollection']
 assert len(forest_markers)==len(m.sites)
 points=sorted(float(c.get_offsets()[0,0]) for c in forest_markers)
 assert np.allclose(points,sorted(m.sites.Slope_Percent))
 original_save(fig,name,*args,**kwargs)
m.save=save;m.main12()
(W/'render_audit.json').write_text(json.dumps(m.manifest,ensure_ascii=False,indent=2))
from PIL import Image
im=Image.open(O/'Figure_2.png');im.thumbnail((1700,2200));im.save(W/'preview.png')
print('Figure saved:',O/'Figure_2.png','height',height,'meta sites',len(m.sites))
