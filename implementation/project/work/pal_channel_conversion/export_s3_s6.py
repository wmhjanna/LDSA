from pathlib import Path
import sys,importlib.util,inspect,json,shutil,zipfile,hashlib
R=Path('@PROJECT_ROOT@');sys.path[:0]=[str(R/'work/ldsa_support/deps'),str(R/'work')]
import pandas as pd
from matplotlib.font_manager import fontManager,FontProperties
from PIL import Image,ImageOps,ImageDraw
W=Path(__file__).parent/'s3_s6';W.mkdir(exist_ok=True)
A=R/'outputs/pal_channel_conversion/analysis'
O=Path('@DELIVERY_ROOT@/PAL_station_diagnostics');O.mkdir(exist_ok=True)
font='@FONT_ROOT@/Times New Roman Bold.ttf'
for name in ['Times New Roman.ttf','Times New Roman Bold.ttf','Times New Roman Italic.ttf','Times New Roman Bold Italic.ttf']:fontManager.addfont('@FONT_ROOT@/'+name)
spec=importlib.util.spec_from_file_location('template',R/'work/workbook_processing/redraw.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.ANA=A;m.WORK=W.parent;m.DEST=O
m.tr=pd.read_csv(A/'trends.csv');m.bp=pd.read_csv(A/'breakpoints.csv');m.bp['Breakpoint']=pd.to_datetime(m.bp.Breakpoint,format='mixed');m.wk=pd.read_pickle(W.parent/'weekly.pkl');m.sp=pd.read_csv(A/'site_eligibility.csv')
assert len(m.sp)==65 and m.sp.Site.is_unique and len(m.bp)==61
oldclear=m.clear_breakpoint_label
def clear(fig,ax):
 for t in fig.findobj(m.mpl.text.Text):t.set_fontproperties(FontProperties(fname=font,size=t.get_fontsize(),weight='bold'))
 oldclear(fig,ax)
m.clear_breakpoint_label=clear
src=inspect.getsource(m.s5)
src=src.replace("    scope=globals().copy();scope['plot_data']", "    body=body.replace('figsize=(15, 7)', 'figsize=(20, 7)').replace(\"[2.4, 1.6]\", \"[3.5, 1.6]\").replace('labelrotation=45', 'labelrotation=90').replace(\"ha='right'\", \"ha='center'\")\n    scope=globals().copy();scope['plot_data']")
exec(compile(src,'s5_65site_template','exec'),m.__dict__)
original_save=m.save
def save(fig,name,*args,**kwargs):
 if name=='Figure_S5':
  fig.axes[0].set_xlim(-.6,len(m.bp)-.4)
  for ax in fig.axes:
   ax.yaxis.set_major_locator(m.mpl.ticker.MultipleLocator(2))
   ax.yaxis.set_major_formatter(m.mpl.ticker.FormatStrFormatter('%d'))
  for t in fig.findobj(m.mpl.text.Text):
   t.set_text(t.get_text().replace('Distribution by Station Type','Distribution by Station Class').replace('Station Type','Station class'))
 for t in fig.findobj(m.mpl.text.Text):t.set_fontproperties(FontProperties(fname=font,size=t.get_fontsize(),weight='bold'))
 original_save(fig,name,*args,**kwargs)
m.save=save
if '--s5-only' in sys.argv:m.s5()
else:m.s3();m.s456()
files=sorted(O.rglob('*.png'));s6=sorted((O/'Figure_S6').glob('*.png'));s4=sorted(O.glob('Figure_S4*.png'))
assert len(s6)==65 and len(s4)==5 and len(files)==72
assert {p.stem.removeprefix('Figure_S6_') for p in s6}==set(m.sp.Site)
for item in m.manifest:
 for t in item.get('text',[]):
  assert '19.6–500' not in t['text'] and '31.65–500' not in t['text'] and 'μm²/cm³' not in t['text']
for p in files:
 im=Image.open(p);im.verify()
for p in [O/'Figure_S3.png',O/'Figure_S5.png',*s4]:
 im=Image.open(p);im.thumbnail((1800,1300));im.save(W/(p.stem+'_preview.png'))
for start in range(0,len(s6),12):
 page=Image.new('RGB',(1560,4*350),'white')
 for n,p in enumerate(s6[start:start+12]):
  im=Image.open(p);im.thumbnail((510,340));page.paste(im,((n%3)*520,(n//3)*350))
 page.save(W/f'contact_{start//12+1}.jpg')
report={'analysis':str(A),'PAL_assumption':'2000–2006 Source_ID13 uses the study assumption of per-bin DeltaN; the provider unit definition is unconfirmed','S3':'existing external context with unchanged LND breakpoint','S4_sites':[p.stem.removeprefix('Figure_S4_') for p in s4],'S5_breakpoints':len(m.bp),'S6_sites':len(s6),'PNG_count':len(files),'files':[{ 'file':str(p.relative_to(O)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]}
(O/'output_notes.md').write_text('For the study protocol and reproduction instructions, see docs/METHODS.md and docs/WORKFLOW.md.\n')
(W/'audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
for zipname,selected,base in [('Figure_S6_station_diagnostics.zip',s6,O/'Figure_S6'),('Fig_S3-S6_allPNG.zip',files,O)]:
 with zipfile.ZipFile(O/zipname,'w',zipfile.ZIP_DEFLATED) as z:
  for p in selected:z.write(p,p.relative_to(base))
print(json.dumps({k:v for k,v in report.items() if k!='files'},ensure_ascii=False,indent=2))
