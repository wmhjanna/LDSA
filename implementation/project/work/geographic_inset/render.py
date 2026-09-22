from pathlib import Path
import sys,json,hashlib
W=Path(__file__).resolve().parent
ROOT=W.parents[1]
sys.path[:0]=[str(W/'deps'),str(ROOT/'work/ldsa_support/deps')]
import shapefile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MPath
from matplotlib.font_manager import FontProperties
BASE=Path('@DELIVERY_ROOT@')
OUT=BASE/'geographic_inset';OUT.mkdir(exist_ok=True)
SRC=BASE/'europe_map/Figure_1_a_hollow_preview.png'
im=mpimg.imread(SRC);h,w=im.shape[:2]
font='@FONT_ROOT@/Times New Roman Bold.ttf'
plt.rcParams.update({'font.family':'Times New Roman','font.weight':'bold','axes.titleweight':'bold','pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(w/150,h/150),facecolor='white',dpi=300)
bg=fig.add_axes([0,0,1,1]);bg.imshow(im,aspect='auto');bg.axis('off')
# Code-native scientific inset, on unused eastern-European background.
ax=fig.add_axes([.69,.52,.29,.19]);ax.set_facecolor('#eaf8fd')
extent=(-82,-66,38,48)
def visible(s):
 a,b,c,d=s.bbox
 return not(c<extent[0] or a>extent[1] or d<extent[2] or b>extent[3])
def shapes(name):
 return shapefile.Reader(str(next((W/name).glob('*.shp')))).shapes()
def polygon(s,fc,ec,lw,z):
 pts=np.array(s.points);parts=list(s.parts)+[len(pts)]
 verts=[];codes=[]
 for i,j in zip(parts[:-1],parts[1:]):
  ring=pts[i:j]
  if len(ring)<3:continue
  verts.extend(ring);codes.extend([MPath.MOVETO]+[MPath.LINETO]*(len(ring)-2)+[MPath.CLOSEPOLY])
 ax.add_patch(PathPatch(MPath(verts,codes),facecolor=fc,edgecolor=ec,lw=lw,zorder=z))
for s in shapes('admin_0_countries'):
 if visible(s):polygon(s,'#eeeeec','#9c9c9c',.55,1)
for s in shapes('admin_1_states_provinces_lines'):
 if visible(s):
  pts=np.array(s.points);parts=list(s.parts)+[len(pts)]
  for i,j in zip(parts[:-1],parts[1:]):ax.plot(pts[i:j,0],pts[i:j,1],color='#b7b7b7',lw=.35,zorder=2)
for s in shapes('lakes'):
 if visible(s):polygon(s,'#eaf8fd','#b8d6e2',.35,3)
ax.set_xlim(extent[:2]);ax.set_ylim(extent[2:]);ax.set_aspect(1/np.cos(np.deg2rad(43)),adjustable='box')
ax.set_xticks([]);ax.set_yticks([])
for sp in ax.spines.values():sp.set_linewidth(.7);sp.set_color('#777777')
ax.set_title('Northeastern US',fontsize=10,pad=5,fontproperties=FontProperties(fname=font,size=10))
lon,lat=-77.55,43.15
# ROC mean = 23.9988: blue 20–30 bin; UB circle; hollow centre.
ax.plot(lon,lat,'o',mfc='none',mec='white',mew=3.5,ms=8,zorder=5)
ax.plot(lon,lat,'o',mfc='none',mec='#1764bc',mew=1.6,ms=8,zorder=6)
ax.annotate('ROC_UB',(lon,lat),xytext=(-74.9,42.5),fontsize=9,fontproperties=FontProperties(fname=font,size=9),va='center',arrowprops={'arrowstyle':'-','color':'#444444','lw':.6},zorder=7)
for label,x,y in [('CANADA',-74.5,46.4),('USA',-76.5,39.5)]:ax.text(x,y,label,ha='center',color='#555555',fontproperties=FontProperties(fname=font,size=7.5))
for fmt in ['png','pdf','svg']:fig.savefig(OUT/f'Figure_1_a_ROC_UB_inset.{fmt}',dpi=300,facecolor='white',pad_inches=0)
fig.savefig(OUT/'preview.png',dpi=150,facecolor='white',pad_inches=0)
plt.close(fig)
(OUT/'provenance.json').write_text(json.dumps({'source_map':str(SRC),'source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'ROC_UB':{'lon':lon,'lat':lat,'coordinate_source':str(ROOT/'outputs/german_RB_revision/site_map_coordinates.csv'),'mean_LDSA':23.99881867777555,'class':'UB','bin':'20–30','statistical_source':str(ROOT/'outputs/spatial_statistics/site_spatial_statistics.csv')},'inset_geometry':'Natural Earth 1:50m countries, state/province boundaries, lakes','extent':extent,'note':'Inset added using geographic coordinates; existing European raster retained. No statistical recalculation; PDF/SVG contain raster European base with vector inset.'},indent=2,ensure_ascii=False))
print(OUT)
