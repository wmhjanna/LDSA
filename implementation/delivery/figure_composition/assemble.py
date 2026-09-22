from pathlib import Path
import json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
D=Path(__file__).parent
S=Path('@ATTACHMENT_ROOT@')
paths=[S/'a666a0de2a4d17bafc1ac1e39a88d00b.png',S/'figure_reference.png']
ims=[mpimg.imread(p) for p in paths]
h=max(im.shape[0] for im in ims)
widths=[round(im.shape[1]*h/im.shape[0]) for im in ims]
gap=24;w=sum(widths)+gap
fig=plt.figure(figsize=(w/300,h/300),dpi=300,facecolor='white')
x=0
for im,width in zip(ims,widths):
 ax=fig.add_axes([x/w,0,width/w,1]);ax.imshow(im,interpolation='lanczos',aspect='auto');ax.set_axis_off();x+=width+gap
# Infer the font scale from a, use the same scale for b, and align b to the 80 label.
yy,xx=np.where(np.all(ims[0][:50,:45,:3]<.3,axis=2))
scale=h/ims[0].shape[0]
top=yy.min()*scale; glyph_height=(yy.max()-yy.min()+1)*scale
font=FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf')
path=TextPath((0,0),'b',size=1,prop=font);bb=path.get_extents()
a_box=TextPath((0,0),'a',size=1,prop=font).get_extents();k=glyph_height/a_box.height
# Source tick 80 occupies this region, excluding the axis line and y-axis title.
ty,tx=np.where(np.all(ims[1][235:305,130:210,:3]<.3,axis=2))
left=widths[0]+gap+(130+tx.min())*(h/ims[1].shape[0])
trans=Affine2D().scale(k).translate(left-bb.x0*k,h-top-bb.y1*k).scale(1/w,1/h)+fig.transFigure
fig.add_artist(PathPatch(path,transform=trans,facecolor='black',edgecolor='none',clip_on=False))
fig.savefig(D/'Figure_1_map_boxplot.png',dpi=300,facecolor='white',pad_inches=0)
fig.savefig(D/'preview.png',dpi=120,facecolor='white',pad_inches=0)
plt.close(fig)
(D/'figure_composition.json').write_text(json.dumps({'layout':'map left, boxplot right; equal height; proportional scaling; 24 px gap; b added at right panel top, aligned to the left edge of tick label 80; same inferred font size as a; a and b glyph tops aligned','pixels':[w,h],'source_images':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths],'note':'Image-only assembly of supplied panels; no data recalculation; map retains supplied preview provenance'},ensure_ascii=False,indent=2))
print(w,h)
