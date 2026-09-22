from pathlib import Path
import sys,os,json,importlib.util,hashlib,shutil,zipfile
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent;O=R/'outputs/integration_sensitivity';D=Path('@DELIVERY_ROOT@/integration_sensitivity_figures');O.mkdir(exist_ok=True);D.mkdir(exist_ok=True);os.environ['MPLCONFIGDIR']=str(W/'mplconfig');sys.path.insert(0,str(R/'work/ldsa_support/deps'))
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager,FontProperties
from scipy.stats import pearsonr
from PIL import Image
for n in ['Times New Roman.ttf','Times New Roman Bold.ttf','Times New Roman Italic.ttf','Times New Roman Bold Italic.ttf']:fontManager.addfont('@FONT_ROOT@/'+n)
font='@FONT_ROOT@/Times New Roman Bold.ttf'
def load(n):
 spec=importlib.util.spec_from_file_location(n,R/'outputs'/f'{n}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
s1=load('redraw_figure_s1');s2=load('redraw_figure_s2')
sens=pd.concat([pd.read_csv(W/'workbook_sensitivity.csv'),pd.read_csv(W/'netcdf_sensitivity.csv')],ignore_index=True);conf=pd.concat([pd.read_csv(W/'workbook_coverage.csv'),pd.read_csv(W/'netcdf_coverage.csv')],ignore_index=True);coverage=conf.groupby('Site')[['covers_10_500','covers_10_800']].all();coverage['cohort']=np.where(coverage.index.isin(pd.read_csv(W/'workbook_sources.csv').Site),'Workbook input','netcdf');assert len(coverage)==65 and not sens.duplicated(['Site','Range','Region']).any();sens.to_csv(O/'size_range_sensitivity_current65.csv',index=False);coverage.to_csv(O/'integration_interval_coverage.csv');conf.to_csv(O/'integration_configuration_coverage.csv',index=False)
checks=[];plotdata=[];stats=[]
def save(fig,name):
 plt.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.bf':'Times New Roman:bold','mathtext.it':'Times New Roman:italic','mathtext.fallback':'stix'})
 for t in fig.findobj(matplotlib.text.Text):
  t.set_text(t.get_text().replace('μm²/cm³',r'μm² cm$^{\mathbf{-3}}$'));t.set_fontproperties(FontProperties(fname=font,size=t.get_fontsize(),weight='bold'))
 if name=='Figure_S1':
  for ax in fig.axes:
   for t in ax.texts:
    if t.get_text() in list('abcdef'):t.set_position((-.060,1.025));t.set_va('bottom')
 fig.canvas.draw()
 if name=='Figure_S1':
  a=fig.axes[3];limits=a.get_ylim();a.set_yticks([v for v in a.get_yticks() if not np.isclose(v,80)]);a.set_ylim(limits)
 fig.canvas.draw();ren=fig.canvas.get_renderer()
 for a in fig.axes:
  for t in [a.xaxis.label,a.yaxis.label,a.title,*a.texts]:
   if t.get_text():
    b=t.get_window_extent(ren);assert b.x0>=-1 and b.y0>=-1 and b.x1<=fig.bbox.width+1 and b.y1<=fig.bbox.height+20,(name,t.get_text(),b.bounds)
  points=np.concatenate([c.get_offsets() for c in a.collections if len(c.get_offsets())],axis=0) if a.collections else np.empty((0,2))
  if len(points):assert (points[:,1]>=a.get_ylim()[0]).all() and (points[:,1]<=a.get_ylim()[1]).all()
 checks.append(dict(figure=name,axes=len(fig.axes),texts=[t.get_text() for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()],labels_within_canvas=True))
 p=D/(name+'.png');fig.savefig(p,dpi=600,bbox_inches='tight',facecolor='white');plt.close(fig);im=Image.open(p);im.thumbnail((1900,1500));im.save(W/(name+'_preview.png'));print('SAVED',name,flush=True)
s1.configure_style();fig,axs=plt.subplots(2,3,figsize=(11.6,6.65))
for row,stat in enumerate(['Mean','Median']):
 wide=sens[sens.Region.eq('Total')].pivot(index='Site',columns='Range',values=stat)
 for j,(a,b) in enumerate([('Primary','10-500'),('Primary','10-800'),('10-500','10-800')]):
  eligible=coverage.index[coverage['covers_10_500' if j==0 else 'covers_10_800']];pair=wide.loc[eligible,[a,b]].dropna();pair.columns=['x','y'];lo=pair.min().min()*.9;hi=pair.max().max()*1.1;letter=chr(97+row*3+j);corr=s1.draw_panel(axs[row,j],pair,'17–500 nm' if a=='Primary' else a.replace('-','–')+' nm',b.replace('-','–')+' nm',stat.lower(),letter,(lo,hi),False)
  # Truncate displayed values only when rounding would produce 1.000.
  if corr < 1 and f'{corr:.3f}' == '1.000':
   for text in axs[row,j].texts:
    if text.get_text().startswith('Pearson r ='):text.set_text(f'Pearson r = {np.trunc(corr*1000)/1000:.3f}')
  stats.append(dict(Figure='S1',Panel=letter,Statistic=stat,Comparison=a+' vs '+b,n=len(pair),Workbook_sites=int(coverage.loc[pair.index].cohort.eq('Workbook input').sum()),NetCDF_sites=int(coverage.loc[pair.index].cohort.eq('netcdf').sum()),Pearson_r=corr,Pearson_p=pearsonr(pair.x,pair.y).pvalue))
  for site,r in pair.iterrows():plotdata.append(dict(Figure='S1',Panel=letter,Site=site,Region='Total',x=r.x,y=r.y,Statistic=stat))
fig.subplots_adjust(left=.075,right=.992,bottom=.095,top=.985,wspace=.24,hspace=.22);save(fig,'Figure_S1')
s2.configure_style();fig,axs=plt.subplots(2,2,figsize=(9.3,6));eligible=coverage.index[coverage.covers_10_800]
for j,stat in enumerate(['Mean','Median']):
 wide=sens[sens.Site.isin(eligible)].pivot(index='Site',columns=['Region','Range'],values=stat)
 vals=[]
 for region in ['HA','TB','ALV','Total']:
  v=(100*(wide[region,'10-800']-wide[region,'10-500'])/wide[region,'10-800']).replace([np.inf,-np.inf],np.nan).dropna();assert len(v)==30;vals.append(v.to_numpy());letter=chr(97+j)
  for site,x in v.items():plotdata.append(dict(Figure='S2',Panel=letter,Site=site,Region=region,y=x,Statistic=stat))
  stats.append(dict(Figure='S2',Panel=letter,Statistic=stat,Region=region,n=len(v),Mean=v.mean(),Median=v.median(),Minimum=v.min(),Maximum=v.max()))
 extent=np.concatenate(vals);limits=(min(-1,extent.min()-1),np.ceil((extent.max()+1)/5)*5);s2.draw_boxplot(axs[0,j],vals,['HA','TB','ALV','Total'],['#1f77b4','#ff7f0e','#2ca02c','#d62728'],chr(97+j),'Site-level '+stat.lower()+'s','Contribution of 500–800 nm (%)',limits,next(x for x in [2,5,10,20,50] if x>=(limits[1]-limits[0])/10),101+j)
 vals=[]
 for region in ['HA','TB','ALV']:
  v=(100*(wide[region,'10-800']/wide['Total','10-800']-wide[region,'10-500']/wide['Total','10-500'])).replace([np.inf,-np.inf],np.nan).dropna();assert len(v)==30;vals.append(v.to_numpy());letter=chr(99+j)
  for site,x in v.items():plotdata.append(dict(Figure='S2',Panel=letter,Site=site,Region=region,y=x,Statistic=stat))
  stats.append(dict(Figure='S2',Panel=letter,Statistic=stat,Region=region,n=len(v),Mean=v.mean(),Median=v.median(),Minimum=v.min(),Maximum=v.max()))
 extent=np.concatenate(vals);limits=(np.floor((extent.min()-.2)*2)/2,np.ceil((extent.max()+.2)*2)/2);s2.draw_boxplot(axs[1,j],vals,['HA','TB','ALV'],['#1f77b4','#ff7f0e','#2ca02c'],chr(99+j),'10–800 nm minus 10–500 nm ('+stat.lower()+'s)','Change in fraction (percentage points)',limits,next(x for x in [.5,1,2,2.5,5,10] if x>=(limits[1]-limits[0])/8),103+j)
fig.subplots_adjust(left=.095,right=.992,bottom=.085,top=.975,wspace=.23,hspace=.24);save(fig,'Figure_S2')
pd.DataFrame(stats).to_csv(O/'integration_sensitivity_statistics.csv',index=False);pd.DataFrame(plotdata).to_csv(O/'integration_sensitivity_plot_data.csv',index=False);(W/'render_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2));print(pd.DataFrame(stats).head(6).to_string(index=False))
for p in O.iterdir():shutil.copy2(p,D/p.name)
