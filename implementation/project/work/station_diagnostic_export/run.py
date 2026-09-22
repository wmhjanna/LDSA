from pathlib import Path
import json,hashlib,shutil,inspect,zipfile
R=Path('@PROJECT_ROOT@');TASK=Path(__file__).parent;TASK.mkdir(exist_ok=True)
B=Path('@DELIVERY_ROOT@')
DELIVERY=B/'station_diagnostic_figures';PNG=DELIVERY/'PNG';PNG.mkdir(parents=True,exist_ok=True)
A=R/'outputs/quality_control/analysis'
before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in A.iterdir() if p.is_file()}
source=(R/'work/quality_control/render_s3_s6.py').read_text().split("if '--s5-only' in sys.argv:",1)[0]
source=source.replace("O=Path('@DELIVERY_ROOT@/provider_quality_control_analysis/PNG')",f'O=Path({str(PNG)!r})')
source=source.replace("W=Path(__file__).parent/'s3_s6'", "W=R/'work/quality_control/s3_s6'")
exec(compile(source,'current_s6_template','exec'),globals())
original_clear=m.clear_breakpoint_label
audit=[]
def transparent(fig,ax):
 original_clear(fig,ax)
 for text in fig.findobj(m.mpl.text.Text):
  box=text.get_bbox_patch()
  if box is not None:box.set_facecolor('none');box.set_edgecolor('none');box.set_alpha(0)
 legend=ax.get_legend()
 if legend is not None:legend.get_frame().set_alpha(0)
 lines=[line for line in ax.lines if line.get_linestyle()=='--' and line.get_color()=='black']
 for line in lines:line.set_ydata([0,1])
 fig.canvas.draw();renderer=fig.canvas.get_renderer();text=ax.texts[0];bbox=text.get_window_extent(renderer)
 crosses=False
 for line in lines:
  x=ax.transData.transform((m.mdates.date2num(line.get_xdata()[0]),0))[0]
  crosses|=bbox.x0<=x<=bbox.x1
 assert all(list(line.get_ydata())==[0,1] for line in lines)
 assert text.get_bbox_patch() is None or text.get_bbox_patch().get_alpha()==0
 audit.append(dict(Site=fig.texts[0].get_text(),annotation=text.get_text(),background_transparent=True,breakpoint_lines=len(lines),full_height=True,line_crosses_annotation=bool(crosses)))
m.clear_breakpoint_label=transparent
src=inspect.getsource(m.s456)
src=src.replace('        if site in [\'LEI_TR\',\'ZUR_UB\',\'MEL_RB\',\'PAR_SUB\',\'PRA_SUB\']:shutil.copy2(target,DEST/f"Figure_S4_{site}.png")','').replace('    s5()','')
exec(compile(src,'s6_only','exec'),m.__dict__);m.s456()
assert len(audit)==65 and sum(r['breakpoint_lines'] for r in audit)==61
assert before=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in A.iterdir() if p.is_file()}
pointer_path=R/'outputs/current_analysis.json';current=json.load(open(pointer_path))
copied=[]
for key,target in [('Figure_2','Figure_2.png'),('Figure_3','Figure_3.png'),('Figure_S5','Figure_S5.png')]:
 origin=Path(current['figure_overrides'][key]['figure']);shutil.copy2(origin,PNG/target);copied.append(dict(figure=key,source=str(origin),target=str(PNG/target),sha256=hashlib.sha256((PNG/target).read_bytes()).hexdigest()))
images=sorted(PNG.rglob('*.png'));assert len(images)==68,len(images)
verification=dict(data_hashes_unchanged=True,total_PNG=len(images),S6_sites=65,S6_breakpoints=61,S6_without_breakpoint=4,annotation_overlap_sites=[r['Site'] for r in audit if r['line_crosses_annotation']],copied=copied,site_audit=audit)
(TASK/'verification.json').write_text(json.dumps(verification,ensure_ascii=False,indent=2));shutil.copy2(TASK/'verification.json',DELIVERY/'figure_validation.json')
(DELIVERY/'figure_inventory.md').write_text('For the study protocol and reproduction instructions, see docs/METHODS.md and docs/WORKFLOW.md.\n')
with zipfile.ZipFile(B/'station_diagnostic_figures.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(DELIVERY.rglob('*')):
  if p.is_file():z.write(p,str(Path(DELIVERY.name)/p.relative_to(DELIVERY)))
current['figure_overrides']['Figure_S6']=dict(version='transparent_annotation',directory=str(PNG/'Figure_S6'),work=str(TASK),scope='All 65 S6 panels: transparent annotation background, full-height breakpoint lines; numerical results unchanged')
current['station_diagnostic_package']=dict(directory=str(DELIVERY),zip=str(B/'station_diagnostic_figures.zip'),PNG_count=68)
pointer_path.write_text(json.dumps(current,ensure_ascii=False,indent=2))
# Small verification sheet: representative panels where the line crosses the text.
samples=[r['Site'] for r in audit if r['line_crosses_annotation']][:4]
canvas=Image.new('RGB',(1640,1040),'white')
for i,site in enumerate(samples):
 im=Image.open(PNG/'Figure_S6'/f'Figure_S6_{site}.png');im.thumbnail((820,520));canvas.paste(im,((i%2)*820,(i//2)*520))
canvas.save(TASK/'contact.png')
print('DONE',len(images),'PNG;',len(verification['annotation_overlap_sites']),'lines crossing text;',DELIVERY,flush=True)
