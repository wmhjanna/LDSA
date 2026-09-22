from pathlib import Path
import json,hashlib,shutil,inspect,zipfile
R=Path('@PROJECT_ROOT@');TASK=Path(__file__).parent;TASK.mkdir(exist_ok=True)
B=Path('@DELIVERY_ROOT@')
DELIVERY=B/'record_exclusion_analysis';PNG=DELIVERY/'PNG';PNG.mkdir(parents=True,exist_ok=True)
A=R/'outputs/record_exclusions/analysis'
before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in A.iterdir() if p.is_file()}
source=(R/'work/record_exclusions/render_s3_s6.py').read_text().split("if '--s5-only' in sys.argv:",1)[0]
source=source.replace("O=Path('@DELIVERY_ROOT@/provider_quality_control_analysis/PNG')",f'O=Path({str(PNG)!r})')
source=source.replace("W=Path(__file__).parent/'s3_s6'", "W=R/'work/record_exclusions/s3_s6'")
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

(TASK/'s6_verification.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
