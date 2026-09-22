from pathlib import Path
import importlib.util,json,hashlib,sys
ROOT=Path('@PROJECT_ROOT@')
W=ROOT/'work/spatial_template'
path=ROOT/'work/workbook_processing/redraw.py'
spec=importlib.util.spec_from_file_location('figure_template',path)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.sp=m.pd.read_csv(W/'source_site_statistics.csv')
m.DEST=Path('@DELIVERY_ROOT@/Figure1_workbook_sites')
m.DEST.mkdir(exist_ok=True)
bold_b_only='--bold-b-only' in sys.argv
original_save=m.save
def save(fig,name,source='',dpi=300):
    if bold_b_only and name!='Figure_1_b':
        m.plt.close(fig)
        return
    # Keep the former export bounds to preserve panel size while removing the letters.
    fig.canvas.draw()
    bounds=fig.get_tightbbox(fig.canvas.get_renderer()).padded(m.plt.rcParams['savefig.pad_inches'])
    for ax in fig.axes:
        for t in list(ax.texts):
            if t.get_text() in ['b','c']:t.remove()
    if name=='Figure_1_b':
        from matplotlib.font_manager import FontProperties
        for ax in fig.axes:
            ax.set_ylabel(r'Site-mean total LDSA ($\mathbf{\mu m^{2}\,cm^{-3}}$)',fontsize=12)
        for t in fig.findobj(m.mpl.text.Text):
            t.set_fontproperties(FontProperties(fname='@FONT_ROOT@/Times New Roman Bold.ttf',size=t.get_fontsize(),weight='bold'))
    labels=[t.get_text() for t in fig.findobj(m.mpl.text.Text) if t.get_visible() and t.get_text()]
    assert not any(t in ['b','c'] for t in labels)
    assert not any('Primary LDSA:' in t or 'physical monitoring sites' in t for t in labels)
    fig.savefig(m.DEST/(name+'.png'),dpi=dpi,bbox_inches=bounds,facecolor='white')
    m.manifest.append({'figure':name,'source':'station_locations.xlsx','sites':34,'PAR_period':'2015-2019','dpi':dpi,'panel_letter_removed':True,'labels':labels})
    m.plt.close(fig);print('SAVED',name)
m.save=save
m.main12(only_fig1=True)
if not bold_b_only:assert 'PAR_SUB' in m.manifest[1]['labels'] and m.manifest[1]['labels'].count('PAR_SUB')==1
(W/('verification_bold_b.json' if bold_b_only else 'verification.json')).write_text(json.dumps(m.manifest,ensure_ascii=False,indent=2))
src=Path('@DELIVERY_ROOT@/station_locations.xlsx')
(W/'provenance.json').write_text(json.dumps({'source':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scope':'Replot Figure1 b and c only using supplied workbook-format table; official analysis datasets unchanged.'},indent=2))
