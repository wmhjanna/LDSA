def align_fig2_layout(fig):
    """Align outer text by rendered bounds, leaving all data artists intact."""
    if len(fig.axes)!=7:return
    a,c,b,study,ci,weight,forest=fig.axes
    a.get_subplotspec().get_gridspec().update(hspace=.18)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    # Balance the visible white space, accounting for b having an x-axis title.
    gap_ab=min(t.get_window_extent(renderer).y0 for t in a.get_xticklabels() if t.get_visible())-b._left_title.get_window_extent(renderer).y1
    gap_bc=b.xaxis.label.get_window_extent(renderer).y0-c._left_title.get_window_extent(renderer).y1
    pos=b.get_position();b.set_position([pos.x0,pos.y0+(gap_ab-gap_bc)/2/fig.bbox.height,pos.width,pos.height])
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    alabel=next(t for t in a.texts if t.get_text()=='a')
    dlabel=next(t for t in study.texts if t.get_text()=='d')
    dtitle=next(t for t in study.texts if t.get_text().startswith('Random-effects synthesis'))
    for text,target in [(dlabel,alabel),(dtitle,a.title)]:
        box=text.get_window_extent(renderer);y=target.get_window_extent(renderer).y0
        x=fig.transFigure.inverted().transform((box.x0,y))[0]
        text.set_transform(fig.transFigure);text.set_position((x,y/fig.bbox.height));text.set_va('bottom');text.set_ha('left')
    footer=next(t for t in study.texts if t.get_text().startswith('Subgroup difference'))
    clabel=c.xaxis.label
    cb=clabel.get_window_extent(renderer);fb=footer.get_window_extent(renderer)
    target_bottom=(cb.y0+fb.y0)/2
    y=c.transAxes.inverted().transform((0,target_bottom))[1]
    c.xaxis.set_label_coords(.5,y);clabel.set_va('bottom')
    footer.set_transform(fig.transFigure);footer.set_position((fb.x0/fig.bbox.width,target_bottom/fig.bbox.height));footer.set_va('bottom');footer.set_ha('left')
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    assert abs(dlabel.get_window_extent(renderer).y0-alabel.get_window_extent(renderer).y0)<1
    assert abs(dtitle.get_window_extent(renderer).y0-a.title.get_window_extent(renderer).y0)<1
    assert abs(clabel.get_window_extent(renderer).y0-footer.get_window_extent(renderer).y0)<1
    report={'a_d_letters_bottom_delta_px':float(dlabel.get_window_extent(renderer).y0-alabel.get_window_extent(renderer).y0),'a_d_titles_bottom_delta_px':float(dtitle.get_window_extent(renderer).y0-a.title.get_window_extent(renderer).y0),'c_d_bottom_text_delta_px':float(clabel.get_window_extent(renderer).y0-footer.get_window_extent(renderer).y0),'hspace':.18}
    (Path(__file__).parent/'alignment_verification.json').write_text(json.dumps(report,indent=2))
