"""Word table typography helpers."""
def rprop(p):
 rr=p.find('.//w:rPr',ns);r=deepcopy(rr) if rr is not None else E.Element(qn('rPr'))
 for v in r.findall('w:vertAlign',ns):r.remove(v)
 return r

def set_p(p,text,bold=None):
 rp=rprop(p)
 if bold is not None:
  for tag in ['b','bCs']:
   x=rp.find('w:'+tag,ns)
   if x is not None:rp.remove(x)
   if bold:E.SubElement(rp,qn(tag))
 for child in list(p):
  if child.tag!=qn('pPr'):p.remove(child)
 r=E.SubElement(p,qn('r'));r.append(rp);t=E.SubElement(r,qn('t'));t.set('{http://www.w3.org/XML/1998/namespace}space','preserve');t.text=text

def set_cell(cell,text,chemical=False):
 ps=cell.findall('w:p',ns);p=ps[0];rp=rprop(p);set_p(p,str(text))
 for extra in ps[1:]:cell.remove(extra)
 if chemical:
  for child in list(p):
   if child.tag!=qn('pPr'):p.remove(child)
  pattern=r'(PM)(2\.5|10|1)|(NO|SO)(2)|(O)(3)';pos=0
  chunks=[]
  for m in re.finditer(pattern,str(text)):
   if m.start()>pos:chunks.append((str(text)[pos:m.start()],False))
   chunks.extend([(m.group(1) or m.group(3) or m.group(5),False),(m.group(2) or m.group(4) or m.group(6),True)]);pos=m.end()
  chunks.append((str(text)[pos:],False))
  for text,sub in chunks:
   if not text:continue
   r=E.SubElement(p,qn('r'));pr=deepcopy(rp)
   if sub:E.SubElement(pr,qn('vertAlign')).set(qn('val'),'subscript')
   r.append(pr);t=E.SubElement(r,qn('t'));t.set('{http://www.w3.org/XML/1998/namespace}space','preserve');t.text=text

def fmt(x,n=2,sign=False):
 if not np.isfinite(x):return '—'
 v=f'{x:+.{n}f}' if sign else f'{x:.{n}f}';return v.replace('-','−')

def caption(p,title,explanation):
 orig_runs=p.findall('w:r',ns);boldrp=deepcopy(orig_runs[0].find('w:rPr',ns));reg=next((r for r in reversed(orig_runs) if r.find('w:rPr/w:b',ns) is None),orig_runs[-1]);regrp=deepcopy(reg.find('w:rPr',ns))
 for rp in [boldrp,regrp]:
  if rp is not None:
   for x in rp.findall('w:vertAlign',ns):rp.remove(x)
 if regrp is not None:
  for tag in ['b','bCs']:
   for x in regrp.findall('w:'+tag,ns):regrp.remove(x)
 for child in list(p):
  if child.tag!=qn('pPr'):p.remove(child)
 for text,rp in [(title,boldrp),(' '+explanation,regrp)]:
  if not text.strip():continue
  r=E.SubElement(p,qn('r'))
  if rp is not None:r.append(deepcopy(rp))
  t=E.SubElement(r,qn('t'));t.set('{http://www.w3.org/XML/1998/namespace}space','preserve');t.text=text
 pp=p.find('w:pPr',ns)
 if pp is None:pp=E.SubElement(p,qn('pPr'))
 E.SubElement(pp,qn('keepNext'))
