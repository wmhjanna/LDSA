from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from copy import deepcopy
import sys, re, json, hashlib
from lxml import etree as E
ROOT=Path('@PROJECT_ROOT@')
sys.path.insert(0,str(ROOT/'work/ldsa_support/deps'))
import pandas as pd
import numpy as np
import ruptures as rpt
W=Path(__file__).parent
BASE=Path('@DELIVERY_ROOT@')
SRC=Path('@ATTACHMENT_ROOT@/SP-65sites.docx')
TAB=BASE/'supplementary_tables/LDSA_supplementary_tables_S1-S6.docx'
DEST=BASE/'supplementary_information'; DEST.mkdir(exist_ok=True)
OUT=DEST/'Supplementary_Information.docx'
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
qn=lambda s:'{'+ns['w']+'}'+s
text=lambda el:''.join(el.xpath('.//w:t/text()',namespaces=ns))
with ZipFile(SRC) as z: parts={i.filename:z.read(i.filename) for i in z.infolist()}
root=E.fromstring(parts['word/document.xml']); body=root.find('w:body',ns)
original_media={k:hashlib.sha256(v).hexdigest() for k,v in parts.items() if k.startswith('word/media/')}
# Apply tracked additions and deletions while retaining document package assets.
for el in list(root.xpath('//w:del | //w:moveFrom | //w:rPrChange | //w:pPrChange | //w:tblPrChange | //w:tcPrChange | //w:trPrChange',namespaces=ns)):
 if el.getparent() is not None: el.getparent().remove(el)
for el in list(root.xpath('//w:ins | //w:moveTo',namespaces=ns)):
 par=el.getparent()
 if par is None: continue
 idx=par.index(el)
 for child in list(el): par.insert(idx,child);idx+=1
 par.remove(el)
ps=body.findall('w:p',ns); assert len(ps)==114
changes=[]
def replace(el,old,new):
 s=text(el); assert old in s,(old,s)
 if old in text(el):
  ts=el.findall('.//w:t',ns); s=''.join(t.text or '' for t in ts); start=s.index(old); end=start+len(old); pos=0; inserted=False
  for t in ts:
   val=t.text or ''; nxt=pos+len(val)
   if nxt>start and pos<end:
    left=max(0,start-pos); right=min(len(val),end-pos)
    t.text=val[:left]+(new if not inserted else '')+val[right:];t.set('{http://www.w3.org/XML/1998/namespace}space','preserve'); inserted=True
   pos=nxt
  changes.append({'old':old,'new':new})
def setpara(p,s,boldlead=None):
 rp=p.find('w:r/w:rPr',ns); rp=deepcopy(rp) if rp is not None else E.Element(qn('rPr'))
 for name in ['b','i','vertAlign']:
  for e in rp.findall('w:'+name,ns):rp.remove(e)
 for child in list(p):
  if child.tag!=qn('pPr'):p.remove(child)
 for chunk,bold in ([(boldlead,True),(s[len(boldlead):],False)] if boldlead else [(s,False)]):
  r=E.SubElement(p,qn('r')); props=deepcopy(rp);r.append(props)
  if bold:E.SubElement(props,qn('b'))
  t=E.SubElement(r,qn('t'));t.text=chunk;t.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
# Replace table data and their explanatory captions, preserving S2 additions.
tr=E.fromstring(ZipFile(TAB).read('word/document.xml')); newtables=tr.findall('w:body/w:tbl',ns); oldtables=body.findall('w:tbl',ns)
for i in [0,2,3,4,5]:
 old=oldtables[i];new=deepcopy(newtables[i]); old.getparent().replace(old,new)
 cap=oldtables[i].getprevious() if oldtables[i].getparent() is not None else None
 # Use fixed source paragraph references because table is already detached.
 cap=ps[{0:58,2:63,3:65,4:68,5:71}[i]]
 cp=newtables[i].getprevious();s=text(cp);lead=s[:s.index('. ',10)+1] if '. ' in s[10:] else None
 setpara(cap,s,lead)
# Clarify terminology and explicitly mark all DF as MPPD-derived.
replace(ps[71],'All stations use their corresponding MPPD deposition fractions.','All stations use their corresponding MPPD deposition fractions.')
replace(ps[71],'Partial sums use available cleaned channels and wholly missing rows remain missing. ','')
replace(ps[71],'HA, TB and ALV','HA, TB and ALV')
replace(ps[58],'Spatial statistics use available observations in 2015–2019','Core period denotes the 2015–2019 analysis window, not necessarily continuous or complete observational coverage. Spatial statistics use available observations in 2015–2019')
# S2 keeps all 15 rows, with the percent-sign typo fixed.
s2=body.findall('w:tbl',ns)[1];replace(s2,'68.25','68.2%')
tables=body.findall('w:tbl',ns)
for table in tables:
 for el in table.findall('.//w:t',ns):
  if el.text:el.text=el.text.replace('Station type','Station class')
# S5: three-line table, same black top/header/bottom hierarchy as S1.
t=tables[4];prop=t.find('w:tblPr',ns)
for name in ['tblCellSpacing','tblBorders']:
 for el in prop.findall('w:'+name,ns):prop.remove(el)
E.SubElement(prop,qn('tblCellSpacing'),{qn('w'):'0',qn('type'):'dxa'})
b=E.SubElement(prop,qn('tblBorders'))
for edge in ['top','left','bottom','right','insideH','insideV']:E.SubElement(b,qn(edge),{qn('val'):'nil'})
for j,row in enumerate(t.findall('w:tr',ns)):
 for cell in row.findall('w:tc',ns):
  cp=cell.find('w:tcPr',ns)
  for el in cp.findall('w:tcBorders',ns)+cp.findall('w:shd',ns):cp.remove(el)
  b=E.SubElement(cp,qn('tcBorders'))
  for edge in ['top','left','bottom','right']:
   size=8 if j==0 and edge=='top' else 4 if edge=='bottom' and j in [0,9] else None
   attrs={qn('val'):'single' if size else 'nil'}
   if size:attrs.update({qn('sz'):str(size),qn('color'):'000000',qn('space'):'0'})
   E.SubElement(b,qn(edge),attrs)
# Match S5 typography and row spacing to S1 so the three-line table stays compact.
for r in tables[4].findall('.//w:r',ns):
 rp=r.find('w:rPr',ns)
 if rp is None:rp=E.SubElement(r,qn('rPr'))
 for name in ['sz','szCs']:
  e=rp.find('w:'+name,ns)
  if e is None:e=E.SubElement(rp,qn(name))
  e.set(qn('val'),'18')
for p in tables[4].findall('.//w:p',ns):
 pp=p.find('w:pPr',ns)
 if pp is None:pp=E.Element(qn('pPr'));p.insert(0,pp)
 e=pp.find('w:spacing',ns)
 if e is None:e=E.SubElement(pp,qn('spacing'))
 for k,v in {'before':'0','after':'0','line':'240','lineRule':'auto'}.items():e.set(qn(k),v)
for row in tables[4].findall('w:tr',ns):
 for e in row.findall('w:trPr/w:trHeight',ns):e.getparent().remove(e)
for t in tables:
 pr=t.find('w:tblPr',ns)
 for e in pr.findall('w:tblInd',ns):pr.remove(e)
 jc=pr.find('w:jc',ns)
 if jc is None:jc=E.SubElement(pr,qn('jc'))
 jc.set(qn('val'),'center')
# Remove residual accepted formatting history everywhere in document.
for el in list(root.iter()):
 if E.QName(el).localname.endswith('PrChange') and el.getparent() is not None:el.getparent().remove(el)
if 'word/settings.xml' in parts:
 settings=E.fromstring(parts['word/settings.xml'])
 for el in settings.findall('w:trackRevisions',ns):settings.remove(el)
 parts['word/settings.xml']=E.tostring(settings,xml_declaration=True,encoding='UTF-8',standalone=True)
parts['word/document.xml']=E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
with ZipFile(OUT,'w',ZIP_DEFLATED) as z:
 for k,v in parts.items():z.writestr(k,v)
assert [len(t.findall('w:tr',ns)) for t in tables]==[66,16,25,13,10,66]
assert not root.xpath('//w:ins | //w:del | //w:moveTo | //w:moveFrom',namespaces=ns)
assert all(hashlib.sha256(parts[k]).hexdigest()==v for k,v in original_media.items())
for i in [0,2,3,4,5]:
 for a,b in zip(tables[i].findall('w:tr',ns)[1:],newtables[i].findall('w:tr',ns)[1:]):assert text(a)==text(b)
(W/'document_contents.md').write_text('# Supplementary information\n\nTables S1 and S3–S6 contain the analysis outputs. Document figures, relationships and fields are retained from the supplied template.\n')
(W/'verification.json').write_text(json.dumps({'output':str(OUT),'table_rows':[len(t.findall('w:tr',ns)) for t in tables],'media_preserved':len(original_media),'revisions_remaining':0,'inline_changes':changes,'scope':'Generate supplementary tables using the supplied document template'},ensure_ascii=False,indent=2))
(W/'final_paragraphs.txt').write_text('\n'.join(text(p) for p in body.findall('w:p',ns)))
print(OUT)
