from pathlib import Path
import json,xml.etree.ElementTree as E,concurrent.futures
from download_all import fetch,ddscheck
W=Path('work/ebas_download');O=Path('@DELIVERY_ROOT@/EBAS_observations')
d={r['file']:r for r in map(json.loads,(O/'download_log.jsonl').read_text().splitlines())};ns={'t':'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0'};catalog=E.parse(W/'catalog_refreshed.xml');items={x.get('name'):x for x in catalog.findall('.//t:dataset',ns) if x.get('urlPath')}
replacements=[];checks=[]
for r in d.values():
 if r['status']=='downloaded' or '404' not in r.get('error',''):continue
 old=r['file'];p=old.split('.');key=p[:2]+p[3:];matches=[n for n in items if n.split('.')[:2]+n.split('.')[3:]==key];assert len(matches)==1,(old,matches)
 name=matches[0];assert name!=old
 el=items[name];sz=el.find('t:dataSize',ns);size=float(sz.text)*{'bytes':1,'Kbytes':1000,'Mbytes':1000000,'Gbytes':1000000000}[sz.get('units')]
 x={'name':name,'urlPath':el.get('urlPath'),'size_bytes':size};check=ddscheck(x);checks.append(check);assert check['status']=='metadata_checked' and check['variables']
 replacements.append({'old_file':old,'new_file':name,'reason':'EBAS daily NRT revision replaced the catalog filename during download; refreshed catalog 2026-09-18','original_url':r['url'],'new_url':'https://thredds.nilu.no/thredds/fileServer/'+x['urlPath']})
 result=fetch(x,r['category'],check['variables']);assert result['status']=='downloaded',result
(W/'replacement_metadata_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
(O/'source_replacements.json').write_text(json.dumps(replacements,ensure_ascii=False,indent=2))
print('REPLACED',len(replacements),flush=True)
