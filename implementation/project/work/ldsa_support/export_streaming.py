"""Fallback after Artifact Tool exceeds 8 GiB exporting a 53,385-row sheet."""
from pathlib import Path
import gzip,csv,json,datetime,sys,zipfile,re
import xlsxwriter
ROOT=Path('@PROJECT_ROOT@')
for d in sorted((ROOT/'outputs/ldsa_cleaned_reanalysis/sites').iterdir()):
    if not (d/'complete.json').exists():continue
    p=d/f'{d.name}_LDSA_standardized.xlsx'
    if p.exists():
        try:
            with zipfile.ZipFile(p) as z:dim=re.search(rb'<dimension ref="[^"]*?([0-9]+)"',z.read('xl/worksheets/sheet1.xml')[:1000])
            if dim and int(dim.group(1))==json.loads((d/'complete.json').read_text())['rows']+1:continue
        except zipfile.BadZipFile:pass
    wb=xlsxwriter.Workbook(p,{'constant_memory':True});ws=wb.add_worksheet('Data')
    h=wb.add_format({'font_name':'Arial','font_size':10,'bold':True,'font_color':'white','bg_color':'#1F4E78','text_wrap':True,'valign':'vcenter'})
    num=wb.add_format({'font_name':'Arial','font_size':10,'num_format':'0.000'});dt=wb.add_format({'font_name':'Arial','font_size':10,'num_format':'yyyy-mm-dd hh:mm:ss'})
    with gzip.open(d/f'{d.name}_LDSA_standardized.csv.gz','rt') as f:
        reader=csv.reader(f);headers=next(reader);ws.write_row(0,0,headers,h)
        for i,row in enumerate(reader,1):
            values=[None if v.strip()=='' else datetime.datetime.fromisoformat(v) if j in [2,3] else float(v) if j>=4 else v for j,v in enumerate(row)]
            ws.write_row(i,0,values)
    ws.set_column(0,0,16);ws.set_column(1,1,25);ws.set_column(2,3,22,dt);ws.set_column(4,len(headers)-1,19,num);ws.set_row(0,44);ws.freeze_panes(1,4);ws.autofilter(0,0,i,len(headers)-1);ws.hide_gridlines(2);wb.close()
    print(d.name,i,flush=True)
