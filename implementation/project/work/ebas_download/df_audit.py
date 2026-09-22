from pathlib import Path
import json,math,hashlib
import openpyxl
W=Path('work/ebas_download')
base=Path('@INPUT_ROOT@/data_sources')
paths=[base/'processed_excel_data/LDSA/DF.xlsx',base/'processed_excel_data/processed_excel_data/DF.xlsx']
blocks=[]
for path in paths[:1]:
 w=openpyxl.load_workbook(path,read_only=True,data_only=True)
 for sheet in w:
  rr=list(sheet.iter_rows(values_only=True));h=rr[0]
  for j in range(0,len(h),5):
   label=str(h[j]); points=[]
   for i,row in enumerate(rr[1:],2):
    a=row[j:j+5]
    if len(a)==5 and all(isinstance(v,(int,float)) and math.isfinite(v) for v in a) and a[0]>0 and all(0<=v<=1 for v in a[1:]):points.append([a[0]*1000,*a[1:]])
   if points:blocks.append({'label':label,'column':j+1,'headers':h[j:j+5],'n':len(points),'diameters':[p[0] for p in points],'points':points})
 w.close()
selected=json.load(open('@DELIVERY_ROOT@/EBAS_input_selection/selected_netcdf_records.json'))
detail=json.load(open('@DELIVERY_ROOT@/EBAS_station_inventory/station_inventory_details.json'))
result=[]
for st in selected:
 code=st['code']; own=[b for b in blocks if b['label'].split('_')[0]==code]
 grids={r['grid']:r['diameters_nm'] for r in detail['files'] if r['code']==code and r['channels']>=3}
 tests=[]
 for grid,dp in grids.items():
  primary=[d for d in dp if 17<=d<=500]
  best=None
  for b in own:
   distances=[min(abs(x-y) for y in b['diameters']) for x in primary]
   score=max(distances) if distances else 1e9
   val={'label':b['label'],'max_delta_nm':score,'exact_count':sum(v<1e-6 for v in distances),'within_0.05nm_count':sum(v<=.050001 for v in distances),'n_target':len(primary)}
   if best is None or score<best['max_delta_nm']:best=val
  tests.append({'grid':grid,'n_channels':len(dp),'low':min(dp),'high':max(dp),'best_same_code':best})
 result.append({'code':code,'name':st['name'],'df_blocks':len(own),'labels':[b['label'] for b in own],'grids':tests})
out={'source':str(paths[0]),'same_sha256_copy':hashlib.sha256(paths[0].read_bytes()).hexdigest()==hashlib.sha256(paths[1].read_bytes()).hexdigest(),'blocks':blocks,'stations':result}
(W/'europe_df_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
print('blocks',len(blocks),'duplicate_copy',out['same_sha256_copy'])
for r in result:
 exact=sum(t['best_same_code'] is not None and t['best_same_code']['max_delta_nm']<1e-6 for t in r['grids'])
 near=sum(t['best_same_code'] is not None and t['best_same_code']['max_delta_nm']<=.050001 for t in r['grids'])
 print(r['code'],r['name'],'blocks',r['df_blocks'],'grids',len(r['grids']),'exact',exact,'within0.05nm',near)
print('sites with own code',sum(bool(r['df_blocks']) for r in result))
print('all grids exact',sum(all(t['best_same_code'] is not None and t['best_same_code']['max_delta_nm']<1e-6 for t in r['grids']) for r in result))
print('all grids rounded',sum(all(t['best_same_code'] is not None and t['best_same_code']['max_delta_nm']<=.050001 for t in r['grids']) for r in result))
