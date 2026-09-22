"""Read-only, cross-platform validation of deposition-fraction inputs. Python 3.10+, standard library only."""
from pathlib import Path
import csv,hashlib,json,math,re
ROOT=Path(__file__).resolve().parents[1]
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def require(condition,message):
    if not condition:raise ValueError(message)
def close(a,b):return abs(float(a)-float(b))<=1e-13
def main():
    count=0
    for line in (ROOT/'SHA256SUMS.txt').read_text(encoding='utf8').splitlines():
        expected,name=line.split('  ',1);p=(ROOT/name).resolve()
        require(p.is_relative_to(ROOT),'Manifest path outside package')
        require(p.is_file(),'Missing file: '+name)
        require(hashlib.sha256(p.read_bytes()).hexdigest()==expected,'SHA256 mismatch: '+name);count+=1
    fields=['DF_head','DF_TB','DF_ALV','DF_total']
    native=read(ROOT/'raw_exports/direct/engine_values.csv')
    targets=read(ROOT/'inputs/target_diameters.csv')
    require(len(native)==len(targets)==3561,'Wrong unique-diameter count')
    by={r['diameter_id']:r for r in native};require(len(by)==3561,'Duplicate diameter IDs')
    for r,t in zip(native,targets):
        require(r['diameter_id']==t['diameter_id'],'Target ID mismatch')
        require(r['diameter_nm']==t['diameter_nm'] and r['diameter_um']==t['diameter_um'],'Exact input diameter mismatch')
        require(abs(float(r['diameter_um'])*1000-float(r['diameter_nm']))<1e-9,'Diameter unit mismatch')
        vals=[float(r[k]) for k in fields]
        require(all(math.isfinite(v) and 0<=v<=1 for v in vals),'DF outside [0,1]')
        require(abs(sum(vals[:3])-vals[3])<1e-12,'Regional sum mismatch')
    mappings={}
    for scope,expected in [('direct_europe31',6663),('direct_europe31_17_500',5196),('direct_izana',118),('direct_all',6781)]:
        rows=read(ROOT/f'results/{scope}/DF_by_station_grid.csv');require(len(rows)==expected,'Wrong channel count: '+scope)
        for r in rows:
            original=by[r['diameter_id']]
            require(r['calculation_method']=='MPPD_direct','Unexpected interpolation or missing DF')
            require(all(close(r[k],original[k]) for k in fields),'Mapped DF mismatch')
            require(close(r['diameter_nm'],original['diameter_nm']),'Mapped diameter mismatch')
        mappings[scope]=len(rows)
    comparison=read(ROOT/'manual_check/manual_comparison.csv')
    fresh=read(ROOT/'raw_exports/manual_check/direct/engine_values.csv')
    require(len(comparison)==len(fresh)==42,'Manual comparison count mismatch')
    region_matches=0;all_matches=0;sum_rounding=0
    for old,new in zip(comparison,fresh):
        reference=by[old['diameter_id']]
        require(close(old['exact_diameter_nm'],new['diameter_nm']),'Manual check input mismatch')
        require(all(close(old['recomputed_'+k],new[k]) and close(new[k],reference[k]) for k in fields),'Fresh calculation mismatch')
        diffs=[abs(float(old['manual_'+k])-float(new[k])) for k in fields]
        regions=all(v<=0.000050000001 for v in diffs[:3]);allfour=all(v<=0.000050000001 for v in diffs)
        region_matches+=regions;all_matches+=allfour
        if regions and not allfour:
            require(abs(float(old['manual_DF_total'])-sum(round(float(new[k]),4) for k in fields[:3]))<1e-10,'Total difference not explained by rounding')
            sum_rounding+=1
    require((region_matches,all_matches,sum_rounding)==(40,22,18),'Manual discrepancy summary mismatch')
    for um in ['0.02527','0.1122']:
        text=(ROOT/f'raw_exports/manual_check/gui_{um}.log').read_text(encoding='utf-8-sig')
        vals=re.search(r'RESULT head=(\S+) tb=(\S+) pul=(\S+) total=(\S+)',text)
        require(vals is not None,'Missing native Run result')
        r=min(fresh,key=lambda r:abs(float(r['diameter_um'])-float(um)))
        require(close(r['diameter_um'],um),'Missing GUI diameter')
        require(all(close(a,r[k]) for a,k in zip(vals.groups(),fields)),'Native Run comparison mismatch')
    print(json.dumps({'status':'PASS','sha256_files_verified':count,'unique_diameters':len(native),'channel_rows':mappings,
       'manual_rows':42,'regional_match_rows':region_matches,'all_four_match_rows':all_matches,
       'total_rounding_rows':sum_rounding,'regional_discrepancy_rows':2,'anomalous_rows_native_Run_check':'PASS'},indent=2))
if __name__=='__main__':
    try:main()
    except (ValueError,FileNotFoundError,KeyError) as e:raise SystemExit('FAIL: '+str(e))
