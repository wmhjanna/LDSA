"""Match verified MPPD numeric exports to native PNSD diameters. Python 3.10+, standard library only."""
import argparse, bisect, csv, hashlib, json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ['diameter_nm', 'DF_head', 'DF_TB', 'DF_ALV', 'DF_total']

def read_csv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def write_csv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_curve(path, sum_tolerance=0.00021):
    rows=read_csv(path)
    if not rows: raise ValueError('Curve has no rows. Do not use the empty template as model output.')
    result={}
    for n,r in enumerate(rows,2):
        p=[float(r[k]) for k in FIELDS]
        if not all(math.isfinite(x) for x in p): raise ValueError(f'Non-finite value in row {n}')
        if p[0]<=0 or any(x<0 or x>1 for x in p[1:]): raise ValueError(f'Invalid diameter or fraction in row {n}; use nm and fractions, not percent')
        if abs(sum(p[1:4])-p[4])>sum_tolerance: raise ValueError(f'Head + TB + ALV != Total in row {n}; inspect exported definitions and precision')
        if p[0] in result and p[1:]!=result[p[0]]: raise ValueError(f'Conflicting DF rows at {p[0]} nm')
        result[p[0]]=p[1:]
    return sorted(result),result

def lookup(d, sizes, values, mode='exact'):
    k=bisect.bisect_left(sizes,d)
    candidates=[i for i in (k-1,k) if 0<=i<len(sizes)]
    i=min(candidates,key=lambda i:abs(sizes[i]-d))
    if abs(sizes[i]-d)<=1e-6:
        return values[sizes[i]],'MPPD_direct',sizes[i],sizes[i]
    if mode=='exact' or k==0 or k==len(sizes):return None,'missing',None,None
    lo,hi=sizes[k-1],sizes[k]
    w=(math.log(d)-math.log(lo))/(math.log(hi)-math.log(lo))
    return [a+w*(b-a) for a,b in zip(values[lo],values[hi])],'log_diameter_linear_interpolation',lo,hi

def validate_inputs(root):
    expected=json.loads((root/'inputs'/'input_summary.json').read_text(encoding='utf8'))
    channels=read_csv(root/'inputs'/'station_grid_channels.csv')
    targets=read_csv(root/'inputs'/'target_diameters.csv')
    mapping={r['diameter_id']:float(r['diameter_nm']) for r in targets}
    if len(mapping)!=len(targets):raise ValueError('Duplicate diameter IDs')
    keys=set()
    for r in channels:
        key=(r['site_code'],r['grid_id'],r['channel_index'])
        if key in keys:raise ValueError(f'Duplicate channel {key}')
        keys.add(key)
        d=float(r['diameter_nm'])
        if mapping[r['diameter_id']]!=d:raise ValueError('Channel/target mapping mismatch')
        if abs(float(r['diameter_um'])*1000-d)>1e-9:raise ValueError('nm/um mismatch')
        if int(r['in_17_500_nm'])!=int(17<=d<=500):raise ValueError('Analysis window flag mismatch')
    if len(channels)!=expected['channel_rows'] or len(targets)!=expected['unique_diameters_with_izana']:
        raise ValueError('Count mismatch')
    print(json.dumps(expected,ensure_ascii=False,indent=2))
    return expected

def verify_manifest(root):
    for line in (root/'SHA256SUMS.txt').read_text(encoding='utf8').splitlines():
        digest,name=line.split('  ',1)
        if sha(root/name)!=digest:raise ValueError(f'Package file changed or damaged: {name}')
    print('Package checksums verified. Run this before editing model_profile.json.')

def validate_profile(path, root):
    cfg=json.loads(Path(path).read_text(encoding='utf8'))
    if cfg.get('status')!='confirmed':raise ValueError('Model parameters are unconfirmed. Resolve model_profile.json before producing DF results.')
    for key in ['software_name','software_version','parameter_evidence','diameter_interpretation','region_mapping','parameters']:
        if not cfg.get(key):raise ValueError('Missing model profile field: '+key)
    required=['species_age','lung_geometry','FRC_ml','upper_airway_volume_ml','breathing_route','tidal_volume_ml','breathing_frequency_min_1','inspiratory_fraction','pause_fraction','particle_density_g_cm3','shape_assumption','hygroscopicity','inhalability_adjustment','size_distribution','diameter_input_definition']
    for key in required:
        if key not in cfg['parameters'] or cfg['parameters'][key] is None or cfg['parameters'][key]=='':
            raise ValueError('Unresolved model parameter: '+key)
    if cfg.get('curve_diameter_unit')!='nm' or cfg.get('df_unit')!='fraction':raise ValueError('Normalize diameter to nm and DF to 0–1 fraction, with original units documented')
    if cfg.get('df_definition')!='fraction_of_inhaled_particles':raise ValueError('Confirm deposited fraction denominator; regional conditional efficiency is not interchangeable')
    raw=cfg.get('raw_export_files',[])
    if not raw:raise ValueError('Keep actual MPPD .mpp/.rpt or equivalent numeric exports and list them in raw_export_files')
    evidence={}
    for name in raw:
        p=root/name
        if not p.is_file():raise ValueError('Missing raw export: '+name)
        evidence[name]=sha(p)
    return cfg,evidence

def match(args):
    root=args.root.resolve();validate_inputs(root)
    cfg,evidence=validate_profile(args.profile,root)
    sizes,values=load_curve(args.curve)
    targets=read_csv(root/'inputs'/'target_diameters.csv')
    channels=read_csv(root/'inputs'/'station_grid_channels.csv')
    channels=[r for r in channels if (args.scope=='all' or r['scope']=='europe31') and (args.window=='native' or r['in_17_500_nm']=='1')]
    wanted={r['diameter_id'] for r in channels};out=[];missing=[]
    for r in targets:
        if r['diameter_id'] not in wanted:continue
        d=float(r['diameter_nm']);df,method,lo,hi=lookup(d,sizes,values,args.mode)
        item={k:r[k] for k in ['diameter_id','diameter_nm','diameter_um']}
        item.update(dict(zip(FIELDS[1:],df if df is not None else ['']*4)))
        item.update(calculation_method=method,source_lower_nm=lo if lo is not None else '',source_upper_nm=hi if hi is not None else '')
        out.append(item)
        if df is None:missing.append(item)
    directory=args.out.resolve();directory.mkdir(parents=True,exist_ok=True)
    by={r['diameter_id']:r for r in out}
    fields=['diameter_id',*FIELDS[:1],'diameter_um',*FIELDS[1:],'calculation_method','source_lower_nm','source_upper_nm']
    write_csv(directory/'DF_by_diameter.csv',out,fields)
    write_csv(directory/'unmatched_diameters.csv',missing,fields)
    mapped=[]
    for c in channels:
        r=dict(c);r.update({k:v for k,v in by[c['diameter_id']].items() if k not in r});mapped.append(r)
    write_csv(directory/'DF_by_station_grid.csv',mapped,list(mapped[0]) if mapped else list(channels[0]))
    interpolated=sum(r['calculation_method']=='log_diameter_linear_interpolation' for r in out)
    report={'status':'incomplete' if missing else 'interpolated_requires_validation' if interpolated else 'complete_direct',
        'scope':args.scope,'window':args.window,'mode':args.mode,'target_diameters':len(out),'unmatched':len(missing),
        'direct':len(out)-len(missing)-interpolated,'interpolated':interpolated,
        'software':cfg['software_name'],'version':cfg['software_version'],'profile_sha256':sha(args.profile),
        'normalized_curve_sha256':sha(args.curve),'raw_export_sha256':evidence,
        'note':'Mapping and arithmetic checks do not independently validate the model, its scientific assumptions, or the transcription of native exports.'}
    (directory/'validation_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    return 2 if missing else 3 if interpolated else 0

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT)
    sub=p.add_subparsers(dest='cmd',required=True)
    sub.add_parser('verify-package');sub.add_parser('check-inputs')
    m=sub.add_parser('match');m.add_argument('--curve',type=Path,required=True);m.add_argument('--profile',type=Path,required=True)
    m.add_argument('--mode',choices=['exact','interpolate'],default='exact')
    m.add_argument('--scope',choices=['europe31','all'],default='europe31')
    m.add_argument('--window',choices=['native','17-500'],default='native')
    m.add_argument('--out',type=Path,default=ROOT/'results')
    a=p.parse_args()
    if a.cmd=='check-inputs':validate_inputs(a.root);return 0
    if a.cmd=='verify-package':verify_manifest(a.root);return 0
    return match(a)

if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,KeyError,FileNotFoundError) as e:raise SystemExit('STOP: '+str(e))
