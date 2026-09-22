#!/usr/bin/env python3
"""Materialize source scripts with configured portable path tokens; never execute them."""
from pathlib import Path
import argparse, json, re, shutil, ast
ROOT=Path(__file__).resolve().parents[1]
TEXT={'.py','.r','.json','.csv','.js','.sh','.ps1'}
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--destination',type=Path,required=True)
    args=ap.parse_args();config=args.config.resolve()
    values=json.loads(config.read_text())
    for key,value in values.items():
        if key not in {'PYTHON','RSCRIPT'}:
            p=Path(value).expanduser();values[key]=str((config.parent/p).resolve()) if not p.is_absolute() else str(p)
    if args.destination.exists():raise SystemExit('Destination exists; choose an empty destination.')
    source=ROOT/'implementation'
    unresolved=[];count=0
    for p in source.rglob('*'):
        if not p.is_file():continue
        out=args.destination/p.relative_to(source);out.parent.mkdir(parents=True,exist_ok=True)
        if p.suffix.lower() in TEXT:
            s=p.read_text(encoding='utf-8-sig')
            for key,value in values.items():s=s.replace('@'+key+'@',value.replace('\\','/'))
            tokens=re.findall(r'@[A-Z_]+@',s)
            if tokens:unresolved.append({'file':str(p.relative_to(source)),'tokens':sorted(set(tokens))})
            if p.suffix=='.py':ast.parse(s)
            out.write_text(s,encoding='utf-8')
        else:shutil.copyfile(p,out)
        count+=1
    (args.destination/'staging_report.json').write_text(json.dumps({'files':count,'unresolved_tokens':unresolved,'data_supplied':False,'pipeline_executed':False},indent=2))
    if unresolved:raise SystemExit('Unresolved tokens remain; see staging_report.json.')
    print(f'Staged {count} source and resource files. No data downloaded and no pipeline executed.')
if __name__=='__main__':main()
