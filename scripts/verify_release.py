#!/usr/bin/env python3
"""Verify every release file against SHA256SUMS."""
from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parents[1]
def main():
    bad=[];count=0
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        digest,name=line.split('  ',1);p=ROOT/name
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:bad.append(name)
        count+=1
    if bad:raise SystemExit('Verification failed: '+', '.join(bad))
    print(f'PASS: {count} release files verified.')
if __name__=='__main__':main()
