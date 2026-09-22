from pathlib import Path
import hashlib,json,shutil,subprocess,sys,tempfile,zipfile
P=Path('@DELIVERY_ROOT@/MPPD_deposition_fractions')
shutil.copy2('work/ebas_download/test_df_windows.py',P/'scripts/test_pipeline.py')
(P/'results/README.md').write_text('Place validated MPPD deposition-fraction exports in this directory. Use the input checks and numerical tests before mapping fractions to native diameter channels.\n',encoding='utf8')
(P/'validation_scope.json').write_text(json.dumps({'checks':['input schemas','diameter mapping','numerical unit tests'],'test_data':'Synthetic values are used only in temporary unit-test inputs; they are not study deposition fractions.'},indent=2),encoding='utf8')
files=sorted(p for p in P.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256SUMS.txt')
(P/'SHA256SUMS.txt').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(P).as_posix()+'\n' for p in files),encoding='utf8')
files.append(P/'SHA256SUMS.txt')
dest=P.with_suffix('.zip')
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,P.name+'/'+p.relative_to(P).as_posix())
with zipfile.ZipFile(dest) as z:
 assert z.testzip() is None
 with tempfile.TemporaryDirectory(prefix='mppd_validation_') as tmp:
  z.extractall(tmp);root=Path(tmp)/P.name
  for args in [['verify-package'],['check-inputs']]:
   r=subprocess.run([sys.executable,str(root/'scripts/df_pipeline.py'),*args],cwd=root,capture_output=True,text=True)
   assert r.returncode==0,r.stderr+r.stdout
  r=subprocess.run([sys.executable,str(root/'scripts/test_pipeline.py')],cwd=root,capture_output=True,text=True)
  assert r.returncode==0,r.stderr+r.stdout
print(json.dumps({'zip':str(dest),'files':len(files),'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'zip_crc':'passed','extracted_standalone_inputs_and_tests':'passed'},ensure_ascii=False,indent=2))
