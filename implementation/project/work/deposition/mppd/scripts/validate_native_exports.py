"""Validate unmodified MPPD exports, normalize native full-precision values, and record evidence."""
from pathlib import Path
import csv, json, hashlib, re, math, shutil
ROOT=Path(__file__).resolve().parents[1]
def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rows,fields):
    with Path(p).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
raw=ROOT/'raw_exports/direct'
targets=read(ROOT/'inputs/target_diameters.csv');rows=read(raw/'engine_values.csv')
assert len(rows)==len(targets)==3561
assert [r['diameter_id'] for r in rows]==[r['diameter_id'] for r in targets]
fields=['diameter_nm','DF_head','DF_TB','DF_ALV','DF_total']
labels=['Total head deposition fraction','Total TB deposition fraction','Total pulmonary deposition fraction','Total deposition fraction']
max_rounding=0;max_sum=0;hashes={}
for r,t in zip(rows,targets):
    assert r['diameter_nm']==t['diameter_nm'] and r['diameter_um']==t['diameter_um']
    values=[float(r[k]) for k in fields[1:]]
    assert all(math.isfinite(v) and 0<=v<=1 for v in values)
    max_sum=max(max_sum,abs(sum(values[:3])-values[3]));assert max_sum<1e-12
    p=raw/r['raw_report'];text=p.read_text(encoding='gb18030',errors='replace')
    assert 'Lung Geometry: Humansymmetric' in text and 'Breathing route: nasal' in text
    assert 'FRC = 3300.00' in text and 'Tidal volume: 625.0' in text and 'Breathing Frequency: 12.0' in text
    assert 'Mass density: 1.0' in text and 'Sigma_g (GSD): 1.000E+0' in text
    assert 'Inhalable fraction: 1.0000' in text
    diameter_token=re.search(r'^Diameter:\s+(\S+)',text,re.M)[1]
    exported_um=float(diameter_token)
    mantissa,exponent=diameter_token.upper().split('E')
    report_diameter_rounding=0.5*10**(int(exponent)-len(mantissa.split('.')[1]))
    assert abs(exported_um-float(t['diameter_um'])) <= report_diameter_rounding+1e-14
    for label,v in zip(labels,values):
        v_report=float(re.search(r'^'+re.escape(label)+r':\s+(\S+)',text,re.M)[1])
        delta=abs(v-v_report);max_rounding=max(max_rounding,delta);assert delta<=0.0000500001
    par=raw/r['raw_parameters'];par_lines=par.read_text().splitlines()
    assert math.isclose(float(par_lines[7]),float(t['diameter_um']),rel_tol=1e-14)
    for p in [p,par]:hashes[str(p.relative_to(ROOT)).replace('\\','/')]=sha(p)
write(ROOT/'raw_exports/curve.csv',[{k:r[k] for k in fields} for r in rows],fields)
# Compare the standard GUI action/worker route with direct native engine calls.
refs=read(ROOT/'raw_exports/validation/direct/engine_values.csv');checks=[]
for r in refs[:6]:
    um=str(float(r['diameter_um']));um='1' if um=='1.0' else um
    log=ROOT/f'raw_exports/validation/gui_{um}.log'
    text=log.read_text(encoding='utf-8-sig')
    vals=[float(x) for x in re.search(r'RESULT head=(\S+) tb=(\S+) pul=(\S+) total=(\S+)',text).groups()]
    delta=max(abs(a-float(r[k])) for a,k in zip(vals,fields[1:]));assert delta<=1e-12
    checks.append({'diameter_nm':float(r['diameter_nm']),'maximum_absolute_difference':delta})
assert all(refs[1][k]==refs[6][k] for k in fields[1:])
main_by_nm={float(r['diameter_nm']):r for r in rows}
for r in refs[:6]:
    if float(r['diameter_nm']) in main_by_nm:
        assert all(abs(float(r[k])-float(main_by_nm[float(r['diameter_nm'])][k]))<=1e-12 for k in fields[1:])
pilot=read(ROOT/'raw_exports/validation/pilot/engine_values.csv')
assert all(all(a[k]==b[k] for k in fields) for a,b in zip(pilot,rows[:10]))
write(ROOT/'results/validation_reference_points.csv',[{k:r[k] for k in fields} for r in refs[:6]],fields)
report={'status':'passed','unique_diameters':len(rows),'direct':len(rows),'interpolated':0,'missing':0,
 'raw_reports_checked':len(rows),'raw_parameter_files_checked':len(rows),'max_region_sum_absolute_error':max_sum,
 'max_native_report_rounding_error':max_rounding,'gui_vs_engine_checks':checks,
 'repeat_17nm_after_1000nm':'exact_match','first_ten_repeat':'exact_match',
 'precision_note':'curve.csv retains native engine doubles; original .rpt values are rounded by MPPD, usually to 4 decimals.',
 'scientific_scope':'Verification of specified scenario, software execution, export, and mapping; not external validation of human lung deposition predictions.'}
(ROOT/'results/native_validation_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'raw_exports/raw_file_sha256.json').write_text(json.dumps(hashes,indent=2),encoding='utf8')
install=Path('C:/Program Files (x86)/Applied Research Associates, Inc/MPPD')
software={str(p.relative_to(install)).replace('\\','/'):sha(p) for p in [install/'MPPD3.04.jar',install/'MPPD3.04.exe',install/'ModelData/hum_a.dat',install/'ModelData/hum_s.dat',install/'Jama-1.0.2.jar',install/'jxl.jar']}
(ROOT/'references/local_software_sha256.json').write_text(json.dumps(software,indent=2),encoding='utf8')
cfg=json.loads((ROOT/'model_profile.json').read_text(encoding='utf8'))
cfg.update(status='confirmed',software_name='ARA MPPD',software_version='3.04',
 parameter_evidence='MPPD 3.04 configuration with each native diameter supplied as the particle diameter; additional settings are listed below.',
 diameter_interpretation='Native PNSD mobility diameter is treated as spherical physical diameter at density 1 g/cm3; nm divided by 1000 for CMD input. GSD=1 single-size calculation.',
 region_mapping={'DF_head':'MPPD Head, inhalation plus exhalation','DF_TB':'MPPD TB','DF_ALV':'MPPD Pulmonary/P (project ALV/P label)','DF_total':'MPPD Total = Head + TB + P'},
 raw_export_files=['raw_exports/direct/engine_values.csv','raw_exports/raw_file_sha256.json'],
 execution_method='Unmodified local MPPD3.04.jar native SetParameters/MakeTree/test classes, with native ReportCommon/ReportSingle/SaveParameters exports. Custom Java runner; not an official CLI. Six diameters cross-checked against MPPDframe Run action and LongTask worker.',
 notes='Deposition fractions are calculated directly without interpolation. This profile describes the specified calculation and does not establish identical settings for other published calculations.')
cfg['parameters'].update(species_age='Human adult (Yeh/Schum standard; no specific chronological age entered)',lung_geometry='Yeh / Schum Symmetric',FRC_ml=3300,upper_airway_volume_ml=50,breathing_route='Nasal',tidal_volume_ml=625,breathing_frequency_min_1=12,inspiratory_fraction=0.5,pause_fraction=0,particle_density_g_cm3=1,shape_assumption='Spherical, aspect ratio 1',hygroscopicity='No growth; installed MPPD Particle.getHygroscopic() false',inhalability_adjustment=False,size_distribution='Single, monodisperse; GSD diameter 1, GSD length 1, correlation 0',diameter_input_definition='CMD in micrometres')
cfg['additional_settings']={'gravity_cm_s2':981,'body_orientation':'Upright','aerosol_concentration_mg_m3':1,'equivalent_diameter_model':False,'flow_model':0,'nanoparticle_specific_axial_diffusion_switch':False,'hidden_switch_evidence':'MPPD 3.04 session parameter: false.','calculation':'deposition only, constant exposure'}
(ROOT/'model_profile.json').write_text(json.dumps(cfg,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(report,ensure_ascii=True,indent=2))
