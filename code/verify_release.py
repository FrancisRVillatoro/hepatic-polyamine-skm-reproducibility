#!/usr/bin/env python3
"""Static integrity checks for the v1.0.1 code-and-data release object."""
from pathlib import Path
import hashlib
import sys

ROOT = Path(__file__).resolve().parent.parent
forbidden = {'.tex', '.pdf', '.bib', '.doc', '.docx'}
bad = [p.relative_to(ROOT) for p in ROOT.rglob('*')
       if p.is_file() and p.suffix.lower() in forbidden]
if bad:
    print('ERROR: manuscript/submission-like files found:')
    for p in bad: print(' ', p)
    sys.exit(1)

expected = [
    'code/BIOMD0000000190_url.xml','code/BIOMD0000000674_url.xml','code/BIOMD0000000450_url.xml',
    'requirements-lock.txt','.zenodo.json','CITATION.cff','V1_0_1_ADDITIONS_MANIFEST.tsv','SHA256SUMS_v1.0.1.txt',
    'code/compatibility_spectral_projection.py','code/cusp_nondegeneracy_check.py','code/mu_inf_precision_audit.py',
    'code/native_meal_cycle_control.py','code/reproduce_v10_2_audits.py',
    'data/compatibility_spectral_projection_tail.csv','data/compatibility_spectral_projection_summary.json',
    'data/cusp_nondegeneracy_check.json','data/mu_inf_precision_estimates_v10_2.csv',
    'data/mu_inf_precision_summary_v10_2.json','data/native_meal_cycle_day26_keypoints.csv',
    'data/native_meal_cycle_summary_v10_2.json']
missing = [x for x in expected if not (ROOT/x).is_file()]
if missing:
    print('ERROR: missing required release files:',*missing,sep='\n  '); sys.exit(1)
version=(ROOT/'VERSION').read_text(encoding='utf-8').strip()
if version!='1.0.1': print('ERROR: VERSION is',version,'expected 1.0.1'); sys.exit(1)
checks={}
for line in (ROOT/'SHA256SUMS_v1.0.1.txt').read_text(encoding='utf-8').splitlines():
    if not line.strip(): continue
    h,name=line.split(None,1); checks[name.strip().removeprefix('./')]=h
errors=[]
for rel,h in checks.items():
    p=ROOT/rel
    if not p.is_file(): errors.append(f'missing {rel}'); continue
    if hashlib.sha256(p.read_bytes()).hexdigest()!=h: errors.append(f'checksum mismatch {rel}')
if errors: print('ERROR:',*errors,sep='\n  '); sys.exit(1)
print('release verification OK')
print('version:',version)
print('v1.0.1 checksummed files:',len(checks))
print('python files:',len(list((ROOT/'code').glob('*.py'))))
print('SBML files:',len(list((ROOT/'code').glob('*.xml'))))
print('CSV files:',len(list((ROOT/'data').glob('*.csv'))))
print('JSON data files:',len(list((ROOT/'data').glob('*.json'))))
print('manuscript/submission files: 0')
