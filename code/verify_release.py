#!/usr/bin/env python3
"""Static checks for the code-and-data release object."""
from pathlib import Path
import hashlib, sys
ROOT=Path(__file__).resolve().parent.parent
forbidden={'.tex','.pdf','.bib','.doc','.docx'}
bad=[p.relative_to(ROOT) for p in ROOT.rglob('*') if p.is_file() and p.suffix.lower() in forbidden]
if bad:
    print('ERROR: manuscript/submission-like files found:')
    for p in bad: print(' ',p)
    sys.exit(1)
expected=['code/BIOMD0000000190_url.xml','code/BIOMD0000000674_url.xml','code/BIOMD0000000450_url.xml',
          'requirements-lock.txt','.zenodo.json','CITATION.cff','SHA256SUMS.txt']
missing=[x for x in expected if not (ROOT/x).is_file()]
if missing:
    print('ERROR: missing required release files:',*missing,sep='\n  '); sys.exit(1)
# Validate checksums, excluding the checksum file itself.
checks={}
for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
    if not line.strip(): continue
    h,name=line.split(None,1); checks[name.strip().removeprefix('./')]=h
errors=[]
for rel,h in checks.items():
    p=ROOT/rel
    if not p.is_file(): errors.append(f'missing {rel}'); continue
    got=hashlib.sha256(p.read_bytes()).hexdigest()
    if got!=h: errors.append(f'checksum mismatch {rel}')
if errors:
    print('ERROR:',*errors,sep='\n  '); sys.exit(1)
print('release verification OK')
print('python files:',len(list((ROOT/'code').glob('*.py'))))
print('SBML files:',len(list((ROOT/'code').glob('*.xml'))))
print('CSV files:',len(list((ROOT/'data').glob('*.csv'))))
print('manuscript/submission files: 0')
