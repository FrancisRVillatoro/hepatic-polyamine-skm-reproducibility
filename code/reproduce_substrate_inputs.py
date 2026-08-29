#!/usr/bin/env python3
"""Optional full reproduction driver for the substrate-input extension.

This workflow is intentionally separate from ``reproduce_closures.py --full``
because the pseudo-arclength branches approaching the high-SAM asymptote can
be substantially more expensive. Outputs are written under
``data/recomputed_substrate`` and never overwrite archived canonical data.
"""
from pathlib import Path
import argparse, subprocess, sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OUT=ROOT/'data'/'recomputed_substrate'
OUT.mkdir(parents=True,exist_ok=True)

ap=argparse.ArgumentParser(); ap.add_argument('--smoke',action='store_true')
a=ap.parse_args()

def run(*args):
    cmd=[sys.executable,str(HERE/args[0]),*map(str,args[1:])]
    print('RUN',' '.join(map(str,args)),flush=True); subprocess.run(cmd,check=True)

run('substrate_input_gains.py','--output-dir',OUT)
if a.smoke:
    run('substrate_input_palc.py','--protocol','methionine','--rho','0','--max-steps','8','--p-stop','1.2','--output',OUT/'methionine_phys_smoke.csv')
    sys.exit(0)

jobs=[
 ('methionine',0,'up',None,1e9), ('methionine',0,'down',0.45,None),
 ('methionine',1,'up',7.16,None), ('methionine',1,'down',0.476,None),
 ('common-aa',0,'up',None,1e9), ('common-aa',0,'down',0.19,None),
 ('common-aa',1,'up',22.77,None), ('common-aa',1,'down',0.198,None),
]
for protocol,rho,direction,pstop,samstop in jobs:
    name=f"{protocol.replace('-','_')}_rho{rho}_{direction}.csv"
    args=['substrate_input_palc.py','--protocol',protocol,'--rho',str(rho),'--direction',direction,'--max-steps','700','--output',OUT/name]
    if pstop is not None: args += ['--p-stop',str(pstop)]
    if samstop is not None: args += ['--sam-stop',str(samstop)]
    run(*args)
print('DONE',OUT)
