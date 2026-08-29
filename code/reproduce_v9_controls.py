#!/usr/bin/env python3
"""Regenerate the v9-specific MAT controls and asymptotic audits.

Fast mode regenerates the three gain sweeps, SBML diff/conserved-pool audit,
and the simple-pole asymptotic reduction. With --full it additionally
regenerates the three rho=1 methionine PALC branches.
"""
from pathlib import Path
import argparse, os, subprocess, sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DATA=ROOT/'data'

def run(args):
    env=os.environ.copy(); env.setdefault('OPENBLAS_NUM_THREADS','1'); env.setdefault('OMP_NUM_THREADS','1')
    print('+',' '.join(map(str,args)),flush=True)
    subprocess.run(args,cwd=ROOT,env=env,check=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--full',action='store_true',help='also regenerate all three methionine PALC branches')
    ap.add_argument('--figures',action='store_true',help='regenerate v9 control figures into figures_v9/')
    a=ap.parse_args()
    py=sys.executable
    for pol in ('published','constant-total','flux-matched'):
        run([py,'code/mat_control_gain_sweeps.py','--policy',pol,'--nrho','21','--outdir','data'])
    run([py,'code/sbml_v1_v2_diff.py','--outdir','data'])
    run([py,'code/sam_asymptotic_boundary.py','--outdir','data'])
    if a.full:
        for pol in ('published','constant-total','flux-matched'):
            out=f'data/mat_control_methionine_palc_{pol}.csv'
            run([py,'code/mat_control_palc.py','--policy',pol,'--rho','1','--p-stop','15','--output',out])
    if a.figures:
        run([py,'code/make_v9_control_figures.py','--outdir','figures_v9'])
    print('v9 control reproduction complete')

if __name__=='__main__': main()
