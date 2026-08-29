#!/usr/bin/env python3
"""Recompute basal logarithmic SAM gains across the proliferative coordinate."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import substrate_input_core as sic

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"


def run(nrho=21):
    rows=[]
    y=sic.basal_state(0.0)
    for rho in np.linspace(0.0,1.0,nrho):
        y,res,_=sic.solve_equilibrium(y,float(rho),1.0,"methionine")
        gm=sic.logarithmic_sam_gain_fast(y,float(rho),1.0,"methionine")
        ga=sic.logarithmic_sam_gain_fast(y,float(rho),1.0,"common-aa")
        lam=sic.dominant_nonconservation_eigenvalue(y,float(rho),1.0,"methionine")
        rows.append(dict(rho=rho,SAM=y[sic.IDX['sam']],met=y[sic.IDX['met']],
                         gain_dlnSAM_dlnmu=gm,alpha_real=lam.real,alpha_imag=lam.imag,resid=res,
                         gain_common_aa=ga))
    df=pd.DataFrame(rows)
    met=df[['rho','SAM','met','gain_dlnSAM_dlnmu','alpha_real','alpha_imag','resid']]
    cmp=pd.DataFrame({'rho':df.rho,'gain_common_aa':df.gain_common_aa,
                      'gain_methionine_only':df.gain_dlnSAM_dlnmu,
                      'common_vs_met_ratio':df.gain_common_aa/df.gain_dlnSAM_dlnmu})
    return met,cmp

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output-dir',type=Path,default=DATA); ap.add_argument('--nrho',type=int,default=21)
    a=ap.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
    met,cmp=run(a.nrho)
    met.to_csv(a.output_dir/'methionine_gain_vs_rho_recomputed.csv',index=False)
    cmp.to_csv(a.output_dir/'gain_commonAA_vs_methionine_recomputed.csv',index=False)
    print(met[['rho','gain_dlnSAM_dlnmu']].iloc[[0,-1]].to_string(index=False))
    print(cmp[['rho','gain_common_aa']].iloc[[0,-1]].to_string(index=False))
