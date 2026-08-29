#!/usr/bin/env python3
"""Pseudo-arclength continuation for frozen substrate-input protocols.

Continuation variables are (log x, log p), where p is either the methionine-only
multiplier mu or the common-amino-acid multiplier eta.  New files are written by
default, so archived canonical CSV data are never overwritten accidentally.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import substrate_input_core as sic

HERE=Path(__file__).resolve().parent.parent; DATA=HERE/'data'; N=sic.N


def Fscaled(w,rho,protocol,scale):
    y=np.exp(w[:N]); p=float(np.exp(w[-1]))
    r=sic.rhs_input(y,rho,p,protocol)/scale
    for s,k in sic.core.REPL.items():
        r[sic.IDX[s]]=(sum(y[sic.IDX[n]] for n in sic.core.INV[k])-sic.core.target[k])/sic.core.target[k]
    return r


def Jscaled(w,rho,protocol,scale):
    y=np.exp(w[:N]); p=float(np.exp(w[-1]))
    Jx,Jp=sic.equilibrium_log_jacobian(y,rho,p,protocol,scale=scale)
    return np.c_[Jx,Jp]


def tangent(w,rho,protocol,orient=None):
    scale=np.maximum(np.exp(w[:N]),1.0)
    _,_,vh=np.linalg.svd(Jscaled(w,rho,protocol,scale),full_matrices=True)
    t=vh[-1]; t/=np.linalg.norm(t)
    if orient is None:
        if t[-1]<0: t=-t
    elif np.dot(t,orient)<0: t=-t
    return t


def correct(wpred,tprev,rho,protocol,tol=2e-11,max_nfev=80):
    scale=np.maximum(np.exp(wpred[:N]),1.0)
    def fun(w): return np.r_[Fscaled(w,rho,protocol,scale),np.dot(w-wpred,tprev)]
    def jac(w): return np.vstack([Jscaled(w,rho,protocol,scale),tprev])
    q=least_squares(fun,wpred,jac=jac,xtol=tol,ftol=tol,gtol=tol,max_nfev=max_nfev,x_scale='jac')
    rr=Fscaled(q.x,rho,protocol,np.maximum(np.exp(q.x[:N]),1.0))
    return q.x,float(np.max(np.abs(rr))),int(q.nfev),q.success


def record(w,rho,protocol,s,ds,nit,t,with_gain=False):
    y=np.exp(w[:N]); p=float(np.exp(w[-1])); lam=sic.dominant_nonconservation_eigenvalue(y,rho,p,protocol)
    rr=Fscaled(w,rho,protocol,np.maximum(y,1.0))
    out=dict(rho_program=rho,s=s,SAM=y[sic.IDX['sam']],met=y[sic.IDX['met']],alpha_real=lam.real,
             alpha_imag=lam.imag,resid=float(np.max(np.abs(rr))),ds=ds,nit=nit)
    if protocol=='methionine':
        out.update(mu=p,b_met=30*p,tlogmu=t[-1])
        if with_gain:
            try: out['gain_dlnSAM_dlnmu']=sic.logarithmic_sam_gain_fast(y,rho,p,protocol)
            except np.linalg.LinAlgError: out['gain_dlnSAM_dlnmu']=np.nan
    else: out.update(eta=p,tlogeta=t[-1])
    return out


def continue_branch(rho,protocol,direction='up',max_steps=500,ds0=0.05,p_stop=None,sam_stop=None,with_gain=False):
    y=sic.basal_state(rho,protocol); w=np.r_[np.log(y),0.0]; t=tangent(w,rho,protocol)
    if direction=='down': t=-t
    rows=[record(w,rho,protocol,0.0,0.0,0,t,with_gain)]; ds=ds0; arc=0.0
    attempts=0
    while len(rows)-1<max_steps:
        attempts+=1
        if attempts>max_steps*8: raise RuntimeError('too many rejected PALC attempts')
        wp=w+ds*t; wc,res,nit,ok=correct(wp,t,rho,protocol)
        if (not ok) or res>2e-7:
            ds*=0.5
            if ds<2e-4: raise RuntimeError(f'PALC corrector failed, residual={res:g}')
            continue
        arc+=float(np.linalg.norm(wc-w)); tn=tangent(wc,rho,protocol,t); w,t=wc,tn
        rows.append(record(w,rho,protocol,arc,ds,nit,t,with_gain))
        if nit<=5: ds=min(ds*1.25,0.35)
        elif nit>=15: ds=max(ds*0.7,5e-4)
        p=float(np.exp(w[-1])); sam=float(np.exp(w[sic.IDX['sam']]))
        if p_stop is not None and ((direction=='up' and p>=p_stop) or (direction=='down' and p<=p_stop)): break
        if sam_stop is not None and sam>=sam_stop: break
    return pd.DataFrame(rows)


def default_name(protocol,rho,direction):
    state='physiological' if rho==0 else 'proliferative'; stem='methionine_palc' if protocol=='methionine' else 'aainput_palc'
    if direction=='down': stem+=f'_{state}_down'
    else: stem+=f'_{state}'
    return stem+'_recomputed.csv'

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--protocol',choices=['methionine','common-aa'],required=True)
    ap.add_argument('--rho',type=float,choices=[0.0,1.0],required=True); ap.add_argument('--direction',choices=['up','down'],default='up')
    ap.add_argument('--max-steps',type=int,default=500); ap.add_argument('--ds',type=float,default=0.05); ap.add_argument('--p-stop',type=float); ap.add_argument('--sam-stop',type=float)
    ap.add_argument('--with-gain',action='store_true'); ap.add_argument('--output',type=Path)
    a=ap.parse_args(); df=continue_branch(a.rho,a.protocol,a.direction,a.max_steps,a.ds,a.p_stop,a.sam_stop,a.with_gain)
    out=a.output or (DATA/default_name(a.protocol,a.rho,a.direction)); out.parent.mkdir(parents=True,exist_ok=True); df.to_csv(out,index=False)
    print(f'wrote {len(df)} points to {out}'); print(df.tail(3).to_string(index=False))
