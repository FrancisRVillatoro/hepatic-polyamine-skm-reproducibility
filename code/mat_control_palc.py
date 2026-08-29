#!/usr/bin/env python3
"""Pseudo-arclength continuation in methionine input for MAT control policies."""
from pathlib import Path
import argparse, numpy as np, pandas as pd
from scipy.optimize import root
import mat_control_gain_sweeps as mc
import step3_core as s3
HERE=Path(__file__).resolve().parent; DATA=HERE.parent/'data'; N=mc.N; idx=mc.idx

def F(w,rho,policy,protocol='methionine'):
    z=w[:N]; p=np.exp(w[-1]); return mc.fun(z,rho,p,protocol,policy)

def J(w,rho,policy,protocol='methionine',h=1e-30):
    z=w[:N]; y=np.exp(z); p=float(np.exp(w[-1])); Jx=mc.jaclog(z,rho,p,protocol,policy)
    pc=complex(p,p*h); rp=mc.rhs(y.astype(complex),rho,pc,protocol,policy); Jp=np.imag(rp)/h/mc.row_scale
    for ss in s3.REPL:Jp[idx[ss]]=0
    return np.c_[Jx,Jp]

def tangent(w,rho,policy,orient=None):
    _,_,vh=np.linalg.svd(J(w,rho,policy),full_matrices=True);t=vh[-1];t/=np.linalg.norm(t)
    if orient is None:
        if t[-1]<0:t=-t
    elif np.dot(t,orient)<0:t=-t
    return t

def correct(wpred,tprev,rho,policy):
    def ff(w):return np.r_[F(w,rho,policy),np.dot(w-wpred,tprev)]
    def jj(w):return np.vstack([J(w,rho,policy),tprev])
    q=root(ff,wpred,jac=jj,method='hybr',options={'xtol':2e-9,'maxfev':120});res=np.max(np.abs(F(q.x,rho,policy)))
    if res>2e-7:
        q=root(ff,q.x,jac=jj,method='lm',options={'ftol':1e-10,'xtol':1e-10,'gtol':1e-10,'maxiter':180});res=np.max(np.abs(F(q.x,rho,policy)))
    return q.x,float(res),bool(q.success)

def initial_state(policy,rho=1.0):
    f=DATA/f'mat_control_gains_{policy}_fullstate.csv';d=pd.read_csv(f);r=d.iloc[(d.rho-rho).abs().argmin()];y=np.array([r['x_'+n] for n in mc.M.state],float);return y

def rec(w,rho,policy,s,ds,t):
    y=np.exp(w[:N]);p=float(np.exp(w[-1]));aval=mc.alpha(y,rho,policy)
    return dict(policy=policy,rho=rho,s=s,mu=p,b_met=30*p,SAM=y[idx['sam']],met=y[idx['met']],tlogmu=t[-1],alpha=aval,residual=float(np.max(np.abs(F(w,rho,policy)))),ds=ds)

def run(policy,rho=1.0,max_steps=240,ds0=0.05,p_stop=15.0,sam_stop=1e9):
    y=initial_state(policy,rho);w=np.r_[np.log(y),0.0];t=tangent(w,rho,policy);rows=[rec(w,rho,policy,0,0,t)];ds=ds0;s=0.
    for _ in range(max_steps):
        wp=w+ds*t;wc,res,ok=correct(wp,t,rho,policy)
        if (not ok) or res>2e-7:
            ds*=.5
            if ds<2e-4:break
            continue
        s+=float(np.linalg.norm(wc-w));tn=tangent(wc,rho,policy,t);w,t=wc,tn;rows.append(rec(w,rho,policy,s,ds,t))
        if res<1e-9: ds=min(ds*1.18,.25)
        p=np.exp(w[-1]);sam=np.exp(w[idx['sam']])
        print(policy,len(rows)-1,f'mu={p:.8g}',f'SAM={sam:.6g}',f'tmu={t[-1]:.3e}',f'a={rows[-1]["alpha"]:.3e}',f'res={res:.1e}',flush=True)
        if p>=p_stop or sam>=sam_stop:break
    return pd.DataFrame(rows)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--policy',choices=['published','constant-total','flux-matched'],required=True);ap.add_argument('--rho',type=float,default=1.0);ap.add_argument('--max-steps',type=int,default=240);ap.add_argument('--ds',type=float,default=.05);ap.add_argument('--p-stop',type=float,default=15.0);ap.add_argument('--sam-stop',type=float,default=1e9);ap.add_argument('--output',type=Path)
    a=ap.parse_args();d=run(a.policy,a.rho,a.max_steps,a.ds,a.p_stop,a.sam_stop);out=a.output or DATA/f'mat_control_methionine_palc_{a.policy}.csv';d.to_csv(out,index=False);print('wrote',out,'rows',len(d));print(d.tail().to_string(index=False))
