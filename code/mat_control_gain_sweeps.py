#!/usr/bin/env python3
"""Gain sweeps for published, constant-total and flux-matched MAT controls."""
from pathlib import Path
import argparse, json, time
import numpy as np, pandas as pd
from scipy.optimize import root
import step3_core as s3

HERE=Path(__file__).resolve().parent; DATA=HERE.parent/'data'
M=s3.M; idx=s3.idx; N=len(M.state); row_scale=np.maximum(np.asarray(s3.yref,float),1.0)

def phys_flux_matched_vmax():
    y=np.asarray(s3.yref,float)
    ov={'Vm_MAT1':260.,'Vm_MAT3':220.,'parameter_19':0.,'H2O2':0.01,
        'fasting':1.,'breakfast':1.,'lunch':1.,'dinner':1.,'b_met_basal':30.}
    fl=M.flux(y,0,ov); target=fl['V_MATI']+fl['V_MATIII']
    ov2=dict(ov);ov2.update({'Vm_MAT1':0.,'Vm_MAT3':0.,'parameter_19':1.})
    per=M.flux(y,0,ov2)['reaction_13']
    return float(target/per),float(target),float(per)

def v2cap(policy):
    if policy=='published': return 220.0
    if policy=='constant-total': return 480.0
    if policy=='flux-matched': return phys_flux_matched_vmax()[0]
    raise ValueError(policy)

def overrides(rho,p,protocol,policy):
    v2=v2cap(policy)
    if protocol=='methionine': fs=br=lu=di=1.; b=30*p
    else: fs=br=lu=di=p; b=30.
    return {'Vm_MAT1':260*(1-rho),'Vm_MAT3':220*(1-rho),'parameter_19':v2*rho,
            'H2O2':0.01*(1+0.5*rho),'fasting':fs,'breakfast':br,'lunch':lu,'dinner':di,'b_met_basal':b}

def rhs(y,rho,p,protocol,policy):
    ov=overrides(rho,p,protocol,policy);e=M.env(y,0,ov);fl=M.flux(y,0,ov);d={n:0.0 for n in M.state}
    for rx in M.rx:
        for sp,nu in rx['st'].items():
            if sp in d:d[sp]+=nu*fl[rx['id']]/M.comp[M.spec[sp]['comp']]
    D,S=e['species_4'],e['species_3'];free=1/(1+e['parameter_5']*(D+S));g=(1-rho)+rho*65.06/e['sam']
    d['parameter_1']=60*g*e['parameter_7']*free-e['parameter_6']*e['parameter_4']*e['parameter_1']
    d['parameter_2']=60*e['parameter_9']*(1-free)-e['parameter_8']*free*e['parameter_2']
    d['parameter_3']=60*g*e['parameter_11']*free-e['parameter_10']*e['parameter_3']
    d['parameter_4']=e['parameter_13']*(1-1/(1+e['parameter_5']*.01*(D+S)))-e['parameter_12']*e['parameter_4']
    return np.array([d[n] for n in M.state])

def jac(y,rho,p,protocol,policy,h=1e-30):
    J=np.empty((N,N))
    for j in range(N):
        z=y.astype(complex);z[j]+=1j*h;J[:,j]=np.imag(rhs(z,rho,p,protocol,policy))/h
    return J

def fun(z,rho,p,protocol,policy):
    y=np.exp(z);rr=rhs(y,rho,p,protocol,policy)/row_scale
    for ss,k in s3.REPL.items():rr[idx[ss]]=(sum(y[idx[n]] for n in s3.INV[k])-s3.target[k])/s3.target[k]
    return rr

def jaclog(z,rho,p,protocol,policy):
    y=np.exp(z);J=jac(y,rho,p,protocol,policy)*y[np.newaxis,:]/row_scale[:,None]
    for ss,k in s3.REPL.items():
        i=idx[ss];J[i,:]=0
        for n in s3.INV[k]:J[i,idx[n]]=y[idx[n]]/s3.target[k]
    return J

def solve(z0,rho,p,protocol,policy):
    args=(rho,p,protocol,policy)
    sol=root(lambda z:fun(z,*args),z0,jac=lambda z:jaclog(z,*args),method='hybr',options={'xtol':2e-10,'maxfev':100})
    z=sol.x;res=float(np.max(np.abs(fun(z,*args))))
    if res>2e-7:
        sol=root(lambda z:fun(z,*args),z,jac=lambda z:jaclog(z,*args),method='lm',options={'ftol':1e-11,'xtol':1e-11,'gtol':1e-11,'maxiter':200})
        z=sol.x;res=float(np.max(np.abs(fun(z,*args))))
    return z,res,bool(sol.success)

def loggain(y,rho,protocol,policy,p=1.0,h=1e-30):
    z=np.log(y);J=jaclog(z,rho,p,protocol,policy)
    pc=complex(p,p*h); rp=rhs(y.astype(complex),rho,pc,protocol,policy)
    Jp=np.imag(rp)/h/row_scale
    for ss in s3.REPL:Jp[idx[ss]]=0
    dz=np.linalg.solve(J,-Jp);return float(dz[idx['sam']])

def alpha(y,rho,policy):
    ev=np.linalg.eigvals(jac(y,rho,1.,'methionine',policy));ev=ev[np.abs(ev)>1e-8];return float(ev[np.argmax(ev.real)].real)

def run(policy,nrho=21):
    z=np.log(np.asarray(s3.yref,float));rows=[];full=[]
    for rho in np.linspace(0,1,nrho):
        z,res,ok=solve(z,float(rho),1.,'methionine',policy);y=np.exp(z)
        gm=loggain(y,float(rho),'methionine',policy);ga=loggain(y,float(rho),'common-aa',policy)
        a=alpha(y,float(rho),policy)
        rows.append(dict(policy=policy,rho=rho,MATII_Vmax_cap=v2cap(policy),SAM=y[idx['sam']],met=y[idx['met']],Put=y[idx['species_2']],Spd=y[idx['species_4']],Spm=y[idx['species_3']],G_Met=gm,G_AA=ga,alpha=a,residual=res,solver_ok=ok))
        fr=dict(policy=policy,rho=rho,MATII_Vmax_cap=v2cap(policy),G_Met=gm,G_AA=ga,alpha=a,residual=res)
        fr.update({f'x_{n}':y[idx[n]] for n in M.state});full.append(fr)
        print(policy,f'rho={rho:.2f}',f'SAM={y[idx["sam"]]:.6g}',f'GM={gm:.6g}',f'GAA={ga:.6g}',f'res={res:.2e}',flush=True)
    return pd.DataFrame(rows),pd.DataFrame(full)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--policy',choices=['published','constant-total','flux-matched'],required=True);ap.add_argument('--nrho',type=int,default=21);ap.add_argument('--outdir',type=Path,default=DATA)
    a=ap.parse_args();t=time.time();d,f=run(a.policy,a.nrho);a.outdir.mkdir(parents=True,exist_ok=True);d.to_csv(a.outdir/f'mat_control_gains_{a.policy}.csv',index=False);f.to_csv(a.outdir/f'mat_control_gains_{a.policy}_fullstate.csv',index=False)
    vm,flux,per=phys_flux_matched_vmax();meta={'flux_matched_MATII_Vmax':vm,'physiological_MATI_plus_MATIII_flux':flux,'MATII_flux_per_unit_Vmax_at_physiological_state':per};(a.outdir/'mat_control_flux_match_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('elapsed',time.time()-t);print(d.iloc[[0,-1]].to_string(index=False))
