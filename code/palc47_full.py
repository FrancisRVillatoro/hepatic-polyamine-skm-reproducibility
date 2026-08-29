#!/usr/bin/env python3
"""
Full 47-dimensional pseudo-arclength continuation of the data-informed illustrative
polyamine oncogenic extension.

State:
  46 autonomous metabolic variables from the Reyes-Palomares 2012 common
  network + one slow variable Y.

Continuation variable:
  a, the self-reinforcement amplitude.

Three redundant metabolic balances are replaced by the exact conserved pools
(cytosolic folate, mitochondrial folate, CoA+acetyl-CoA), leaving 47
independent equilibrium equations. A pseudo-arclength constraint closes the
48-dimensional corrector for (x, Y, a).

Primary calibration:
  q=n=3,
  m = abs(log(2.5)/log(0.26)),
  tau_Y = 72 h (tau_Y does not affect equilibrium locations).
"""
from pathlib import Path
import importlib.util, math
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.interpolate import PchipInterpolator, CubicSpline
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
spec = importlib.util.spec_from_file_location("s3", HERE/"step3_core.py")
s3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s3)

idx = s3.idx
N = len(s3.M.state)
branch = pd.read_csv(DATA/"rho_branch.csv")
Sfun = PchipInterpolator(branch.rho.values, branch.SAM_uM.values)

SAMH = float(Sfun(0.0))
m = abs(np.log(2.5)/np.log(0.26))
q = n = 3
tau = 72.0
rho_low = 0.05
YL = (rho_low/(1-rho_low))**(1/q)
A = YL/(SAMH/float(Sfun(rho_low)))**m

cons_rows = {
    idx["c_thf"]:"cfol",
    idx["m_thf"]:"mfol",
    idx["species_9"]:"coa",
}
row_scale = np.maximum(np.asarray(s3.yref,float),1.0)

def rhoY(Y):
    return Y**q/(1+Y**q)

def HY(Y):
    return Y**n/(1+Y**n)

def equilibrium_residual(u):
    """47 independent equations in u=(log x, log Y, log a)."""
    x = np.exp(u[:N])
    Y = np.exp(u[N])
    a = np.exp(u[N+1])
    r = rhoY(Y)
    rr = s3.rhs(x,r,r,r)
    out = np.empty(N+1,dtype=np.result_type(u,rr))
    for i in range(N):
        if i in cons_rows:
            key = cons_rows[i]
            out[i] = (
                sum(x[idx[name]] for name in s3.INV[key]) - s3.target[key]
            )/s3.target[key]
        else:
            out[i] = rr[i]/row_scale[i]
    out[N] = A*(SAMH/x[idx["sam"]])**m*(1+a*HY(Y))-Y
    return out

def jacobian_equilibrium(u,h=1e-30):
    J = np.empty((N+1,N+2))
    for j in range(N+2):
        z = u.astype(complex)
        z[j] += 1j*h
        J[:,j] = np.imag(equilibrium_residual(z))/h
    return J

def tangent(u, previous=None):
    J = jacobian_equilibrium(u)
    _,_,vh = np.linalg.svd(J,full_matrices=True)
    t = vh[-1]
    t /= np.linalg.norm(t)
    if previous is not None and np.dot(t,previous)<0:
        t = -t
    return t

def fixed_a_equilibrium(a,xguess=None,Yguess=None):
    if xguess is None:
        xguess = s3.yref.copy()
    if Yguess is None:
        Yguess = YL
    b = math.log(a)
    def fun(v):
        return equilibrium_residual(np.r_[v,b])
    v0 = np.r_[np.log(xguess),math.log(Yguess)]
    sol = least_squares(
        fun,v0,jac="2-point",xtol=2e-12,ftol=2e-12,gtol=2e-12,
        max_nfev=120,x_scale="jac"
    )
    return np.r_[sol.x,b], sol

def corrector(upred,t):
    def H(u):
        return np.r_[equilibrium_residual(u),np.dot(t,u-upred)]
    return least_squares(
        H,upred,jac="2-point",xtol=2e-11,ftol=2e-11,gtol=2e-11,
        max_nfev=30,x_scale="jac"
    )

def physical_jacobian(x,Y,a,h=1e-30):
    def rhs47(z):
        xx=z[:N]
        YY=z[N]
        r=rhoY(YY)
        dx=s3.rhs(xx,r,r,r)
        dY=(A*(SAMH/xx[idx["sam"]])**m*(1+a*HY(YY))-YY)/tau
        return np.r_[dx,dY]
    z=np.r_[x,Y]
    J=np.empty((N+1,N+1))
    for j in range(N+1):
        zz=z.astype(complex)
        zz[j]+=1j*h
        J[:,j]=np.imag(rhs47(zz))/h
    return J

def diagnostics(u):
    x=np.exp(u[:N])
    Y=float(np.exp(u[N]))
    a=float(np.exp(u[N+1]))
    r=float(rhoY(Y))
    ev=np.linalg.eigvals(physical_jacobian(x,Y,a))
    order=np.argsort(np.abs(ev))
    keep=np.ones(len(ev),bool)
    keep[order[:3]]=False   # exact conservation modes
    ev=ev[keep]
    dom=ev[np.argmax(ev.real)]
    return dict(
        a=a,Y=Y,rho=r,SAM=float(x[idx["sam"]]),
        Put=float(x[idx["species_2"]]),
        Spd=float(x[idx["species_4"]]),
        Spm=float(x[idx["species_3"]]),
        dom_real=float(dom.real),dom_imag=float(dom.imag),
        resnorm=float(np.linalg.norm(equilibrium_residual(u),np.inf)),
    )

u0,s0=fixed_a_equilibrium(2.3)
u1,s1=fixed_a_equilibrium(2.4,np.exp(u0[:N]),np.exp(u0[N]))

tprev=(u1-u0)/np.linalg.norm(u1-u0)
if tprev[-1]<0:
    tprev=-tprev

rows=[]
sarc=0.0
for u,sol in [(u0,s0),(u1,s1)]:
    d=diagnostics(u)
    rows.append({**d,"s":sarc,"tb":np.nan,"ds":np.nan,"nfev":sol.nfev})
    if u is u0:
        sarc=np.linalg.norm(u1-u0)

u=u1.copy()
ds=0.035
for k in range(260):
    t=tangent(u,tprev)
    sol=corrector(u+ds*t,t)
    uc=sol.x
    if (not sol.success) or np.linalg.norm(equilibrium_residual(uc),np.inf)>2e-7:
        ds*=0.5
        if ds<0.002:
            raise RuntimeError("Continuation step failed")
        continue

    sarc += np.linalg.norm(uc-u)
    tnew=tangent(uc,t)
    d=diagnostics(uc)
    rows.append({**d,"s":sarc,"tb":tnew[-1],"ds":ds,"nfev":sol.nfev})
    u=uc
    tprev=tnew

    if sol.nfev<=5:
        ds=min(1.15*ds,0.06)
    elif sol.nfev>=12:
        ds=max(0.7*ds,0.008)

    # after both folds, continue into the stable high branch
    signs=np.sign([r["tb"] for r in rows if np.isfinite(r["tb"])])
    nturn=np.sum(signs[:-1]*signs[1:]<0)
    if nturn>=2 and d["a"]>3.5 and tnew[-1]>0:
        break

df=pd.DataFrame(rows)
df["stability"]=np.where(df.dom_real<0,"stable","unstable")
df.to_csv(DATA/"palc47_branch.csv",index=False)

# Refine fold positions by cubic interpolation in pseudo-arclength.
turns=[]
tb=df.tb.to_numpy()
for i in range(len(df)-1):
    if np.isfinite(tb[i]) and np.isfinite(tb[i+1]) and tb[i]*tb[i+1]<0:
        turns.append(i)

fold_rows=[]
for k,i in enumerate(turns):
    sub=df.iloc[max(0,i-3):min(len(df),i+4)]
    s=sub.s.to_numpy()
    csa=CubicSpline(s,sub.a.to_numpy())
    sf=brentq(csa.derivative(),df.s.iloc[i],df.s.iloc[i+1])
    row={"fold_index":k+1,"s":sf}
    for col in ["a","Y","rho","SAM","Put","Spd","Spm","dom_real"]:
        row[col]=float(CubicSpline(s,sub[col].to_numpy())(sf))
    fold_rows.append(row)

pd.DataFrame(fold_rows).to_csv(DATA/"palc47_folds.csv",index=False)
print(pd.DataFrame(fold_rows).to_string(index=False))
