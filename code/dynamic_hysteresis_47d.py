#!/usr/bin/env python3
"""Full 47D dynamic hysteresis sweeps for the data-informed illustrative model."""
from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location("s3", HERE/"step3_core.py")
s3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s3)

idx=s3.idx
N=len(s3.M.state)
rho_table=pd.read_csv(DATA/"rho_branch.csv")
Sfun=PchipInterpolator(rho_table.rho,rho_table.SAM_uM)

SAMH=float(Sfun(0.0))
m=abs(np.log(2.5)/np.log(0.26))
q=n=3
tauY=72.0
rho0=.05
Y0=(rho0/(1-rho0))**(1/q)
A=Y0/(SAMH/float(Sfun(rho0)))**m

amin,amax=2.3,3.5
fold_up=3.108646018
fold_down=2.723623028

def rhoY(Y): return Y**3/(1+Y**3)
def Hn(Y): return Y**3/(1+Y**3)

def rhs47(t,z,a):
    x=z[:N]; Y=z[N]; r=rhoY(Y)
    dx=s3.rhs(x,r,r,r)
    dY=(A*(SAMH/x[idx["sam"]])**m*(1+a*Hn(Y))-Y)/tauY
    return np.r_[dx,dY]

# Endpoint equilibria.
lo=solve_ivp(lambda t,z:rhs47(t,z,amin),(0,2500),np.r_[s3.yref,.45],
             method="BDF",rtol=1e-8,atol=1e-10,max_step=20)
zlo=lo.y[:,-1]

x=s3.yref.copy()
for r in np.linspace(0,1,21)[1:]:
    x=s3.polish(x,(r,r,r))
hi=solve_ivp(lambda t,z:rhs47(t,z,amax),(0,3500),np.r_[x,2.2],
             method="BDF",rtol=1e-8,atol=1e-10,max_step=20)
zhi=hi.y[:,-1]

def sweep(T,direction,z0):
    def aa(t):
        if direction=="up":
            return amin+(amax-amin)*t/T
        return amax-(amax-amin)*t/T
    def f(t,z):
        return rhs47(t,z,aa(t))
    sol=solve_ivp(f,(0,T),z0,method="BDF",rtol=2e-7,atol=1e-9,
                  max_step=max(20,T/2500),dense_output=True)
    tt=np.linspace(0,T,5001)
    zz=sol.sol(tt)
    a=np.array([aa(t) for t in tt])
    Y=zz[N]; SAM=zz[idx["sam"]]
    dY=(A*(SAMH/SAM)**m*(1+a*Hn(Y))-Y)/tauY

    if direction=="up":
        hit=np.where((Y[:-1]<1)&(Y[1:]>=1))[0]
    else:
        hit=np.where((Y[:-1]>1)&(Y[1:]<=1))[0]
    ac=np.nan
    if len(hit):
        i=int(hit[0])
        w=(1-Y[i])/(Y[i+1]-Y[i])
        ac=a[i]+w*(a[i+1]-a[i])
    ip=int(np.argmax(np.abs(dY)))

    tr=pd.DataFrame({
        "t_h":tt,"a":a,"Y":Y,"rho":rhoY(Y),
        "SAM_uM":SAM,"Put_uM":zz[idx["species_2"]],
        "dYdt_per_h":dY,
    })
    return ac,a[ip],abs(dY[ip]),tr,len(sol.t)

rows=[]
slow_up=slow_down=None
for T in [5000.,10000.,20000.,40000.]:
    uc,up,ud,tu,nu=sweep(T,"up",zlo)
    dc,dp,dd,td,nd=sweep(T,"down",zhi)
    rows.append({
        "T_h":T,"rate_per_h":(amax-amin)/T,
        "a_up_Y1":uc,"a_up_peak":up,"up_delay_from_fold":uc-fold_up,
        "a_down_Y1":dc,"a_down_peak":dp,"down_delay_from_fold":fold_down-dc,
        "dynamic_loop_width_Y1":uc-dc,
        "up_solver_steps":nu,"down_solver_steps":nd,
    })
    if T==40000:
        slow_up,slow_down=tu,td

conv=pd.DataFrame(rows)
conv.to_csv(DATA/"dynamic_switching_convergence.csv",index=False)
slow_up.to_csv(DATA/"dynamic_hysteresis_up_T40000.csv",index=False)
slow_down.to_csv(DATA/"dynamic_hysteresis_down_T40000.csv",index=False)

fig,ax=plt.subplots(figsize=(6.6,4.5))
ax.plot(slow_up.a,slow_up.Y,label="increasing $a$")
ax.plot(slow_down.a,slow_down.Y,label="decreasing $a$")
ax.axvline(fold_up,linewidth=.8); ax.axvline(fold_down,linewidth=.8)
ax.set_xlabel("Self-reinforcement amplitude $a$")
ax.set_ylabel("Oncogenic/epigenetic state $Y$")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(FIG/"fig_dynamic_hysteresis_Y.pdf",bbox_inches="tight")
plt.close(fig)

fig,ax=plt.subplots(figsize=(6.6,4.5))
x=1/conv.T_h
ax.plot(x,conv.a_up_Y1,marker="o",label="up-sweep switch")
ax.plot(x,conv.a_down_Y1,marker="o",label="down-sweep switch")
ax.axhline(fold_up,linewidth=.8); ax.axhline(fold_down,linewidth=.8)
ax.set_xlabel("Inverse ramp time $1/T$ (h$^{-1}$)")
ax.set_ylabel("Dynamic switching value of $a$")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(FIG/"fig_dynamic_delay_convergence.pdf",bbox_inches="tight")
plt.close(fig)

print(conv.to_string(index=False))
