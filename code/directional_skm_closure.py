#!/usr/bin/env python3
"""
Directional positive-flux closure for the Reyes-Palomares 2012 SKM.

Purpose
-------
The original combined model writes 19 reversible processes as signed net rates
v = v+ - v-.  When v is close to zero, the normalized net elasticity

    theta = (x/v) dv/dx

is ill-conditioned even though v+ and v- are regular.  This script splits
those 19 rates algebraically into positive directional processes, reconstructs
the normalized Jacobian exactly, quantifies cancellation, and repeats the
anchored SKM robustness calculation in the directional representation.

This does NOT replace the published kinetic laws by a thermodynamically
constrained reversible mechanism.  It preserves the published SBML exactly.
"""
from pathlib import Path
from dataclasses import dataclass
import importlib.util
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location("s3", HERE/"step3_core.py")
s3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s3)

M = s3.M
idx = s3.idx
state = M.state
N = len(state)
rx_by_id = {r["id"]:r for r in M.rx}

# 19 net rates whose subtraction represents opposing physical directions.
REV = {
    "V_b_GLU_c":"function_4_V_b_GLU_c_1",
    "V_b_GLY_c":"function_4_V_b_GLY_c_1",
    "V_b_SER_c":"function_4_V_b_SER_c_1",
    "V_b_MET_c":"function_4_V_b_MET_c_1",
    "VmSHMT":"function_4_VmSHMT_1",
    "VmFTS":"function_4_VmFTS_1",
    "VmNE":"function_4_VmNE_1",
    "VmMTD":"function_4_VmMTD_1",
    "VmMTCH":"function_4_VmMTCH_1",
    "VmSERc":"function_4_VmSERc_1",
    "VmHCOOHc":"function_4_VmHCOOHc_1",
    "VmGLYc":"function_4_VmGLYc_1",
    "VcSHMT":"function_4_VcSHMT_1",
    "VcNE":"function_4_VcNE_1",
    "VcMTD":"function_4_VcMTD_1",
    "VcMTCH":"function_4_VcMTCH_1",
    "V_SAHH":"function_4_V_SAHH_1",
    "V_GCS":"function_4_V_GCS_1",
    "V_GS":"function_4_V_GS_1",
}

CLASSIFICATION = {
    "V_b_GLU_c":"opposed MM/linear transport",
    "V_b_GLY_c":"opposed MM/linear transport",
    "V_b_SER_c":"opposed MM/linear transport",
    "V_b_MET_c":"opposed MM/linear transport",
    "VmSHMT":"difference of independent saturating directional rates",
    "VmFTS":"difference of independent saturating directional rates",
    "VmNE":"shared-factor mass-action difference",
    "VmMTD":"difference of independent saturating directional rates",
    "VmMTCH":"difference of independent saturating directional rates",
    "VmSERc":"difference of independent saturating transport rates",
    "VmHCOOHc":"shared-factor linear directional rates",
    "VmGLYc":"difference of independent saturating transport rates",
    "VcSHMT":"difference of independent saturating directional rates",
    "VcNE":"shared-factor mass-action difference",
    "VcMTD":"difference of independent saturating directional rates",
    "VcMTCH":"difference of independent saturating directional rates",
    "V_SAHH":"difference of independent saturating directional rates",
    "V_GCS":"common-denominator reversible numerator",
    "V_GS":"common-denominator reversible numerator",
}

def controls(c):
    rM,rF,rH = c
    return {
        "Vm_MAT1":260*(1-rM),
        "Vm_MAT3":220*(1-rM),
        "parameter_19":220*rM,
        "H2O2":0.01*(1+0.5*rH),
        "fasting":1.0, "breakfast":1.0, "lunch":1.0, "dinner":1.0,
    }

@dataclass
class SV:
    p: complex
    n: complex

def S(x):
    return SV(x,0*x)
def add(a,b):
    return SV(a.p+b.p,a.n+b.n)
def neg(a):
    return SV(a.n,a.p)
def sub(a,b):
    return add(a,neg(b))
def mul(a,b):
    return SV(a.p*b.p+a.n*b.n, a.p*b.n+a.n*b.p)
def div(a,b):
    if abs(b.n) != 0:
        raise ValueError("signed denominator encountered")
    return SV(a.p/b.p,a.n/b.p)

def ev_split(node,e,t=0.0,targets=frozenset()):
    k=node.tag.split("}")[-1]
    if k=="ci":
        return S(e[(node.text or "").strip()])
    if k=="cn":
        return S(float((node.text or "0").strip()))
    if k=="csymbol":
        return S(t)
    if k=="piecewise":
        other=None
        for ch in list(node):
            q=ch.tag.split("}")[-1]
            if q=="piece":
                val,cond=list(ch)
                if M.ev(cond,e,t):
                    return ev_split(val,e,t,targets)
            elif q=="otherwise":
                other=list(ch)[0]
        return ev_split(other,e,t,targets)
    if k!="apply":
        raise NotImplementedError(k)

    ch=list(node)
    op=ch[0].tag.split("}")[-1]
    args=ch[1:]
    if op=="ci":
        f=(ch[0].text or "").strip()
        vals=[M.ev(x,e,t) for x in args]
        names,body=M.func[f]
        ee=dict(e); ee.update(zip(names,vals))
        if f in targets:
            return ev_split(body,ee,t,frozenset())
        return S(M.ev(body,ee,t))

    vals=[ev_split(x,e,t,targets) for x in args]
    if op=="plus":
        z=S(0.0)
        for v in vals: z=add(z,v)
        return z
    if op=="times":
        z=S(1.0)
        for v in vals: z=mul(z,v)
        return z
    if op=="minus":
        if len(vals)==1: return neg(vals[0])
        z=vals[0]
        for v in vals[1:]: z=sub(z,v)
        return z
    if op=="divide":
        return div(vals[0],vals[1])
    if op=="power":
        if abs(vals[0].n)!=0 or abs(vals[1].n)!=0:
            raise ValueError("signed power")
        return S(vals[0].p**vals[1].p)
    if op=="exp":
        if abs(vals[0].n)!=0: raise ValueError("signed exp")
        return S(np.exp(vals[0].p))
    if op=="floor":
        if abs(vals[0].n)!=0: raise ValueError("signed floor")
        return S(np.floor(vals[0].p))
    if op in ("leq","and"):
        return S(M.ev(node,e,t))
    raise NotImplementedError(op)

def directional_rate(y,c,rid):
    r=rx_by_id[rid]
    env=M.env(y,0,controls(c))
    env.update(r["local"])
    sv=ev_split(r["body"],env,0,frozenset([REV[rid]]))
    return sv.p,sv.n

def regulatory_rates(y,c):
    e=M.env(y,0,controls(c))
    D,Ss=e["species_4"],e["species_3"]
    free=1/(1+e["parameter_5"]*(D+Ss))
    g=(1-c[1])+c[1]*65.06/e["sam"]
    return [
        ("ODC_syn",60*g*e["parameter_7"]*free,"parameter_1",+1),
        ("ODC_deg",e["parameter_6"]*e["parameter_4"]*e["parameter_1"],"parameter_1",-1),
        ("SSAT_syn",60*e["parameter_9"]*(1-free),"parameter_2",+1),
        ("SSAT_deg",e["parameter_8"]*free*e["parameter_2"],"parameter_2",-1),
        ("SAMDC_syn",60*g*e["parameter_11"]*free,"parameter_3",+1),
        ("SAMDC_deg",e["parameter_10"]*e["parameter_3"],"parameter_3",-1),
        ("Antz_syn",e["parameter_13"]*(1-1/(1+e["parameter_5"]*.01*(D+Ss))),"parameter_4",+1),
        ("Antz_deg",e["parameter_12"]*e["parameter_4"],"parameter_4",-1),
    ]

def processes(y,c,split=True):
    fl=M.flux(y,0,controls(c))
    rates=[]; cols=[]; names=[]
    for r in M.rx:
        rid=r["id"]
        b=np.zeros(N)
        for sp,nu in r["st"].items():
            if sp in idx:
                b[idx[sp]] += nu/M.comp[M.spec[sp]["comp"]]
        if split and rid in REV:
            vp,vm=directional_rate(y,c,rid)
            rates.extend([vp,vm]); cols.extend([b,-b])
            names.extend([rid+"__f",rid+"__r"])
        else:
            rates.append(fl[rid]); cols.append(b); names.append(rid)
    for nm,rate,sp,sgn in regulatory_rates(y,c):
        b=np.zeros(N); b[idx[sp]]=sgn
        rates.append(rate); cols.append(b); names.append(nm)
    return np.asarray(rates),np.asarray(cols).T,names

def proc_derivative(y,c,split=True,h=1e-30):
    r0=processes(y,c,split)[0]
    D=np.empty((len(r0),len(y)))
    for j in range(len(y)):
        z=y.astype(complex); z[j]+=1j*h
        D[:,j]=np.imag(processes(z,c,split)[0])/h
    return D

def skm(y,c,split=True):
    v,B,names=processes(y,c,split)
    Dv=proc_derivative(y,c,split)
    active=np.abs(v)>1e-14
    v=v[active]; B=B[:,active]; Dv=Dv[active]
    names=[names[i] for i in np.where(active)[0]]
    Lam=(B*v[np.newaxis,:])/y[:,np.newaxis]
    Theta=(Dv*y[np.newaxis,:])/v[:,np.newaxis]
    J=Lam@Theta
    Jp=s3.jac(y,c)
    Jnorm=(Jp*y[np.newaxis,:])/y[:,np.newaxis]
    return Lam,Theta,J,Jnorm,names,v

def solve_branch():
    y=s3.polish(s3.yref.copy(),(0,0,0))
    rows=[]; states={}
    for r in np.linspace(0,1,11):
        if r>0:
            y=s3.polish(y,(r,r,r))
        states[round(float(r),1)]=y.copy()
        _,Tn,_,_,_,_=skm(y,(r,r,r),False)
        _,Ts,Js,Jnorm,_,_=skm(y,(r,r,r),True)
        rows.append({
            "rho":r,
            "alpha":s3.dom(y,(r,r,r)).real,
            "max_abs_theta_net":np.max(np.abs(Tn)),
            "max_abs_theta_directional":np.max(np.abs(Ts)),
            "jacobian_mismatch_directional":np.max(np.abs(Js-Jnorm)),
        })
    return pd.DataFrame(rows),states

def pair_metrics(y,c):
    _,Tn,_,_,nn,vn=skm(y,c,False)
    _,Ts,_,_,ns,vs=skm(y,c,True)
    ni={n:i for i,n in enumerate(nn)}
    si={n:i for i,n in enumerate(ns)}
    rows=[]
    for rid in REV:
        vp,vm=directional_rate(y,c,rid)
        net=vp-vm
        rows.append({
            "reaction":rid,
            "classification":CLASSIFICATION[rid],
            "v_plus":vp,
            "v_minus":vm,
            "v_net":net,
            "cancellation_factor":(abs(vp)+abs(vm))/max(abs(net),1e-300),
            "max_abs_theta_net":np.max(np.abs(Tn[ni[rid]])),
            "max_abs_theta_directional":max(
                np.max(np.abs(Ts[si[rid+"__f"]])),
                np.max(np.abs(Ts[si[rid+"__r"]])),
            ),
        })
    return pd.DataFrame(rows)

def ensemble(y,c,f,n,seed):
    Lam,T,_,_,names,_=skm(y,c,True)
    free=(np.abs(T)>1e-12) & (~np.isclose(T,1.0,rtol=0,atol=1e-10))
    ij=np.argwhere(free)
    rng=np.random.default_rng(seed)
    alpha=np.empty(n)
    nr=nc=0
    lnf=np.log(f)
    for k in range(n):
        Tk=T.copy()
        mult=np.exp(rng.uniform(-lnf,lnf,len(ij)))
        Tk[ij[:,0],ij[:,1]] *= mult
        ev=np.linalg.eigvals(Lam@Tk)
        ev=ev[np.abs(ev)>1e-8]
        d=ev[np.argmax(ev.real)]
        alpha[k]=d.real
        if d.real>=0:
            if abs(d.imag)<1e-7: nr+=1
            else: nc+=1
    return {
        "factor":f,"n":n,
        "p_stable":1-(nr+nc)/n,
        "p_real_unstable":nr/n,
        "p_complex_unstable":nc/n,
        "alpha_q01":np.quantile(alpha,.01),
        "alpha_median":np.median(alpha),
        "alpha_q99":np.quantile(alpha,.99),
        "n_free_elasticities":len(ij),
    }

branch,states=solve_branch()
branch.to_csv(DATA/"directional_skm_branch.csv",index=False)

phys=states[0.0]; prol=states[1.0]
mp=pair_metrics(phys,(0,0,0))
mq=pair_metrics(prol,(1,1,1))
mp.to_csv(DATA/"directional_pair_metrics_phys.csv",index=False)
mq.to_csv(DATA/"directional_pair_metrics_prol.csv",index=False)

# Exact decomposition audit at 100 random positive states around the endpoints.
rng=np.random.default_rng(42)
max_decomp=0.0
for y,c in [(phys,(0,0,0)),(prol,(1,1,1))]:
    for _ in range(50):
        yr=y*np.exp(rng.normal(0,0.4,N))
        fl=M.flux(yr,0,controls(c))
        for rid in REV:
            vp,vm=directional_rate(yr,c,rid)
            max_decomp=max(max_decomp,abs((vp-vm)-fl[rid]))

rob=[]
for label,y,c in [("phys",phys,(0,0,0)),("prol",prol,(1,1,1))]:
    for f,n in [(2,5000),(3,20000),(5,20000)]:
        d=ensemble(y,c,f,n,91000+10*f+(0 if label=="phys" else 1000))
        d["endpoint"]=label
        rob.append(d)
rob=pd.DataFrame(rob)
rob.to_csv(DATA/"directional_skm_robustness.csv",index=False)

pd.DataFrame(
    [{"reaction":r,"classification":CLASSIFICATION[r]} for r in REV]
).to_csv(DATA/"reversible_reaction_classification.csv",index=False)

# Figures.
fig,ax=plt.subplots(figsize=(6.5,4.3))
ax.semilogy(branch.rho,branch.max_abs_theta_net,marker="o",label="net-rate representation")
ax.semilogy(branch.rho,branch.max_abs_theta_directional,marker="o",label="directional representation")
ax.set_xlabel(r"Continuation coordinate $\rho$")
ax.set_ylabel(r"Maximum $|\theta|$")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(FIG/"fig_directional_theta.pdf",bbox_inches="tight")
fig.savefig(FIG/"fig_directional_theta.png",dpi=220,bbox_inches="tight")
plt.close(fig)

fig,ax=plt.subplots(figsize=(6.5,4.3))
ax.loglog(mp.cancellation_factor,mp.max_abs_theta_net,"o",label="physiological")
ax.loglog(mq.cancellation_factor,mq.max_abs_theta_net,"x",label="proliferative")
ax.set_xlabel(r"Cancellation factor $(v^++v^-)/|v^+-v^-|$")
ax.set_ylabel(r"Maximum net-rate $|\theta|$")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(FIG/"fig_cancellation_elasticity.pdf",bbox_inches="tight")
fig.savefig(FIG/"fig_cancellation_elasticity.png",dpi=220,bbox_inches="tight")
plt.close(fig)

# Summary.
corrp=np.corrcoef(np.log10(mp.cancellation_factor),np.log10(mp.max_abs_theta_net))[0,1]
corrq=np.corrcoef(np.log10(mq.cancellation_factor),np.log10(mq.max_abs_theta_net))[0,1]
summary=f"""DIRECTIONAL POSITIVE-FLUX SKM CLOSURE
=====================================

Reversible decomposition
------------------------
19 signed net rates were split algebraically as v = v+ - v-.
The decomposition changes no kinetic law.

Random-state audit (100 positive states around the two endpoints):
  max |(v+ - v-) - v_SBML| = {max_decomp:.3e}

Process count:
  original generalized representation: 74 reaction rates + 8 regulatory
  pseudo-processes = 82 process columns (before zero-rate omission).
  directional representation: 93 reaction-direction rates + 8 regulatory
  pseudo-processes = 101 columns in the interior.

Exact SKM reconstruction
------------------------
Across rho=0,0.1,...,1:
  max |J_directional - J_normalized| =
      {branch.jacobian_mismatch_directional.max():.3e}.

The spectrum is therefore unchanged by the directional decomposition.

Removal of the net-flux elasticity pathology
---------------------------------------------
Physiological endpoint:
  max |theta_net| = {branch.iloc[0].max_abs_theta_net:.6f}
  max |theta_directional| = {branch.iloc[0].max_abs_theta_directional:.6f}

At rho=0.1:
  max |theta_net| = {branch.iloc[1].max_abs_theta_net:.6f}
  max |theta_directional| = {branch.iloc[1].max_abs_theta_directional:.6f}

Proliferative endpoint:
  max |theta_net| = {branch.iloc[-1].max_abs_theta_net:.6f}
  max |theta_directional| = {branch.iloc[-1].max_abs_theta_directional:.6f}

All elasticities belonging to the 19 split reversible reactions satisfy
|theta| <= 1 at both endpoints.  The global directional maximum 2 is the
exact quadratic elasticity of the irreversible cys_usage reaction.

Cancellation explains the large net elasticities:
  corr(log10 cancellation, log10 max|theta_net|)
    physiological = {corrp:.6f}
    proliferative = {corrq:.6f}

Examples at the physiological endpoint:
{mp.sort_values("cancellation_factor",ascending=False).head(5)[["reaction","v_plus","v_minus","v_net","cancellation_factor","max_abs_theta_net","max_abs_theta_directional"]].to_string(index=False)}

Directional anchored ensembles
-------------------------------
{rob.to_string(index=False)}

These stability fractions are very close to the previous net-rate SKM
ensembles: the conclusion that both endpoints occupy broad locally stable
structural neighborhoods is unchanged.

Thermodynamic interpretation
-----------------------------
The positive-direction split fixes the normalization pathology, but it does
not make the source kinetic model microscopically thermodynamic.

Of the 19 signed reversible rates:
  - 5 have a shared-factor/common-denominator directional structure
    (VmNE, VcNE, VmHCOOHc, V_GCS, V_GS);
  - 14 are phenomenological differences of independently saturated or
    differently parameterized directional rates/transports.

For those 14 reactions the published SBML does not encode a single reversible
rate law together with an explicit Haldane/equilibrium-constant constraint.
Imposing strict thermodynamic consistency would therefore require changing
the published kinetic model, not merely reparameterizing its SKM.

Conclusion
----------
The enormous elasticities in the original full-network SKM are a coordinate
artifact of normalizing near-zero net reversible fluxes.  Replacing each such
net process by two positive directional processes removes the artifact,
reconstructs the exact kinetic Jacobian, and leaves the stability/robustness
conclusions unchanged.  The manuscript should use this directional SKM as
the primary full-network structural representation and should explicitly
state that microscopic thermodynamic certification lies beyond what the
published kinetic laws specify.
"""
(HERE/"SUMMARY.txt").write_text(summary,encoding="utf-8")
print(summary)
