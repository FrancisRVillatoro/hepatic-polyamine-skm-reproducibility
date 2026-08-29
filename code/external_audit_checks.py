#!/usr/bin/env python3
"""Reproduce numerical checks added after external manuscript audit.

Prerequisites: run audit_control_sweeps.py first; canonical rho_branch.csv and
 directional_skm_robustness.csv must be present in ../data.
"""
from pathlib import Path
import importlib.util, hashlib, xml.etree.ElementTree as ET, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, root
from scipy.stats import beta

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent; DATA=ROOT/'data'
spec=importlib.util.spec_from_file_location('s3audit',HERE/'step3_core.py')
s3=importlib.util.module_from_spec(spec); spec.loader.exec_module(s3)
idx=s3.idx; N=len(s3.M.state)

# Fixed-H2O2 regulatory variant.
ctrl=pd.read_csv(DATA/'audit_control_sweeps.csv')
fh=ctrl[ctrl.name.eq('MAT_feedback_Hfixed0')].sort_values('r')
Sfh=PchipInterpolator(fh.r,fh.SAM); dSfh=Sfh.derivative()
m=abs(np.log(2.5)/np.log(.26)); q=n=3; tau=72.; rho0=.05
SAMH=float(Sfh(0)); Y0=(rho0/(1-rho0))**(1/q); A=Y0/(SAMH/float(Sfh(rho0)))**m
rhoY=lambda Y:Y**3/(1+Y**3); Hn=lambda Y:Y**3/(1+Y**3)
def fold_eval(Y):
    r=rhoY(Y); S=float(Sfh(r)); H=Hn(Y); K=Y/(A*(SAMH/S)**m)
    a=(K-1)/H; em=-m*(float(dSfh(r))/S)*3*r*(1-r); es=3*(K-1)*(1-H)/K
    return em+es-1,a,em,es,r,S
grid=np.logspace(-3,2,10000); gg=np.array([fold_eval(y)[0] for y in grid])
folds=[]
for i in np.where(gg[:-1]*gg[1:]<0)[0]:
    Y=brentq(lambda z:fold_eval(z)[0],grid[i],grid[i+1])
    _,a,em,es,r,S=fold_eval(Y)
    if a>0: folds.append(dict(a=a,Y=Y,rho=r,SAM=S,E_met=em,E_self=es))
pd.DataFrame(sorted(folds,key=lambda d:d['a'])).to_csv(DATA/'fixedH_reduced_folds.csv',index=False)

def scalarF(Y,a):
    r=rhoY(Y); return A*(SAMH/float(Sfh(r)))**m*(1+a*Hn(Y))-Y
vv=np.array([scalarF(y,3.) for y in grid])
roots=[brentq(lambda z:scalarF(z,3.),grid[i],grid[i+1]) for i in np.where(vv[:-1]*vv[1:]<0)[0]]
full=[]
# Fast constrained Newton solve in log variables for metabolic equilibria.
row_scale=np.maximum(np.asarray(s3.yref,float),1.0)
def log_jac(z,c):
    y=np.exp(z); J=s3.jac(y,c)*y[np.newaxis,:]/row_scale[:,None]
    for ss,k in s3.REPL.items():
        i=idx[ss]; J[i,:]=0
        for nm in s3.INV[k]: J[i,idx[nm]]=y[idx[nm]]/s3.target[k]
    return J
def solve_metabolic(r):
    z=np.log(s3.yref.copy())
    for rr in np.linspace(0,r,max(2,int(np.ceil(r*10))+1))[1:]:
        c=(float(rr),float(rr),0.0)
        sol=root(lambda zz:s3.residual(zz,*c),z,jac=lambda zz:log_jac(zz,c),method='hybr',options={'xtol':5e-10,'maxfev':100})
        z=sol.x
    return np.exp(z)

for Y in roots:
    r=rhoY(Y)
    x=solve_metabolic(r)
    def f47(z):
        xx=z[:N]; YY=z[N]; rr=rhoY(YY)
        return np.r_[s3.rhs(xx,rr,rr,0),(A*(SAMH/xx[idx['sam']])**m*(1+3*Hn(YY))-YY)/tau]
    z=np.r_[x,Y]; J=np.empty((N+1,N+1)); h=1e-30
    for j in range(N+1):
        zz=z.astype(complex); zz[j]+=1j*h; J[:,j]=np.imag(f47(zz))/h
    ev=np.linalg.eigvals(J); keep=np.ones(len(ev),bool); keep[np.argsort(np.abs(ev))[:3]]=False; ev=ev[keep]
    dom=ev[np.argmax(ev.real)]
    full.append(dict(a=3.,Y=Y,rho=r,SAM=x[idx['sam']],Put=x[idx['species_2']],Spd=x[idx['species_4']],Spm=x[idx['species_3']],dominant_real=dom.real,dominant_imag=dom.imag,n_unstable=int(np.sum(ev.real>1e-8))))
pd.DataFrame(full).to_csv(DATA/'fixedH_full47_a3.csv',index=False)

# Endpoint states and identity of tracked polyamine mode.
states={0.:s3.polish(s3.yref.copy(),(0,0,0))}
# Follow the canonical diagonal with the same constrained-Newton formulation.
z=np.log(states[0.])
for rr in np.linspace(0,1,11)[1:]:
    c=(float(rr),float(rr),float(rr))
    sol=root(lambda zz:s3.residual(zz,*c),z,jac=lambda zz:log_jac(zz,c),method='hybr',options={'xtol':5e-10,'maxfev':100})
    z=sol.x
states[1.]=np.exp(z)
rb=pd.read_csv(DATA/'rho_branch.csv'); parts=[]
for r in [0.,1.]:
    y=states[r]; J=s3.jac(y,(r,r,r)); row=rb.iloc[0 if r==0 else -1]
    target=complex(row.polyamine_pair_real_per_h,row.polyamine_pair_imag_per_h)
    vals,R=np.linalg.eig(J); k=np.argmin(np.abs(vals-target)); lam=vals[k]; rv=R[:,k]
    valsL,L=np.linalg.eig(J.T); l=np.argmin(np.abs(valsL-lam)); lv=L[:,l]
    pp=np.abs(lv*rv); pp/=pp.sum()
    for name,p in zip(s3.M.state,pp): parts.append(dict(rho=r,eigen_real=lam.real,eigen_imag=lam.imag,state=name,participation=p))
pd.DataFrame(parts).to_csv(DATA/'polyamine_mode_participation_2012.csv',index=False)

# Complex-step vs centered finite differences.
rows=[]
for r,y in states.items():
    Jc=s3.jac(y,(r,r,r)); Jf=np.empty_like(Jc)
    for j in range(N):
        hj=2e-6*max(abs(y[j]),1.); yp=y.copy(); ym=y.copy(); yp[j]+=hj; ym[j]-=hj
        Jf[:,j]=(s3.rhs(yp,r,r,r)-s3.rhs(ym,r,r,r))/(2*hj)
    rows.append(dict(rho=r,frobenius_relative=np.linalg.norm(Jc-Jf)/np.linalg.norm(Jc),max_scaled_difference=np.max(np.abs(Jc-Jf)/np.maximum(1,np.abs(Jc))),max_absolute_difference=np.max(np.abs(Jc-Jf))))
pd.DataFrame(rows).to_csv(DATA/'jacobian_complexstep_vs_finite_difference.csv',index=False)

# Proper XML tag inventory for non-holomorphic constructs.
from collections import Counter
nh=[]
for fn in ['BIOMD0000000674_url.xml','BIOMD0000000450_url.xml']:
    c=Counter(el.tag.split('}')[-1] for el in ET.parse(HERE/fn).getroot().iter())
    nh.append(dict(file=fn,floor=c['floor'],piecewise=c['piecewise'],abs=c['abs'],min=c['min'],max=c['max']))
pd.DataFrame(nh).to_csv(DATA/'mathml_nonholomorphic_constructs.csv',index=False)

# SBML identities and hashes.
mi=[]
for fn in ['BIOMD0000000190_url.xml','BIOMD0000000674_url.xml','BIOMD0000000450_url.xml']:
    p=HERE/fn; txt=p.read_text(errors='ignore'); root=ET.parse(p).getroot(); uri=root.tag.split('}')[0].strip('{'); model=root.find(f'{{{uri}}}model')
    created=re.findall(r'<dcterms:created[^>]*>.*?<dcterms:W3CDTF>(.*?)</dcterms:W3CDTF>',txt,re.S)
    modified=re.findall(r'<dcterms:modified[^>]*>.*?<dcterms:W3CDTF>(.*?)</dcterms:W3CDTF>',txt,re.S)
    mi.append(dict(file=fn,model_id=model.attrib.get('id',''),model_name=model.attrib.get('name',''),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),metadata_created=';'.join(created[:3]),metadata_modified=';'.join(modified[:3])))
pd.DataFrame(mi).to_csv(DATA/'sbml_identity_checksums.csv',index=False)

# Exact binomial intervals for directional structural ensembles.
rob=pd.read_csv(DATA/'directional_skm_robustness.csv'); ci=[]
for _,r in rob.iterrows():
    nn=int(r.n); kk=int(round(r.p_stable*nn)); alpha=.05
    lo=0 if kk==0 else beta.ppf(alpha/2,kk,nn-kk+1); hi=1 if kk==nn else beta.ppf(1-alpha/2,kk+1,nn-kk)
    ci.append(dict(endpoint=r.endpoint,factor=r.factor,n=nn,stable_count=kk,unstable_count=nn-kk,stable_fraction=kk/nn,cp95_low=lo,cp95_high=hi))
pd.DataFrame(ci).to_csv(DATA/'directional_skm_robustness_counts_ci.csv',index=False)

# Audit the N1C hypothesis: count kinetic dependencies from species/states
# with zero stoichiometric coefficient in the corresponding reaction.
y0=states[0.]
ov={'Vm_MAT1':260,'Vm_MAT3':220,'parameter_19':0,'H2O2':.01,'fasting':1.,'breakfast':1.,'lunch':1.,'dinner':1.}
n1c=[]; h=1e-30
for rx in s3.M.rx:
    stoich={nm for nm in rx['st'] if nm in idx}
    for i,nm in enumerate(s3.M.state):
        z=y0.astype(complex); z[i]+=1j*h
        dv=np.imag(s3.M.flux(z,0,ov)[rx['id']])/h
        if abs(dv)>1e-10 and nm not in stoich:
            n1c.append(dict(reaction=rx['id'],state=nm,dv_dx=dv))
pd.DataFrame(n1c).to_csv(DATA/'nonstoichiometric_regulatory_dependencies.csv',index=False)

# Figures for the capacity and timescale controls.
FIG=ROOT/'figures'
FIG.mkdir(parents=True, exist_ok=True)
orig=ctrl[ctrl.name.eq('original_MAT_only')].sort_values('r')
ct=ctrl[ctrl.name.eq('constant_total_MAT_only')].sort_values('r')
fig,ax=plt.subplots(figsize=(6.4,4.3))
ax.plot(orig.r,orig.SAM,label='published MAT reparameterization')
ax.plot(ct.r,ct.SAM,label='constant total nominal MAT capacity')
ax.set_xlabel(r'MAT program coordinate $r_M$'); ax.set_ylabel(r'SAM ($\mu$M)'); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(FIG/'fig_mat_capacity_control.pdf',bbox_inches='tight'); plt.close(fig)
if (DATA/'tauY_branch_summary.csv').exists():
    td=pd.read_csv(DATA/'tauY_branch_summary.csv')
    fig,ax=plt.subplots(figsize=(6.4,4.3))
    ax.plot(td.tau,-td.max_dom_stable,marker='o',label=r'stable: $-\Re\lambda_{\max}$')
    ax.plot(td.tau,td.min_dom_unstable,marker='o',label=r'middle: $\Re\lambda_{\max}$')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel(r'$\tau_Y$ (h)')
    ax.set_ylabel(r'Critical eigenvalue distance from zero (h$^{-1}$)'); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(FIG/'fig_tauY_stability.pdf',bbox_inches='tight'); plt.close(fig)
print('external audit checks complete')
