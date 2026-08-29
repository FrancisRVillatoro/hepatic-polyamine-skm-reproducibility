#!/usr/bin/env python3
"""Reduced slow-feedback sensitivity analysis used in Supplementary Section S9."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, root
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; DATA=ROOT/'data'; FIG=ROOT/'figures'
FIG.mkdir(parents=True, exist_ok=True)
df=pd.read_csv(DATA/'rho_branch.csv'); Sfun=PchipInterpolator(df.rho,df.SAM_uM); dS=Sfun.derivative(); SAMH=float(Sfun(0))
def saddles(m=.680208,q=3,n=3,rho0=.05):
 Y0=(rho0/(1-rho0))**(1/q); A=Y0/(SAMH/float(Sfun(rho0)))**m
 def ev(Y):
  Y=np.asarray(Y); r=Y**q/(1+Y**q); S=Sfun(r); Em=-m*(dS(r)/S)*q*r*(1-r); H=Y**n/(1+Y**n); K=Y/(A*(SAMH/S)**m); a=(K-1)/H; Es=n*(K-1)*(1-H)/K
  return Em+Es-1,a,Em,Es,r,S
 ys=np.logspace(-4,6,5000); g=ev(ys)[0]; ii=np.where(np.isfinite(g[:-1])&np.isfinite(g[1:])&(g[:-1]*g[1:]<0))[0]; out=[]
 for i in ii:
  y=brentq(lambda z:float(ev(z)[0]),ys[i],ys[i+1]); _,a,Em,Es,r,S=ev(y)
  if a>0: out.append(dict(Y=y,a=float(a),rho=float(r),SAM=float(S),E_met=float(Em),E_self=float(Es),A=A))
 return sorted(out,key=lambda z:z['a'])
def row(m=.680208,q=3,n=3,rho0=.05):
 z=saddles(m,q,n,rho0); d=dict(m=m,q=q,n=n,rho0=rho0,n_saddles=len(z),status='two_positive_saddles' if len(z)>=2 else 'no_positive_saddle_pair')
 if len(z)>=2: d.update(A=z[0]['A'],a_SN_on=z[0]['a'],a_SN_off=z[-1]['a'],Y_SN_on=z[0]['Y'],Y_SN_off=z[-1]['Y'],rho_SN_on=z[0]['rho'],rho_SN_off=z[-1]['rho'],SAM_SN_on=z[0]['SAM'],SAM_SN_off=z[-1]['SAM'],E_met_on=z[0]['E_met'],E_self_crit=z[0]['E_self'],E_met_off=z[-1]['E_met'],E_self_off=z[-1]['E_self'],bistable_width=z[-1]['a']-z[0]['a'])
 return d
m=pd.DataFrame([row(m=float(x)) for x in np.linspace(.3,1.2,46)]); m.to_csv(DATA/'m_sensitivity_q3n3.csv',index=False)
sel=pd.DataFrame([row(m=x) for x in [.5,.6,.680208,.8,1.0]]); sel.to_csv(DATA/'m_sensitivity_selected.csv',index=False)
r0=pd.DataFrame([row(rho0=float(x)) for x in np.arange(.005,.091,.001)]); r0.to_csv(DATA/'rho0_sensitivity_q3n3.csv',index=False)
qn=pd.DataFrame([row(q=q,n=n) for q in range(1,6) for n in range(2,6)]); qn.to_csv(DATA/'qn_complete_grid.csv',index=False)
# Cusp where the two reduced folds merge as rho0 is varied.
def gain(Y,rho0,m0=.680208,q=3,n=3):
 Y0=(rho0/(1-rho0))**(1/q); A=Y0/(SAMH/float(Sfun(rho0)))**m0; r=Y**q/(1+Y**q); S=float(Sfun(r)); Em=-m0*(float(dS(r))/S)*q*r*(1-r); H=Y**n/(1+Y**n); K=Y/(A*(SAMH/S)**m0); Es=n*(K-1)*(1-H)/K
 return Em+Es-1
def dg(Y,rho0):
 h=1e-5*max(1,abs(Y)); return (gain(Y+h,rho0)-gain(Y-h,rho0))/(2*h)
sol=root(lambda z:[gain(z[0],z[1]),dg(z[0],z[1])],[.88,.0818]); Yc,r0c=sol.x
Y0=(r0c/(1-r0c))**(1/3); A=Y0/(SAMH/float(Sfun(r0c)))**.680208; r=Yc**3/(1+Yc**3); S=float(Sfun(r)); H=Yc**3/(1+Yc**3); K=Yc/(A*(SAMH/S)**.680208); a=(K-1)/H; Em=-.680208*(float(dS(r))/S)*3*r*(1-r); Es=3*(K-1)*(1-H)/K
pd.DataFrame([dict(rho0_cusp=r0c,A_cusp=A,Y_cusp=Yc,a_cusp=a,rho_cusp=r,SAM_cusp=S,E_met_cusp=Em,E_self_cusp=Es)]).to_csv(DATA/'rho0_cusp_q3n3.csv',index=False)
fig,ax=plt.subplots(figsize=(6.4,4.3)); ax.plot(m.m,m.E_self_crit); ax.scatter(sel.m,sel.E_self_crit); ax.axvline(.680208,linewidth=.8); ax.set_xlabel(r'SAM-response exponent $m$'); ax.set_ylabel(r'Critical self-elasticity $E_{\rm self}^{\rm crit}$'); fig.tight_layout(); fig.savefig(FIG/'fig_m_sensitivity.pdf',bbox_inches='tight'); plt.close(fig)
ok=r0[r0.n_saddles>=2]; fig,ax=plt.subplots(figsize=(6.4,4.3)); ax.plot(ok.rho0,ok.E_self_crit); ax.axvline(.05,linewidth=.8); ax.axvline(r0c,linewidth=.8); ax.set_xlabel(r'Basal program target $\rho_0$ at $a=0$'); ax.set_ylabel(r'Critical self-elasticity $E_{\rm self}^{\rm crit}$'); fig.tight_layout(); fig.savefig(FIG/'fig_rho0_sensitivity.pdf',bbox_inches='tight'); plt.close(fig)
print(sel[['m','a_SN_on','a_SN_off','E_self_crit']].to_string(index=False)); print('rho0_cusp',r0c,'a_cusp',a,'E_self',Es)
