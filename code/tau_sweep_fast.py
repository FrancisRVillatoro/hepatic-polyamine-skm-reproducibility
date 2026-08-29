from pathlib import Path
import sys, numpy as np, pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
import step3_core as s3
from scipy.interpolate import PchipInterpolator
HERE=Path(__file__).resolve().parent.parent
branch=pd.read_csv(HERE/'data/palc47_branch.csv')
rhodata=pd.read_csv(HERE/'data/rho_branch.csv')
Sfun=PchipInterpolator(rhodata.rho,rhodata.SAM_uM)
SAMH=float(Sfun(0.0)); m=abs(np.log(2.5)/np.log(0.26)); q=n=3; rho0=.05
Y0=(rho0/(1-rho0))**(1/q); A=Y0/(SAMH/float(Sfun(rho0)))**m
idx=s3.idx; N=len(s3.M.state)
# reconstruct x by ascending rho continuation
y=s3.yref.copy(); xs={}
for k in np.argsort(branch.rho.values):
    r=float(branch.rho.iloc[k]); y=s3.polish(y,(r,r,r)); xs[k]=y.copy()
# precompute metabolic blocks
blocks={}
h=1e-30
for k,row in branch.iterrows():
    x=xs[k]; r=float(row.rho); Y=float(row.Y)
    Jxx=s3.jac(x,(r,r,r))
    # derivative along diagonal control wrt r, then rY derivative
    rr=r+1j*h
    dr_diag=np.imag(s3.rhs(x.astype(complex),rr,rr,rr))/h
    drdY=q*r*(1-r)/Y
    JxY=dr_diag*drdY
    blocks[k]=(x,Jxx,JxY)

taus=[6.,12.,24.,48.,72.,96.,168.,336.,720.,1440.]
rows=[]; summary=[]
for tau in taus:
    hopf_cross=False; max_stable=-1e9; min_unstable=1e9; max_imag_near=0; stable_n=unstable_n=0
    worst_st=None; closest_un=None
    for k,row in branch.iterrows():
        x,Jxx,JxY=blocks[k]; Y=float(row.Y); a=float(row.a); SAM=float(x[idx['sam']])
        H=Y**n/(1+Y**n); Hp=n*Y**(n-1)/(1+Y**n)**2
        JyX=np.zeros(N); JyX[idx['sam']]=-m*Y/(tau*SAM)
        base=A*(SAMH/SAM)**m
        JyY=(base*a*Hp-1)/tau
        J=np.zeros((N+1,N+1)); J[:N,:N]=Jxx; J[:N,N]=JxY; J[N,:N]=JyX; J[N,N]=JyY
        ev=np.linalg.eigvals(J)
        # remove conservation zeros
        keep=np.ones(len(ev),bool); keep[np.argsort(np.abs(ev))[:3]]=False; ev=ev[keep]
        dom=ev[np.argmax(ev.real)]; npos=int(np.sum(ev.real>1e-8))
        rows.append(dict(tau=tau,k=k,a=a,Y=Y,rho=row.rho,dom_real=dom.real,dom_imag=dom.imag,npos=npos))
        expected=0 if row.stability=='stable' else 1
        if expected==0:
            stable_n+=1
            if dom.real>max_stable: max_stable=dom.real; worst_st=(a,Y,row.rho,dom.real,dom.imag)
        else:
            unstable_n+=1
            if dom.real<min_unstable: min_unstable=dom.real; closest_un=(a,Y,row.rho,dom.real,dom.imag)
        if abs(dom.imag)>max_imag_near: max_imag_near=abs(dom.imag)
    summary.append(dict(tau=tau,max_dom_stable=max_stable,min_dom_unstable=min_unstable,
                        worst_stable_a=worst_st[0],worst_stable_imag=worst_st[4],
                        closest_unstable_a=closest_un[0],closest_unstable_imag=closest_un[4],
                        max_abs_dom_imag=max_imag_near,stable_points=stable_n,unstable_points=unstable_n))
    print(summary[-1])
pd.DataFrame(rows).to_csv(HERE/'data/tauY_branch_sweep.csv',index=False)
pd.DataFrame(summary).to_csv(HERE/'data/tauY_branch_summary.csv',index=False)
