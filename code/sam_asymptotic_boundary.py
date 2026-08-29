#!/usr/bin/env python3
"""Conditional-equilibrium audit of the physiological high-SAM boundary."""
from pathlib import Path
import argparse, numpy as np, pandas as pd
from scipy.optimize import root, brentq, curve_fit
import substrate_input_core as sic
import step3_core as core
HERE=Path(__file__).resolve().parent;DATA=HERE.parent/'data';N=sic.N;IDX=sic.IDX;IS=IDX['sam']
KEEP_ROWS=np.array([i for i in range(N) if i!=IS]); KEEP_COLS=KEEP_ROWS.copy()

def archive_df(): return pd.read_csv(DATA/'methionine_palc_physiological_fullstate.csv')
def archive_state_near_S(S):
    d=archive_df(); r=d.iloc[(np.log(d.SAM)-np.log(S)).abs().argmin()]
    return np.array([r['x_'+n] for n in core.M.state],float), float(r.mu)

def residual_other(logu,S,mu):
    y=np.empty(N);y[IS]=S;y[KEEP_COLS]=np.exp(logu)
    scale=np.maximum(y,1.0); rr=sic.rhs_input(y,0.0,mu,'methionine')/scale
    for ss,k in core.REPL.items():
        i=IDX[ss];rr[i]=(sum(y[IDX[n]] for n in core.INV[k])-core.target[k])/core.target[k]
    return rr[KEEP_ROWS]

def jac_other(logu,S,mu):
    y=np.empty(N);y[IS]=S;y[KEEP_COLS]=np.exp(logu);scale=np.maximum(y,1.0)
    Jp=sic.physical_jacobian(y,0.0,mu,'methionine')/scale[:,None];Jlog=Jp*y[None,:]
    for ss,k in core.REPL.items():
        i=IDX[ss];Jlog[i,:]=0
        for n in core.INV[k]:Jlog[i,IDX[n]]=y[IDX[n]]/core.target[k]
    return Jlog[np.ix_(KEEP_ROWS,KEEP_COLS)]

def solve_conditional(y0,S,mu):
    u0=np.log(np.asarray(y0)[KEEP_COLS]);args=(S,mu)
    q=root(lambda u:residual_other(u,*args),u0,jac=lambda u:jac_other(u,*args),method='hybr',options={'xtol':2e-10,'maxfev':100})
    u=q.x;res=np.max(np.abs(residual_other(u,*args)))
    if res>2e-7:
        q=root(lambda u:residual_other(u,*args),u,jac=lambda u:jac_other(u,*args),method='lm',options={'ftol':1e-11,'xtol':1e-11,'gtol':1e-11,'maxiter':200});u=q.x;res=np.max(np.abs(residual_other(u,*args)))
    y=np.empty(N);y[IS]=S;y[KEEP_COLS]=np.exp(u)
    return y,float(res),bool(q.success)

def sam_balance(y,mu): return float(sic.rhs_input(y,0.0,mu,'methionine')[IS])
def flux_balance(y,mu):
    ov=sic.input_overrides(0.0,mu,'methionine');fl=core.M.flux(y,0,ov)
    ids=['V_MATI','V_MATIII','reaction_13','V_GNMT','V_DNMT','reaction_2'];return {k:float(fl[k]) for k in ids}

def branch_mu_at_S(S,mu0=None,y0=None):
    if y0 is None or mu0 is None:y0,mu0=archive_state_near_S(S)
    base=y0.copy()
    def eval_mu(mu):
        y,res,_=solve_conditional(base,S,float(mu));
        if res>2e-7: raise RuntimeError((S,mu,res))
        return sam_balance(y,float(mu)),y
    width=max(5e-5,0.002);lo=max(.01,mu0-width);hi=mu0+width;flo,_=eval_mu(lo);fhi,_=eval_mu(hi)
    for _ in range(8):
        if flo*fhi<=0:break
        width*=2;lo=max(.01,mu0-width);hi=mu0+width;flo,_=eval_mu(lo);fhi,_=eval_mu(hi)
    if flo*fhi>0: raise RuntimeError(('no bracket',S,mu0,flo,fhi))
    mu=brentq(lambda m:eval_mu(m)[0],lo,hi,xtol=2e-13,rtol=2e-13,maxiter=60)
    y,res,_=solve_conditional(base,S,mu);return mu,y,res,sam_balance(y,mu)

def run(outdir=DATA):
    rows=[]
    for S in [1e5,3e5,1e6,3e6,1e7,3e7,1e8,3e8,1e9,3e9,1e10]:
        y0,mu0=archive_state_near_S(S);mu,y,res,B=branch_mu_at_S(S,mu0,y0);fl=flux_balance(y,mu)
        r=dict(SAM=S,mu=mu,delta_mu_times_S=np.nan,residual_other=res,sam_balance=B,met=y[IDX['met']],sah=y[IDX['sah']],c_5mf=y[IDX['c_5mf']]);r.update(fl);rows.append(r);print(S,mu,'K?',res,flush=True)
    df=pd.DataFrame(rows);fits=[]
    for Smin in [1e5,1e6,1e7,1e8,1e9]:
        q=df[df.SAM>=Smin];x=1/q.SAM.values;A=np.c_[np.ones(len(q)),x];coef=np.linalg.lstsq(A,q.mu.values,rcond=None)[0];muinf=coef[0];K=-coef[1];pred=A@coef;fits.append(dict(Smin=Smin,n=len(q),mu_inf=muinf,K=K,max_abs_resid=float(np.max(np.abs(pred-q.mu.values)))))
    fits=pd.DataFrame(fits)
    qfit=df[(df.SAM>=1e7)&(df.SAM<=1e8)]
    def model(S,muinf,C,p): return muinf-C*S**(-p)
    popt,_=curve_fit(model,qfit.SAM.values,qfit.mu.values,p0=[2.15180254,132.0,1.0],maxfev=10000)
    muinf=float(popt[0]); Cfit=float(popt[1]); pfit=float(popt[2]);df['delta_mu_times_S']=(muinf-df.mu)*df.SAM
    Sref=1e7; yref,_=archive_state_near_S(Sref); yref,_,_=solve_conditional(yref,Sref,muinf); Bref=sam_balance(yref,muinf); B1eff=Bref*Sref
    dm=1e-5; vals=[]
    for dmu in [-2*dm,-dm,dm,2*dm]:
        yy,_,_=solve_conditional(yref,Sref,muinf+dmu); vals.append((dmu,sam_balance(yy,muinf+dmu)))
    xx=np.array([v[0] for v in vals]); yyv=np.array([v[1] for v in vals]); B0prime=float(np.dot(xx,yyv)/np.dot(xx,xx)); Kcoeff=float(B1eff/B0prime)
    S=float(df.iloc[-1].SAM);mu=float(df.iloc[-1].mu);y0,_=archive_state_near_S(S);y,_,_=solve_conditional(y0,S,mu)
    scale=np.maximum(y,1.0); Jraw=sic.physical_jacobian(y,0.0,mu,'methionine')/scale[:,None]
    for ss,k in core.REPL.items():
        i=IDX[ss]; Jraw[i,:]=0.0
        for n in core.INV[k]: Jraw[i,IDX[n]]=1.0/core.target[k]
    sv=np.linalg.svd(Jraw[np.ix_(KEEP_ROWS,KEEP_COLS)],compute_uv=False)
    last=flux_balance(y,mu);out=sum(last[k] for k in ['V_GNMT','V_DNMT','reaction_2']);inp=last['V_MATI']+last['V_MATIII']+last['reaction_13']
    summary={'mu_inf_fit':muinf,'C_fit':Cfit,'p_fit':pfit,'B0prime_mu':B0prime,'B1_effective':B1eff,'K_from_balance_coefficients':Kcoeff,'min_singular_value_reduced_raw_scaled_jacobian':float(sv[-1]),'condition_number_reduced_raw_scaled_jacobian':float(sv[0]/sv[-1]),'input_flux':inp,'output_flux':out,'DNMT_fraction_output':last['V_DNMT']/out,'GNMT_fraction_output':last['V_GNMT']/out,'SAMDC_fraction_output':last['reaction_2']/out,'MATIII_fraction_input':last['V_MATIII']/inp,'MATI_fraction_input':last['V_MATI']/inp}
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True);df.to_csv(outdir/'sam_asymptotic_conditional_branch.csv',index=False);fits.to_csv(outdir/'sam_asymptotic_simple_pole_fits.csv',index=False);pd.DataFrame([summary]).to_csv(outdir/'sam_asymptotic_summary.csv',index=False)
    print(fits.to_string(index=False));print(summary);return df,fits,summary
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=DATA);a=ap.parse_args();run(a.outdir)
