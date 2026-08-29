from pathlib import Path
import importlib.util, itertools, numpy as np, pandas as pd, sys, time
from scipy.optimize import root
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; DATA=ROOT/'data'; spec=importlib.util.spec_from_file_location('s3',HERE/'step3_core.py'); s3=importlib.util.module_from_spec(spec); spec.loader.exec_module(s3)
row_scale=np.maximum(np.asarray(s3.yref,float),1.0); ref=pd.read_csv(DATA/'rho_branch.csv').iloc[0]; lam_pa0=complex(ref.polyamine_pair_real_per_h,ref.polyamine_pair_imag_per_h)
def fun(z,c):
 y=np.exp(z); rr=s3.rhs(y,*c)/row_scale
 for ss,k in s3.REPL.items(): rr[s3.idx[ss]]=(sum(y[s3.idx[n]] for n in s3.INV[k])-s3.target[k])/s3.target[k]
 return rr
def jaclog(z,c):
 y=np.exp(z); J=s3.jac(y,c)*y[np.newaxis,:]/row_scale[:,None]
 for ss,k in s3.REPL.items():
  i=s3.idx[ss]; J[i,:]=0
  for n in s3.INV[k]: J[i,s3.idx[n]]=y[s3.idx[n]]/s3.target[k]
 return J
paths=[('diagonal',1.,1.,1.)]+[(f'async_{p[0]:g}_{p[1]:g}_{p[2]:g}',*p) for p in itertools.product([0.5,2.],repeat=3)]
lo=int(sys.argv[1]) if len(sys.argv)>1 else 0; hi=int(sys.argv[2]) if len(sys.argv)>2 else len(paths); rows=[]; t0=time.time()
for name,pM,pF,pH in paths[lo:hi]:
 y=s3.polish(s3.yref,(0,0,0)); z=np.log(y); prev_pa=lam_pa0
 for ss in np.linspace(0,1,31):
  c=(0 if ss==0 else ss**pM,0 if ss==0 else ss**pF,0 if ss==0 else ss**pH)
  if ss>0:
   sol=root(lambda zz:fun(zz,c),z,jac=lambda zz:jaclog(zz,c),method='hybr',options={'xtol':5e-10,'maxfev':80}); z=sol.x
   if np.max(np.abs(fun(z,c)))>2e-6:
    sol=root(lambda zz:fun(zz,c),z,jac=lambda zz:jaclog(zz,c),method='lm',options={'ftol':1e-10,'xtol':1e-10,'gtol':1e-10,'maxiter':100}); z=sol.x
  y=np.exp(z); res=float(np.max(np.abs(fun(z,c))))
  ev=np.linalg.eigvals(s3.jac(y,c)); order=np.argsort(np.abs(ev)); keep=np.ones(len(ev),bool); keep[order[:3]]=False; e=ev[keep]; dom=e[np.argmax(e.real)]
  pos=e[e.imag>1e-7]; pa=pos[np.argmin(np.abs(pos-prev_pa))] if len(pos) else np.nan+1j*np.nan
  if len(pos): prev_pa=pa
  rows.append((name,pM,pF,pH,ss,*c,res,dom.real,dom.imag,pa.real,pa.imag,y[s3.idx['sam']],y[s3.idx['species_2']],y[s3.idx['species_4']],y[s3.idx['species_3']]))
 print(name,'done',time.time()-t0,flush=True)
cols=['path','pM','pF','pH','s','rM','rF','rH','residual','alpha','alpha_imag','pa_real','pa_imag','SAM','Put','Spd','Spm']; df=pd.DataFrame(rows,columns=cols); out=DATA/('homotopy_9paths.csv' if lo==0 and hi==len(paths) else f'homotopy_subset_{lo}_{hi}.csv'); df.to_csv(out,index=False); summary=df.loc[df.groupby('path').alpha.idxmax(),['path','pM','pF','pH','s','rM','rF','rH','alpha','pa_real']].rename(columns={'alpha':'max_alpha','pa_real':'max_pa_real'}); summary.to_csv(DATA/'homotopy_9paths_summary.csv',index=False) if lo==0 and hi==len(paths) else None
