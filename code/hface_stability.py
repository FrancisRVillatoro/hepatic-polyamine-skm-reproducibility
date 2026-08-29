from pathlib import Path
import importlib.util, numpy as np, pandas as pd, time
from scipy.optimize import root
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; DATA=ROOT/'data'; spec=importlib.util.spec_from_file_location('s3',HERE/'step3_core.py'); s3=importlib.util.module_from_spec(spec); spec.loader.exec_module(s3)
row_scale=np.maximum(np.asarray(s3.yref,float),1.0)
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
# First reach H-only corner smoothly.
y=s3.polish(s3.yref,(0,0,0)); z=np.log(y)
for h in np.linspace(.1,1,10):
 c=(0.,0.,h); sol=root(lambda zz:fun(zz,c),z,jac=lambda zz:jaclog(zz,c),method='hybr',options={'xtol':5e-10,'maxfev':80}); z=sol.x
# Traverse H=1 face from (M,F)=(0,0) in snake order.
levels=np.linspace(0,1,9); rows=[]; t0=time.time(); seg=0
for i,rM in enumerate(levels):
 fs=levels if i%2==0 else levels[::-1]
 for rF in fs:
  c=(rM,float(rF),1.)
  sol=root(lambda zz:fun(zz,c),z,jac=lambda zz:jaclog(zz,c),method='hybr',options={'xtol':5e-10,'maxfev':80}); z=sol.x
  if np.max(abs(fun(z,c)))>2e-6:
   sol=root(lambda zz:fun(zz,c),z,jac=lambda zz:jaclog(zz,c),method='lm',options={'ftol':1e-10,'xtol':1e-10,'gtol':1e-10,'maxiter':100}); z=sol.x
  y=np.exp(z); res=float(np.max(abs(fun(z,c))))
  ev=np.linalg.eigvals(s3.jac(y,c)); order=np.argsort(np.abs(ev)); keep=np.ones(len(ev),bool); keep[order[:3]]=False; e=ev[keep]; dom=e[np.argmax(e.real)]
  rows.append((rM,rF,res,dom.real,dom.imag,y[s3.idx['sam']],y[s3.idx['species_2']]))
df=pd.DataFrame(rows,columns=['rM','rF','residual','alpha','alpha_imag','SAM','Put']); df.to_csv(DATA/'hface_stability_9x9.csv',index=False)
print('elapsed',time.time()-t0); print('max alpha'); print(df.loc[df.alpha.idxmax()].to_string()); print('positive',int((df.alpha>=0).sum()),'maxres',df.residual.max())
