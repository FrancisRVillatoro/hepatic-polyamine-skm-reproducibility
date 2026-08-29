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
levels=np.linspace(0,1,5); points=[]; seg=0
for i,rM in enumerate(levels):
 fs=levels if i%2==0 else levels[::-1]
 for rF in fs:
  hs=levels if seg%2==0 else levels[::-1]; seg+=1
  for rH in hs: points.append((float(rM),float(rF),float(rH)))
# grid-neighbor path
assert max(sum(abs(a-b) for a,b in zip(points[t],points[t+1])) for t in range(len(points)-1))<0.251
y=s3.polish(s3.yref,(0,0,0)); z=np.log(y); rows=[]; t0=time.time()
for n,c in enumerate(points):
 if n>0:
  sol=root(lambda zz:fun(zz,c),z,jac=lambda zz:jaclog(zz,c),method='hybr',options={'xtol':5e-10,'maxfev':80}); z=sol.x
  if np.max(abs(fun(z,c)))>2e-6:
   sol=root(lambda zz:fun(zz,c),z,jac=lambda zz:jaclog(zz,c),method='lm',options={'ftol':1e-10,'xtol':1e-10,'gtol':1e-10,'maxiter':120}); z=sol.x
 y=np.exp(z); res=float(np.max(abs(fun(z,c))))
 ev=np.linalg.eigvals(s3.jac(y,c)); order=np.argsort(np.abs(ev)); keep=np.ones(len(ev),bool); keep[order[:3]]=False; e=ev[keep]; dom=e[np.argmax(e.real)]
 rows.append((*c,res,dom.real,dom.imag,y[s3.idx['sam']],y[s3.idx['species_2']],y[s3.idx['species_4']],y[s3.idx['species_3']]))
df=pd.DataFrame(rows,columns=['rM','rF','rH','residual','alpha','alpha_imag','SAM','Put','Spd','Spm']); df.to_csv(DATA/'control_cube_5x5x5.csv',index=False)
print('elapsed',time.time()-t0); print('max alpha'); print(df.loc[df.alpha.idxmax()].to_string()); print('min alpha'); print(df.loc[df.alpha.idxmin()].to_string()); print('positive',int((df.alpha>=0).sum()),'maxres',df.residual.max())
