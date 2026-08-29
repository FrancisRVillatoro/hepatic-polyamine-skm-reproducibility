from pathlib import Path
import importlib.util, numpy as np, pandas as pd, time
from scipy.optimize import root
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; DATA=ROOT/'data'
spec=importlib.util.spec_from_file_location('s3',HERE/'step3_core.py'); s3=importlib.util.module_from_spec(spec); spec.loader.exec_module(s3)
M=s3.M; idx=s3.idx; N=len(M.state); row_scale=np.maximum(np.asarray(s3.yref,float),1.0)

def controls(rM,rH,const_total=False):
 if const_total:
  vm1=260*(1-rM); vm3=220*(1-rM); vm2=480*rM
 else:
  vm1=260*(1-rM); vm3=220*(1-rM); vm2=220*rM
 return {'Vm_MAT1':vm1,'Vm_MAT3':vm3,'parameter_19':vm2,'H2O2':0.01*(1+0.5*rH),'fasting':1.,'breakfast':1.,'lunch':1.,'dinner':1.}

def rhs(y,rM,rF,rH,const_total=False,fbconst=65.06):
 ov=controls(rM,rH,const_total); e=M.env(y,0,ov); fl=M.flux(y,0,ov); d={n:0.0 for n in M.state}
 for rx in M.rx:
  for sp,nu in rx['st'].items():
   if sp in d: d[sp]+=nu*fl[rx['id']]/M.comp[M.spec[sp]['comp']]
 D,S=e['species_4'],e['species_3']; free=1/(1+e['parameter_5']*(D+S)); g=(1-rF)+rF*fbconst/e['sam']
 d['parameter_1']=60*g*e['parameter_7']*free-e['parameter_6']*e['parameter_4']*e['parameter_1']
 d['parameter_2']=60*e['parameter_9']*(1-free)-e['parameter_8']*free*e['parameter_2']
 d['parameter_3']=60*g*e['parameter_11']*free-e['parameter_10']*e['parameter_3']
 d['parameter_4']=e['parameter_13']*(1-1/(1+e['parameter_5']*.01*(D+S)))-e['parameter_12']*e['parameter_4']
 return np.array([d[n] for n in M.state])

def jac(y,c,const_total=False,fbconst=65.06,h=1e-30):
 J=np.empty((N,N))
 for j in range(N):
  z=y.astype(complex); z[j]+=1j*h; J[:,j]=np.imag(rhs(z,*c,const_total=const_total,fbconst=fbconst))/h
 return J

def fun(z,c,const_total=False,fbconst=65.06):
 y=np.exp(z); rr=rhs(y,*c,const_total=const_total,fbconst=fbconst)/row_scale
 for ss,k in s3.REPL.items(): rr[idx[ss]]=(sum(y[idx[n]] for n in s3.INV[k])-s3.target[k])/s3.target[k]
 return rr

def jaclog(z,c,const_total=False,fbconst=65.06):
 y=np.exp(z); J=jac(y,c,const_total,fbconst)*y[np.newaxis,:]/row_scale[:,None]
 for ss,k in s3.REPL.items():
  i=idx[ss];J[i,:]=0
  for n in s3.INV[k]:J[i,idx[n]]=y[idx[n]]/s3.target[k]
 return J

def solve_path(name,const_total=False,full=False,fbconst=65.06,rH_fixed=None):
 z=np.log(s3.polish(s3.yref,(0,0,0))); rows=[]
 for r in np.linspace(0,1,21):
  c=(float(r),float(r if full else 0),float(r if rH_fixed is None and full else (rH_fixed or 0)))
  if r>0:
   sol=root(lambda zz:fun(zz,c,const_total,fbconst),z,jac=lambda zz:jaclog(zz,c,const_total,fbconst),method='hybr',options={'xtol':5e-10,'maxfev':80});z=sol.x
   if np.max(abs(fun(z,c,const_total,fbconst)))>2e-6:
    sol=root(lambda zz:fun(zz,c,const_total,fbconst),z,jac=lambda zz:jaclog(zz,c,const_total,fbconst),method='lm',options={'ftol':1e-10,'xtol':1e-10,'gtol':1e-10,'maxiter':120});z=sol.x
  y=np.exp(z); res=np.max(abs(fun(z,c,const_total,fbconst))); ev=np.linalg.eigvals(jac(y,c,const_total,fbconst)); keep=np.ones(len(ev),bool);keep[np.argsort(abs(ev))[:3]]=False;ee=ev[keep];dom=ee[np.argmax(ee.real)]
  ov=controls(c[0],c[2],const_total)
  rows.append(dict(name=name,r=r,rM=c[0],rF=c[1],rH=c[2],MAT_total=ov['Vm_MAT1']+ov['Vm_MAT3']+ov['parameter_19'],SAM=y[idx['sam']],Put=y[idx['species_2']],Spd=y[idx['species_4']],Spm=y[idx['species_3']],alpha=dom.real,residual=res))
 return rows

allrows=[];t=time.time()
allrows += solve_path('original_MAT_only',False,False)
allrows += solve_path('constant_total_MAT_only',True,False)
allrows += solve_path('original_full_65p06',False,True,65.06)
allrows += solve_path('original_full_66p5',False,True,66.5)
allrows += solve_path('MAT_feedback_Hfixed0',False,True,65.06,rH_fixed=0.0)
df=pd.DataFrame(allrows);df.to_csv(DATA/'audit_control_sweeps.csv',index=False)
print('elapsed',time.time()-t)
for name,g in df.groupby('name'):
 print('\n',name);print(g.iloc[[0,-1]][['r','MAT_total','SAM','Put','Spd','Spm','alpha','residual']].to_string(index=False))
