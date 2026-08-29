#!/usr/bin/env python3
"""Step 3 homotopy for Reyes-Palomares et al. 2012 BioModels.

Requires in the same directory:
  BIOMD0000000674_url.xml  (physiological)
  BIOMD0000000450_url.xml  (proliferative)

The stability analysis is autonomous: the explicit meal cycle is normalized to
 aa_input=1 in both models. The homotopy uses the SBML-encoded changes:
 MATI/III -> MATII, inverse-SAM feedback in ODC/SAMDC synthesis, and +50% H2O2.
"""
from pathlib import Path
import xml.etree.ElementTree as ET
import collections, csv
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
P1=HERE/'BIOMD0000000674_url.xml'
P2=HERE/'BIOMD0000000450_url.xml'

class SBML:
    def __init__(self,path):
        root=ET.parse(path).getroot(); uri=root.tag.split('}')[0].strip('{')
        self.ns={'s':uri,'m':'http://www.w3.org/1998/Math/MathML'}
        self.model=root.find('s:model',self.ns)
        self.comp={e.attrib['id']:float(e.attrib.get('size','1')) for e in self.model.findall('./s:listOfCompartments/s:compartment',self.ns)}
        self.spec={}
        for e in self.model.findall('./s:listOfSpecies/s:species',self.ns):
            self.spec[e.attrib['id']]={'name':e.attrib.get('name',''),'init':float(e.attrib.get('initialConcentration',e.attrib.get('initialAmount','0'))),'boundary':e.attrib.get('boundaryCondition','false')=='true','constant':e.attrib.get('constant','false')=='true','comp':e.attrib['compartment']}
        self.par={e.attrib['id']:float(e.attrib.get('value','0')) for e in self.model.findall('./s:listOfParameters/s:parameter',self.ns)}
        self.func={}
        for fd in self.model.findall('./s:listOfFunctionDefinitions/s:functionDefinition',self.ns):
            lam=list(fd.find('m:math',self.ns))[0]; args=[]; body=None
            for ch in list(lam):
                if ch.tag.split('}')[-1]=='bvar': args.append((list(ch)[0].text or '').strip())
                else: body=ch
            self.func[fd.attrib['id']]=(args,body)
        self.assign={}; self.rate={}
        for r in self.model.findall('./s:listOfRules/*',self.ns):
            body=list(r.find('m:math',self.ns))[0]; typ=r.tag.split('}')[-1]
            (self.assign if typ=='assignmentRule' else self.rate)[r.attrib['variable']]=body
        self.rx=[]
        for r in self.model.findall('./s:listOfReactions/s:reaction',self.ns):
            local={p.attrib['id']:float(p.attrib.get('value','0')) for p in r.findall('./s:kineticLaw/s:listOfParameters/s:parameter',self.ns)}
            st=collections.defaultdict(float)
            for sr in r.findall('./s:listOfReactants/s:speciesReference',self.ns): st[sr.attrib['species']]-=float(sr.attrib.get('stoichiometry','1'))
            for sr in r.findall('./s:listOfProducts/s:speciesReference',self.ns): st[sr.attrib['species']]+=float(sr.attrib.get('stoichiometry','1'))
            self.rx.append({'id':r.attrib['id'],'body':list(r.find('./s:kineticLaw/m:math',self.ns))[0],'local':local,'st':dict(st)})
        self.dynsp=[s for s,d in self.spec.items() if not d['boundary'] and not d['constant'] and s not in self.assign]
        self.state=self.dynsp+list(self.rate)
    def ev(self,n,e,t=0.0):
        k=n.tag.split('}')[-1]
        if k=='ci': return e[(n.text or '').strip()]
        if k=='cn': return float((n.text or '0').strip())
        if k=='csymbol': return t
        if k=='piecewise':
            other=None
            for c in list(n):
                q=c.tag.split('}')[-1]
                if q=='piece':
                    val,cond=list(c)
                    if self.ev(cond,e,t): return self.ev(val,e,t)
                elif q=='otherwise': other=list(c)[0]
            return self.ev(other,e,t)
        if k=='apply':
            c=list(n); op=c[0].tag.split('}')[-1]; a=c[1:]
            if op=='ci':
                f=(c[0].text or '').strip(); vals=[self.ev(x,e,t) for x in a]; names,body=self.func[f]; ee=dict(e); ee.update(zip(names,vals)); return self.ev(body,ee,t)
            v=[self.ev(x,e,t) for x in a]
            if op=='plus': return sum(v)
            if op=='times':
                z=1
                for x in v: z*=x
                return z
            if op=='minus': return -v[0] if len(v)==1 else v[0]-sum(v[1:])
            if op=='divide': return v[0]/v[1]
            if op=='power': return v[0]**v[1]
            if op=='exp': return np.exp(v[0])
            if op=='floor': return np.floor(v[0])
            if op=='leq': return v[0]<=v[1]
            if op=='and': return all(v)
            raise NotImplementedError(op)
        raise NotImplementedError(k)
    def env(self,y,t=0.0,ov=None):
        e={}; e.update(self.comp); e.update({s:d['init'] for s,d in self.spec.items()}); e.update(self.par)
        e.update(zip(self.state,y))
        if ov: e.update(ov)
        for _ in range(3):
            for q,n in self.assign.items(): e[q]=self.ev(n,e,t)
        return e
    def flux(self,y,t=0.0,ov=None):
        e=self.env(y,t,ov); z={}
        for r in self.rx:
            ee=dict(e); ee.update(r['local']); z[r['id']]=self.ev(r['body'],ee,t)
        return z
    def y0(self): return np.array([self.spec[n]['init'] if n in self.spec else self.par[n] for n in self.state],float)

M=SBML(P2); M1=SBML(P1); idx={n:i for i,n in enumerate(M.state)}
assert M.state==M1.state

def rhs(y,rM,rF,rH):
    ov={'Vm_MAT1':260*(1-rM),'Vm_MAT3':220*(1-rM),'parameter_19':220*rM,'H2O2':0.01*(1+0.5*rH),'fasting':1.,'breakfast':1.,'lunch':1.,'dinner':1.}
    e=M.env(y,0,ov); fl=M.flux(y,0,ov); d={n:0.0 for n in M.state}
    for r in M.rx:
        for s,nu in r['st'].items():
            if s in d: d[s]+=nu*fl[r['id']]/M.comp[M.spec[s]['comp']]
    D,S=e['species_4'],e['species_3']; free=1/(1+e['parameter_5']*(D+S)); g=(1-rF)+rF*65.06/e['sam']
    d['parameter_1']=60*g*e['parameter_7']*free-e['parameter_6']*e['parameter_4']*e['parameter_1']
    d['parameter_2']=60*e['parameter_9']*(1-free)-e['parameter_8']*free*e['parameter_2']
    d['parameter_3']=60*g*e['parameter_11']*free-e['parameter_10']*e['parameter_3']
    d['parameter_4']=e['parameter_13']*(1-1/(1+e['parameter_5']*.01*(D+S)))-e['parameter_12']*e['parameter_4']
    return np.array([d[n] for n in M.state])

# Three exact conserved pools.
INV={'cfol':['c_thf','c_5mf','c_2cf','c_1cf','c_10f','c_dhf'],'mfol':['m_thf','m_2cf','m_1cf','m_10f'],'coa':['species_8','species_9']}
REPL={'c_thf':'cfol','m_thf':'mfol','species_9':'coa'}
# Physiological autonomous state: version-1 initial state is stationary at aa_input=1.
yref=M1.y0(); target={k:sum(yref[idx[n]] for n in ns) for k,ns in INV.items()}

def residual(logy,rM,rF,rH):
    y=np.exp(logy); r=rhs(y,rM,rF,rH)/np.maximum(np.abs(y),1.)
    for s,k in REPL.items(): r[idx[s]]=(sum(y[idx[n]] for n in INV[k])-target[k])/target[k]
    return r

def polish(y,c):
    q=least_squares(lambda z:residual(z,*c),np.log(y),xtol=3e-11,ftol=3e-11,gtol=3e-11,max_nfev=80,x_scale='jac')
    return np.exp(q.x)

def jac(y,c,h=1e-30):
    J=np.empty((len(y),len(y)))
    for j in range(len(y)):
        z=y.astype(complex); z[j]+=1j*h; J[:,j]=np.imag(rhs(z,*c))/h
    return J

def dom(y,c):
    ev=np.linalg.eigvals(jac(y,c)); ev=ev[np.abs(ev)>1e-8]; return ev[np.argmax(ev.real)]

