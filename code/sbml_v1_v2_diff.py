#!/usr/bin/env python3
"""Automated structural diff of the two Reyes-Palomares 2012 SBML records."""
from pathlib import Path
import xml.etree.ElementTree as ET
import pandas as pd
HERE=Path(__file__).resolve().parent;DATA=HERE.parent/'data'
A=HERE/'BIOMD0000000674_url.xml';B=HERE/'BIOMD0000000450_url.xml'

def root_model(path):
    r=ET.parse(path).getroot();uri=r.tag.split('}')[0].strip('{');ns={'s':uri,'m':'http://www.w3.org/1998/Math/MathML'};return r.find('s:model',ns),ns

def canon(el):
    if el is None:return ''
    def rec(x):
        tag=x.tag.split('}')[-1];attrs=tuple(sorted((k.split('}')[-1],v) for k,v in x.attrib.items() if k.split('}')[-1] not in {'metaid'}));txt=' '.join((x.text or '').split())
        return (tag,attrs,txt,tuple(rec(c) for c in list(x)))
    return repr(rec(el))

def maps(path):
    m,ns=root_model(path);out={}
    out['species']={e.attrib['id']:(e.attrib.get('initialConcentration',e.attrib.get('initialAmount','')),e.attrib.get('boundaryCondition','false'),e.attrib.get('constant','false')) for e in m.findall('./s:listOfSpecies/s:species',ns)}
    out['parameters']={e.attrib['id']:e.attrib.get('value','') for e in m.findall('./s:listOfParameters/s:parameter',ns)}
    out['functions']={e.attrib['id']:canon(e.find('m:math',ns)) for e in m.findall('./s:listOfFunctionDefinitions/s:functionDefinition',ns)}
    out['rules']={e.attrib.get('variable',''):(e.tag.split('}')[-1],canon(e.find('m:math',ns))) for e in m.findall('./s:listOfRules/*',ns)}
    rx={}
    for e in m.findall('./s:listOfReactions/s:reaction',ns):
        local={p.attrib['id']:p.attrib.get('value','') for p in e.findall('./s:kineticLaw/s:listOfParameters/s:parameter',ns)}
        rx[e.attrib['id']]=(e.attrib.get('name',''),canon(e.find('./s:kineticLaw/m:math',ns)),local)
    out['reactions']=rx
    return out

def classify(kind,key,a,b):
    if kind=='species' and key=='c_gly': return 'rounding-only initial value'
    if key=='H2O2': return 'oxidative-stress control'
    if key in {'Vm_MAT1','Vm_MAT3','parameter_19'} or key=='function_4_V_MATII' or key=='reaction_13': return 'MAT reparameterization'
    if kind=='rules' and key in {'parameter_1','parameter_3'}: return 'inverse-SAM ODC/SAMDC feedback'
    if key in {'b_met','b_ser','V_oGly_b','V_oCys_b','V_oGlu_b','aa_input','fasting','breakfast','lunch','dinner'}: return 'meal-forcing/input bookkeeping'
    if key in {'Constant_flux__reversible','function_1'} or (kind=='reactions' and key in {'b_cys_import','b_glu_import','b_gly_import'}): return 'syntactic identity-function rename'
    return 'other; inspect'

def autonomous_effect(classification):
    if classification in {'meal-forcing/input bookkeeping','syntactic identity-function rename','rounding-only initial value'}: return 'removed/neutralized by autonomous basal-input normalization or numerically null'
    return 'retained as an explicit control in the common autonomous model'

def main(outdir=DATA):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    a=maps(A);b=maps(B);rows=[]
    for kind in a:
        keys=sorted(set(a[kind])|set(b[kind]))
        for k in keys:
            av=a[kind].get(k,'<missing>');bv=b[kind].get(k,'<missing>')
            if av!=bv:
                cl=classify(kind,k,av,bv); rows.append(dict(kind=kind,key=k,value_v1=str(av),value_v2=str(bv),classification=cl,autonomous_effect=autonomous_effect(cl)))
    df=pd.DataFrame(rows);df.to_csv(outdir/'sbml_v1_v2_structural_diff.csv',index=False)
    print(df.groupby(['kind','classification']).size().to_string());print('\nDETAIL');print(df[['kind','key','classification']].to_string(index=False))
    from step3_core import M1,M,INV
    rr=[]
    for pool,names in INV.items():
        v1=sum(M1.spec[n]['init'] for n in names);v2=sum(M.spec[n]['init'] for n in names);rr.append(dict(pool=pool,total_v1=v1,total_v2=v2,difference=v2-v1,components=';'.join(names)))
    pd.DataFrame(rr).to_csv(outdir/'sbml_conserved_pool_totals_v1_v2.csv',index=False);print('\nPOOLS');print(pd.DataFrame(rr).to_string(index=False))
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=DATA);a=ap.parse_args();main(a.outdir)
