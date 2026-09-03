#!/usr/bin/env python3
"""Nondegeneracy audit of the reduced q=n=3 regulatory cusp."""
from pathlib import Path
import json, math
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
M = 0.680208; Q = N = 3.0


def rho_y(y): return y**Q/(1.0+y**Q)
def hn(y): return y**N/(1.0+y**N)


def poly_derivative(fun,x0,order,h,radius=4):
    dx=np.arange(-radius,radius+1,dtype=float)*h
    vals=np.array([fun(x0+d) for d in dx]); c=np.polynomial.polynomial.polyfit(dx,vals,2*radius)
    return float(math.factorial(order)*c[order])


def run():
    br=pd.read_csv(DATA/"rho_branch.csv"); cp=pd.read_csv(DATA/"rho0_cusp_q3n3.csv").iloc[0]
    p=PchipInterpolator(br.rho.to_numpy(float),np.log(br.SAM_uM.to_numpy(float)))
    samh=float(np.exp(p(0.0)))
    def sam(r): return float(np.exp(p(r)))
    def amp(r0):
        y0=(r0/(1.0-r0))**(1.0/Q); return y0*(sam(r0)/samh)**M
    def F(y,a,r0): return amp(r0)*(samh/sam(rho_y(y)))**M*(1.0+a*hn(y))-y
    yc=float(cp.Y_cusp); ac=float(cp.a_cusp); rc=float(cp.rho0_cusp)
    hs=[4e-4,6e-4,8e-4,1e-3,1.5e-3]
    fy=float(np.median([poly_derivative(lambda y:F(y,ac,rc),yc,1,h) for h in hs]))
    fyy=float(np.median([poly_derivative(lambda y:F(y,ac,rc),yc,2,h) for h in hs]))
    fyyy=float(np.median([poly_derivative(lambda y:F(y,ac,rc),yc,3,h) for h in hs]))
    ha,hr,hy=1e-5,1e-6,5e-4
    Fa=(F(yc,ac+ha,rc)-F(yc,ac-ha,rc))/(2*ha)
    Fr=(F(yc,ac,rc+hr)-F(yc,ac,rc-hr))/(2*hr)
    def Fy(a,r): return poly_derivative(lambda y:F(y,a,r),yc,1,hy)
    Fya=(Fy(ac+ha,rc)-Fy(ac-ha,rc))/(2*ha); Fyr=(Fy(ac,rc+hr)-Fy(ac,rc-hr))/(2*hr)
    out={"rho0_cusp":rc,"Y_cusp":yc,"a_cusp":ac,"F":float(F(yc,ac,rc)),
         "F_Y":fy,"F_YY":fyy,"F_YYY":fyyy,
         "unfolding_determinant_a_rho0":float(Fa*Fyr-Fr*Fya),
         "reference_F_YYY":-4.804,"reference_unfolding_determinant":-0.812}
    (DATA/"cusp_nondegeneracy_check.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2)); return out


if __name__ == "__main__": run()
