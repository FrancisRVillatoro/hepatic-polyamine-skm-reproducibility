#!/usr/bin/env python3
"""
Structural Kinetic Modeling (SKM) analysis of
BioModels BIOMD0000000190 (Rodriguez-Caso et al., JBC 2006).

The 4 SBML rate rules are split into positive synthesis and negative
degradation processes, giving an 11-state / 19-process generalized SKM.

Outputs:
  skm_exact_elasticities.csv
  skm_sensitivity.csv
  skm_robustness.csv
  skm_summary.txt
"""

from pathlib import Path
import importlib.util
import csv
import numpy as np
from numpy.linalg import eigvals, eig
from scipy.integrate import solve_ivp
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
XML = HERE / "BIOMD0000000190_url.xml"
AUDIT = HERE / "audit_BIOMD0000000190.py"

spec = importlib.util.spec_from_file_location("audit", AUDIT)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
m = audit.SBMLModel(XML)

sol = solve_ivp(m.rhs, (0, 2e5), m.initial_state(),
                method="BDF", rtol=1e-12, atol=1e-14)
if not sol.success:
    raise RuntimeError(sol.message)
xstar = sol.y[:, -1]

states = m.state_names
si = {s:i for i,s in enumerate(states)}

met = ["ODC","SAMdc","SSAT_for_S","SSAT_for_D",
       "PAO_for_aD","PAO_for_aS","SpdS","SpmS","MAT",
       "P_efflux","aD_efflux"]
reg = ["ODC_syn","ODC_deg","SSAT_syn","SSAT_deg",
       "SAMdc_syn","SAMdc_deg","Antz_syn","Antz_deg"]
procs = met + reg

N = np.zeros((len(states), len(procs)))
def col(name, changes):
    j = procs.index(name)
    for s, nu in changes.items():
        N[si[s], j] = nu

col("ODC", {"P":1})
col("SAMdc", {"SAM":-1,"A":1})
col("SSAT_for_S", {"S":-1,"aS":1})
col("SSAT_for_D", {"D":-1,"aD":1})
col("PAO_for_aD", {"aD":-1,"P":1})
col("PAO_for_aS", {"aS":-1,"D":1})
col("SpdS", {"A":-1,"P":-1,"D":1})
col("SpmS", {"A":-1,"D":-1,"S":1})
col("MAT", {"SAM":1})
col("P_efflux", {"P":-1})
col("aD_efflux", {"aD":-1})
col("ODC_syn", {"Vmaxodc":1})
col("ODC_deg", {"Vmaxodc":-1})
col("SSAT_syn", {"Vmaxssat":1})
col("SSAT_deg", {"Vmaxssat":-1})
col("SAMdc_syn", {"Vmaxsamdc":1})
col("SAMdc_deg", {"Vmaxsamdc":-1})
col("Antz_syn", {"Antz":1})
col("Antz_deg", {"Antz":-1})

def process_fluxes(y):
    env = m.environment(y)
    f = m.fluxes(y)
    D, S = env["D"], env["S"]
    out = [f[n] for n in met]
    out += [
        env["Ksodc"]/(1+env["Keq"]*(D+S)),
        env["Kdodc"]*env["Antz"]*env["Vmaxodc"],
        env["Ksssat"]*(1-1/(1+env["Keq"]*(D+S))),
        env["Kdssat"]/(1+env["Keq"]*(D+S))*env["Vmaxssat"],
        env["Kssamdc"]/(1+env["Keq"]*(D+S)),
        env["Kdsamdc"]*env["Vmaxsamdc"],
        env["Ksantz"]*(1-1/(1+env["Keq"]*0.01*(D+S))),
        env["Kdantz"]*env["Antz"],
    ]
    return np.asarray(out)

vstar = process_fluxes(xstar)

def elasticity_matrix(y, h=1e-30):
    v0 = np.asarray(process_fluxes(y), dtype=float)
    T = np.zeros((len(v0), len(y)))
    for k in range(len(y)):
        z = y.astype(complex)
        z[k] += 1j*h
        dv = np.imag(process_fluxes(z))/h
        T[:,k] = y[k]/v0 * dv
    return T

Theta = elasticity_matrix(xstar)
Lambda = np.diag(1/xstar) @ N @ np.diag(vstar)

Jraw = audit.jacobian_complex_step(m.rhs, xstar)
Jnormalized = np.diag(1/xstar) @ Jraw @ np.diag(xstar)
Jskm = Lambda @ Theta

errJ = np.max(np.abs(Jskm-Jnormalized))
e0 = eigvals(Jskm)
alpha0 = np.max(e0.real)
dom = e0[np.argmax(e0.real)]

# Fixed exact linear/mass-action elasticities.
fixed_pairs = [
    ("ODC","Vmaxodc"), ("SAMdc","Vmaxsamdc"),
    ("SSAT_for_S","Vmaxssat"), ("SSAT_for_D","Vmaxssat"),
    ("P_efflux","P"), ("aD_efflux","aD"),
    ("ODC_deg","Antz"), ("ODC_deg","Vmaxodc"),
    ("SSAT_deg","Vmaxssat"), ("SAMdc_deg","Vmaxsamdc"),
    ("Antz_deg","Antz"),
]
fixed = {(procs.index(p),states.index(s)) for p,s in fixed_pairs}
free = []
for j in range(len(procs)):
    for k in range(len(states)):
        if abs(Theta[j,k]) > 1e-12 and (j,k) not in fixed:
            free.append((j,k,1 if Theta[j,k] > 0 else -1))

# Exact elasticities.
with open(DATA/"skm_exact_elasticities.csv","w",newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["process","state","theta_exact","fixed_linear"])
    for j,p in enumerate(procs):
        for k,s in enumerate(states):
            if abs(Theta[j,k]) > 1e-12:
                w.writerow([p,s,f"{Theta[j,k]:.16g}",int((j,k) in fixed)])

# Dominant-mode participation.
vals, vr = eig(Jskm)
idom = np.argmax(vals.real)
lam = vals[idom]
v = vr[:,idom]
valsL, vl = eig(Jskm.T)
ileft = np.argmin(np.abs(valsL-lam))
wleft = vl[:,ileft]
wleft /= np.dot(wleft, v)
participation = np.abs(wleft*v)
participation /= participation.sum()

# Local sensitivity of the dominant branch to each free elasticity.
def branch(T):
    ee = eigvals(Lambda@T)
    return ee[np.argmin(np.abs(ee-lam))]

sens = []
h = 1e-6
for j,k,sgn in free:
    Tp, Tm = Theta.copy(), Theta.copy()
    Tp[j,k] += h
    Tm[j,k] -= h
    dlam = (branch(Tp)-branch(Tm))/(2*h)
    sens.append((procs[j],states[k],Theta[j,k],dlam.real,dlam.imag))

with open(DATA/"skm_sensitivity.csv","w",newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["process","state","theta_exact",
                "d_Re_lambda_d_theta","d_Im_lambda_d_theta"])
    for row in sorted(sens,key=lambda z:abs(z[3]),reverse=True):
        w.writerow(row)

# Monte Carlo ensembles.
def classify(T):
    ee = eigvals(Lambda@T)
    q = np.argmax(ee.real)
    a = ee[q].real
    im = abs(ee[q].imag)
    return a, im

robust = []

# Broad sign-only ensemble.
rng = np.random.default_rng(20260801)
n_broad = 200_000
alphas = np.empty(n_broad)
ims = np.empty(n_broad)
X = rng.random((n_broad,len(free)))
for i in range(n_broad):
    T = np.zeros_like(Theta)
    for j,k in fixed:
        T[j,k] = Theta[j,k]
    for q,(j,k,sgn) in enumerate(free):
        T[j,k] = sgn*X[i,q]
    alphas[i],ims[i] = classify(T)

stable = alphas < 0
robust.append([
    "sign_only", "", n_broad, stable.mean(),
    np.mean((alphas>=0)&(ims<1e-8)),
    np.mean((alphas>=0)&(ims>=1e-8)),
    *np.quantile(alphas,[0.01,0.50,0.99])
])

# Spearman ranking in broad ensemble.
rank = []
for q,(j,k,sgn) in enumerate(free):
    actual = sgn*X[:,q]
    rho,_ = spearmanr(actual,alphas)
    rank.append((abs(rho),rho,procs[j],states[k],Theta[j,k]))
rank.sort(reverse=True)

# Kinetically anchored multiplicative envelopes.
for z,factor in enumerate([1.5,2,3,5,10]):
    rng = np.random.default_rng(20260811+z)
    n = 100_000
    aa = np.empty(n)
    ii = np.empty(n)
    L = np.log(factor)
    for i in range(n):
        T = np.zeros_like(Theta)
        for j,k in fixed:
            T[j,k] = Theta[j,k]
        mult = np.exp(rng.uniform(-L,L,len(free)))
        for q,(j,k,sgn) in enumerate(free):
            mag = min(abs(Theta[j,k])*mult[q],1.0)
            T[j,k] = sgn*mag
        aa[i],ii[i] = classify(T)
    robust.append([
        "anchored", factor, n, np.mean(aa<0),
        np.mean((aa>=0)&(ii<1e-8)),
        np.mean((aa>=0)&(ii>=1e-8)),
        *np.quantile(aa,[0.01,0.50,0.99])
    ])

with open(DATA/"skm_robustness.csv","w",newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["ensemble","factor","n","p_stable",
                "p_real_instability","p_complex_instability",
                "alpha_q01","alpha_median","alpha_q99"])
    w.writerows(robust)

turnover = (np.maximum(N,0)@vstar)/xstar

with open(DATA/"skm_summary.txt","w") as fh:
    print("BIOMD0000000190 generalized SKM",file=fh)
    print(f"states={len(states)} processes={len(procs)} free_elasticities={len(free)}",file=fh)
    print(f"max |J_SKM-J_normalized| = {errJ:.6e}",file=fh)
    print(f"dominant eigenvalue = {dom.real:+.12e} {dom.imag:+.12e} i 1/min",file=fh)
    print("\nDominant-mode participation:",file=fh)
    for q in np.argsort(participation)[::-1]:
        print(f"{states[q]:12s} {participation[q]:.8f}",file=fh)
    print("\nTurnover rates (1/min):",file=fh)
    for s,t in zip(states,turnover):
        print(f"{s:12s} {t:.12g}",file=fh)
    print("\nBroad sign-only stability:",file=fh)
    print(f"P(stable)={stable.mean():.8f}",file=fh)
    print("\nTop broad-ensemble Spearman correlations with spectral abscissa:",file=fh)
    for _,rho,p,s,th in rank[:12]:
        print(f"{p:16s} {s:12s} rho={rho:+.6f} theta*={th:+.6f}",file=fh)

print(f"states={len(states)}, processes={len(procs)}, free elasticities={len(free)}")
print(f"max SKM Jacobian error: {errJ:.3e}")
print(f"dominant eigenvalue: {dom.real:+.9e} {dom.imag:+.9e} i 1/min")
print("robustness:")
for r in robust:
    print(r)
