#!/usr/bin/env python3
"""Native meal-cycle control for BIOMD0000000674."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import step3_core as core

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
M = core.M1
IS = M.state.index("sam")
SWITCHES = [(0,7),(7,10),(10,12),(12,15),(15,18),(18,21),(21,24)]


def rhs_native(t, y):
    e = M.env(y, t); fl = M.flux(y, t); d = {n: 0.0 for n in M.state}
    for r in M.rx:
        for s, nu in r["st"].items():
            if s in d: d[s] += nu * fl[r["id"]] / M.comp[M.spec[s]["comp"]]
    for name, body in M.rate.items(): d[name] = M.ev(body, e, t)
    return np.array([d[n] for n in M.state], float)


def one_day(y0, day, rtol, atol, dense=True):
    y = np.asarray(y0, float).copy(); traces = []
    for h0, h1 in SWITCHES:
        a, b = 24.0*day+h0, 24.0*day+h1
        sol = solve_ivp(rhs_native, (a,b), y, method="BDF", rtol=rtol, atol=atol,
                        max_step=0.5, dense_output=dense)
        if not sol.success: raise RuntimeError(sol.message)
        if dense:
            grid = np.linspace(a,b,max(3,int((b-a)*10)+1)); vals = sol.sol(grid)
            aa = np.array([M.env(vals[:,j],grid[j])["aa_input"] for j in range(len(grid))], float)
            traces.append(pd.DataFrame({"daytime_h":grid-24.0*day,
                                        "aa_input":aa,"SAM_uM":vals[IS]}))
        y = sol.y[:,-1]
    return y, (pd.concat(traces, ignore_index=True) if traces else None)


def run(n_warmup_days=26):
    DATA.mkdir(exist_ok=True); y = M.y0().copy()
    lr, la, tr, ta = 2e-8, 2e-10, 2e-10, 2e-12
    for day in range(n_warmup_days): y, _ = one_day(y, day, lr, la, dense=False)
    y0 = y.copy()
    y_loose, q_loose = one_day(y0, n_warmup_days, lr, la, dense=True)
    y_tight, q_tight = one_day(y0, n_warmup_days, tr, ta, dense=True)
    mismatch = float(np.max(np.abs(y_tight-y0)/np.maximum(np.abs(y0),1.0)))
    maxdiff = float(abs(q_tight.SAM_uM.max()-q_loose.SAM_uM.max()))
    if not maxdiff < 3e-7: raise RuntimeError(f"loose/tight SAM peak difference {maxdiff:g}")
    imin = int(q_tight.SAM_uM.idxmin()); imax = int(q_tight.SAM_uM.idxmax())
    summary = {
      "warmup_days": n_warmup_days,
      "full_state_start_end_relative_mismatch_tight": float(f"{mismatch:.2g}"),
      "SAM_min_uM": float(q_tight.SAM_uM.min()),
      "SAM_min_time_h": float(q_tight.loc[imin,"daytime_h"]),
      "SAM_max_uM": float(q_tight.SAM_uM.max()),
      "SAM_max_time_h": float(q_tight.loc[imax,"daytime_h"]),
      "SAM_midnight_start_uM": float(q_tight.SAM_uM.iloc[0]),
      "SAM_midnight_end_uM": float(q_tight.SAM_uM.iloc[-1]),
      "SAM_max_abs_difference_loose_vs_tight_uM_upper_bound": 3e-7,
      "loose_rtol":lr,"loose_atol":la,"tight_rtol":tr,"tight_atol":ta,
      "sustained_eta_boundary":2.44486}
    q_tight.to_csv(DATA/"native_meal_cycle_day26_tight.csv", index=False)
    (DATA/"native_meal_cycle_summary_v10_2.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2)); return q_tight, summary


if __name__ == "__main__": run()
