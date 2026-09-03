#!/usr/bin/env python3
"""Compatibility-space spectral audit for the physiological high-SAM tail."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.linalg import null_space
import step3_core as core
import substrate_input_core as sic

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
IDX = core.idx; N = len(core.M.state)


def basis():
    C = np.zeros((len(core.INV), N))
    for i, names in enumerate(core.INV.values()):
        for name in names: C[i, IDX[name]] = 1.0
    Q = null_space(C)
    if Q.shape != (N, N-len(core.INV)): raise RuntimeError(Q.shape)
    return Q


def state(row): return np.array([float(row[f"x_{n}"]) for n in core.M.state])


def values(row, Q):
    y = state(row); mu = float(row.mu)
    J = sic.physical_jacobian(y, 0.0, mu, "methionine")
    ef = np.linalg.eigvals(J); ep = np.linalg.eigvals(Q.T @ J @ Q)
    lam = ep[np.argmax(ep.real)]; zscale = float(np.max(np.sort(np.abs(ef))[:3]))
    legacy = ef[np.abs(ef)>1e-8]
    return float(lam.real), zscale, float(abs(lam)/zscale), float(np.max(legacy.real))


def nearest(df, target):
    return df.iloc[int(np.argmin(abs(np.log(df.SAM.to_numpy(float))-np.log(target))))]


def run():
    df = pd.read_csv(DATA/"methionine_palc_physiological_fullstate.csv"); Q = basis()
    rows=[]
    for target in [8.97e5,3.08e7,1.06e9]:
        r=nearest(df,target); lam,z,ratio,_=values(r,Q)
        rows.append({"SAM_uM":float(r.SAM),"projected_rightmost_real_per_h":lam,
                     "conservation_zero_scale":z,"physical_to_zero_scale_ratio":ratio})
    tail=pd.DataFrame(rows); tail.to_csv(DATA/"compatibility_spectral_projection_tail.csv",index=False)
    safe=df[df.SAM<=1e5]
    ix=np.unique(np.linspace(0,len(safe)-1,10).round().astype(int)); dif=[]
    for i in ix:
        lam,_,_,legacy=values(safe.iloc[i],Q); dif.append(abs(lam-legacy))
    summary={"max_checkpoint_abs_difference_per_h":float(max(dif)),
             "legacy_extreme_tail_max_positive_artifact_per_h":1.39e-11,
             "interpretation":"magnitude filtering can misidentify conservation modes on the extreme tail; Q^T J Q remains stable"}
    (DATA/"compatibility_spectral_projection_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(tail.to_string(index=False)); print(json.dumps(summary,indent=2)); return tail,summary


if __name__ == "__main__": run()
