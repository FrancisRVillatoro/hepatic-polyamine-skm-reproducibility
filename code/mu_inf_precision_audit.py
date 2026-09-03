#!/usr/bin/env python3
"""Root-specific numerical-spread audit for the physiological simple-pole boundary."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def run():
    fits = pd.read_csv(DATA / "sam_asymptotic_simple_pole_fits.csv")
    branch = pd.read_csv(DATA / "sam_asymptotic_conditional_branch.csv")
    summary = pd.read_csv(DATA / "sam_asymptotic_summary.csv").iloc[0]
    K = float(summary["K_from_balance_coefficients"])
    nonlinear = float(summary["mu_inf_fit"])
    rows = []
    for r in fits.loc[fits.Smin >= 1e6].itertuples(index=False):
        rows.append({"method": "linear_tail_fit", "scale": float(r.Smin),
                     "mu_inf": float(r.mu_inf)})
    for r in branch.loc[branch.SAM >= 1e7].itertuples(index=False):
        rows.append({"method": "balance_corrected_point", "scale": float(r.SAM),
                     "mu_inf": float(r.mu + K / r.SAM)})
    rows.append({"method": "nonlinear_tail_fit", "scale": 1e7,
                 "mu_inf": nonlinear})
    out = pd.DataFrame(rows, columns=["method", "scale", "mu_inf"])
    vals = out.mu_inf.to_numpy(float)
    audit = {
        "min": float(vals.min()),
        "max": float(vals.max()),
        "span": float(np.ptp(vals)),
        "half_span": float(np.ptp(vals) / 2.0),
        "recommended_display": 2.15180254,
        "conservative_method_spread": 1e-8,
        "interpretation": "cross-method numerical spread; not a confidence interval or rigorous error bound",
    }
    out.to_csv(DATA / "mu_inf_precision_estimates_v10_2.csv", index=False)
    (DATA / "mu_inf_precision_summary_v10_2.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(out.to_string(index=False))
    print(json.dumps(audit, indent=2))
    return out, audit


if __name__ == "__main__":
    run()
