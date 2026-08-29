#!/usr/bin/env python3
"""Shared autonomous substrate-input controls for the 2012 hepatic model.

This module extends ``step3_core`` without changing its kinetic equations.
Two frozen input protocols are implemented:

* ``methionine``: vary ``b_met_basal = 30 * p`` while all meal multipliers are 1;
* ``common-aa``: set fasting/breakfast/lunch/dinner all equal to ``p`` while
  keeping ``b_met_basal = 30``.

The proliferative coordinate rho simultaneously activates the published MAT,
inverse-SAM ODC/SAMDC, and H2O2 changes, exactly as in ``step3_core``.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import least_squares
import step3_core as core

N = len(core.M.state)
IDX = core.idx


def input_overrides(rho: float, p: float, protocol: str) -> dict:
    if protocol == "methionine":
        fasting = breakfast = lunch = dinner = 1.0
        b_met_basal = 30.0 * p
    elif protocol in {"common-aa", "common_aa", "aa"}:
        fasting = breakfast = lunch = dinner = p
        b_met_basal = 30.0
    else:
        raise ValueError(f"unknown protocol {protocol!r}")
    return {
        "Vm_MAT1": 260.0 * (1.0 - rho),
        "Vm_MAT3": 220.0 * (1.0 - rho),
        "parameter_19": 220.0 * rho,
        "H2O2": 0.01 * (1.0 + 0.5 * rho),
        "fasting": fasting,
        "breakfast": breakfast,
        "lunch": lunch,
        "dinner": dinner,
        "b_met_basal": b_met_basal,
    }


def rhs_input(y: np.ndarray, rho: float, p: float, protocol: str) -> np.ndarray:
    ov = input_overrides(rho, p, protocol)
    e = core.M.env(y, 0.0, ov)
    fl = core.M.flux(y, 0.0, ov)
    d = {n: 0.0 for n in core.M.state}
    for r in core.M.rx:
        for s, nu in r["st"].items():
            if s in d:
                d[s] += nu * fl[r["id"]] / core.M.comp[core.M.spec[s]["comp"]]

    D, S = e["species_4"], e["species_3"]
    free = 1.0 / (1.0 + e["parameter_5"] * (D + S))
    g = (1.0 - rho) + rho * 65.06 / e["sam"]
    d["parameter_1"] = 60.0 * g * e["parameter_7"] * free - e["parameter_6"] * e["parameter_4"] * e["parameter_1"]
    d["parameter_2"] = 60.0 * e["parameter_9"] * (1.0 - free) - e["parameter_8"] * free * e["parameter_2"]
    d["parameter_3"] = 60.0 * g * e["parameter_11"] * free - e["parameter_10"] * e["parameter_3"]
    d["parameter_4"] = e["parameter_13"] * (1.0 - 1.0 / (1.0 + e["parameter_5"] * 0.01 * (D + S))) - e["parameter_12"] * e["parameter_4"]
    return np.array([d[n] for n in core.M.state])


def constrained_residual(logy: np.ndarray, rho: float, p: float, protocol: str) -> np.ndarray:
    y = np.exp(logy)
    r = rhs_input(y, rho, p, protocol) / np.maximum(y, 1.0)
    for s, k in core.REPL.items():
        r[IDX[s]] = (sum(y[IDX[n]] for n in core.INV[k]) - core.target[k]) / core.target[k]
    return r


def solve_equilibrium(y0: np.ndarray, rho: float, p: float, protocol: str,
                      tol: float = 3e-11, max_nfev: int = 200):
    q = least_squares(
        lambda z: constrained_residual(z, rho, p, protocol),
        np.log(np.asarray(y0, float)), xtol=tol, ftol=tol, gtol=tol,
        max_nfev=max_nfev, x_scale="jac"
    )
    y = np.exp(q.x)
    return y, float(np.max(np.abs(q.fun))), int(q.nfev)


def basal_state(rho: float, protocol: str = "methionine") -> np.ndarray:
    """Track the p=1 equilibrium from the physiological state to requested rho."""
    y = core.M1.y0().copy()
    if rho == 0:
        y, _, _ = solve_equilibrium(y, 0.0, 1.0, protocol)
        return y
    nstep = max(2, int(np.ceil(abs(rho) / 0.05)))
    for rr in np.linspace(0.0, rho, nstep + 1):
        y, res, _ = solve_equilibrium(y, float(rr), 1.0, protocol)
        if res > 2e-7:
            raise RuntimeError(f"equilibrium residual {res:g} at rho={rr:g}")
    return y


def fd_log_jacobian(y: np.ndarray, rho: float, p: float, protocol: str, eps: float = 1e-6):
    z = np.log(y)
    J = np.empty((N, N))
    for j in range(N):
        zp = z.copy(); zm = z.copy()
        zp[j] += eps; zm[j] -= eps
        J[:, j] = (constrained_residual(zp, rho, p, protocol) -
                   constrained_residual(zm, rho, p, protocol)) / (2.0 * eps)
    pp = p * np.exp(eps); pm = p * np.exp(-eps)
    Jp = (constrained_residual(z, rho, pp, protocol) -
          constrained_residual(z, rho, pm, protocol)) / (2.0 * eps)
    return J, Jp


def logarithmic_sam_gain(y: np.ndarray, rho: float, p: float, protocol: str) -> float:
    J, Jp = fd_log_jacobian(y, rho, p, protocol)
    dz = np.linalg.solve(J, -Jp)
    return float(dz[IDX["sam"]])


def physical_jacobian(y: np.ndarray, rho: float, p: float, protocol: str, h: float = 1e-30):
    J = np.empty((N, N))
    for j in range(N):
        z = np.asarray(y, complex).copy()
        z[j] += 1j * h
        J[:, j] = np.imag(rhs_input(z, rho, p, protocol)) / h
    return J


def dominant_nonconservation_eigenvalue(y: np.ndarray, rho: float, p: float, protocol: str):
    ev = np.linalg.eigvals(physical_jacobian(y, rho, p, protocol))
    ev = ev[np.abs(ev) > 1e-8]
    return ev[np.argmax(ev.real)]


def equilibrium_log_jacobian(y: np.ndarray, rho: float, p: float, protocol: str,
                             scale: np.ndarray | None = None, h: float = 1e-30):
    """Jacobian of a fixed-scaled constrained equilibrium system in log variables.

    At an equilibrium, replacing the dynamic-row scaling by any fixed positive
    scale does not alter the nullspace or implicit logarithmic gain.  This
    implementation uses complex-step derivatives of the physical vector field
    and exact derivatives of the three conservation constraints.
    """
    y=np.asarray(y,float)
    if scale is None: scale=np.maximum(y,1.0)
    scale=np.asarray(scale,float)
    Jphys=physical_jacobian(y,rho,p,protocol,h=h)
    Jx=(Jphys/scale[:,None])*y[None,:]
    # log-p derivative by complex step
    pc=complex(p, p*h)
    rp=rhs_input(y.astype(complex),rho,pc,protocol)
    Jp=np.imag(rp)/h/scale
    for s,k in core.REPL.items():
        i=IDX[s]; Jx[i,:]=0.0; Jp[i]=0.0
        for n in core.INV[k]: Jx[i,IDX[n]]=y[IDX[n]]/core.target[k]
    return Jx,Jp


def logarithmic_sam_gain_fast(y: np.ndarray, rho: float, p: float, protocol: str) -> float:
    J,Jp=equilibrium_log_jacobian(y,rho,p,protocol)
    dz=np.linalg.solve(J,-Jp)
    return float(dz[IDX['sam']])
