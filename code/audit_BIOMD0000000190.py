#!/usr/bin/env python3
"""
Exact reproduction/audit of BioModels BIOMD0000000190
Rodriguez-Caso2006_Polyamine_Metabolism.

Reads the SBML itself (no libSBML required), evaluates the MathML kinetic laws
and rate rules, integrates from the exact SBML initial conditions, reports the
corrected BioModels steady state, fluxes, residual, and local Jacobian spectrum.

Usage:
    python audit_BIOMD0000000190.py BIOMD0000000190_url.xml
"""

import sys
import xml.etree.ElementTree as ET
import numpy as np
from scipy.integrate import solve_ivp

SBML_NS = "http://www.sbml.org/sbml/level2/version3"
MATH_NS = "http://www.w3.org/1998/Math/MathML"
NS = {"s": SBML_NS, "m": MATH_NS}

TARGET_NAMES = [
    "P", "D", "S", "SAM", "A", "aD", "aS",
    "Antz", "Vmaxodc", "Vmaxssat", "Vmaxsamdc",
]
TARGET = np.array([
    104.681,
    76.7492,
    58.0135,
    52.327,
    0.0101962,
    0.832236,
    0.0245375,
    0.574038,
    1.28315,
    0.673814,
    0.36829,
], dtype=float)


def tag(x):
    return x.tag.split("}")[-1]


def mathml_to_python(node):
    t = tag(node)
    if t == "math":
        return mathml_to_python(list(node)[0])
    if t in ("ci", "cn"):
        return (node.text or "").strip()
    if t != "apply":
        raise NotImplementedError(f"Unsupported MathML element: {t}")

    children = list(node)
    op = tag(children[0])
    a = [mathml_to_python(x) for x in children[1:]]

    if op == "plus":
        return "(" + " + ".join(a) + ")"
    if op == "times":
        return "(" + " * ".join(a) + ")"
    if op == "minus":
        if len(a) == 1:
            return f"(-({a[0]}))"
        return "(" + " - ".join(a) + ")"
    if op == "divide":
        return f"(({a[0]}) / ({a[1]}))"
    if op == "power":
        return f"(({a[0]}) ** ({a[1]}))"
    if op == "exp":
        return f"np.exp({a[0]})"
    if op == "ln":
        return f"np.log({a[0]})"

    raise NotImplementedError(f"Unsupported MathML operator: {op}")


class SBMLModel:
    def __init__(self, filename):
        tree = ET.parse(filename)
        self.model = tree.getroot().find("s:model", NS)

        self.species = {}
        for s in self.model.findall("./s:listOfSpecies/s:species", NS):
            self.species[s.attrib["id"]] = {
                "init": float(s.attrib.get("initialConcentration", "0")),
                "boundary": s.attrib.get("boundaryCondition", "false") == "true",
                "constant": s.attrib.get("constant", "false") == "true",
            }

        self.parameters = {}
        for p in self.model.findall("./s:listOfParameters/s:parameter", NS):
            self.parameters[p.attrib["id"]] = float(p.attrib.get("value", "0"))

        self.assignment_rules = {}
        self.rate_rules = {}
        for rule in self.model.findall("./s:listOfRules/*", NS):
            expression = mathml_to_python(rule.find("m:math", NS))
            if tag(rule) == "assignmentRule":
                self.assignment_rules[rule.attrib["variable"]] = expression
            elif tag(rule) == "rateRule":
                self.rate_rules[rule.attrib["variable"]] = expression

        self.reactions = {}
        for r in self.model.findall("./s:listOfReactions/s:reaction", NS):
            rid = r.attrib["id"]
            local = {}
            for p in r.findall(
                "./s:kineticLaw/s:listOfParameters/s:parameter", NS
            ):
                local[p.attrib["id"]] = float(p.attrib["value"])

            stoich = {}
            for sr in r.findall(
                "./s:listOfReactants/s:speciesReference", NS
            ):
                sid = sr.attrib["species"]
                stoich[sid] = stoich.get(sid, 0.0) - float(
                    sr.attrib.get("stoichiometry", "1")
                )
            for sr in r.findall(
                "./s:listOfProducts/s:speciesReference", NS
            ):
                sid = sr.attrib["species"]
                stoich[sid] = stoich.get(sid, 0.0) + float(
                    sr.attrib.get("stoichiometry", "1")
                )

            self.reactions[rid] = {
                "expr": mathml_to_python(r.find("./s:kineticLaw/m:math", NS)),
                "local": local,
                "stoich": stoich,
            }

        dynamic_species = [
            sid for sid, info in self.species.items()
            if not info["boundary"] and not info["constant"]
        ]

        # Use the scientifically convenient order used in the paper/audit.
        preferred = ["P", "D", "S", "SAM", "A", "aD", "aS"]
        self.dynamic_species = [x for x in preferred if x in dynamic_species]
        self.dynamic_species += [
            x for x in dynamic_species if x not in self.dynamic_species
        ]

        preferred_rules = ["Antz", "Vmaxodc", "Vmaxssat", "Vmaxsamdc"]
        self.dynamic_rules = [x for x in preferred_rules if x in self.rate_rules]
        self.dynamic_rules += [
            x for x in self.rate_rules if x not in self.dynamic_rules
        ]

        self.state_names = self.dynamic_species + self.dynamic_rules

    def environment(self, y):
        env = {"np": np, "cytosol": 1.0}
        env.update({sid: info["init"] for sid, info in self.species.items()})
        env.update(self.parameters)

        for name, value in zip(self.state_names, y):
            env[name] = value

        # Assignment rules in this model are nonrecursive.
        for variable, expr in self.assignment_rules.items():
            env[variable] = eval(expr, {"np": np}, env)

        return env

    def fluxes(self, y):
        env = self.environment(y)
        values = {}
        for rid, reaction in self.reactions.items():
            renv = dict(env)
            renv.update(reaction["local"])
            values[rid] = eval(reaction["expr"], {"np": np}, renv)
        return values

    def rhs(self, t, y):
        env = self.environment(y)
        flux = self.fluxes(y)
        dy = {name: 0.0 * y[0] for name in self.state_names}

        for rid, reaction in self.reactions.items():
            v = flux[rid]
            for sid, nu in reaction["stoich"].items():
                if sid in dy:
                    dy[sid] += nu * v

        for variable, expr in self.rate_rules.items():
            dy[variable] = eval(expr, {"np": np}, env)

        return np.array([dy[name] for name in self.state_names])

    def initial_state(self):
        y0 = []
        for name in self.dynamic_species:
            y0.append(self.species[name]["init"])
        for name in self.dynamic_rules:
            y0.append(self.parameters[name])
        return np.array(y0, dtype=float)


def jacobian_complex_step(fun, y, h=1e-30):
    y = np.asarray(y, dtype=float)
    n = y.size
    J = np.empty((n, n), dtype=float)
    for j in range(n):
        z = y.astype(complex)
        z[j] += 1j * h
        J[:, j] = np.imag(fun(0.0, z)) / h
    return J


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: python audit_BIOMD0000000190.py BIOMD0000000190_url.xml"
        )

    m = SBMLModel(sys.argv[1])
    print("Model state:", ", ".join(m.state_names))
    print("Number of SBML reactions:", len(m.reactions))
    print("Number of rate-rule variables:", len(m.rate_rules))

    y0 = m.initial_state()
    print("\nExact SBML initial state:")
    for n, v in zip(m.state_names, y0):
        print(f"  {n:12s} {v:.15g}")

    sol = solve_ivp(
        m.rhs,
        (0.0, 2.0e5),
        y0,
        method="BDF",
        rtol=1e-12,
        atol=1e-14,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    yss = sol.y[:, -1]
    resid = np.max(np.abs(m.rhs(sol.t[-1], yss)))

    # Reorder against BioModels' corrected reported values.
    idx = {name: i for i, name in enumerate(m.state_names)}
    ytarget_order = np.array([yss[idx[n]] for n in TARGET_NAMES])
    rel = np.abs((ytarget_order - TARGET) / TARGET)

    print("\nCorrected basal steady state:")
    print(f"{'variable':12s} {'computed':>16s} {'BioModels':>16s} {'rel.err.':>14s}")
    for name, x, ref, err in zip(TARGET_NAMES, ytarget_order, TARGET, rel):
        print(f"{name:12s} {x:16.9g} {ref:16.9g} {err:14.3e}")
    print(f"\nmax relative error = {rel.max():.6e}")
    print(f"max ODE residual    = {resid:.6e}")

    print("\nSteady-state reaction fluxes:")
    for rid, v in m.fluxes(yss).items():
        print(f"  {rid:16s} {float(np.real(v)):.12g}")

    J = jacobian_complex_step(m.rhs, yss)
    eig = np.linalg.eigvals(J)
    eig = sorted(eig, key=lambda z: z.real, reverse=True)

    print("\nLocal kinetic Jacobian eigenvalues (1/min):")
    for z in eig:
        print(f"  {z.real:+.12e} {z.imag:+.12e}i")
    print(f"\nspectral abscissa = {eig[0].real:+.12e} 1/min")


if __name__ == "__main__":
    main()
