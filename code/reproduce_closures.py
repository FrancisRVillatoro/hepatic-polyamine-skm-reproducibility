#!/usr/bin/env python3
"""Run the paper reproducibility workflows.

Default mode runs the fast deterministic audits/sensitivity calculations.
Use --full to rerun Monte Carlo, continuation, path/cube, timescale, and
finite-rate hysteresis calculations.
"""
from pathlib import Path
import subprocess,sys
HERE=Path(__file__).resolve().parent
quick=[
    ["audit_BIOMD0000000190.py","BIOMD0000000190_url.xml"],
    ["slow_feedback_sensitivity.py"],
    ["substrate_input_gains.py"],
]
full=[
    ["skm_BIOMD0000000190.py"],
    ["directional_skm_closure.py"],
    ["audit_control_sweeps.py"],
    ["homotopy_subset.py"],
    ["cube_smooth_5.py"],
    ["hface_stability.py"],
    ["palc47_full.py"],
    ["tau_sweep_fast.py"],
    ["external_audit_checks.py"],
    ["dynamic_hysteresis_47d.py"],
]
for args in quick+(full if "--full" in sys.argv else []):
    cmd=[sys.executable,str(HERE/args[0])]+[str(HERE/a) if a.endswith('.xml') else a for a in args[1:]]
    print('RUN',*args,flush=True); subprocess.run(cmd,check=True)
print('DONE')
