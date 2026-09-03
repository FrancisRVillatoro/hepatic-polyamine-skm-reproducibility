#!/usr/bin/env python3
"""Reproduce post-v1.0.0 audits incorporated through manuscript v10.2."""
import argparse
import compatibility_spectral_projection as spec
import cusp_nondegeneracy_check as cusp
import mu_inf_precision_audit as mu
import native_meal_cycle_control as meal


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--skip-meal',action='store_true',help='skip the slower 26-day native meal-cycle integration')
    a=ap.parse_args()
    print('== compatibility-space spectral projection =='); spec.run()
    print('\n== reduced-cusp nondegeneracy =='); cusp.run()
    print('\n== mu_inf root-specific numerical spread =='); mu.run()
    if not a.skip_meal:
        print('\n== native physiological meal cycle =='); meal.run()


if __name__=='__main__': main()
