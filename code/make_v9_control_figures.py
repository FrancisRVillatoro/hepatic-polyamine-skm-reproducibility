#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
DATA=HERE.parent/'data'
LABELS={
    'published':'published program',
    'constant-total':'constant nominal MAT capacity',
    'flux-matched':'baseline-flux-matched MATII',
}

def gains(which,out):
    fig,ax=plt.subplots(figsize=(6.2,4.5))
    for pol in LABELS:
        d=pd.read_csv(DATA/f'mat_control_gains_{pol}.csv')
        ax.plot(d.rho,d[which],marker='o',markersize=2.8,label=LABELS[pol])
    ylab=(r'$G_{\rm Met}=d\ln\mathrm{SAM}/d\ln\mu$' if which=='G_Met' else r'$G_{\rm AA}=d\ln\mathrm{SAM}/d\ln\eta$')
    ax.set_xlabel(r'Proliferative coordinate $\rho$')
    ax.set_ylabel(ylab)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True,alpha=.25)
    fig.tight_layout(); fig.savefig(out,bbox_inches='tight'); plt.close(fig)

def palc(out):
    fig,ax=plt.subplots(figsize=(6.2,4.5))
    for pol in LABELS:
        d=pd.read_csv(DATA/f'mat_control_methionine_palc_{pol}.csv')
        ax.plot(d.mu,d.SAM,marker='o',markersize=2.5,label=LABELS[pol])
    ax.set_yscale('log')
    ax.set_xlabel(r'Methionine-input multiplier $\mu$')
    ax.set_ylabel(r'Steady-state SAM ($\mu$M)')
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True,which='both',alpha=.25)
    fig.tight_layout(); fig.savefig(out,bbox_inches='tight'); plt.close(fig)

def simple_pole(out):
    d=pd.read_csv(DATA/'sam_asymptotic_conditional_branch.csv')
    s=pd.read_csv(DATA/'sam_asymptotic_summary.csv').iloc[0]
    q=d[(d.SAM>=1e6)&(d.SAM<=1e8)].copy()
    x=1/q.SAM.to_numpy(); y=float(s.mu_inf_fit)-q.mu.to_numpy()
    xx=np.logspace(np.log10(x.min()),np.log10(x.max()),200)
    yy=float(s.C_fit)*xx**float(s.p_fit)
    fig,ax=plt.subplots(figsize=(6.2,4.5))
    ax.loglog(x,y,'o',label='conditional equilibria')
    ax.loglog(xx,yy,label=fr'fit: slope $p={float(s.p_fit):.4f}$')
    ax.set_xlabel(r'$1/\mathrm{SAM}$ ($\mu$M$^{-1}$)')
    ax.set_ylabel(r'$\mu_\infty-\mu$')
    ax.legend(frameon=False)
    ax.grid(True,which='both',alpha=.25)
    fig.tight_layout(); fig.savefig(out,bbox_inches='tight'); plt.close(fig)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--outdir',default=str(HERE.parent/'figures_v9'))
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    gains('G_Met',out/'fig_mat_control_gain_met.pdf')
    gains('G_AA',out/'fig_mat_control_gain_aa.pdf')
    palc(out/'fig_mat_control_palc.pdf')
    simple_pole(out/'fig_sam_simple_pole.pdf')
    print('wrote',out)

if __name__=='__main__': main()
