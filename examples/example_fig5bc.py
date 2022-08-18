#!/usr/bin/env python
# -*- coding: utf-8 -*-

# --------------------------
# Plot Fig. 5bc (see Dahm and Hainzl, 2022)
# --------------------------

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

EXAMPLE_DIR = Path(__file__).parent
REPO_ROOT = EXAMPLE_DIR.parent.absolute()
DATA_DIR = REPO_ROOT / "data"
sys.path.insert(0, str(REPO_ROOT))

from tdsr import TDSR1, CFM, RSD, RSD1  # noqa: E402
from tdsr.loading import CustomLoading  # noqa: E402
from tdsr.utils import Zvalues, X0uniform, X0steady, X0gaussian

figname = REPO_ROOT / "plots/fig5bc"

tdsr = TDSR1()
cfm = CFM()
rsd = RSD()
rsd1 = RSD1()

# Asig = dsig in the RS model
dsig = 1.0
depthS = -dsig
# [MPa] gap for CF and RS_subcrit
dS0 = 8.5
# [MPa] gap for TDSR with uniform subcritical stress state
dS0_TDSR = 22.0
# [MPa] mean stress for TDSR with gaussian subcritical stress state
dS0mean = 24.5
# [MPa] standard deviation " " "
dS0std = 1.0
# [years] relaxation time ta=dsig/dotsigc
ta = 3e5
# [years] relaxation time in the subcritical RS model
taRSsub = 60.0
dotsigc = dsig / ta
strend = dotsigc
dotsigc_sub = dsig / taRSsub
# [years] mean failure time of a source at critical stress in the TDRS model
t0 = 1e-4

# discretizing stress axis for integration
# increment do discretize Coulomb stress axis
deltaS = -depthS / 60.0
# maximum depth on Coulomb axis (limit of integral)
sigma_max = 3000.0 * deltaS

# start time for model simulation
tstart = 1960
# end time for model simulation
tend = 2022.0
# magnitude cutoff for the fitted earthquakes
Mcut = 1.45
# start time period for fit
T1 = 1991
# end time period for fit
T2 = tend

# earthquake data
data = np.loadtxt(DATA_DIR / "groningenEQ.dat")
teq = data[:, 0]
# pressure data
data = np.loadtxt(DATA_DIR / "groningenP.dat")
tgrid = data[:, 0]
p = data[:, 1]

# coulomb stress function
S = -p
dS = np.ediff1d(S, to_begin=0.0)
S = np.cumsum(dS)

NEQ = len(teq[((teq >= T1) & (teq <= T2))])
dt = tgrid[1] - tgrid[0]
nt = int(np.ceil((tend - tstart) / dt))
t00 = np.min(tgrid[((S >= dS0) & (tgrid <= teq[0]))])
print(
    "time when EQ start, t00=",
    t00,
    tgrid[0],
    dotsigc * (t00 - tgrid[0]),
    dotsigc_sub * (t00 - tgrid[0]),
)

# loading assuming that faults under critical stress exists
common_loading = dict(
    data=np.transpose([tgrid, S]),
    scal_t=1.0,
    scal_cf=1.0,
    c_tstart=tstart,
    tstart=tstart,
    tend=tend,
    deltat=dt,
)
loading = CustomLoading(strend=dotsigc, **common_loading)

common = dict(
    chi0=1.0,
    t0=t0,
    depthS=depthS,
    deltat=dt,
    tstart=tstart,
    tend=tend,
    taxis_log=False,
    deltaS=deltaS,
    sigma_max=sigma_max,
)

# CF model with critical initial stress state:
t, chiz, cf, R_CF, xn = cfm(
    loading=loading,
    Sshadow=0.0,
    **common,
)
# renormalize to observed EQ-number within [T1, T2]:
facCFc = NEQ / np.sum(R_CF[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_CF *= facCFc

# CF model with subcritical initial stress state:
t, chiz, cf, R_CFsub, xn = cfm(
    loading=loading,
    Sshadow=dS0,
    **common,
)
facCFsub = NEQ / np.sum(R_CFsub[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_CFsub *= facCFsub

# RS model for critical initial stress state:
t, chiz, cf, R_RS, xn = rsd1(
    loading=loading,
    Sshadow=0.0,
    **common,
)
fac = NEQ / np.sum(R_RS[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_RS *= fac

# TDSR model for critical initial stress state:
t, chiz, cf, R_TDSR, xn = tdsr(
    loading=loading,
    Sshadow=0.0,
    **common,
)
facS = NEQ / np.sum(R_TDSR[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_TDSR *= facS

# assuming that no faults exists under critical stress
# (needed for TDSR and RS)
sub_loading = CustomLoading(strend=dotsigc_sub, **common_loading)

# RS model for subcritical initial stress state:
t, chiz, cf, R_RSsub, xn = rsd1(
    loading=sub_loading,
    Sshadow=dS0,
    **common,
)
fac = NEQ / np.sum(R_RSsub[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_RSsub *= fac

# TDSR model with uniform subcritical stress:
t, chiz, cf, R_TDSRsub_uniform, xn = tdsr(
    loading=sub_loading,
    Sshadow=dS0_TDSR,
    iX0="uniform",
    **common,
)
facSsub_uniform = NEQ / np.sum(R_TDSRsub_uniform[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_TDSRsub_uniform *= facSsub_uniform

# TDSR model with gaussian subcritical prestress:
t, chiz, cf, R_TDSRsub_gauss, xn = tdsr(
    loading=sub_loading,
    Sshadow=0.0,
    iX0="gaussian",
    Zmean=dS0mean,
    Zstd=dS0std,
    **common,
)
fac = NEQ / np.sum(R_TDSRsub_gauss[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_TDSRsub_gauss *= fac

# plot results
plt.rc("font", family="sans-serif")
plt.rc("font", size=12)
plt.rc("legend", fontsize=8)
fig, ax = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 2]}, figsize=(7, 6.5))
Z = Zvalues(S, 0.0, 0.0, dsig)
ax[0].plot(Z, facCFc*X0uniform(Z, 0.0, 1.0), c='k', ls='dashed', alpha=0.3, label=r'CF')
ax[0].plot(Z, facCFsub*X0uniform(Z, dS0, 1.0), c='k', alpha=0.3, label=r'CF$_\mathrm{subcrit}$')
ax[0].plot(Z, facS*X0steady(Z,1.0*dotsigc,t0,dsig,dotsigc), c='b', lw=3, alpha=0.3, label=r'TDRS$_\mathrm{crit}$')
ax[0].plot(Z, facSsub_uniform*X0uniform(Z, dS0_TDSR, 1.0), 'r--', label=r'TDSR$_\mathrm{uniform}$')
ax[0].plot(Z, fac*X0gaussian(Z, dS0mean, dS0std, 1.0), 'r-', label=r'TDSR$_\mathrm{gauss}$')
ax[0].set_xlim(-2, 30)
ax[0].set_ylim(0, 180)
ax[0].set_xlabel(r"$\zeta$    [MPa]", labelpad=1)
ax[0].set_ylabel(r"$\chi$   [#/MPa]")
ax[0].legend()

tbin = np.arange(1960, tend + 0.1, 1.0)
nobs, tbin1 = np.histogram(teq, bins=tbin)
ax[1].scatter(
    0.5 * (tbin[:-1] + tbin[1:]),
    nobs,
    c="k",
    marker="o",
    s=60,
    lw=3,
    alpha=0.5,
    label="observed",
)
ax[1].plot(t, R_CF, c="k", ls="dashed", alpha=0.3, label=r"CF")
ax[1].plot(t, R_CFsub, c="k", alpha=0.3, label=r"CF$_\mathrm{subcrit}$")
ax[1].plot(t, R_RS, c="b", lw=3, alpha=0.3, label=r"RS")
ax[1].plot(t, R_RSsub, "g--", label=r"RS$_\mathrm{subcrit}$")
ax[1].plot(t, R_TDSR, "b--", label=r"TDSR$_\mathrm{crit}$")
ax[1].plot(t, R_TDSRsub_uniform, "r--", label=r"TDSR$_\mathrm{uniform}$")
ax[1].plot(t, R_TDSRsub_gauss, "r-", label=r"TDSR$_\mathrm{gauss}$")
ax[1].set_xlim(1980, tend)
ax[1].set_xlabel("Time   [year]")
ax[1].set_ylabel("Rate   [#/year]")
ax[1].legend()
ax[0].text(
    -0.12,
    0.93,
    "(b)",
    fontsize=14,
    color="k",
    horizontalalignment="left",
    transform=ax[0].transAxes,
)
ax[1].text(
    -0.12,
    0.95,
    "(c)",
    fontsize=14,
    color="k",
    horizontalalignment="left",
    transform=ax[1].transAxes,
)

plt.show()
fig.savefig(str(figname) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
fig.savefig(str(figname) + ".png", format="png", dpi=300, bbox_inches="tight")
print(f"\n\t OUTPUT: {str(figname) + '.pdf'}\n")
