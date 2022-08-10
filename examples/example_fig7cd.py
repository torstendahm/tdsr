#!/usr/bin/env python
# -*- coding: utf-8 -*-

# --------------------------
# Example Fig 7cd in Dahm and Hainzl - Morsleben induced seismicity case
# --------------------------

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

EXAMPLE_DIR = Path(__file__).parent
REPO_ROOT = EXAMPLE_DIR.parent.absolute()
sys.path.insert(0, str(REPO_ROOT))

from tdsr import TDSR1, Traditional  # noqa: E402
from tdsr.loading import CustomLoading  # noqa: E402

data_cfs = REPO_ROOT / "data/stresschange_morsleben.dat"
fn_obsrate = REPO_ROOT / "data/morsleben_EQ.dat"
figname = REPO_ROOT / "plots/fig7cd"

tdsr = TDSR1()
trad = Traditional()

km = 1000.0
days = 1.0

tstart = 0.0 * days
tend = 100.0 * days
deltat = 0.01 * days
nt = np.floor((tend - tstart) / deltat).astype(int)
tgrid = np.linspace(tstart, tend, nt)
dt = np.ediff1d(tgrid, to_end=tgrid[-1] - tgrid[-2])

t0 = 0.01
# no tectonic loading
strend = 0.0
# susceptibility to trigger earthquakes by unit stress increase
chi0 = 1.0

# increment do discretize Coulomb stress axis
deltaS = 0.3 / 500.0
# maximum depth on Coulomb axis (limit of integral)
sigma_max = 10000.0 * deltaS

# uniform distribution (and strend=0)
iX0 = "uniform"

# data provided in Pa, to be scaled to  MPa
scal_cf = 1.0e-6
# time provilded in units of days, to be scaled to days
scal_t = 1.0
c_tstart = 0.0

# read observed data again to get some paramter for plotting and scaling
data = np.loadtxt(fn_obsrate)
data = data[((data[:, 1] >= tstart) & (data[:, 1] <= tend))]
teq = data[:, 1]
Neq = data[:, 0]
NEQ = np.sum(Neq)
# Sebastion used the full EQ catalog
# I constrained in my input only to the Southern region
# NEQ = 70446.
dteq = np.ediff1d(teq, to_end=teq[-1] - teq[-2])
Req = Neq / dteq

# read morsleben data
morsleben_data = np.loadtxt(data_cfs)

# calculate earthquake rates with tdsm, lcm and rsm
ns = 3
cfs = np.zeros(nt)
r_tdsr = np.zeros((ns, nt))

for i in range(3):
    if i == 0:
        depthS = -0.3
        Sshadow = 3.0
    if i == 1:
        depthS = -0.4
        Sshadow = 3.5
    if i == 2:
        depthS = -0.5
        Sshadow = 4.5

    loading = CustomLoading(
        data=morsleben_data[:, [1, 0]],
        scal_t=scal_t,
        scal_cf=scal_cf,
        strend=strend,
        c_tstart=c_tstart,
        tstart=tstart,
        tend=tend,
        deltat=deltat,
    )
    t, chiz, cfs, r, xn = tdsr(
        loading=loading,
        chi0=chi0,
        t0=t0,
        depthS=depthS,
        Sshadow=Sshadow,
        deltaS=deltaS,
        sigma_max=sigma_max,
        iX0=iX0,
        deltat=deltat,
        taxis_log=0,
        tstart=tstart,
        tend=tend,
    )
    X0 = NEQ / np.sum(r * dt)
    r *= X0  # X0 set to match the observed number
    r_tdsr[i, :] = r[:]

# calculate stress shadow ad hoc  by using the linear coulomb failure class
t, cf_shad, cf, r, xn = trad(
    loading=loading, chi0=chi0, deltat=deltat, tstart=tstart, tend=tend
)

# plot results

cb = ["r", "b", "m", "r", "b", "m", "r", "b", "m"]
fac = 1.0 / 24.0  # scaling rate from days to hours
plt.rc("font", family="sans-serif")
plt.rc("font", size=12)
plt.rc("legend", fontsize=8)

fig, ax = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 1]}, figsize=(7, 8))
plt.subplots_adjust(hspace=0, wspace=0.2)

ax[0].plot(t, cfs - cfs[0], c="k", ls="solid", lw=1, alpha=0.9)
ax[0].plot(t, cf_shad, c="k", ls="dotted", lw=1, alpha=0.3)
ax[0].set_xlim(tstart, tend)
ax[0].set_ylim(0, 6.5)
ax[0].set_ylabel(r"$\sigma_c$    [MPa]", labelpad=10)
ax[0].set_xticks([])

ax[1].fill_between(teq, 0, fac * Req, color="k", lw=0, alpha=0.2)
for i in range(3):
    if i == 0:
        dsig = 0.3
        t0 = 0.01
        dS0 = 3.0
    if i == 1:
        dsig = 0.4
        t0 = 0.01
        dS0 = 3.5
    if i == 2:
        dsig = 0.5
        t0 = 0.01
        dS0 = 4.5
    ax[1].plot(
        t[:],
        fac * r_tdsr[i, :],
        c=cb[i],
        lw=1,
        label=r"$\delta\sigma=%.1f$ $\zeta_\mathrm{min}=%.1f$" % (dsig, dS0),
    )

ax[1].set_xlim(tstart, tend)
ax[1].set_ylim(
    0,
)
ax[1].set_xlabel("Time   [days]")
ax[1].set_ylabel("Rate   [#/hour]")
ax[1].legend()
ax[0].text(
    -0.14,
    0.93,
    "(c)",
    fontsize=14,
    color="k",
    horizontalalignment="left",
    transform=ax[0].transAxes,
)
ax[1].text(
    -0.14,
    0.95,
    "(d)",
    fontsize=14,
    color="k",
    horizontalalignment="left",
    transform=ax[1].transAxes,
)


plt.show()
fig.savefig(str(figname) + ".pdf", dpi=300, format="pdf", bbox_inches="tight")
fig.savefig(str(figname) + ".png", dpi=300, format="png", bbox_inches="tight")
print(f"\n\t OUTPUT: {str(figname) + '.pdf'}\n")
