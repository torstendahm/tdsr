#!/usr/bin/env python
# -*- coding: utf-8 -*-

# --------------------------
# Plot Fig. 4c and d (see Dahm and Hainzl, 2022)
# --------------------------

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

EXAMPLE_DIR = Path(__file__).parent
REPO_ROOT = EXAMPLE_DIR.parent.absolute()
sys.path.insert(0, str(REPO_ROOT))

from tdsr import TDSR1, CFM, RSD  # noqa: E402
from tdsr.loading import CyclicLoading  # noqa: E402

figname = REPO_ROOT / "plots/fig4cd"

tdsr = TDSR1()
cfm = CFM()
rsd = RSD()

# time unit
hours = 1.0
# (=-dsig) skin depth in MPa
# must be negativ, will be changed with next release
depthS = -1.0
# mean failure time in sec for critically stressed source
t0 = 1.0
# (=dotsigc) in MPa/timeunit: tectonic stressing rate before loading
strend = 1.0
# susceptibility to trigger critically stressed sources
# if a unit step increas is applied (r0 = chi0*strend)
chi0 = 1.0

# Amplitude of sinusoidal loading
ampsin = 5.0e0
Tsin = ampsin / strend
Tmax = 3.0 * Tsin

tstart = 0.0
tend = Tmax

deltat = 0.02 * hours
t0 = deltat
# =NT = len(t)
nt = int(np.ceil((tend - tstart) / deltat))
# linear time axis discretization
taxis_log = False

# discretizing stress axis for integration
# increment do discretize Coulomb stress axis
deltaS = -depthS / 60.0
# maximum depth on Coulomb axis (limit of integral)
sigma_max = 3000.0 * deltaS

# calculate earthquake rates for cyclic loading
loading = CyclicLoading(
    strend=strend,
    ampsin=ampsin,
    Tsin=Tsin,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)

common = dict(
    loading=loading,
    chi0=chi0,
    depthS=depthS,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
    taxis_log=False,
    deltaS=deltaS,
    sigma_max=sigma_max,
    t0=t0,
    iX0="equilibrium",
)

t, chiz, cf, r, xn = tdsr(**common)
cfs = cf
r_tdsr = r

t, chiz, cf, r, xn = cfm(**common)
cf_shad = chiz[0:nt]
r_cfm = r

t, chiz, cf, r, xn = rsd(**common)
r_rsd = r

# plot results

# make a bigger global font size and sans-serif style
plt.rc("font", family="sans-serif")
plt.rc("font", size=12)
plt.rc("legend", fontsize=10)
fig, ax = plt.subplots(2, 1, figsize=(7, 6))
plt.subplots_adjust(hspace=0.0)

ax[0].plot(t, cf_shad, c="k", ls="dotted", lw=1, alpha=0.3)
ax[0].plot(t, cf, c="k")
ax[0].set_xlim(0, Tmax)
ax[0].set_ylim(1.1 * np.min(cf), 1.1 * np.max(cf))
ax[0].set_ylabel(r"$\sigma_c(t) - \sigma_c(0)$")
ax[0].set_xticks([])
ax[0].text(
    -0.18,
    0.93,
    "(c)",
    fontsize=14,
    color="k",
    horizontalalignment="left",
    transform=ax[0].transAxes,
)

ax[1].plot(t, r_tdsr, c="b", lw=3, alpha=0.6, label="TDSR")
ax[1].plot(t, r_rsd, c="r", ls="dashed", lw=1, label="RS")
ax[1].plot(t, r_cfm, c="g", ls="dashed", lw=1, label="CF")

ax[1].set_xlim(0, Tmax)
ax[1].set_ylim(0, 1.07 * np.max((np.max(r_tdsr), np.max(r_rsd), np.max(r_cfm))))
ax[1].set_xlabel("Time")
ax[1].set_ylabel("Seismicity rate  $R$ / $r_0$")
ax[1].legend()
ax[1].text(
    -0.18,
    0.93,
    "(d)",
    fontsize=14,
    color="k",
    horizontalalignment="left",
    transform=ax[1].transAxes,
)

plt.show()
fig.savefig(str(figname) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
fig.savefig(str(figname) + ".png", format="png", dpi=300, bbox_inches="tight")
print(f"\n\t OUTPUT: {str(figname) + '.pdf'}\n")
