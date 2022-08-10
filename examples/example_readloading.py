#!/usr/bin/env python
# -*- coding: utf-8 -*-

# --------------------------
# Example how to read in a Coulomb stress loading function
# see supplemenrtary figure in Dahm and Hainzl (2022) for complex CSF evolution
# --------------------------

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

EXAMPLE_DIR = Path(__file__).parent
REPO_ROOT = EXAMPLE_DIR.parent.absolute()
sys.path.insert(0, str(REPO_ROOT))

from tdsr import TDSR, TDSR1, RSD1  # noqa: E402
from tdsr.loading import CustomLoading  # noqa: E402

figname = REPO_ROOT / "plots/fig_custom_loading"
data_file = REPO_ROOT / "data/CFSloading_synthetic.dat"

tdsm = TDSR()
tdsr = TDSR1()
rsd1 = RSD1()

hours = 1.0
tstart = 0.0 * hours
tend = 10 * hours
deltat = 0.001 * hours
t0 = 1.0 * deltat
nt = np.floor((tend - tstart) / deltat).astype(int)

strend = 1.0
# susceptibility to trigger earthquakes by unit stress increase
chi0 = 1.0
# skin depth in MPa (must be negative)
depthS = -1.0

# increment do discretize Coulomb stress axis
deltaS = -depthS / 500.0
# maximum depth on Coulomb axis (limit of integral)
sigma_max = 10000.0 * deltaS

# linear time axis discretisation
taxis_log = False

# steady state distribution
iX0 = "equilibrium"

# data provided in Pa, to be scaled to  MPa
scal_cf = 1.0
# time provilded in units of days, to be scaled to seconds
scal_t = 1.0
c_tstart = 0.0

# calculate stress loading function
dsig = -depthS
t = np.linspace(tstart, tend, nt)
dS = np.zeros(nt)

N1 = int(0.1 * nt)
N2 = int(0.3 * nt)
N3 = int(0.5 * nt)
N4 = int(0.6 * nt)
N5 = int(0.8 * nt)

dS[:N1] = strend * deltat
dS[N2:N3] = 3 * strend * deltat
dS[N4:] = 0.5 * strend * deltat
dS[N1] += 5 * dsig
dS[N2] -= 5 * dsig
dS[N4] += 2 * dsig
S = np.cumsum(dS)

# calculate earthquake rates with tdsm, lcm and rsm
cfs = np.zeros(nt)
r_tdsr = np.zeros(nt)
r_rsd = np.zeros(nt)

loading = CustomLoading(
    data=np.transpose([t, S]),
    scal_t=scal_t,
    scal_cf=scal_cf,
    strend=strend,
    c_tstart=c_tstart,
    tstart=tstart,
    tend=tend,
    deltat=deltat,
)
common = dict(
    loading=loading,
    chi0=chi0,
    depthS=depthS,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
    taxis_log=taxis_log,
    iX0=iX0,
    deltaS=deltaS,
    sigma_max=sigma_max,
    t0=t0,
)
t, chiz, cf, r_tdsr, xn = tdsr(**common)
t, chiz, cf, r_rsd, xn = rsd1(**common)

# plot results
plt.rc("font", family="sans-serif")
plt.rc("font", size=12)
plt.rc("legend", fontsize=10)
fig, ax = plt.subplots(2, 1, figsize=(7, 6))
plt.subplots_adjust(hspace=0.0)

ax[0].plot(t, S, c="k")
ax[0].set_xlim(tstart, tend)
ax[0].set_ylim(1.1 * np.min(S), 1.1 * np.max(S))
ax[0].set_ylabel(r"$\sigma_c(t) - \sigma_c(0)$")
ax[0].set_xticks([])
ax[1].plot(t, r_tdsr, c="b", lw=3, alpha=0.4, label="TDSR")
ax[1].plot(t, r_rsd, c="r", ls="dashed", lw=1, label="RS")
ax[1].set_xlim(tstart, tend)
ax[1].set_yscale("log")
ax[1].set_xlabel("Time")
ax[1].set_ylabel("Seismicity rate  $R$ / $r_0$")
ax[1].legend()

plt.show()
fig.savefig(str(figname) + ".pdf", dpi=300, format="pdf", bbox_inches="tight")
fig.savefig(str(figname) + ".png", dpi=300, format="png", bbox_inches="tight")
print(f"\n\t OUTPUT: {str(figname) + '.pdf'}\n")
