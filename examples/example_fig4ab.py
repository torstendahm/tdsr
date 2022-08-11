#!/usr/bin/env python
# -*- coding: utf-8 -*-

# --------------------------
# Plot Fig. 4a and b (see Dahm and Hainzl, 2022)
# --------------------------

import sys
from pathlib import Path
import numpy as np
import math
import matplotlib.pyplot as plt

EXAMPLE_DIR = Path(__file__).parent
REPO_ROOT = EXAMPLE_DIR.parent.absolute()
sys.path.insert(0, str(REPO_ROOT))

from tdsr import TDSR1, CFM, RSD1  # noqa: E402
from tdsr.loading import StepLoading  # noqa: E402
from tdsr.utils import Eq7  # noqa: E402

figname = REPO_ROOT / "plots/fig4ab"

tdsr = TDSR1()
cfm = CFM()
rsd = RSD1()

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
# stress steps applied at time tstep, in MPa
# formerly Sstepvalues = (-2.0, -4.0)
sstep = [
    -2.0,
    -4.0,
]
tstep = 0.0 * hours
tstep = 0.01 * hours

# simulation start
tstart = -1.0 * hours
# simulation end
tend = 10.0 * hours

deltat = 0.02 * hours
t0 = deltat
# =NT = len(t)
nt = int(math.ceil((tend - tstart) / deltat))
# linear time axis discretization
taxis_log = False

# discretizing stress axis for integration
# increment do discretize Coulomb stress axis
deltaS = -depthS / 60.0
# maximum depth on Coulomb axis (limit of integral)
sigma_max = 3000.0 * deltaS

# calculare earthquake rates
cfs = np.zeros((2, nt))
cf_shad = np.zeros((2, nt))
r_tdsr = np.zeros((2, nt))
r_lcm = np.zeros((2, nt))
r_rsd = np.zeros((2, nt))
R_TDSR_theory = np.zeros((2, nt))

for i, Sstep in enumerate(sstep):
    print("i=", i, " Sstep=", Sstep, " sstep[i]=", sstep[i], " tstep=", tstep)
    print("deltat=", deltat, " tstart=", tstart, " tend=", tend)
    loading = StepLoading(
        strend=strend,
        sstep=Sstep,
        tstep=tstep,
        taxis_log=taxis_log,
        deltat=deltat,
        tstart=tstart,
        tend=tend,
    )
    common = dict(
        loading=loading,
        chi0=chi0,
        t0=t0,
        depthS=depthS,
        deltaS=deltaS,
        sigma_max=sigma_max,
        iX0="equilibrium",
        deltat=deltat,
        taxis_log=taxis_log,
        tstart=tstart,
        tend=tend,
    )

    t, chiz, cf, r, xn = tdsr(**common)
    cfs[i, :] = cf[:]
    r_tdsr[i, :] = r[:]

    t, cf_shadow, cf, r, xn = rsd(**common)
    r_rsd[i, :] = r

    t, cf_shadow, cf, r, xn = cfm(**common)
    cf_shad[i, :] = cf_shadow[:]
    r_lcm[i, :] = r

    R_TDSR_theory[i, :] = Eq7(t, sstep[i], chi0 * strend, -depthS, strend, strend)

# plot results

# make a bigger global font size and sans-serif style
lstype = ["solid", "dashed"]
plt.rc("font", family="sans-serif")
plt.rc("font", size=12)
plt.rc("legend", fontsize=10)
fig, ax = plt.subplots(2, 1, figsize=(7, 6))
plt.subplots_adjust(hspace=0.0)

al = (0.8, 0.2)
good = t > 0.0
for i, Sstep in enumerate(sstep):
    ax[0].plot(t, cf_shad[i, :], c="k", ls="dotted", lw=1, alpha=al[i])
    ax[0].plot(
        t,
        cfs[i, :],
        c="k",
        ls="solid",
        lw=2,
        alpha=al[i],
        label=r"$\Delta\sigma_c / \delta\sigma$=%.1f" % (Sstep),
    )

ax[0].axhline(0, c="gray", ls="dotted", lw=1)
ax[0].set_xticks([])
ax[0].set_xlim(tstart, tend)
ax[0].set_ylim(-5, 1.1 * (strend * tend))
ax[0].set_ylabel(r"$\sigma_c(t) - \sigma_c(0)$")
ax[0].legend()
ax[0].text(
    -0.18,
    0.93,
    "(a)",
    fontsize=14,
    color="k",
    horizontalalignment="left",
    transform=ax[0].transAxes,
)

for i, Sstep in enumerate(sstep):
    if i == 0:
        ax[1].plot(
            t, r_tdsr[i, :], c="b", ls="solid", lw=3, alpha=al[i], label="TDSR, RS"
        )
        ax[1].plot(t, r_lcm[i, :], c="g", ls="dashed", lw=1, alpha=al[i], label="CF")
    else:
        ax[1].plot(t, r_tdsr[i, :], c="b", ls="solid", lw=3, alpha=al[i])
        ax[1].plot(t, r_lcm[i, :], c="g", ls="dashed", lw=1, alpha=al[i])
    ax[1].plot(t[good], R_TDSR_theory[i, good], c="r", ls="dashed", lw=1)

ax[1].set_xlim(tstart, tend)
ax[1].set_ylim(0, 1.1 * chi0 * strend)
ax[1].set_xlabel(r"Time  $t$ / ($\delta\sigma / \dot\sigma_c$)")
ax[1].set_ylabel("Seismicity rate  $R$ / $r_0$")
ax[1].legend()
ax[1].text(
    -0.18,
    0.93,
    "(b)",
    fontsize=14,
    color="k",
    horizontalalignment="left",
    transform=ax[1].transAxes,
)

plt.show()
fig.savefig(str(figname) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
fig.savefig(str(figname) + ".png", format="png", dpi=300, bbox_inches="tight")
print(f"\n\t OUTPUT: {str(figname) + '.pdf'}\n")
