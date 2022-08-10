#!/usr/bin/env python
# -*- coding: utf-8 -*-

# --------------------------
# Plot Fig. 3b (see Dahm and Hainzl, 2022)
# --------------------------

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

EXAMPLE_DIR = Path(__file__).parent
REPO_ROOT = EXAMPLE_DIR.parent.absolute()
sys.path.insert(0, str(REPO_ROOT))

from tdsr import TDSR1  # noqa: E402
from tdsr.loading import BackgroundLoading  # noqa: E402
from tdsr.utils import Eq7  # noqa: E402

figname = REPO_ROOT / "plots/fig3b"

tdsr = TDSR1()

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

# stress steps in e.g. MPa
sstep = np.asarray([2.0, 4.0, 6.0, 8.0, 10.0])
# time of step must be zero if logarithmic time axis sampling
tstep = 0.0 * hours

# if logarithmic time axis, tstart is overwritten by tstartl
tstart = 0.0
tend = 100.0 * hours

# obsolete if logarithmic time axis
deltat = 0.02 * hours
# =NT = len(t) - overwritten if logarithmic time axis
nt = int(np.ceil((tend - tstart) / deltat))
# logarithmic time axis discretization
taxis_log = True
ntlog = 10000
tstartl = 1.0e-6 * hours

# discretizing stress axis for integration
# increment do discretize Coulomb stress axis
deltaS = -depthS / 60.0
# maximum depth on Coulomb axis (limit of integral)
sigma_max = 3000.0 * deltaS

# steady state distribution
iX0 = "equilibrium"

# calculate earthquake rates for cyclic loading
ns = len(sstep)
r_tdsr = np.zeros((ns, ntlog))

for k in range(ns):
    loading = BackgroundLoading(
        strend=strend,
        sstep=sstep[k],
        taxis_log=taxis_log,
        ntlog=ntlog,
        tstart=tstartl,
        tend=tend,
    )
    t, chiz, cf, r, xn = tdsr(
        loading=loading,
        chi0=chi0,
        t0=t0,
        depthS=depthS,
        deltaS=deltaS,
        sigma_max=sigma_max,
        iX0=iX0,
        Sshadow=sstep[k],
        taxis_log=taxis_log,
        ntlog=ntlog,
        tstart=tstartl,
        tend=tend,
    )
    r_tdsr[k, :] = r

# plot results
Tmin = 2e-5
Tmax = 10
al = [0.1, 0.3, 0.5, 0.7, 0.9]

# make a bigger global font size and sans-serif style
plt.rc("font", family="sans-serif")
plt.rc("font", size=12)
plt.rc("legend", fontsize=10)
fig, ax = plt.subplots(1, 1, figsize=(7, 5))

for k, Sstep in enumerate(sstep):
    ax.plot(
        t,
        r_tdsr[k, :],
        c="b",
        ls="solid",
        lw=3,
        alpha=al[k],
        label=r"$\Delta\sigma_c / \delta\sigma$=%.1f" % (Sstep),
    )
    r_theo = Eq7(t, sstep[k], chi0 * strend, -depthS, strend, strend)
    ax.plot(t, r_theo, c="r", ls="dashed", lw=1, alpha=1, zorder=10)

ax.axhline(1, c="k", ls="dotted", lw=1, alpha=0.1)
ax.set_xlim(Tmin, Tmax)
ax.set_ylim(0.5, 1.1 * np.max(chi0 * strend * np.exp(sstep / -depthS)))
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"Time  $t$ / ($\delta\sigma / \dot\sigma_c$)")
ax.set_ylabel(r"Rate $R$ / $r_0$", labelpad=-3)
ax.legend()

plt.show()
fig.savefig(str(figname) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
fig.savefig(str(figname) + ".png", format="png", dpi=300, bbox_inches="tight")
print(f"\n\t OUTPUT: {str(figname) + '.pdf'}\n")
