#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test figure 3a"""

import inspect
import numpy as np
import matplotlib.pyplot as plt

from .utils import dump_values, load_values, PLOT_DIR
from tdsr import CFM, RSD1, TDSR1
from tdsr.loading import BackgroundLoading
from tdsr.utils import Eq5


def plot_fig3a(out, t, t0, Sshadow, strend, r_tdsr, r_cfm, r_rsd, chi0, depthS):
    Tmin = 0
    Tmax = 70
    cb = ["gray", "b", "g"]
    # make a bigger global font size and sans-serif style
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=12)
    plt.rc("legend", fontsize=10)
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))

    r0 = chi0 * strend
    for k, Zmin in enumerate(Sshadow):
        r_theo = Eq5(t, Sshadow[k], chi0, t0, -depthS, strend)
        ax.plot(t, r_cfm[k, :] / r0, c=cb[k], lw=1, alpha=0.3, ls="dashed", label=r"CF")
        ax.plot(
            t,
            r_rsd[k, :],
            c=cb[k],
            lw=1,
            alpha=0.9,
            ls="dotted",
            label=r"RS$_\mathrm{subcrit}$",
        )
        ax.plot(
            t,
            r_tdsr[k, :] / r0,
            c=cb[k],
            lw=3,
            alpha=0.5,
            zorder=2,
            label=r"$\zeta_\mathrm{min}=%.0f$" % (Zmin),
        )
        ax.plot(t, r_theo / r0, c=cb[k], lw=1, alpha=1.0, zorder=2, label="Eq.(5)")

    ax.set_xlim(Tmin, Tmax)
    ax.set_ylim(-0.2, 5.2)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"Rate  $R \, / \, (\chi_0\, \dot\sigma_{c})$")
    ax.legend()

    fig.savefig(str(out) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".png", format="png", dpi=300, bbox_inches="tight")


def test_fig3a():
    tdsr = TDSR1()
    cfm = CFM()
    rsd = RSD1()

    # time unit
    hours = 1.0
    # (=-dsig) skin depth in MPa
    # must be negative, will be changed with next release
    depthS = -5.0
    # mean failure time in sec for critically stressed source
    t0 = 1.0
    # (=dotsigc) in MPa/timeunit: tectonic stressing rate before loading
    strend = 1.0
    # susceptibility to trigger critically stressed sources
    # if a unit step increas is applied (r0 = chi0*strend)
    chi0 = 1.0

    # define Sshadow (Zmin)
    Sshadow = [0.0, 30.0, 60.0]

    tstart = 1.0e-4 * hours
    tend = 100.0 * hours

    # obsolete if logarithmic time axis
    deltat = 0.1 * hours
    # =NT = len(t) - overwritten if logarithmic time axis
    nt = int(np.ceil((tend - tstart) / deltat))
    # linear time axis discretization
    taxis_log = False
    ntlog = 1000
    tstartl = 1.0e-4 * hours

    # discretizing stress axis for integration
    # increment do discretize Coulomb stress axis
    deltaS = -depthS / 60.0
    # maximum depth on Coulomb axis (limit of integral)
    sigma_max = 3000.0 * deltaS

    # uniform distribution
    iX0 = "uniform"

    # calculate earthquake rates for cyclic loading
    ns = len(Sshadow)
    r_tdsr = np.zeros((ns, nt))
    r_cfm = np.zeros((ns, nt))
    r_rsd = np.zeros((ns, nt))

    for k, Zmin in enumerate(Sshadow):
        loading = BackgroundLoading(
            strend=strend,
            taxis_log=taxis_log,
            ntlog=ntlog,
            deltat=deltat,
            tstart=tstartl,
            tend=tend,
        )

        common = dict(
            loading=loading,
            chi0=chi0,
            depthS=depthS,
            deltaS=deltaS,
            Sshadow=Sshadow[k],
            taxis_log=taxis_log,
            ntlog=ntlog,
            deltat=deltat,
            tstart=tstartl,
            tend=tend,
            t0=t0,
            sigma_max=sigma_max,
            iX0=iX0,
        )

        t, chiz, cf, r, xn = tdsr(**common)
        r_tdsr[k, :] = r

        t, chiz, cf, r, xn = cfm(**common)
        r_cfm[k, :] = r

        t, chiz, cf, r, xn = rsd(**common)
        r_rsd[k, :] = r

    # plot values
    test_name = inspect.currentframe().f_code.co_name
    plot_fig3a(
        out=PLOT_DIR / test_name,
        t=t,
        t0=t0,
        Sshadow=Sshadow,
        strend=strend,
        r_tdsr=r_tdsr,
        r_cfm=r_cfm,
        r_rsd=r_rsd,
        chi0=chi0,
        depthS=depthS,
    )

    # dump values
    values = dict(
        r_tdsr=r_tdsr,
        r_cfm=r_cfm,
        r_rsd=r_rsd,
    )
    dump_values(key=test_name, values=values, overwrite=False)

    # assert correct values
    values = load_values(key=test_name)
    assert np.allclose(values["r_tdsr"], r_tdsr)
    assert np.allclose(values["r_cfm"], r_cfm)
    assert np.allclose(values["r_rsd"], r_rsd)
