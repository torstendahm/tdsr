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
    print("set-up the tdsr, linear Coulomb failure and Rate and State models")
    tdsr = TDSR1()
    cfm = CFM()
    rsd = RSD1()

    # ----------------------------------------
    print("define model parameter for plotting, if different from config file")
    # ----------------------------------------
    hours = 1.0  # time unit
    # (=-dsig) skin depth in MPa (must be negativ, will be changed with next release)
    depthS = -5.0
    t0 = 1.0  # mean failure time in sec for critically stressed source
    strend = 1.0  # (=dotsigc) in MPa/timeunit: tectonic stressing rate before loading
    chi0 = 1.0  # susceptibility to trigger critically stressed sources if a unit step increas is applied (r0 = chi0*strend)

    # ---- define Sshadow (Zmin) --------
    Sshadow = [0.0, 30.0, 60.0]

    tstart = 1.0e-4 * hours
    tend = 100.0 * hours

    deltat = 0.1 * hours  # obsolet if logarithmic time axis
    # =NT = len(t) - overwritten if logarithmic time axis
    nt = int(np.ceil((tend - tstart) / deltat))
    taxis_log = 0  # linear time axis discretization = default in config.toml
    ntlog = 1000
    tstartl = 1.0e-4 * hours

    # ---- discretizing stress axis for integration with
    deltaS = -depthS / 60.0  # increment do discretize Coulomb stress axis
    sigma_max = 3000.0 * deltaS  # maximum depth on Coulomb axis (limit of integral)

    iX0switch = 1  # uniform distribution

    # ----------------------------------------
    # Calculate earthquake rates for cyclic loading
    # ----------------------------------------
    ns = len(Sshadow)
    r_tdsr = np.zeros((ns, nt))
    r_cfm = np.zeros((ns, nt))
    r_rsd = np.zeros((ns, nt))

    for k, Zmin in enumerate(Sshadow):
        loading = BackgroundLoading(
            _config=tdsr.config,
            strend=strend,
            taxis_log=taxis_log,
            ntlog=ntlog,
            deltat=deltat,
            tstart=tstartl,
            tend=tend,
        )
        config, t, chiz, cf, r, xn = tdsr(
            loading=loading,
            chi0=chi0,
            t0=t0,
            depthS=depthS,
            deltaS=deltaS,
            sigma_max=sigma_max,
            iX0switch=iX0switch,
            Sshadow=Sshadow[k],
            taxis_log=taxis_log,
            ntlog=ntlog,
            deltat=deltat,
            tstart=tstartl,
            tend=tend,
        )
        r_tdsr[k, :] = r

        loading = BackgroundLoading(
            _config=cfm.config,
            strend=strend,
            taxis_log=taxis_log,
            ntlog=ntlog,
            deltat=deltat,
            tstart=tstartl,
            tend=tend,
        )
        config, t, chiz, cf, r, xn = cfm(
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
        )
        r_cfm[k, :] = r

        loading = BackgroundLoading(
            _config=rsd.config,
            strend=strend,
            taxis_log=taxis_log,
            ntlog=ntlog,
            deltat=deltat,
            tstart=tstartl,
            tend=tend,
        )
        config, t, chiz, cf, r, xn = rsd(
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
        )
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
