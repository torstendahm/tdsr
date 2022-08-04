#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test figure 3b"""

import inspect
import numpy as np
import matplotlib.pyplot as plt

from .utils import dump_values, load_values, PLOT_DIR
from tdsr import Config, TDSR1
from tdsr.loading import BackgroundLoading
from tdsr.utils import Eq7


def plot_fig3b(out, t, sstep, strend, r_tdsr, chi0, depthS):
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

    fig.savefig(str(out) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".png", format="png", dpi=300, bbox_inches="tight")


def test_fig3b() -> None:
    print("set-up the tdsr, linear Coulomb failure and Rate and State models")
    tdsr = TDSR1()

    # ----------------------------------------
    print("define model parameter for plotting, if different from config file")
    # ----------------------------------------
    hours = 1.0  # time unit
    depthS = (
        -1.0
    )  # (=-dsig) skin depth in MPa (must be negative, will be changed with next release)
    t0 = 1.0  # mean failure time in sec for critically stressed source
    strend = 1.0  # (=dotsigc) in MPa/timeunit: tectonic stressing rate before loading
    chi0 = 1.0  # susceptibility to trigger critically stressed sources if a unit step increas is applied (r0 = chi0*strend)

    # ---- define stress steps --------
    sstep = np.asarray([2.0, 4.0, 6.0, 8.0, 10.0])  # stress steps in e.g. MPa

    tend = 100.0 * hours

    taxis_log = 1  # linear time axis discretization = default in config.toml
    ntlog = 10000
    tstartl = 1.0e-6 * hours

    # ---- discretizing stress axis for integration with
    deltaS = -depthS / 60.0  # increment do discretize Coulomb stress axis
    sigma_max = 3000.0 * deltaS  # maximum depth on Coulomb axis (limit of integral)
    # precision = 18   # wird bisher nicht beneotigt in tdsr

    iX0switch = 0  # steady state distribution

    # ----------------------------------------
    # Calculate earthquake rates for cyclic loading
    # ----------------------------------------
    ns = len(sstep)
    r_tdsr = np.zeros((ns, ntlog))

    for k in range(ns):
        loading = BackgroundLoading(
            _config=tdsr.config,
            strend=strend,
            sstep=sstep[k],
            taxis_log=taxis_log,
            ntlog=ntlog,
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
            Sshadow=sstep[k],
            taxis_log=taxis_log,
            ntlog=ntlog,
            tstart=tstartl,
            tend=tend,
        )
        r_tdsr[k, :] = r

    # plot values
    test_name = inspect.currentframe().f_code.co_name
    plot_fig3b(
        out=PLOT_DIR / test_name,
        t=t,
        sstep=sstep,
        strend=strend,
        r_tdsr=r_tdsr,
        chi0=chi0,
        depthS=depthS,
    )

    # dump values
    values = dict(
        r_tdsr=r_tdsr,
    )
    dump_values(key=test_name, values=values, overwrite=False)

    # assert correct values
    values = load_values(key=test_name)
    assert np.allclose(values["r_tdsr"], r_tdsr)
