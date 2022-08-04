#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test figure 4cd"""

import inspect
import numpy as np
import matplotlib.pyplot as plt

from .utils import dump_values, load_values, PLOT_DIR
from tdsr import CFM, RSD, TDSR1
from tdsr.loading import CyclicLoading


def plot_fig4cd(out, t, Tmax, cf, cf_shad, r_tdsr, r_rsd, r_cfm):
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

    fig.savefig(str(out) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".png", format="png", dpi=300, bbox_inches="tight")


def test_fig4cd():
    print("set-up the tdsr, linear Coulomb failure and Rate and State models")
    tdsr = TDSR1()
    cfm = CFM()
    rsd = RSD()

    # ----------------------------------------
    print("define model parameter for plotting, if different from config file")
    # ----------------------------------------
    hours = 1.0  # time unit
    depthS = (
        -1.0
    )  # (=-dsig) skin depth in MPa (must be negativ, will be changed with next release)
    t0 = 1.0  # mean failure time in sec for critically stressed source
    strend = 1.0  # (=dotsigc) in MPa/timeunit: tectonic stressing rate before loading
    chi0 = 1.0  # susceptibility to trigger critically stressed sources if a unit step increas is applied (r0 = chi0*strend)

    ampsin = 5.0e0  # Amplitude of sinusoidal loading
    Tsin = ampsin / strend
    Tmax = 3.0 * Tsin

    tstart = 0.0
    tend = Tmax

    deltat = 0.02 * hours
    t0 = deltat
    nt = int(np.ceil((tend - tstart) / deltat))  # =NT = len(t)

    # ---- discretizing stress axis for integration with
    deltaS = -depthS / 60.0  # increment do discretize Coulomb stress axis
    sigma_max = 3000.0 * deltaS  # maximum depth on Coulomb axis (limit of integral)

    # ----------------------------------------
    # Calculate earthquake rates for cyclic loading
    # ----------------------------------------
    loading = CyclicLoading(
        _config=tdsr.config,
        strend=strend,
        ampsin=ampsin,
        Tsin=Tsin,
        deltat=deltat,
        tstart=tstart,
        tend=tend,
    )
    config, t, chiz, cf, r, xn = tdsr(
        loading=loading,
        chi0=chi0,
        t0=t0,
        depthS=depthS,
        deltaS=deltaS,
        sigma_max=sigma_max,
        iX0switch=0,
        deltat=deltat,
        taxis_log=0,
        tstart=0,
        tend=tend,
    )
    cfs = cf
    r_tdsr = r

    loading = CyclicLoading(
        _config=cfm.config,
        strend=strend,
        ampsin=ampsin,
        Tsin=Tsin,
        deltat=deltat,
        tstart=0,
        tend=tend,
    )
    config, t, chiz, cf, r, xn = cfm(
        loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend
    )
    cf_shad = chiz[0:nt]
    r_cfm = r

    loading = CyclicLoading(
        _config=rsd.config,
        strend=strend,
        ampsin=ampsin,
        Tsin=Tsin,
        deltat=deltat,
        tstart=tstart,
        tend=tend,
    )
    config, t, chiz, cf, r, xn = rsd(
        loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend
    )
    r_rsd = r

    test_name = inspect.currentframe().f_code.co_name

    # plot values
    plot_fig4cd(
        out=PLOT_DIR / test_name,
        t=t,
        Tmax=Tmax,
        cf=cf,
        cf_shad=cf_shad,
        r_tdsr=r_tdsr,
        r_rsd=r_rsd,
        r_cfm=r_cfm,
    )

    # dump values
    values = dict(
        cfs=cfs,
        cf_shad=cf_shad,
        r_tdsr=r_tdsr,
        r_rsd=r_rsd,
        r_cfm=r_cfm,
    )
    dump_values(key=test_name, values=values, overwrite=False)

    # assert correct values
    values = load_values(key=test_name)
    assert np.allclose(values["cfs"], cfs)
    assert np.allclose(values["cf_shad"], cf_shad)
    assert np.allclose(values["r_tdsr"], r_tdsr)
    assert np.allclose(values["r_rsd"], r_rsd)
    assert np.allclose(values["r_cfm"], r_cfm)
