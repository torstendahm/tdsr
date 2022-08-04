#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test figure 7cd"""

import inspect
import numpy as np
import matplotlib.pyplot as plt

from .utils import dump_values, load_values, PLOT_DIR, DATA_DIR
from tdsr import Traditional, TDSR1
from tdsr.loading import ExternalFileLoading


def plot_fig7cd(out, t, tstart, tend, teq, cfs, cf_shad, Req, r_tdsr):
    # make a bigger global font size and sans-serif style
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
            dS0 = 3.0
        if i == 1:
            dsig = 0.4
            dS0 = 3.5
        if i == 2:
            dsig = 0.5
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

    fig.savefig(str(out) + ".pdf", dpi=300, format="pdf", bbox_inches="tight")
    fig.savefig(str(out) + ".png", dpi=300, format="png", bbox_inches="tight")


def test_fig7cd():
    data_cfs = DATA_DIR / "stresschange_morsleben.dat"
    fn_obsrate = DATA_DIR / "morsleben_EQ.dat"

    print("set-up the tdsm, lcm and rsm class environments")
    tdsr = TDSR1()
    trad = Traditional()

    print("define model parameter for plotting, if different from config file")
    days = 1.0

    tstart = 0.0 * days
    tend = 100.0 * days
    deltat = 0.01 * days
    nt = np.floor((tend - tstart) / deltat).astype(int)
    tgrid = np.linspace(tstart, tend, nt)
    dt = np.ediff1d(tgrid, to_end=tgrid[-1] - tgrid[-2])

    t0 = 0.01
    strend = 0.0  # no tectonic loading
    chi0 = 1.0  # susceptibility to trigger earthquakes by unit stress increase

    deltaS = 0.3 / 500.0  # increment do discretize Coulomb stress axis
    sigma_max = 10000.0 * deltaS  # maximum depth on Coulomb axis (limit of integral)

    iX0switch = 1  # uniform distribution (and strend=0)

    scal_cf = 1.0e-6  # data provided in Pa, to be scaled to  MPa
    scal_t = 1.0  # time provilded in units of days, to be scaled to days
    c_tstart = 0.0

    # ------------------------------------
    # read observed data again to get some paramter for plotting and scaling
    # ------------------------------------
    data = np.loadtxt(fn_obsrate)
    data = data[((data[:, 1] >= tstart) & (data[:, 1] <= tend))]
    teq = data[:, 1]
    Neq = data[:, 0]
    NEQ = np.sum(Neq)
    dteq = np.ediff1d(teq, to_end=teq[-1] - teq[-2])
    Req = Neq / dteq

    # -----------------------------------------------------
    print("Calculate earthquake rates with tdsm, lcm and rsm")
    # -----------------------------------------------------
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

        loading = ExternalFileLoading(
            _config=tdsr.config,
            iascii=True,
            infile=data_cfs,
            scal_t=scal_t,
            scal_cf=scal_cf,
            strend=strend,
            c_tstart=c_tstart,
            tstart=tstart,
            tend=tend,
            deltat=deltat,
        )
        config, t, chiz, cfs, r, xn = tdsr(
            loading=loading,
            chi0=chi0,
            t0=t0,
            depthS=depthS,
            Sshadow=Sshadow,
            deltaS=deltaS,
            sigma_max=sigma_max,
            iX0switch=iX0switch,
            deltat=deltat,
            taxis_log=0,
            tstart=tstart,
            tend=tend,
        )
        X0 = NEQ / np.sum(r * dt)
        r *= X0  # X0 set to match the observed number
        r_tdsr[i, :] = r[:]

    # calculate stress shadow ad hoc  by using the linear coulomb failure class
    config, t, cf_shad, cf, r, xn = trad(
        loading=loading, chi0=chi0, deltat=deltat, tstart=tstart, tend=tend
    )

    # plot values
    test_name = inspect.currentframe().f_code.co_name
    plot_fig7cd(
        out=PLOT_DIR / test_name,
        t=t,
        tstart=tstart,
        tend=tend,
        teq=teq,
        cfs=cfs,
        cf_shad=cf_shad,
        Req=Req,
        r_tdsr=r_tdsr,
    )

    # dump values
    values = dict(
        cfs=cfs,
        cf_shad=cf_shad,
        Req=Req,
        r_tdsr=r_tdsr,
    )
    dump_values(key=test_name, values=values, overwrite=False)

    # assert correct values
    values = load_values(key=test_name)
    assert np.allclose(values["cfs"], cfs)
    assert np.allclose(values["cf_shad"], cf_shad)
    assert np.allclose(values["Req"], Req)
    assert np.allclose(values["r_tdsr"], r_tdsr)
