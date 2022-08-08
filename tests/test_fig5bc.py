#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test figure 5bc"""

import inspect
import numpy as np
import matplotlib.pyplot as plt

from .utils import dump_values, load_values, PLOT_DIR, DATA_DIR
from tdsr import CFM, RSD1, TDSR1
from tdsr.loading import CustomLoading


def compute_groningen_stress_loading(tstart, tend):
    # earthquake data
    data = np.loadtxt(DATA_DIR / "groningenEQ.dat")
    teq = data[:, 0]
    # pressure data
    data = np.loadtxt(DATA_DIR / "groningenP.dat")
    tgrid = data[:, 0]
    p = data[:, 1]

    # coulomb stress function
    S = -p
    dS = np.ediff1d(S, to_begin=0.0)
    S = np.cumsum(dS)
    return teq, np.transpose([tgrid, S])


def plot_fig5bc(
    out,
    t,
    tend,
    teq,
    R_CF,
    R_CFsub,
    R_RS,
    R_RSsub,
    R_TDSR,
    R_TDSRsub_uniform,
    R_TDSRsub_gauss,
):
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=12)
    plt.rc("legend", fontsize=8)
    fig, ax = plt.subplots(
        2, 1, gridspec_kw={"height_ratios": [1, 2]}, figsize=(7, 6.5)
    )

    ax[0].set_xlim(-2, 30)
    ax[0].set_ylim(0, 180)
    ax[0].set_xlabel(r"$\zeta$    [MPa]", labelpad=1)
    ax[0].set_ylabel(r"$\chi$   [#/MPa]")
    ax[0].legend()

    tbin = np.arange(1960, tend + 0.1, 1.0)
    nobs, tbin1 = np.histogram(teq, bins=tbin)
    ax[1].scatter(
        0.5 * (tbin[:-1] + tbin[1:]),
        nobs,
        c="k",
        marker="o",
        s=60,
        lw=3,
        alpha=0.5,
        label="observed",
    )
    ax[1].plot(t, R_CF, c="k", ls="dashed", alpha=0.3, label=r"CF")
    ax[1].plot(t, R_CFsub, c="k", alpha=0.3, label=r"CF$_\mathrm{subcrit}$")
    ax[1].plot(t, R_RS, c="b", lw=3, alpha=0.3, label=r"RS")
    ax[1].plot(t, R_RSsub, "g--", label=r"RS$_\mathrm{subcrit}$")
    ax[1].plot(t, R_TDSR, "b--", label=r"TDSR$_\mathrm{crit}$")
    ax[1].plot(t, R_TDSRsub_uniform, "r--", label=r"TDSR$_\mathrm{uniform}$")
    ax[1].plot(t, R_TDSRsub_gauss, "r-", label=r"TDSR$_\mathrm{gauss}$")
    ax[1].set_xlim(1980, tend)
    ax[1].set_xlabel("Time   [year]")
    ax[1].set_ylabel("Rate   [#/year]")
    ax[1].legend()
    ax[0].text(
        -0.12,
        0.93,
        "(b)",
        fontsize=14,
        color="k",
        horizontalalignment="left",
        transform=ax[0].transAxes,
    )
    ax[1].text(
        -0.12,
        0.95,
        "(c)",
        fontsize=14,
        color="k",
        horizontalalignment="left",
        transform=ax[1].transAxes,
    )

    fig.savefig(str(out) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".png", format="png", dpi=300, bbox_inches="tight")


def test_fig5bc():
    tdsr = TDSR1()
    cfm = CFM()
    rsd1 = RSD1()

    # Asig = dsig in the RS model
    dsig = 1.0
    depthS = -dsig
    # [MPa] gap for CF and RS_subcrit
    dS0 = 8.5
    # [MPa] gap for TDSR with uniform subcritical stress state
    dS0_TDSR = 22.0
    # [MPa] mean stress for TDSR with gaussian subcritical stress state
    dS0mean = 24.5
    # [MPa] standard deviation " " "
    dS0std = 1.0
    # [years] relaxation time ta=dsig/dotsigc
    ta = 3e5
    # [years] relaxation time in the subcritical RS model
    taRSsub = 60.0
    dotsigc = dsig / ta
    dotsigc_sub = dsig / taRSsub

    # [years] mean failure time of a source at critical stress
    # in the TDRS model
    t0 = 1e-4

    # discretizing stress axis for integration
    # increment do discretize Coulomb stress axis
    deltaS = -depthS / 60.0
    # maximum depth on Coulomb axis (limit of integral)
    sigma_max = 3000.0 * deltaS

    tstart = 1960  # start time for model simulation
    tend = 2022.0  # end time for model simulation
    T1 = 1991  # start time period for fit
    T2 = tend  # end time period for fit

    teq, data = compute_groningen_stress_loading(tstart=tstart, tend=tend)
    tgrid = data[:, 0]
    S = data[:, 1]
    NEQ = len(teq[((teq >= T1) & (teq <= T2))])
    dt = tgrid[1] - tgrid[0]
    nt = int(np.ceil((tend - tstart) / dt))
    t00 = np.min(tgrid[((S >= dS0) & (tgrid <= teq[0]))])
    print(
        "time when EQ start, t00=",
        t00,
        tgrid[0],
        dotsigc * (t00 - tgrid[0]),
        dotsigc_sub * (t00 - tgrid[0]),
    )

    # loading assuming that faults under critical stress exists
    common_loading = dict(
        data=data,
        scal_t=1.0,
        scal_cf=1.0,
        c_tstart=tstart,
        tstart=tstart,
        tend=tend,
        deltat=dt,
    )
    loading = CustomLoading(strend=dotsigc, **common_loading)

    common = dict(
        chi0=1.0,
        t0=t0,
        depthS=depthS,
        deltat=dt,
        tstart=tstart,
        tend=tend,
        taxis_log=False,
        deltaS=deltaS,
        sigma_max=sigma_max,
    )

    # CF model with critical initial stress state:
    t, chiz, cf, R_CF, xn = cfm(
        loading=loading,
        Sshadow=0.0,
        **common,
    )
    # renormalize to observed EQ-number within [T1, T2]:
    fac = NEQ / np.sum(R_CF[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
    R_CF *= fac

    # CF model with subcritical initial stress state:
    t, chiz, cf, R_CFsub, xn = cfm(
        loading=loading,
        Sshadow=dS0,
        **common,
    )
    fac = NEQ / np.sum(R_CFsub[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
    R_CFsub *= fac

    # RS model for critical initial stress state:
    t, chiz, cf, R_RS, xn = rsd1(
        loading=loading,
        Sshadow=0.0,
        **common,
    )
    fac = NEQ / np.sum(R_RS[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
    R_RS *= fac

    # TDSR model for critical initial stress state:
    t, chiz, cf, R_TDSR, xn = tdsr(
        loading=loading,
        Sshadow=0.0,
        **common,
    )
    fac = NEQ / np.sum(R_TDSR[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
    R_TDSR *= fac

    # assuming that no faults exists under critical stress
    # (needed for TDSR and RS)
    sub_loading = CustomLoading(strend=dotsigc_sub, **common_loading)

    # RS model for subcritical initial stress state:
    t, chiz, cf, R_RSsub, xn = rsd1(
        loading=sub_loading,
        Sshadow=(dS0 + dotsigc_sub * (t00 - tgrid[0])),
        **common,
    )
    fac = NEQ / np.sum(R_RSsub[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
    R_RSsub *= fac

    # TDSR model with uniform subcritical stress:
    t, chiz, cf, R_TDSRsub_uniform, xn = tdsr(
        loading=sub_loading,
        Sshadow=dS0_TDSR,
        iX0="uniform",
        **common,
    )
    fac = NEQ / np.sum(
        R_TDSRsub_uniform[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt
    )
    R_TDSRsub_uniform *= fac

    # TDSR model with gaussian subcritical prestress:
    t, chiz, cf, R_TDSRsub_gauss, xn = tdsr(
        loading=sub_loading,
        Sshadow=0.0,
        iX0="gaussian",
        Zmean=dS0mean,
        Zstd=dS0std,
        **common,
    )
    fac = NEQ / np.sum(
        R_TDSRsub_gauss[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt
    )
    R_TDSRsub_gauss *= fac

    # plot values
    test_name = inspect.currentframe().f_code.co_name
    plot_fig5bc(
        out=PLOT_DIR / test_name,
        t=t,
        tend=tend,
        teq=teq,
        R_CF=R_CF,
        R_CFsub=R_CFsub,
        R_RS=R_RS,
        R_RSsub=R_RSsub,
        R_TDSR=R_TDSR,
        R_TDSRsub_uniform=R_TDSRsub_uniform,
        R_TDSRsub_gauss=R_TDSRsub_gauss,
    )

    # dump values
    values = dict(
        R_CF=R_CF,
        R_CFsub=R_CFsub,
        R_RS=R_RS,
        R_RSsub=R_RSsub,
        R_TDSR=R_TDSR,
        R_TDSRsub_uniform=R_TDSRsub_uniform,
        R_TDSRsub_gauss=R_TDSRsub_gauss,
    )
    dump_values(key=test_name, values=values, overwrite=False)

    # assert correct values
    values = load_values(key=test_name)
    assert np.allclose(values["R_CF"], R_CF)
    assert np.allclose(values["R_CFsub"], R_CFsub)
    assert np.allclose(values["R_RS"], R_RS)
    assert np.allclose(values["R_RSsub"], R_RSsub)
    assert np.allclose(values["R_TDSR"], R_TDSR)
    assert np.allclose(values["R_TDSRsub_uniform"], R_TDSRsub_uniform)
    assert np.allclose(values["R_TDSRsub_gauss"], R_TDSRsub_gauss)
