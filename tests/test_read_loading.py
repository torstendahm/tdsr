#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test read loading"""

import inspect
import numpy as np
import matplotlib.pyplot as plt

from .utils import dump_values, load_values, PLOT_DIR, DATA_DIR, TEST_DIR
from tdsr import Config, RSD1, TDSR1
from tdsr.loading import ExternalFileLoading


def plot_read_loading(out, t, tstart, tend, S, r_tdsr, r_rsd):
    # make a bigger global font size and sans-serif style
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

    fig.savefig(str(out) + ".pdf", dpi=300, format="pdf", bbox_inches="tight")
    fig.savefig(str(out) + ".png", dpi=300, format="png", bbox_inches="tight")


def test_read_loading():
    config_file = TEST_DIR / "config.toml"
    data_file = DATA_DIR / "CFSloading_synthetic.dat"

    print("set-up the tdsm, lcm and rsm class environments")
    tdsr = TDSR1(config=Config.open(config_file))
    rsd1 = RSD1(config=Config.open(config_file))

    print("define model parameter for plotting, if different from config file")
    hours = 1.0
    tstart = 0.0 * hours
    tend = 10 * hours
    deltat = 0.001 * hours
    t0 = 1.0 * deltat
    nt = np.floor((tend - tstart) / deltat).astype(int)

    strend = 1.0
    chi0 = 1.0  # susceptibility to trigger earthquakes by unit stress increase
    depthS = -1.0  # skin depth in MPa (must be defined negativ)

    deltaS = -depthS / 500.0  # increment do discretize Coulomb stress axis
    sigma_max = 10000.0 * deltaS  # maximum depth on Coulomb axis (limit of integral)

    iX0switch = 0  # steady state distribution

    scal_cf = 1.0  # data provided in Pa, to be scaled to  MPa
    scal_t = 1.0  # time provilded in units of days, to be scaled to seconds
    c_tstart = 0.0

    # -----------------------------------------------------
    print("Calculate stress loading function and save in directory ", data_file)
    # -----------------------------------------------------
    dsig = -depthS
    t = np.linspace(tstart, tend, nt)
    dS = np.zeros(nt)

    N1 = int(0.1 * nt)
    N2 = int(0.3 * nt)
    N3 = int(0.5 * nt)
    N4 = int(0.6 * nt)

    dS[:N1] = strend * deltat
    dS[N2:N3] = 3 * strend * deltat
    dS[N4:] = 0.5 * strend * deltat
    dS[N1] += 5 * dsig
    dS[N2] -= 5 * dsig
    dS[N4] += 2 * dsig
    S = np.cumsum(dS)

    # External loading file is expected to have two headin lines (see tdsm/loading.py mthod ExternalFileLoading
    # first column is Coulomb stress , secind column is time (note that scaling factors can be provided during call)
    np.savetxt(
        data_file,
        np.transpose([S, t]),
        header="\n".join(
            [
                "Ascii loading file is expected to have two header lines",
                "1st column: Coulomb stress, 2nd column: time",
            ]
        ),
    )

    # -----------------------------------------------------
    print("Calculate earthquake rates with tdsr and rsm")
    # -----------------------------------------------------
    r_tdsr = np.zeros(nt)
    r_rsd = np.zeros(nt)

    loading = ExternalFileLoading(
        _config=tdsr.config,
        iascii=True,
        infile=data_file,
        scal_t=scal_t,
        scal_cf=scal_cf,
        strend=strend,
        c_tstart=c_tstart,
        tstart=tstart,
        tend=tend,
        deltat=deltat,
    )
    config, t, chiz, cf, r_tdsr, xn = tdsr(
        loading=loading,
        chi0=chi0,
        t0=t0,
        depthS=depthS,
        deltaS=deltaS,
        sigma_max=sigma_max,
        iX0switch=iX0switch,
        deltat=deltat,
        taxis_log=0,
        tstart=tstart,
        tend=tend,
    )

    loading = ExternalFileLoading(
        _config=rsd1.config,
        iascii=True,
        infile=data_file,
        scal_t=scal_t,
        scal_cf=scal_cf,
        strend=strend,
        c_tstart=c_tstart,
        tstart=tstart,
        tend=tend,
        deltat=deltat,
    )
    config, t, chiz, cf, r_rsd, xn = rsd1(
        loading=loading,
        chi0=chi0,
        depthS=depthS,
        deltat=deltat,
        tstart=tstart,
        tend=tend,
    )

    # plot values
    test_name = inspect.currentframe().f_code.co_name
    plot_read_loading(
        out=PLOT_DIR / test_name,
        t=t,
        tstart=tstart,
        tend=tend,
        S=S,
        r_tdsr=r_tdsr,
        r_rsd=r_rsd,
    )

    # dump values
    values = dict(
        r_tdsr=r_tdsr,
        r_rsd=r_rsd,
    )
    dump_values(key=test_name, values=values, overwrite=False)

    # assert correct values
    values = load_values(key=test_name)
    assert np.allclose(values["r_tdsr"], r_tdsr)
    assert np.allclose(values["r_rsd"], r_rsd)
