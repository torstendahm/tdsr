#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test figure 4ab"""

import inspect
import numpy as np
import matplotlib.pyplot as plt

from .utils import dump_values, load_values, PLOT_DIR
from tdsr import CFM, RSD, TDSR1
from tdsr.loading import StepLoading
from tdsr.utils import Eq7


def plot_fig4ab(
    out,
    t,
    t0,
    tstart,
    tend,
    sstep,
    strend,
    cfs,
    cf_shad,
    r_tdsr,
    r_lcm,
    r_rsd,
    R_TDSR_theory,
    chi0,
    depthS,
):
    # make a bigger global font size and sans-serif style
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=12)
    plt.rc("legend", fontsize=10)
    fig, ax = plt.subplots(2, 1, figsize=(7, 6))
    plt.subplots_adjust(hspace=0.0)

    al = (0.8, 0.2)
    good = t > 0.0
    for i, Sstep in enumerate(sstep):
        ax[0].plot(t, cf_shad[i, :], c="k", ls="solid", lw=1.5, alpha=al[i])
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
            ax[1].plot(
                t, r_lcm[i, :], c="g", ls="dashed", lw=1, alpha=al[i], label="CF"
            )
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

    fig.savefig(str(out) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".png", format="png", dpi=300, bbox_inches="tight")


def test_fig4ab():
    tdsr = TDSR1()
    cfm = CFM()
    rsd = RSD()

    # time unit
    hours = 1.0
    # (=-dsig) skin depth in MPa
    # must be negative, will be changed with next release
    depthS = -1.0
    # mean failure time in sec for critically stressed source t0 = 1.0
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
    nt = int(np.ceil((tend - tstart) / deltat))

    # linear time axis
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

    # calculate equilibrium distribution of tectonic loading,
    # to be used later to avoid initial oscillations
    # loading = BackgroundLoading(strend=strend)
    # t, chiz_background, cf, r, xn = tdsm(
    #   loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS,
    #   sigma_max=sigma_max, deltat=deltat, tstart=0, tend=tend)

    for i, Sstep in enumerate(sstep):
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

        # BUG: cf_shadow is written here, but cf_shad is plotted (zeros only)!
        t, cf_shadow, cf, r, xn = rsd(**common)
        r_rsd[i, :] = r

        t, chiz, cf, r, xn = cfm(**common)
        r_lcm[i, :] = r

        R_TDSR_theory[i, :] = Eq7(t, sstep[i], chi0 * strend, -depthS, strend, strend)

    test_name = inspect.currentframe().f_code.co_name

    # plot values
    plot_fig4ab(
        out=PLOT_DIR / test_name,
        t=t,
        t0=t0,
        tstart=tstart,
        tend=tend,
        sstep=sstep,
        strend=strend,
        cfs=cfs,
        cf_shad=cf_shad,
        r_tdsr=r_tdsr,
        r_lcm=r_lcm,
        r_rsd=r_rsd,
        R_TDSR_theory=R_TDSR_theory,
        chi0=chi0,
        depthS=depthS,
    )

    # dump values
    values = dict(
        cfs=cfs,
        cf_shad=cf_shad,
        r_tdsr=r_tdsr,
        r_lcm=r_lcm,
        r_rsd=r_rsd,
        R_TDSR_theory=R_TDSR_theory,
    )
    dump_values(key=test_name, values=values, overwrite=False)

    # assert correct values
    values = load_values(key=test_name)
    assert np.allclose(values["cfs"], cfs)
    assert np.allclose(values["cf_shad"], cf_shad)
    assert np.allclose(values["r_tdsr"], r_tdsr)
    assert np.allclose(values["r_lcm"], r_lcm)
    assert np.allclose(values["r_rsd"], r_rsd)
    assert np.allclose(values["R_TDSR_theory"], R_TDSR_theory)
