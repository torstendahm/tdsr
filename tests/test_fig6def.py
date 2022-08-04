#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test figure 6 def"""

from pathlib import Path
from .utils import DATA_DIR, TEST_DIR, PLOT_DIR
import matplotlib.pyplot as plt
import matplotlib.ticker
import math
import numpy as np
from scipy import special

from tdsr import Config, TDSR, TDSR1, LCM, Traditional, CFM, RSM, RSD
from tdsr.plotting import plot
from tdsr.loading import StepLoading, BackgroundLoading, TrendchangeLoading, RampLoading


day2sec = 24 * 60 * 60.0


def read_flowdata():
    data = np.loadtxt(DATA_DIR / "ktb_flowrates.dat", skiprows=1, usecols=(0, 1))
    tq = data[:, 0]
    # calculate changes of flow rate:
    dq = np.zeros(len(tq))
    dq[0] = data[0, 1]
    dq[1:] = data[1:, 1] - data[:-1, 1]
    tflow = data[0, 0]
    qflow = data[0, 1]
    for i in range(1, len(tq)):
        tflow = np.append(tflow, tq[i])
        qflow = np.append(qflow, np.sum(dq[:i]))
        tflow = np.append(tflow, tq[i])
        qflow = np.append(qflow, np.sum(dq[: i + 1]))
    # transform units from [l/min] to [m^3/s]:
    dq /= 1000.0 * 60.0
    return tflow, qflow, tq, dq


def read_pressuredata():
    data = np.loadtxt(DATA_DIR / "ktb_HBpressure.dat", skiprows=2, delimiter=",")
    tp = data[:, 0]
    p = data[:, 1]  # [MPa]
    return tp, p


def calculateP(Dm2s, tgrid, rdist, tq, dq):
    """
    Calculation of pore-pressure according to Rudnicki (1986); Eq.25
    Input:
    D:       [m2/s] hydraulic diffusivity
    tq:      [days] time of changes in flow rate
    dq:      [m^3/s] flow rate change
    tgrid:   [days] time grid values
    rdist:   [m] distance
    Return:
    P:       [MPa] pore pressure
    """
    ########################################
    vd = 0.25  # drained Poisson's ratio
    vu = 0.30  # undrained Poisson's ratio
    mu = 1e9  # [Pa] shear modulus
    alpha = 0.1  # Biot coefficient
    lambda_d = 2.0 * mu * vd / (1.0 - 2.0 * vd)  # Segall & Lu 2015, below Eq.(2)
    lambda_u = 2.0 * mu * vu / (1.0 - 2.0 * vu)
    ketafac = np.square(alpha) * (
        lambda_u + 2.0 * mu
    )  # Segall & Lu 2015, Eq.(4): k/eta = D * ketafac
    ketafac /= (lambda_u - lambda_d) * (lambda_d + 2.0 * mu)
    kappa = Dm2s * ketafac
    fac = 1.0 / (4 * np.pi * kappa)
    P = np.zeros(len(tgrid))
    for i in range(len(tgrid)):
        t0 = tgrid[i]
        dqi = dq[(tq < t0)]
        dti = (t0 - tq[(tq < t0)]) * day2sec
        psi = rdist / np.sqrt(Dm2s * dti)
        P[i] = (fac / (Dm2s * rdist)) * np.sum(dqi * special.erfc(0.5 * psi))  # [Pa]
    return P / 1e6


def plot_fig6def(
    tpobs, Pobs, tgrid, Pcalibration, T1, T2, tflow, qflow, teqall, meqall
):
    # -------------------------------------------------
    # Plot results
    # -------------------------------------------------
    # General Plot Properties
    # make a bigger global font size and sans-serif style
    cb = ["b", "r", "g"]
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=12)
    plt.rc("legend", fontsize=12)
    fig, ax = plt.subplots(
        3, 1, gridspec_kw={"height_ratios": [1, 1, 2]}, figsize=(7, 8)
    )
    plt.subplots_adjust(hspace=0, wspace=0.2)

    f = 2.5
    ax[0].plot(tpobs, Pobs, c="r", ls="solid", lw=2, label="observed")
    ax[0].plot(tgrid, Pcalibration, c="r", ls="dashed", lw=2, label="modeled")
    ax[0].set_xlim(T1, T2)
    ax[0].set_ylabel("Pressure [MPa]")
    y1 = 0.8
    ax[0].set_ylim(-y1, f * y1)
    ax[0].set_xticks([])
    ax[0].legend()
    a = ax[0].twinx()
    a.fill_between(tflow, qflow, 0, alpha=0.3)
    a.axhline(0, c="k", ls="dotted", lw=1)
    y2 = 100
    a.set_ylim(-y2, f * y2)
    a.set_ylabel("Flow rate [l/min]")
    # a.set_ylim(0, 6.5)

    ax[1].scatter(teqall, meqall, c="k", s=5, alpha=0.2)
    ax[1].axhline(Mcut, c="k", ls="dotted", lw=1)
    ax[1].set_ylabel("Magnitude")
    ax[1].set_xlim(T1, T2)
    ax[1].set_ylim(
        -3.5,
    )
    ax[1].set_xticks([])

    DT = 5.0
    tbins = np.arange(0, T2 + 1, DT)
    ax[2].hist(teq, bins=tbins, color="k", alpha=0.3, label="observed")
    for i, dsig in enumerate(dsigvalues):
        R = np.zeros(len(tgrid))
        for k in range(len(Ri)):
            dS = dSr[k, :] + strend * dt
            S = np.cumsum(dS)
            # Z = functions.Zvalues(S, 0, t0, dsig)
            # Xs = functions.X0steady(Z, 1.0, t0, dsig, strend)
            # RR = functions.TDSR(tgrid, S, Z, Xs, t0, dsig)
            # RR /= np.sum(RR*dt)
            # R += RR * len(teq)/len(Ri)
        # ax[2].plot(tgrid, R*DT, c=cb[i], lw=2, zorder=10, label=r'$\delta\sigma/\mu=%.1f$ MPa' % (dsig))
    ax[2].set_xlim(T1, T2)
    ax[2].set_ylim(
        0,
    )
    ax[2].set_xlabel("Days from 1/1/2002")
    ax[2].set_ylabel("Rate   [# / %.0f days]" % (DT))
    ax[2].legend()
    ax[0].text(
        -0.14,
        0.93,
        "(d)",
        fontsize=14,
        color="k",
        horizontalalignment="left",
        transform=ax[0].transAxes,
    )
    ax[1].text(
        -0.14,
        0.95,
        "(e)",
        fontsize=14,
        color="k",
        horizontalalignment="left",
        transform=ax[1].transAxes,
    )
    ax[2].text(
        -0.14,
        0.95,
        "(f)",
        fontsize=14,
        color="k",
        horizontalalignment="left",
        transform=ax[2].transAxes,
    )

    plt.show()
    figname.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(figname) + ".pdf", dpi=300, format="pdf", bbox_inches="tight")
    fig.savefig(str(figname) + ".png", dpi=300, format="png", bbox_inches="tight")
    print("\n\t OUTPUT: %s\n" % (figname))


def skip_fig6def():
    """Test fig 6def"""

    config_file = TEST_DIR / "config.toml"
    data_file = DATA_DIR / "CFSloading_KTB.dat"
    figname = PLOT_DIR / "fig6def.png"

    print("set-up the tdsm, lcm and rsm class environments")
    tdsr = TDSR1(config=Config.open(config_file))

    print("define model parameter for plotting, if different from config file")
    km = 1000.0
    hours = 3600.0
    days = hours * 24.0

    tstart = 0.0 * days
    tend = 100.0 * days
    deltat = 0.01 * days
    nt = np.floor((tend - tstart) / deltat).astype(int)
    tgrid = np.linspace(tstart, tend, nt)
    dt = np.ediff1d(tgrid, to_end=tgrid[-1] - tgrid[-2])

    t0 = 0.01  # [days]
    dottauPayr = 30.0
    strend = 1e-6 * dottauPayr / 365.25  # [MPa/day]
    chi0 = 1.0  # susceptibility to trigger earthquakes by unit stress increase

    dsigvalues = (0.9, 1.0, 1.1)  # values for -depthS

    deltaS = 0.3 / 500.0  # increment do discretize Coulomb stress axis
    sigma_max = 10000.0 * deltaS  # maximum depth on Coulomb axis (limit of integral)
    # precision = 18

    iX0switch = 1  # uniform distribution (and strend=0)

    D = 0.033  # [m^2/s]
    Mcut = -2.3

    iascii = True
    scal_cf = 1.0e-6  # data provided in Pa, to be scaled to  MPa
    scal_t = 1.0  # time provilded in units of days, to be scaled to days
    c_tstart = 0.0

    # ------------------------------------
    # read observed data again to get some paramter for plotting and scaling
    # ------------------------------------
    # Graessle et al 2006: Distance between open hole section of PH and MH:
    # open sections: depth: MH: 5.2-5.6 km  PH: 3.85-4.0 and horizontally: 0.2 km
    r_calibration = np.sqrt(np.square(200.0) + np.square(5400 - 3925.0))
    print("\n\t Calibration-Distance= %f km" % (r_calibration / 1000))

    # mean distance to Pilothole [m]:
    rmean = 475.55  # [m] mean EQ-distance
    rstd = 130.38  # [m] standard deviation of EQ-distances
    r10 = 331.02  # [m] 10% quantile of EQ-distances
    r90 = 607.52  # [m] 90% quantile of EQ-distances

    # INPUT-DATA:
    tflow, qflow, tq, dq = read_flowdata()
    tpobs, Pobs = read_pressuredata()

    # Earthquake data:
    data = np.loadtxt(DATA_DIR / "ktb_EQ.dat")
    teqall = data[:, 0]
    meqall = data[:, 1]
    ind = meqall >= Mcut
    teq = teqall[ind]
    meq = meqall[ind]

    T1 = np.minimum(np.min(tflow), np.min(teq))
    T2 = np.max(teqall)
    tgrid = np.linspace(T1, T2, 10000)
    dt = np.ediff1d(tgrid, to_end=tgrid[-1] - tgrid[-2])

    Pcalibration = calculateP(D, tgrid, r_calibration, tq, dq)

    Ri = np.linspace(r10, r90, 10)
    print("\t P calculation in distance range: [%.0f - %.0f]m" % (Ri[0], Ri[-1]))
    dSr = np.zeros((len(Ri), len(tgrid)))
    for i, R in enumerate(Ri):
        P = calculateP(D, tgrid, R, tq, dq)
        dSr[i, :] = np.ediff1d(P, to_begin=0.0)

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

    plot_fig6def(
        tpobs=tpobs,
        Pobs=Pobs,
        tgrid=tgrid,
        Pcalibration=Pcalibration,
        T1=T1,
        T2=T2,
        tflow=tflow,
        qflow=qflow,
        teqall=teqall,
        meqall=meqall,
        Mcut=meqall,
    )

    # todo: assert here
