#!/usr/bin/env python
# -*- coding: utf-8 -*-

# --------------------------
# Example Fig 6def in Dahm and Hainzl (2022): KTB injection induced seismicity
# --------------------------

import sys
import matplotlib.pyplot as plt
from scipy import special
import numpy as np
from pathlib import Path

EXAMPLE_DIR = Path(__file__).parent
REPO_ROOT = EXAMPLE_DIR.parent.absolute()
sys.path.insert(0, str(REPO_ROOT))

from tdsr import TDSR1  # noqa: E402
from tdsr.loading import CustomLoading  # noqa: E402, F401

data_file = REPO_ROOT / "data/CFSloading_KTB.dat"
figname = REPO_ROOT / "plots/fig6def"

tdsr = TDSR1()

day2sec = 24 * 60 * 60.0
km = 1000.0
hours = 3600.0
days = hours * 24.0

t0 = 0.01  # [days]
dottauPayr = 30.0
strend = 1e-6 * dottauPayr / 365.25  # [MPa/day]
# susceptibility to trigger earthquakes by unit stress increase (backgrond rate set to r0=1)
chi0 = 1.0/strend

# values for -depthS (here interpreted as normalised by friction coefficient)
dsigvalues = (0.9, 1.0, 1.1)

# increment do discretize Coulomb stress axis
deltaS = 0.3 / 500.0
# maximum depth on Coulomb axis (limit of integral)
sigma_max = 10000.0 * deltaS
# precision = 18

D = 0.033  # [m^2/s]
Mcut = -2.3

def read_flowdata():
    data = np.loadtxt(REPO_ROOT / "data/ktb_flowrates.dat", skiprows=1, usecols=(0, 1))
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
    data = np.loadtxt(REPO_ROOT / "data/ktb_HBpressure.dat", skiprows=2, delimiter=",")
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

    # Segall & Lu 2015, below Eq.(2)
    lambda_d = 2.0 * mu * vd / (1.0 - 2.0 * vd)
    lambda_u = 2.0 * mu * vu / (1.0 - 2.0 * vu)

    # Segall & Lu 2015, Eq.(4): k/eta = D * ketafac
    ketafac = np.square(alpha) * (lambda_u + 2.0 * mu)
    ketafac /= (lambda_u - lambda_d) * (lambda_d + 2.0 * mu)
    kappa = D * ketafac
    fac = 1.0 / (4 * np.pi * kappa)
    P = np.zeros(len(tgrid))
    for i in range(len(tgrid)):
        tzero = tgrid[i]
        dqi = dq[(tq < tzero)]
        dti = (tzero - tq[(tq < tzero)]) * day2sec
        psi = rdist / np.sqrt(D * dti)
        P[i] = (fac / (D * rdist)) * np.sum(dqi * special.erfc(0.5 * psi))  # [Pa]
    return P / 1e6


# read observed data again to get some paramter for plotting and scaling
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
data = np.loadtxt(REPO_ROOT / "data/ktb_EQ.dat")
teqall = data[:, 0]
meqall = data[:, 1]
ind = meqall >= Mcut
teq = teqall[ind]
meq = meqall[ind]

T1 = np.minimum(np.min(tflow), np.min(teq))
T2 = np.max(teqall)
tstart = T1    # in days
c_tstart = T1
tend = T2      # in days
deltat = 0.1   # in days
nt = np.floor((tend - tstart) / deltat).astype(int)
nt=nt+1
tgrid = np.linspace(tstart, tend, nt)
dt = np.ediff1d(tgrid, to_end=tgrid[-1] - tgrid[-2])

Pcalibration = calculateP(D, tgrid, r_calibration, tq, dq)  # for plotting pressure model

Ri = np.linspace(r10, r90, 10)
print("\t P calculation in distance range: [%.0f - %.0f]m" % (Ri[0], Ri[-1]))
dSr = np.zeros((len(Ri), len(tgrid)))
for i, R in enumerate(Ri):
    P = calculateP(D, tgrid, R, tq, dq)
    dSr[i, :] = np.ediff1d(P, to_begin=0.0)   # pressure change, as later a tectonic loading is added

common_loading = dict(
    scal_t=1.0,
    scal_cf=1.0,
    c_tstart=c_tstart,
    tstart=tstart,
    tend=tend,
    deltat=deltat,
    strend = strend,
)

common = dict(
    chi0=chi0,
    Sshadow=0.0,
    t0=t0,
    iX0='equilibrium',
    deltat=deltat,
    tstart=tstart,
    tend=tend,
    taxis_log=False,
    deltaS=deltaS,
    sigma_max=sigma_max,
)

# ----- calculate earthquake rates with tdsm, lcm and rsm
ns = len(dsigvalues)
cfs = np.zeros(nt)
r_tdsr = np.zeros((ns, nt))

for i, dsig in enumerate(dsigvalues):
    depthS = -dsig
    R = np.zeros(len(tgrid))
    for k in range(len(Ri)):   # calculate mean rate over different distances to injection point
        dS = dSr[k,:] + strend * dt  
        S = np.cumsum(dS)  # stress model (normalized by friction) is sum of pore pressure and tectonic loading 
        loading = CustomLoading(data=np.transpose([tgrid, S]), **common_loading)
        t, chiz, cf, r, xn = tdsr(loading=loading, depthS=depthS, **common)
        r /= np.sum(r*dt)
        r_tdsr[i,:] += r * len(teq)/len(Ri)

# -------- plot results
cb = ["b", "r", "g"]
plt.rc("font", family="sans-serif")
plt.rc("font", size=12)
plt.rc("legend", fontsize=12)
fig, ax = plt.subplots(3, 1, gridspec_kw={"height_ratios": [1, 1, 2]}, figsize=(7, 8))
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
    ax[2].plot(tgrid, r_tdsr[i,:]*DT, c=cb[i], lw=2, zorder=10, label=r'$\delta\sigma/\mu=%.1f$ MPa' % (dsig))
ax[2].set_xlim(T1, T2)
ax[2].set_ylim( 0,)
#ax[2].set_ylim( 0, 100)
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
fig.savefig(str(figname) + ".pdf", dpi=300, format="pdf", bbox_inches="tight")
fig.savefig(str(figname) + ".png", dpi=300, format="png", bbox_inches="tight")
print(f"\n\t OUTPUT: {str(figname) + '.pdf'}\n")
