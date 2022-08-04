# --------------------------
# Plot Fig. 5bc (see Dahm and Hainzl, 2022)
# --------------------------
import sys
from pathlib import Path
import numpy as np
import math
from scipy.stats import norm
import matplotlib.pyplot as plt

EXAMPLE_DIR = Path(__file__).parent
REPO_ROOT = EXAMPLE_DIR.parent.absolute()
sys.path.insert(0, str(REPO_ROOT))

from tdsr import Config, TDSR, TDSR1, LCM, CFM, Traditional, RSM, RSD, RSD1
from tdsr.plotting import plot
from tdsr.loading import StepLoading, BackgroundLoading, ExternalFileLoading

config_file = EXAMPLE_DIR / "config.toml"
figname = REPO_ROOT / "plots/fig5bc"

print("set-up the tdsr, linear Coulomb failure and Rate and State models")
# tdsm = TDSM(config=Config.open(config_file))
tdsr = TDSR1(config=Config.open(config_file))
# lcm  = LCM(config=Config.open(config_file))
cfm = CFM(config=Config.open(config_file))
# trad = Traditional(config=Config.open(config_file))
# rsm  = RSM(config=Config.open(config_file))
rsd = RSD(config=Config.open(config_file))
rsd1 = RSD1(config=Config.open(config_file))

# ----------------------------------------
print("define model parameter for plotting, if different from config file")
# ----------------------------------------
dsig = 1.0  # Asig = dsig in the RS model
depthS = -dsig
dS0 = 8.5  # [MPa] gap for CF and RS_subcrit
dS0_TDSR = 22.0  # [MPa] gap for TDSR with uniform subcritical stress state
dS0mean = 24.5  # [MPa] mean stress for TDSR with gaussian subcritical stress state
dS0std = 1.0  # [MPa] standard deviation " " "
ta = 3e5  # [years] relaxation time ta=dsig/dotsigc
taRSsub = 60.0  # [years] relaxation time in the subcritical RS model
dotsigc = dsig / ta
strend = dotsigc
dotsigc_sub = dsig / taRSsub
t0 = 1e-4  # [years] mean failure time of a source at critical stress in the TDRS model

# ---- discretizing stress axis for integration with
deltaS = -depthS / 60.0  # increment do discretize Coulomb stress axis
sigma_max = 3000.0 * deltaS  # maximum depth on Coulomb axis (limit of integral)
# precision = 18   # wird bisher nicht beneotigt in tdsr

tstart = 1960  # start time for model simulation
tend = 2022.0  # end time for model simulation
Mcut = 1.45  # magnitude cutoff for the fitted earthquakes
T1 = 1991  # start time period for fit
T2 = tend  # end time period for fit

########## INPUTS:
# Earthquake data:
data = np.loadtxt(REPO_ROOT / "data/groningenEQ.dat")
teq = data[:, 0]
meq = data[:, 1]
NEQ = len(teq[((teq >= T1) & (teq <= T2))])
# Pressure data:
data = np.loadtxt(REPO_ROOT / "data/groningenP.dat")
tgrid = data[:, 0]
dt = tgrid[1] - tgrid[0]
nt = int(math.ceil((tend - tstart) / dt))  # =NT = len(t)
p = data[:, 1]

##### Coulomb stress function, saved in tmp file  #####
S = -p
dS = np.ediff1d(S, to_begin=0.0)
S = np.cumsum(dS)
t00 = np.min(tgrid[((S >= dS0) & (tgrid <= teq[0]))])
print(
    "time when EQ start, t00=",
    t00,
    tgrid[0],
    dotsigc * (t00 - tgrid[0]),
    dotsigc_sub * (t00 - tgrid[0]),
)
data_file = REPO_ROOT / "data/groningenCFS.dat"
np.savetxt(
    data_file,
    np.transpose([S, tgrid]),
    header="\n".join(
        [
            "Ascii loading file is expected to have two header lines",
            "1st column: Coulomb stress, 2nd column: time",
        ]
    ),
)

######## MODEL SIMULATIONS: ######

#  ------ loading assuming that faults under critical stress exists --------
loading = ExternalFileLoading(
    _config=tdsr.config,
    iascii=True,
    infile=data_file,
    scal_t=1.0,
    scal_cf=1.0,
    strend=dotsigc,
    c_tstart=tstart,
    tstart=tstart,
    tend=tend,
    deltat=dt,
)

# CF model with critical initial stress state:
config, t, chiz, cf, R_CF, xn = cfm(
    loading=loading,
    chi0=1.0,
    Sshadow=0.0,
    depthS=depthS,
    deltat=dt,
    tstart=tstart,
    tend=tend,
)
fac = NEQ / np.sum(
    R_CF[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt
)  # renormalize to observed EQ-number within [T1, T2]:
R_CF *= fac

# CF model with subcritical initial stress state:
config, t, chiz, cf, R_CFsub, xn = cfm(
    loading=loading,
    chi0=1.0,
    Sshadow=dS0,
    depthS=depthS,
    deltat=dt,
    tstart=tstart,
    tend=tend,
)
fac = NEQ / np.sum(R_CFsub[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_CFsub *= fac

# RS model for critical initial stress state:
config, t, chiz, cf, R_RS, xn = rsd1(
    loading=loading,
    chi0=1.0,
    Sshadow=0.0,
    depthS=depthS,
    deltat=dt,
    tstart=tstart,
    tend=tend,
)
fac = NEQ / np.sum(R_RS[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_RS *= fac

# TDSR model for critical initial stress state:
config, t, chiz, cf, R_TDSR, xn = tdsr(
    loading=loading,
    chi0=1.0,
    t0=t0,
    depthS=depthS,
    Sshadow=0.0,
    deltaS=deltaS,
    sigma_max=sigma_max,
    iX0switch=0,
    deltat=dt,
    taxis_log=0,
    tstart=tstart,
    tend=tend,
)
fac = NEQ / np.sum(R_TDSR[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_TDSR *= fac

#  ------ assuming that no faults exists under critical stress (needed for TDSR and RS) -----
loading = ExternalFileLoading(
    _config=tdsr.config,
    iascii=True,
    infile=data_file,
    scal_t=1.0,
    scal_cf=1.0,
    strend=dotsigc_sub,
    c_tstart=tstart,
    tstart=tstart,
    tend=tend,
    deltat=dt,
)

# RS model for subcritical initial stress state:
config, t, chiz, cf, R_RSsub, xn = rsd1(
    loading=loading,
    chi0=1.0,
    Sshadow=(dS0 + dotsigc_sub * (t00 - tgrid[0])),
    depthS=depthS,
    deltat=dt,
    tstart=tstart,
    tend=tend,
)
# config, t, chiz, cf, R_RSsub, xn = rsd1(loading=loading, chi0=1.0, Sshadow=dS0, depthS=depthS, deltat=dt, tstart=tstart, tend=tend)
fac = NEQ / np.sum(R_RSsub[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_RSsub *= fac

# TDSR model with uniform subcritical stress:
config, t, chiz, cf, R_TDSRsub_uniform, xn = tdsr(
    loading=loading,
    chi0=1.0,
    t0=t0,
    depthS=depthS,
    Sshadow=dS0_TDSR,
    deltaS=deltaS,
    sigma_max=sigma_max,
    iX0switch=1,
    deltat=dt,
    taxis_log=0,
    tstart=tstart,
    tend=tend,
)
fac = NEQ / np.sum(R_TDSRsub_uniform[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_TDSRsub_uniform *= fac

# TDSR model with gaussian subcritical prestress:
config, t, chiz, cf, R_TDSRsub_gauss, xn = tdsr(
    loading=loading,
    chi0=1.0,
    t0=t0,
    depthS=depthS,
    Sshadow=0.0,
    deltaS=deltaS,
    sigma_max=sigma_max,
    iX0switch=2,
    Zmean=dS0mean,
    Zstd=dS0std,
    deltat=dt,
    taxis_log=0,
    tstart=tstart,
    tend=tend,
)
fac = NEQ / np.sum(R_TDSRsub_gauss[((tgrid[0:nt] >= T1) & (tgrid[0:nt] <= T2))] * dt)
R_TDSRsub_gauss *= fac

####### PLOTTING:
plt.rc("font", family="sans-serif")
plt.rc("font", size=12)
plt.rc("legend", fontsize=8)
fig, ax = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 2]}, figsize=(7, 6.5))

# ax[0].plot(Z, XCFc, c='k', ls='dashed', alpha=0.3, label=r'CF')
# ax[0].plot(Z, XCFsub, c='k', alpha=0.3, label=r'CF$_\mathrm{subcrit}$')
# ax[0].plot(Z, Xs, c='b', lw=3, alpha=0.3, label=r'TDRS$_\mathrm{crit}$')
# ax[0].plot(Z, Xuniform, 'r--', label=r'TDSR$_\mathrm{uniform}$')
# ax[0].plot(Z, Xgauss, 'r-', label=r'TDSR$_\mathrm{gauss}$')
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
# ax[1].set_ylim(-1, 30.5)
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

plt.show()
fig.savefig(str(figname) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
fig.savefig(str(figname) + ".png", format="png", dpi=300, bbox_inches="tight")
print("\n\t OUTPUT: %s\n" % (figname))
