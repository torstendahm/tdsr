# --------------------------
# Plot Fig. 4a and b (see Dahm and Hainzl, 2022)
# --------------------------
import sys
from pathlib import Path
import numpy as np
import math
import matplotlib.pyplot as plt

EXAMPLE_DIR = Path(__file__).parent
REPO_ROOT = EXAMPLE_DIR.parent.absolute()
sys.path.insert(0, str(REPO_ROOT))

from tdsr import Config, TDSR, TDSR1, LCM, CFM, Traditional, RSM, RSD
from tdsr.plotting import plot
from tdsr.loading import StepLoading, FourPointLoading, BackgroundLoading, CyclicLoading
from tdsr.utils import Eq7

config_file = EXAMPLE_DIR / "config.toml"
figname = REPO_ROOT / "plots/fig4ab"

print("set-up the tdsr, linear Coulomb failure and Rate and State models")
# tdsm = TDSM(config=Config.open(config_file))
tdsr = TDSR1(config=Config.open(config_file))
# lcm  = LCM(config=Config.open(config_file))
cfm = CFM(config=Config.open(config_file))
# trad = Traditional(config=Config.open(config_file))
# rsm  = RSM(config=Config.open(config_file))
rsd = RSD(config=Config.open(config_file))

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
sstep = [
    -2.0,
    -4.0,
]  # stress steps applied at time tstep, in MPa ( formerly Sstepvalues = (-2.0, -4.0) )
tstep = 0.0 * hours
tstep = 0.01 * hours

tstart = -1.0 * hours  # T0 = -1.0     # simulation start
tend = 10.0 * hours  # T1 = 10.0     # simulation end

# t = np.append((T0, -1e-6), np.linspace(0, T1, 1000))
# S = dotsigc * t
# Ss = np.copy(S)

# deltat = 0.1*hours
# t0 = 0.1*deltat
deltat = 0.02 * hours
t0 = deltat
# t0 = deltat
nt = int(math.ceil((tend - tstart) / deltat))  # =NT = len(t)
taxis_log = 0  # linear time axis discretization = default in config.toml
# ntlog=100
# tstartl = 0.01*hours

# ---- discretizing stress axis for integration with
deltaS = -depthS / 60.0  # increment do discretize Coulomb stress axis
sigma_max = 3000.0 * deltaS  # maximum depth on Coulomb axis (limit of integral)

# ----------------------------------------
# Plot Properties
# ----------------------------------------
# make a bigger global font size and sans-serif style
lstype = ["solid", "dashed"]
plt.rc("font", family="sans-serif")
plt.rc("font", size=12)
plt.rc("legend", fontsize=10)
fig, ax = plt.subplots(2, 1, figsize=(7, 6))
plt.subplots_adjust(hspace=0.0)

# ----------------------------------------
# Calculare earthquake rates
# ----------------------------------------
cfs = np.zeros((2, nt))
cf_shad = np.zeros((2, nt))
r_tdsr = np.zeros((2, nt))
r_lcm = np.zeros((2, nt))
r_rsd = np.zeros((2, nt))
R_TDSR_theory = np.zeros((2, nt))
# R_TDSR = r0 * np.ones(NT)
# R_TDSR_theory = r0 * np.ones(NT)

# ---- calculate equilibrium distribution of tectonic loading, to be used later to avoid initial oscillations
# loading = BackgroundLoading(_config=tdsm.config, strend=strend)
# config, t, chiz_background, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, deltat=deltat, tstart=0, tend=tend)
for i, Sstep in enumerate(sstep):

    print("i=", i, " Sstep=", Sstep, " sstep[i]=", sstep[i], " tstep=", tstep)
    print("deltat=", deltat, " tstart=", tstart, " tend=", tend)
    loading = StepLoading(
        _config=tdsr.config,
        strend=strend,
        sstep=Sstep,
        tstep=tstep,
        taxis_log=0,
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
        tstart=tstart,
        tend=tend,
    )
    cfs[i, :] = cf[:]
    r_tdsr[i, :] = r[:]
    # Z = functions.Zvalues(S, Sstep, t0, dsig)
    # X = functions.X0steady(Z+Sstep, r0, t0, dsig, dotsigc)
    # R_TDSR[2:] = functions.TDSR(t[2:], S[2:], Z, X, t0, dsig)

    # Ss[2:] = S[2:] + Sstep

    loading = StepLoading(
        _config=rsd.config,
        strend=strend,
        sstep=Sstep,
        tstep=tstep,
        deltat=deltat,
        tstart=tstart,
        tend=tend,
    )
    # BUG: cf_shadow is written here, but cf_shad is plotted (zeros only)!
    config, t, cf_shadow, cf, r, xn = rsd(
        loading=loading,
        chi0=chi0,
        depthS=depthS,
        deltat=deltat,
        tstart=tstart,
        tend=tend,
    )
    r_rsd[i, :] = r

    loading = StepLoading(
        _config=cfm.config,
        strend=strend,
        sstep=Sstep,
        tstep=tstep,
        deltat=deltat,
        tstart=tstart,
        tend=tend,
    )
    config, t, chiz, cf, r, xn = cfm(
        loading=loading, chi0=chi0, deltat=deltat, tstart=tstart, tend=tend
    )
    r_lcm[i, :] = r

    R_TDSR_theory[i, :] = Eq7(t, sstep[i], chi0 * strend, -depthS, strend, strend)

# ----------------------------------------
# plot results
# ----------------------------------------
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
        ax[1].plot(t, r_lcm[i, :], c="g", ls="dashed", lw=1, alpha=al[i], label="CF")
    else:
        ax[1].plot(t, r_tdsr[i, :], c="b", ls="solid", lw=3, alpha=al[i])
        ax[1].plot(t, r_lcm[i, :], c="g", ls="dashed", lw=1, alpha=al[i])
    # ax[1].plot(t, r_rsd[i,:], c='r', ls='dashed', lw=1, alpha=al[i])
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

plt.show()
fig.savefig(str(figname) + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
fig.savefig(str(figname) + ".png", format="png", dpi=300, bbox_inches="tight")
print("\n\t OUTPUT: %s\n" % (figname))
