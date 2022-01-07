################################
# Time Dependent Seismicity Model - Plotting results
# T. Dahm, R. Dahm 26.12.2021
################################
"""
    Add description here

    Input:
    return:  tau_rake, tau_plane, tau_max, sign, sign_max 
"""

import math
import pickle

import matplotlib
import numpy as num
from scipy import special as spl

##matplotlib.use('PS')
##matplotlib.use('TkAgg')
matplotlib.use("Qt4Agg")
import matplotlib as mpl
import matplotlib.cm as mpl_cm
import matplotlib.colors as mpl_colors
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

km = 1000.0
hours = 3600.0
days = hours * 24.0

fn_in = "results/tdsm.out"
out = pickle.load(open(fn_in, "rb"))
# --- read meta data parameter
chi0 = out["meta"]["chi0"]
Sshadow = out["meta"]["Sshadow"]
deltaS = out["meta"]["deltaS"]
sigma_max = out["meta"]["sigma_max"]
depthS = out["meta"]["depthS"]
nt = out["meta"]["nt"]
n1 = out["meta"]["n1"]
deltat = out["meta"]["deltat"]
depthS = out["meta"]["depthS"]
chi0 = out["meta"]["chi0"]
sc0 = out["meta"]["sc0"]
sc1 = out["meta"]["sc1"]
sc2 = out["meta"]["sc2"]
sc3 = out["meta"]["sc3"]
tectstressrate = out["meta"]["tectstressrate"]
# ---- read  time series
(
    t,
    cf,
    ratez,
    neqz,
) = out["data"]

file = "CFMprob_ramp2_paper_chi%2.0f" % (chi0)
file2 = "CFMprob_ramp2_paper_totalnumber_chi%2.0f" % (chi0)

Ta = deltat * (num.exp(1.0) - 1.0)
scal = -depthS * chi0 / deltat
# fig = plt.figure(1, figsize=(16, 12))
fig = plt.figure(1, figsize=(12, 8))
tmin = -5
tmax = 25

# ---- Coulomb stress loading (input and interpolated input
ax1a = fig.add_subplot(211)
ax1a.set_ylabel("$\sigma_c/\dot{\sigma}_c\Delta t$", fontsize=20)
ax1a.tick_params(axis="both", which="major", labelsize=20)
ax1a.tick_params(axis="both", which="minor", labelsize=18)

fa = (sc2 - sc1) / (tectstressrate * deltat)
ax1a.plot(
    (t - t[n1]) / Ta,
    cf / (tectstressrate * deltat),
    linewidth=2.0,
    ls="-",
    color="black",
    label=r"$\Delta\sigma/\dot{\sigma}_c\Delta t=$" + "{:d}".format(num.int(fa)),
)

plt.legend(loc="lower right", fontsize=20)
plt.xlim([tmin, tmax])

# ---- Coulomb stress change (at mid sample positions)
ax1b = fig.add_subplot(212)
fa = -depthS / (tectstressrate * deltat)
fb = 0.0
ax1b.set_xlabel(r"$(t-t_0)/\Delta t(e-1)$", fontsize=20)
ax1b.set_ylabel(r"$r / \chi_0 V \delta \dot{\sigma}$", fontsize=20)
ax1b.tick_params(axis="both", which="major", labelsize=20)
ax1b.tick_params(axis="both", which="minor", labelsize=18)

ax1b.plot(
    (t[40 : nt - 1] - t[n1]) / Ta,
    ratez[40 : nt - 1] / scal,
    linewidth=3.5,
    color="red",
    label=r"$\delta \dot{\sigma}/\dot{\sigma}_c=$" + "{:d}".format(num.int(fa)),
)

plt.legend(loc="upper right", fontsize=20)
plt.figtext(0.02, 0.87, "a)", fontsize=20)
plt.xlim([tmin, tmax])
ax1b.set_yscale("log")

# ----- plot and save figure ----
fig.savefig("plots/" + file + ".pdf", format="pdf", bbox_inches="tight")
fig.savefig(
    "plots/" + file + ".png",
    format="png",
    resolution=100.0,
    oversample=2.0,
    size=None,
    width=None,
    height=None,
    psconvert=True,
    bbox_inches="tight",
)

plt.show()

# ------- plot absolute number of earthquakes for own and Coulomb failure models
fig = plt.figure(1, figsize=(12, 4))
ax1c = fig.add_subplot(111)
ax1c.set_xlabel(r"$(t-t_0)/\Delta t(e-1)$", fontsize=20)
ax1c.set_ylabel(r"$n-n[-1]$", fontsize=20)
ax1c.tick_params(axis="both", which="major", labelsize=20)
ax1c.tick_params(axis="both", which="minor", labelsize=18)

ax1c.plot(
    (t[40 : nt - 1] - t[n1]) / Ta,
    neqz[40 : nt - 1] - neqz[n1 - 1],
    linewidth=2.5,
    color="red",
    label=r"$\delta \dot{\sigma}/\dot{\sigma}_c=$" + "{:d}".format(num.int(fa)),
)

plt.legend(loc="lower right", fontsize=20)
plt.figtext(0.02, 0.87, "a)", fontsize=20)
plt.xlim([tmin, tmax])
plt.show()

# ----- plot and save figure ----
fig.savefig("plots/" + file2 + ".pdf", format="pdf", bbox_inches="tight")
fig.savefig(
    "plots/" + file2 + ".png",
    format="png",
    resolution=250.0,
    oversample=2.0,
    size=None,
    width=None,
    height=None,
    psconvert=True,
    bbox_inches="tight",
)
