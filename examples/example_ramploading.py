#!/usr/bin/env python
# -*- coding: utf-8 -*-

# --------------------------
# Example for ramp loading
# --------------------------

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

EXAMPLE_DIR = Path(__file__).parent
REPO_ROOT = EXAMPLE_DIR.parent.absolute()
sys.path.insert(0, str(REPO_ROOT))

from tdsr import TDSR, TDSR1, LCM, Traditional, CFM, RSM, RSD  # noqa: E402
from tdsr.loading import (  # noqa: E402
    BackgroundLoading,
    TrendchangeLoading,
    RampLoading,
)


config_file = EXAMPLE_DIR / "config.toml"
figname = REPO_ROOT / "plots/fig_ramploading"

tdsm = TDSR()
tdsr = TDSR1()
lcm = LCM()
trad = Traditional()
cfm = CFM()
rsm = RSM()
rsd = RSD()


hours = 3600.0
tstart = 0 * hours
tend = 200 * hours
deltat = 600.0
t0 = 1.0 * deltat
nt = np.floor((tend - tstart) / deltat).astype(int)

strend1 = 2.5e-6
strend2 = 1.0e-5
strend3 = strend1
nsample2 = 100

sstep = 0.0
# time when stress step is acting
tstep = 80.0 * hours
nsample1 = np.floor((tstep - tstart) / deltat).astype(int) + 1

tstep3 = tstep + nsample2 * deltat
nstep = np.floor((tstep - tstart) / deltat).astype(int) + 1

# susceptibility to trigger earthquakes by unit stress increase
chi0 = 1.0e4
# skin depth in MPa (must be negative)
depthS = -0.3

# increment do discretize Coulomb stress axis
deltaS = -depthS / 500.0
# maximum depth on Coulomb axis (limit of integral)
sigma_max = 10000.0 * deltaS
precision = 18

# steady state distribution
iX0 = "equilibrium"

print(
    dict(
        deltat=deltat,
        sec=deltat / hours,
        hours=hours,
        end=tend / hours,
        samples=nt,
    )
)

print(
    dict(
        Ta1=-depthS / strend1,
        Ta2=-depthS / strend2,
        sec=-depthS / (strend1 * hours),
        chi0=chi0,
        skinDepth=-depthS,
        deltaS=-deltaS,
        tectonic_trend1=strend1 * 24.0 * hours,
        tectonic_trend2=strend2 * 24.0 * hours,
        sigma_max=-sigma_max,
        dot_sigma1=-strend1 * deltat / depthS,
        dot_sigma2=-strend2 * deltat / depthS,
    )
)

# calculate earthquake rates with tdsm, lcm and rsm
cfs1 = np.zeros(nt)
cfs1_ramp = np.zeros(nt)
cfs2 = np.zeros(nt)

r_tdsm = np.zeros(nt)
r_tdsr = np.zeros(nt)
r_lcm = np.zeros(nt)
r_rsm = np.zeros(nt)
n_tdsm = np.zeros(nt - 1)
n_tdsd = np.zeros(nt - 1)
n_lcm = np.zeros(nt - 1)
n_rsm = np.zeros(nt - 1)

r_theo1 = np.zeros(nt)
r_theo2 = np.zeros(nt)

r_tdsm_ramp = np.zeros(nt)
r_tdsr_ramp = np.zeros(nt)
r_lcm_ramp = np.zeros(nt)
r_rsm_ramp = np.zeros(nt)
r_rsd_ramp = np.zeros(nt)
n_tdsm_ramp = np.zeros(nt - 1)
n_tdsr_ramp = np.zeros(nt - 1)
n_lcm_ramp = np.zeros(nt - 1)
n_rsm_ramp = np.zeros(nt - 1)
n_rsd_ramp = np.zeros(nt - 1)

r_theo1_ramp = np.zeros(nt)
r_theo2_ramp = np.zeros(nt)
r_theo3_ramp = np.zeros(nt)
r_theo1_ramp2 = np.zeros(nt)
r_theo2_ramp2 = np.zeros(nt)
r_theo3_ramp2 = np.zeros(nt)

r_tdsm2 = np.zeros(nt)
r_lcm2 = np.zeros(nt)
r_rsm2 = np.zeros(nt)
n_tdsm2 = np.zeros(nt - 1)
n_lcm2 = np.zeros(nt - 1)
n_rsm2 = np.zeros(nt - 1)

r_theo12 = np.zeros(nt)
r_theo22 = np.zeros(nt)

# calculate equilibrium distribution of tectonic loading,
# to be used later to avoid initial oscillations
loading = BackgroundLoading(strend=strend1, deltat=deltat, tstart=tstart, tend=tend)
t, chiz_background, cf, r, xn = tdsm(
    loading=loading,
    chi0=chi0,
    depthS=depthS,
    deltaS=deltaS,
    sigma_max=sigma_max,
    precision=precision,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)

rinfty1 = strend1 * chi0
Ta1 = -depthS / strend1
good = t >= (tstep - tstart)
good3 = t >= (tstep3 - tstart)

# --- simulation of trend change  at time tstep
loading = TrendchangeLoading(
    strend=strend1,
    strend2=strend2,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = tdsm(
    loading=loading,
    equilibrium=True,
    chiz=chiz_background,
    chi0=chi0,
    depthS=depthS,
    deltaS=deltaS,
    sigma_max=sigma_max,
    precision=precision,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
cfs1 = cf
r_tdsm = r
n_tdsm = xn

loading = TrendchangeLoading(
    strend=strend1,
    strend2=strend2,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = tdsr(
    loading=loading,
    chi0=chi0,
    depthS=depthS,
    deltaS=deltaS,
    sigma_max=sigma_max,
    precision=precision,
    iX0=iX0,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
r_tdsr = r
n_tdsr = xn

loading = TrendchangeLoading(
    strend=strend1,
    strend2=strend2,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = trad(
    loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend
)
r_lcm = r
n_lcm = xn

loading = TrendchangeLoading(
    strend=strend1,
    strend2=strend2,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = rsm(
    loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend
)
r_rsm = r
n_rsm = xn

# -------------------------------------------------
# Figure c, extensionm to three slopes to represent ramp function
# -------------------------------------------------
loading = RampLoading(
    strend=strend1,
    strend2=strend2,
    strend3=strend3,
    nsample2=nsample2,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = tdsm(
    loading=loading,
    equilibrium=True,
    chiz=chiz_background,
    chi0=chi0,
    depthS=depthS,
    deltaS=deltaS,
    sigma_max=sigma_max,
    precision=precision,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
cfs1_ramp = cf
r_tdsm_ramp = r
n_tdsm_ramp = xn

loading = RampLoading(
    strend=strend1,
    strend2=strend2,
    strend3=strend3,
    nsample2=nsample2,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = tdsr(
    loading=loading,
    chi0=chi0,
    depthS=depthS,
    deltaS=deltaS,
    sigma_max=sigma_max,
    precision=precision,
    iX0=iX0,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
r_tdsr_ramp = r
n_tdsr_ramp = xn

loading = RampLoading(
    strend=strend1,
    strend2=strend2,
    strend3=strend3,
    nsample2=nsample2,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = trad(
    loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend
)
r_lcm_ramp = r
n_lcm_ramp = xn

loading = RampLoading(
    strend=strend1,
    strend2=strend2,
    strend3=strend3,
    nsample2=nsample2,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = rsm(
    loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend
)
r_rsm_ramp = r
n_rsm_ramp = xn

loading = RampLoading(
    strend=strend1,
    strend2=strend2,
    strend3=strend3,
    nsample2=nsample2,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = rsd(
    loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend
)
r_rsd_ramp = r
n_rsd_ramp = xn

# -----------------------------------------------------
print("second case (figure b) , two slopes, large trend before small trend")
# -----------------------------------------------------
loading = BackgroundLoading(strend=strend2, deltat=deltat, tstart=tstart, tend=tend)
t, chiz_background, cf, r, xn = tdsm(
    loading=loading,
    chi0=chi0,
    depthS=depthS,
    deltaS=deltaS,
    sigma_max=sigma_max,
    precision=precision,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)

rinfty12 = strend2 * chi0
Ta12 = -depthS / strend2

# --- simulation of trend change  at time tstep
loading = TrendchangeLoading(
    strend=strend2,
    strend2=strend1,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = tdsm(
    loading=loading,
    equilibrium=True,
    chiz=chiz_background,
    chi0=chi0,
    depthS=depthS,
    deltaS=deltaS,
    sigma_max=sigma_max,
    precision=precision,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
cfs2 = cf
r_tdsm2 = r
n_tdsm2 = xn

loading = TrendchangeLoading(
    strend=strend2,
    strend2=strend1,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = tdsr(
    loading=loading,
    chi0=chi0,
    depthS=depthS,
    deltaS=deltaS,
    sigma_max=sigma_max,
    precision=precision,
    iX0=iX0,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
r_tdsr2 = r
n_tdsr2 = xn

loading = TrendchangeLoading(
    strend=strend2,
    strend2=strend1,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = trad(
    loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend
)
r_lcm2 = r
n_lcm2 = xn

loading = TrendchangeLoading(
    strend=strend2,
    strend2=strend1,
    tstep=tstep,
    deltat=deltat,
    tstart=tstart,
    tend=tend,
)
t, chiz, cf, r, xn = rsm(
    loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend
)
r_rsm2 = r
n_rsm2 = xn

# ----------------------------------------------------
# ------  calculate analytical approximations
# ----------------------------------------------------
Ta1 = -depthS / strend1
Ta2 = -depthS / strend2
Ta3 = -depthS / strend3
Ta0 = 0.0
rinfty1 = strend1 * chi0
rinfty2 = strend2 * chi0
rinfty3 = strend3 * chi0

dS = -depthS
t0 = deltat / 10.0
f = 1.0
fact = 1.0
fact1 = 1.0 - np.exp(-(t + Ta1) * np.exp(-(strend2 * t0 / dS)) / t0)
fact2 = 1.0 - np.exp(-(t + Ta2) * np.exp(-(strend2 * t0 / dS)) / t0)
fact1 = -np.exp(-(t + Ta1) / t0)
fact2 = -np.exp(-(t + Ta2) / t0)
fact2 = 1.0
fact1 = fact2
# ---- final correct equation
fact1 = np.exp(-t / Ta1)
fact2 = np.exp(-t / Ta2)
r_theo1 = rinfty2 / (fact2 * (Ta1 - Ta2) / Ta2 + 1)
r_theo2 = rinfty1 / (fact1 * (Ta2 - Ta1) / Ta1 + 1)


# ---------- ramp scenario very rough approximation - not general
r_theo1_ramp = -chi0 * depthS * np.exp(-t / Ta1) / (deltat + t)
r_theo1_ramp = r_theo1_ramp + rinfty1
Ta = 0.25 * (Ta1 + 3 * Ta2)
Ta = 0.5 * (Ta1 + Ta2)
r_theo2_ramp[good] = (
    -chi0 * depthS * np.exp(-(t[good] - tstep) / Ta) / (Ta + (t[good] - tstep))
)
rmax = np.amax(r_theo2_ramp[good])
r_theo2_ramp[good] = (rinfty2 - rinfty1) * (1.0 - r_theo2_ramp[good] / rmax)

Ta = 0.25 * (Ta2 + 3 * Ta3)
Ta = 0.5 * (Ta2 + Ta3)
r_theo3_ramp[good3] = (
    -chi0 * depthS * np.exp(-(t[good3] - tstep) / Ta) / (Ta + (t[good3] - tstep))
)
rmax3 = np.amax(r_theo3_ramp[good3])
r_theo3_ramp[good3] = (rinfty3 - rinfty2) * (1.0 - r_theo3_ramp[good3] / rmax3)

# ---------- ramp scenario , new analytical solution a la Scholz, 30 Jan 2022
dS = -depthS
t0 = deltat / 10.0
t1 = 100 * tstep  # in der simulation wird vorab ein Gleichgewicht chi bestimmt
Ta1b = Ta1
r_theo1_ramp2 = rinfty1 + chi0 * dS * (
    1.0 - np.exp(-np.exp(-t / Ta1b) * (t1 + t) / t0)
) / (t1 + t)

tstep2b = tstep
t0 = deltat / 10.0
t1 = deltat
# Ta aendern bringt nicht viel
# da fuer kleine lag times offensichtlich wenig Einfluss
# - wird ebenso durch t0 bestimmt
Ta2b = Ta2
r_theo2_ramp2[good] = (
    chi0
    * dS
    * (
        1.0
        - np.exp(-np.exp((t[good] - tstep2b) / Ta2b) * (t1 + t[good] - tstep2b) / t0)
    )
    / (t1 + (t[good] - tstep2b))
)
rmax = np.amax(r_theo2_ramp2[good])
print("r_theo2 max =", rmax / rinfty1)
r_theo2_ramp2[good] = (rinfty2 - rinfty1) * (1.0 - r_theo2_ramp2[good] / rmax)

tstep3b = tstep3
t0 = deltat / 10.0
t1 = deltat
Ta3b = Ta3
r_theo3_ramp2[good3] = (
    chi0
    * dS
    * (
        1.0
        - np.exp(-np.exp((t[good3] - tstep3b) / Ta3b) * (t1 + t[good3] - tstep3b) / t0)
    )
    / (t1 + (t[good3] - tstep3b))
)
rmax3 = np.amax(r_theo3_ramp2[good3])
r_theo3_ramp2[good3] = (rinfty3 - rinfty2) * (1.0 - r_theo3_ramp2[good3] / rmax3)


# ------ second case -  calculate analytical approximations
Ta12 = -depthS / strend2
Ta22 = -depthS / strend1
rinfty12 = strend2 * chi0
rinfty22 = strend1 * chi0

Ta = 0.5 * (Ta12 + Ta22)
r_theo12 = -chi0 * depthS * np.exp(-t / Ta12) / (deltat + t)
r_theo12 = r_theo12 + rinfty12
r_theo22[good] = (
    -chi0 * depthS * np.exp(-(t[good] - tstep) / Ta) / (Ta + (t[good] - tstep))
)
rmax = np.amax(r_theo22[good])
print("r_theo2 max =", rmax / rinfty1)
r_theo22[good] = (rinfty22 - rinfty12) * (1.0 - r_theo22[good] / rmax)

# -----------------------------------------------------
print("Plotting Fig. 4d and 4e Appendix")
# -----------------------------------------------------
xrange_set = False
xrange_set = True

ymaxs = (strend1 * tstep + strend2 * (tstep3 - tstep)) / (Ta1 * strend1) + 1.5
lstyle = "-"

if xrange_set:
    nstart = 0
    if strend1 <= strend2:
        tmin = -0.5
        tmax = 2.5
    else:
        tmin = -0.5
        tmax = 10.0
else:
    nstart = 0

# ---- Coulomb stress loading (input and interpolated input)
fig = plt.figure(1, figsize=(12, 6))
ax1a = fig.add_subplot(131)
ax1a.set_ylabel("$\\sigma_c/\\dot{\\sigma}_{c0} T_{0}$", fontsize=20)
ax1a.set_xlabel(r"$(t-t_1)/T_{0}$", fontsize=20)
plt.figtext(0.09, 0.87, "a)", fontsize=20)
ax1a.tick_params(axis="both", which="major", labelsize=20)
ax1a.tick_params(axis="both", which="minor", labelsize=18)

ax1a.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    cfs1[nstart:] / (Ta1 * strend1),
    linewidth=2.0,
    ls=lstyle,
    color="black",
    label=r"$\dot{\sigma}_{c1}/\dot{\sigma}_{c0}=$"
    + "{:.1f}".format(float(strend2 / strend1)),
)

plt.legend(loc="upper left", fontsize=20)
if xrange_set:
    plt.xlim([-2, 5])
    plt.ylim([0, 15])

# ----------------- zweites Beispiel stress loading
ax12a = fig.add_subplot(132)
ax12a.set_xlabel(r"$(t-t_1)/T_{0}$", fontsize=20)
plt.figtext(0.365, 0.87, "b)", fontsize=20)
ax12a.tick_params(axis="both", which="major", labelsize=20)
ax12a.tick_params(axis="both", which="minor", labelsize=18)

ax12a.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    cfs2[nstart:] / (Ta1 * strend1),
    linewidth=2.0,
    ls=lstyle,
    color="black",
    label=r"$\dot{\sigma}_{c1}/\dot{\sigma}_{c0}=$"
    + "{:.1f}".format(float(strend1 / strend2)),
)

plt.legend(loc="lower right", fontsize=20)
if xrange_set:
    plt.xlim([-2, 5])
    plt.ylim([0, 15])

# ----------------- drittes Beispiel ramp function
ax13a = fig.add_subplot(133)
ax13a.set_xlabel(r"$(t-t_1)/T_{0}$", fontsize=20)
plt.figtext(0.635, 0.87, "c)", fontsize=20)
ax13a.tick_params(axis="both", which="major", labelsize=20)
ax13a.tick_params(axis="both", which="minor", labelsize=18)

ax13a.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    cfs1_ramp[nstart:] / (Ta1 * strend1),
    linewidth=2.0,
    ls=lstyle,
    color="black",
)

if xrange_set:
    plt.xlim([-2, 5])
    plt.ylim([0, 15])

plt.show()

# ----------------------------------------------------------
# ---- earthquake rate ------------------------------------
# ----------------------------------------------------------
if strend1 <= strend2:
    ymin = 0.0 * strend2 / strend1
    ymax = 1.2 * strend2 / strend1
else:
    ymin = 0.8 * strend2 / strend1
    ymax = 1.1
fig = plt.figure(2, figsize=(12, 6))
ax1b = fig.add_subplot(131)
ax1b.set_xlabel(r"$(t-t_1)/T_{0}$", fontsize=20)
ax1b.set_ylabel(r"$r / \chi_0 V \dot{\sigma}_{c0}$", fontsize=20)
ax1b.tick_params(axis="both", which="major", labelsize=20)
ax1b.tick_params(axis="both", which="minor", labelsize=18)

# dotted gray line as reference at y=1
ax1b.plot(
    [(t[nstart] - t[nstep]) / Ta1, (t[nt - 1] - t[nstep]) / Ta1],
    [1, 1],
    linewidth=1.0,
    ls="dotted",
    color="gray",
)

ax1b.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    r_lcm[nstart:] / rinfty1,
    linewidth=1.0,
    ls="--",
    color="gray",
    label=r"LCM",
)
ax1b.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    r_tdsr[nstart:] / rinfty1,
    linewidth=3.5,
    ls="-",
    color="red",
    label=r"TDSR",
)
ax1b.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    r_rsm[nstart:] / rinfty1,
    linewidth=2.5,
    ls="--",
    color="blue",
    label=r"RSM",
)
ax1b.plot(
    (t[nstart:] - 0.0 * t[nstep]) / Ta1,
    r_theo1[nstart:] / rinfty1,
    linewidth=2.0,
    ls="--",
    color="blue",
    label=r"theory",
)

# -------------- first part of response in black ddotted line ---
plt.legend(loc="lower right", fontsize=20)
plt.figtext(0.09, 0.87, "a)", fontsize=20)
plt.figtext(
    0.13, 0.82, r"$T_{0}/T_{1}=$" + "{:.1f}".format(float(Ta1 / Ta2)), fontsize=20
)
if xrange_set:
    plt.xlim([tmin, tmax])
    plt.ylim([ymin, ymax])
else:
    print(" nothing defined")

# ---- second case in right panel ------------------------------------
ax12b = fig.add_subplot(132)
ax12b.set_xlabel(r"$(t-t_1)/T_{0}$", fontsize=20)
ax12b.tick_params(axis="both", which="major", labelsize=20)
ax12b.tick_params(axis="both", which="minor", labelsize=18)

ax12b.plot(
    [(t[nstart] - t[nstep]) / Ta1, (t[nt - 1] - t[nstep]) / Ta12],
    [1, 1],
    linewidth=1.0,
    ls="dotted",
    color="gray",
)


ax12b.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    r_lcm2[nstart:] / rinfty1,
    linewidth=1.0,
    ls="--",
    color="gray",
    label=r"LCM",
)
ax12b.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    r_tdsr2[nstart:] / rinfty1,
    linewidth=3.5,
    ls="-",
    color="red",
    label=r"TDSR",
)
ax12b.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    r_rsm2[nstart:] / rinfty1,
    linewidth=2.5,
    ls="--",
    color="blue",
    label=r"RSM",
)
ax12b.plot(
    (t[nstart:] - 0.0 * t[nstep]) / Ta1,
    r_theo2[nstart:] / rinfty1,
    linewidth=2.0,
    ls="--",
    color="blue",
    label=r"theory",
)

# -------------- second part of response in black dotted line ---
plt.figtext(0.365, 0.87, "b)", fontsize=20)
plt.figtext(
    0.41, 0.82, r"$T_{0}/T_{1}=$" + "{:.1f}".format(float(Ta12 / Ta22)), fontsize=20
)
if xrange_set:
    plt.xlim([tmin, tmax])
    plt.ylim([ymin, ymax])
else:
    print("nothing defined")

# ---- third case, ramp function  ------------------------------------
ax13b = fig.add_subplot(133)
ax13b.set_xlabel(r"$(t-t_1)/T_{0}$", fontsize=20)
ax13b.tick_params(axis="both", which="major", labelsize=20)
ax13b.tick_params(axis="both", which="minor", labelsize=18)

ax13b.plot(
    [(t[nstart] - t[nstep]) / Ta1, (t[nt - 1] - t[nstep]) / Ta12],
    [1, 1],
    linewidth=1.0,
    ls="dotted",
    color="gray",
)

ax13b.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    r_lcm_ramp[nstart:] / rinfty1,
    linewidth=1.0,
    ls="--",
    color="gray",
    label=r"LCM",
)
ax13b.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    r_tdsr_ramp[nstart:] / rinfty1,
    linewidth=3.5,
    ls="-",
    color="red",
    label=r"TDSR",
)
ax13b.plot(
    (t[nstart:] - t[nstep]) / Ta1,
    r_rsm_ramp[nstart:] / rinfty1,
    linewidth=1.0,
    ls="--",
    color="blue",
    label=r"RSM",
)

yval = r_theo1[nstart + nsample2]
if strend1 <= strend2:
    n2nd2 = np.argmax(r_theo2 <= yval)
else:
    n2nd2 = np.argmax(r_theo2 >= yval)
ti = -Ta1 * np.log(1.0 - np.exp(-t[nsample2] / Ta2))
n2nd = np.floor(ti / deltat).astype(int) + 1

print(
    "yval=",
    yval / rinfty1,
    " nsample2=",
    nsample2,
    " n2nd_num=",
    n2nd2,
    " nstep=",
    nstep,
)
print("t(n2nd)==", t[n2nd], " ti=", ti, " n2nd=", n2nd)
ax13b.plot(
    (t[nstart : nstart + nsample2]) / Ta1,
    r_theo1[nstart : nstart + nsample2] / rinfty1,
    linewidth=2.0,
    ls="--",
    color="blue",
    label=r"theory",
)
ax13b.plot(
    (t[n2nd:] + t[nsample2] - t[n2nd]) / Ta1,
    r_theo2[n2nd:] / rinfty1,
    linewidth=2.0,
    ls="--",
    color="blue",
    label=r"theory",
)

plt.figtext(0.635, 0.87, "c)", fontsize=20)
if xrange_set:
    plt.xlim([tmin, tmax])
    plt.ylim([ymin, ymax])
else:
    print("nothing defined")

plt.show()
fig.savefig(str(figname) + ".pdf", format="pdf", bbox_inches="tight")
fig.savefig(str(figname) + ".png", format="png", bbox_inches="tight")
print(f"\n\t OUTPUT: {str(figname) + '.pdf'}\n")
