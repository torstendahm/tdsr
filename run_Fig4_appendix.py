#--------------------------
# Plot Fig. 4d Appendix) simulation of trend changes
#--------------------------
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, BackgroundLoading, TrendchangeLoading, RampLoading
import matplotlib.pyplot as plt
import matplotlib.ticker
import math
import numpy as np

current_dir = Path(__file__).parent
config_file = current_dir / "config.toml"
pdf_file1 = current_dir / "plots/Dahm_fig4d_appendix"

print("set-up the tdsm, lcm and rsm class environments")
tdsm = TDSM(config=Config.open(config_file))
lcm  = LCM(config=Config.open(config_file))
trad = Traditional(config=Config.open(config_file))
rsm  = RSM(config=Config.open(config_file))


print("define model parameter for plotting, if different from config file")
hours = 3600.
tstart = 0*hours
#tend   = 30*hours
tend   = 200*hours
#deltat = 0.2*hours
deltat = 600.0
#nt = int(math.ceil((tend - tstart) / deltat))
nt = np.floor( (tend-tstart) / deltat).astype(int)

strend1  = 2.5E-6
strend2  = 2.5E-5
#strend2  = 2.5E-4
strend2  = 5.0E-6
strend3  = strend1
nsample2 = 100
#nsample2 = 150
#strend2  = 2.5E-4
#nsample2 = 5
#nsample2 = 10

#sstep  =  [0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 3.0]  # stress step in MPa
sstep  =  0.0
tstep  = 0.2*hours         # time when stress step is acting
tstep  = 3.0*hours         # time when stress step is acting
tstep  = 5.0*hours         # time when stress step is acting
#tstep  = 40.0*hours         # time when stress step is acting
tstep  = 80.0*hours         # time when stress step is acting

tstep3 = tstep + nsample2*deltat

nstep = np.floor( (tstep-tstart) / deltat).astype(int)+1


chi0   = 1.E4            # susceptibility to trigger earthquakes by unit stress increase
depthS = -0.3            # skin depth in MPa (must be defined negativ)

#deltaS    = -depthS/60.  # increment do discretize Coulomb stress axis
#sigma_max = 3000.*deltaS # maximum depth on Coulomb axis (limit of integral)
deltaS    = -depthS/500.  # increment do discretize Coulomb stress axis
sigma_max = 10000.*deltaS # maximum depth on Coulomb axis (limit of integral)
precision = 18

print('deltat=',deltat,' sec = ',deltat/hours,' hours, from 0 to ',tend/hours,' hours, samples =',nt)
print('Ta1  =',-depthS/strend1,' sec = ',-depthS/(strend1*hours),' hours')
print('Ta2  =',-depthS/strend2,' sec = ',-depthS/(strend2*hours),' hours')
print('suscept. chi0=',chi0,' #/MPa')
#print('stress step=[',sstep,'] MPa')
print('skin depth =',-depthS,' MPa')
print('tectonic trend1 =',strend1*24.*hours,' MPa/days,  dot(sigma)*deltat/deltaS=',-strend1*deltat/depthS)
print('tectonic trend2 =',strend2*24.*hours,' MPa/days,  dot(sigma)*deltat/deltaS=',-strend2*deltat/depthS)
print('deltaS=',deltaS,' MPa, from 0 to ',sigma_max,' MPa')

#-----------------------------------------------------
print("Calculate earthquake rates with tdsm, lcm and rsm")
#-----------------------------------------------------
#ns   = len(sstep)
cfs1   = np.zeros(nt)
cfs1_ramp   = np.zeros(nt)
cfs2   = np.zeros(nt)

r_tdsm = np.zeros(nt)
r_lcm  = np.zeros(nt)
r_rsm  = np.zeros(nt)
n_tdsm = np.zeros(nt-1)
n_lcm  = np.zeros(nt-1)
n_rsm  = np.zeros(nt-1)

r_theo1 = np.zeros(nt)
r_theo2 = np.zeros(nt)

r_tdsm_ramp = np.zeros(nt)
r_lcm_ramp  = np.zeros(nt)
r_rsm_ramp  = np.zeros(nt)
n_tdsm_ramp = np.zeros(nt-1)
n_lcm_ramp  = np.zeros(nt-1)
n_rsm_ramp  = np.zeros(nt-1)

r_theo1_ramp = np.zeros(nt)
r_theo2_ramp = np.zeros(nt)
r_theo3_ramp = np.zeros(nt)
r_theo1_ramp2 = np.zeros(nt)
r_theo2_ramp2 = np.zeros(nt)
r_theo3_ramp2 = np.zeros(nt)

r_tdsm2 = np.zeros(nt)
r_lcm2  = np.zeros(nt)
r_rsm2  = np.zeros(nt)
n_tdsm2 = np.zeros(nt-1)
n_lcm2  = np.zeros(nt-1)
n_rsm2  = np.zeros(nt-1)

r_theo12 = np.zeros(nt)
r_theo22 = np.zeros(nt)

# ---- calculate equilibrium distribution of tectonic loading, to be used later to avoid initial oscillations
# loading = BackgroundLoading(_config=tdsm.config, strend=strend)
loading = BackgroundLoading(_config=tdsm.config, strend=strend1, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz_background, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)

print('background')
print('len(cf)',len(cf))
print('len(chiz_background)',len(chiz_background))
rinfty1 = strend1*chi0
Ta1   = -depthS/strend1
good = t >= (tstep-tstart)
good3 = t >= (tstep3-tstart)

# --- example for step loading
#loading = StepLoading(_config=tdsm.config, strend=strend1, sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
#config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
#cfs1    = cf
#r_tdsm = r
#n_tdsm = xn
#
#loading = StepLoading(_config=trad.config, strend=strend, sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
#config, t, chiz, cf, r, xn = trad(loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend)
#r_lcm = r
#n_lcm = xn
#
#loading = StepLoading(_config=rsm.config, strend=strend, sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
#config, t, chiz, cf, r, xn = rsm(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
#r_rsm = r
#n_rsm = xn

# --- simulation of trend change  at time tstep
loading = TrendchangeLoading(_config=tdsm.config, strend=strend1, strend2=strend2, tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
#config, t, chiz, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
cfs1    = cf
r_tdsm = r
n_tdsm = xn

loading = TrendchangeLoading(_config=trad.config, strend=strend1, strend2=strend2, tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = trad(loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend)
r_lcm = r
n_lcm = xn

loading = TrendchangeLoading(_config=rsm.config, strend=strend1, strend2=strend2, tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = rsm(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
r_rsm = r
n_rsm = xn

#-------------------------------------------------
# Figure c, extensionm to three slopes to represent ramp function
#-------------------------------------------------
loading = RampLoading(_config=tdsm.config, strend=strend1, strend2=strend2, strend3=strend3, nsample2=nsample2, tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
#config, t, chiz, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
cfs1_ramp    = cf
r_tdsm_ramp = r
n_tdsm_ramp = xn

loading = RampLoading(_config=trad.config, strend=strend1, strend2=strend2, strend3=strend3, nsample2=nsample2, tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = trad(loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend)
r_lcm_ramp = r
n_lcm_ramp = xn

loading = RampLoading(_config=rsm.config, strend=strend1, strend2=strend2, strend3=strend3, nsample2=nsample2, tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = rsm(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
r_rsm_ramp = r
n_rsm_ramp = xn

#-----------------------------------------------------
print("second case (figure b) , two slopes, large trend before small trend")
#-----------------------------------------------------
loading = BackgroundLoading(_config=tdsm.config, strend=strend2, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz_background, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)

rinfty12 = strend2*chi0
Ta12   = -depthS/strend2

# --- simulation of trend change  at time tstep
loading = TrendchangeLoading(_config=tdsm.config, strend=strend2, strend2=strend1, tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
#config, t, chiz, cf2, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
cfs2    = cf
r_tdsm2 = r
n_tdsm2 = xn

loading = TrendchangeLoading(_config=trad.config, strend=strend2, strend2=strend1, tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = trad(loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend)
r_lcm2 = r
n_lcm2 = xn

loading = TrendchangeLoading(_config=rsm.config, strend=strend2, strend2=strend1, tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = rsm(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
r_rsm2 = r
n_rsm2 = xn

#----------------------------------------------------
# ------  calculate analytical approximations
#----------------------------------------------------
Ta1   = -depthS/strend1
Ta2   = -depthS/strend2
Ta3   = -depthS/strend3
Ta0 = 0.0
rinfty1 = strend1*chi0
rinfty2 = strend2*chi0
rinfty3 = strend3*chi0

r_theo1 = -chi0*depthS*np.exp(-t/Ta1)/(deltat + t)
r_theo1 = r_theo1 +rinfty1
Ta = 0.5*(Ta2+Ta3)
#r_theo2[good] = -chi0*depthS*np.exp(-(t[good]-tstep)/Ta2)/(deltat+(t[good]-tstep))
#r_theo2[good] = -chi0*depthS*np.exp(-(t[good]-tstep)/(0.25*(Ta1+3*Ta2)))/(Ta1+(t[good]-tstep))
r_theo2[good] = -chi0*depthS*np.exp(-(t[good]-tstep)/Ta)/(Ta+(t[good]-tstep))
#r_theo2[good] = -chi0*depthS*np.exp(-(t[good]-tstep)/(Ta1))/(Ta1+(t[good]-tstep))
rmax = np.amax(r_theo2[good])
print('r_theo2 max =',rmax/rinfty1)
r_theo2[good] = (rinfty2-rinfty1)*(1.-r_theo2[good]/rmax)

# ---------- ramp scenario very rough approximation - not general 
r_theo1_ramp = -chi0*depthS*np.exp(-t/Ta1)/(deltat + t)
r_theo1_ramp = r_theo1_ramp +rinfty1
Ta = 0.25*(Ta1+3*Ta2)
Ta = 0.5*(Ta1+Ta2)
#r_theo2_ramp[good] = -chi0*depthS*np.exp(-(t[good]-tstep)/(Ta1))/(Ta1+(t[good]-tstep))
#r_theo2_ramp[good] = -chi0*depthS*np.exp(-(t[good]-tstep)/(0.5*(Ta1+Ta2)))/(Ta1+(t[good]-tstep))
#r_theo2_ramp[good] = -chi0*depthS*np.exp(-(t[good]-tstep)/(0.25*(Ta1+3*Ta2)))/(Ta1+(t[good]-tstep))
r_theo2_ramp[good] = -chi0*depthS*np.exp(-(t[good]-tstep)/Ta)/(Ta+(t[good]-tstep))
rmax = np.amax(r_theo2_ramp[good])
print('r_theo2 max =',rmax/rinfty1)
r_theo2_ramp[good] = (rinfty2-rinfty1)*(1.-r_theo2_ramp[good]/rmax)

Ta = 0.25*(Ta2+3*Ta3)
Ta = 0.5*(Ta2+Ta3)
#r_theo3_ramp[good3] = -chi0*depthS*np.exp(-(t[good3]-tstep)/(Ta2))/(Ta2+(t[good3]-tstep))
#r_theo3_ramp[good3] = -chi0*depthS*np.exp(-(t[good3]-tstep3)/(0.5*(Ta2+Ta3)))/(Ta2+(t[good3]-tstep3))
#r_theo3_ramp[good3] = -chi0*depthS*np.exp(-(t[good3]-tstep)/(0.25*(Ta2+3*Ta3)))/(Ta2+(t[good3]-tstep))
r_theo3_ramp[good3] = -chi0*depthS*np.exp(-(t[good3]-tstep)/Ta)/(Ta+(t[good3]-tstep))
rmax3 = np.amax(r_theo3_ramp[good3])
print('r_theo3 max =',rmax3/rinfty1)
r_theo3_ramp[good3] = (rinfty3-rinfty2)*(1.-r_theo3_ramp[good3]/rmax3)

# ---------- ramp scenario , new analytical solution a la Scholz, 30 Jan 2022 
dS = -depthS
t0 = deltat/10.
#t1 = deltat
#t1 = t0
#t1 = tstep
t1 = 100*tstep # in der simulation wird vorab ein Gleichgewicht chi bestimmt
Ta1b = Ta1
r_theo1_ramp2 = rinfty1 +chi0*dS*(1.-np.exp(-np.exp(-t/Ta1b)*(t1+t)/t0))/(t1+t)

tstep2b = tstep
#tstep2b = 0.97*tstep
t0 = deltat/10.
t1 = deltat
#t1 = t0
#t1 = 10*deltat
Ta2b = Ta2   # Ta aendern bringt nicht viel, da fuer kleine lag times offensichtlich wenig Einfluss - wird ebenso durch t0 bestimmt 
r_theo2_ramp2[good] = chi0*dS*(1.-np.exp(-np.exp((t[good]-tstep2b)/Ta2b)*(t1+t[good]-tstep2b)/t0))/(t1+(t[good]-tstep2b))
rmax = np.amax(r_theo2_ramp2[good])
print('r_theo2 max =',rmax/rinfty1)
r_theo2_ramp2[good] = (rinfty2-rinfty1)*(1.-r_theo2_ramp2[good]/rmax)

tstep3b = tstep3
t0 = deltat/10.
t1 = deltat
#t1 = t0
#t1 = 0.1*(tstep3b-tstep2b)
Ta3b = Ta3
r_theo3_ramp2[good3] = chi0*dS*(1.-np.exp(-np.exp((t[good3]-tstep3b)/Ta3b)*(t1+t[good3]-tstep3b)/t0))/(t1+(t[good3]-tstep3b))
rmax3 = np.amax(r_theo3_ramp2[good3])
print('r_theo3 max =',rmax3/rinfty1)
r_theo3_ramp2[good3] = (rinfty3-rinfty2)*(1.-r_theo3_ramp2[good3]/rmax3)


# ------ second case -  calculate analytical approximations
Ta12   = -depthS/strend2
Ta22   = -depthS/strend1
rinfty12 = strend2*chi0
rinfty22 = strend1*chi0

Ta = 0.5*(Ta12+Ta22)
r_theo12 = -chi0*depthS*np.exp(-t/Ta12)/(deltat + t)
r_theo12 = r_theo12 +rinfty12
#r_theo22[good] = -chi0*depthS*np.exp(-(t[good]-tstep)/(0.25*(Ta12+3*Ta22)))/(Ta12+(t[good]-tstep))
r_theo22[good] = -chi0*depthS*np.exp(-(t[good]-tstep)/Ta)/(Ta+(t[good]-tstep))
rmax = np.amax(r_theo22[good])
print('r_theo2 max =',rmax/rinfty1)
r_theo22[good] = (rinfty22-rinfty12)*(1.-r_theo22[good]/rmax)

#-----------------------------------------------------
print("Plotting Fig. 4d and 4e Appendix")
#-----------------------------------------------------
tmin = -0.5
tmax = 1.5
tmax = 2.5
xrange_set = False
xrange_set = True
ymaxs = (strend1*tstep+strend2*(tstep3-tstep))/(Ta1*strend1) + 1.5
lstyle = '-'

if xrange_set:
    #nstep = 0
    #nstart = nstep-1
    nstart = 0
else:
    #nstep = 0
    nstart = 0

# ---- Coulomb stress loading (input and interpolated input)
#fig = plt.figure(1, figsize=(12, 8))
fig = plt.figure(1, figsize=(12, 6))
ax1a = fig.add_subplot(131)
ax1a.set_ylabel('$\sigma_c/\dot{\sigma}_{c1} T_{a1}$', fontsize=20)
ax1a.set_xlabel(r'$(t-t_0)/T_{a1}$', fontsize=20)
plt.figtext(0.09, 0.87, 'a)', fontsize=20)
ax1a.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1a.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax1a.plot((t[nstart:]-t[nstep])/Ta1, cfs1[nstart:]/(Ta1*strend1), linewidth=2.0, ls=lstyle, color='black', label=r'$\dot{\sigma}_{c2}/\dot{\sigma}_{c1}=$'+'{:.1f}'.format(float(strend2/strend1)))

plt.legend(loc='upper left',fontsize=20)
if xrange_set:
    #plt.ylim([1.5, 28])
    plt.ylim([1.5, ymaxs])
    plt.xlim([tmin, tmax])

# ----------------- zweites Beispiel stress loading
ax12a = fig.add_subplot(132)
#ax12a.set_ylabel('$\sigma_c/\dot{\sigma}_{c1} T_{a1}$', fontsize=20)
ax12a.set_xlabel(r'$(t-t_0)/T_{a1}$', fontsize=20)
#plt.figtext(0.50, 0.87, 'b)', fontsize=20)
plt.figtext(0.365, 0.87, 'b)', fontsize=20)
ax12a.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax12a.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax12a.plot((t[nstart:]-t[nstep])/Ta1, cfs2[nstart:]/(Ta1*strend1), linewidth=2.0, ls=lstyle, color='black', label=r'$\dot{\sigma}_{c2}/\dot{\sigma}_{c1}=$'+'{:.1f}'.format(float(strend1/strend2)))

plt.legend(loc='lower right',fontsize=20)
if xrange_set:
    #plt.ylim([0, 18])
    plt.ylim([1.5, ymaxs])
    plt.xlim([tmin, tmax])

# ----------------- drittes Beispiel ramp function
ax13a = fig.add_subplot(133)
ax13a.set_xlabel(r'$(t-t_0)/T_{a1}$', fontsize=20)
plt.figtext(0.635, 0.87, 'c)', fontsize=20)
ax13a.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax13a.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax13a.plot((t[nstart:]-t[nstep])/Ta1, cfs1_ramp[nstart:]/(Ta1*strend1), linewidth=2.0, ls=lstyle, color='black')

if xrange_set:
    #plt.ylim([0, 18])
    plt.ylim([1.5, ymaxs])
    plt.xlim([tmin, tmax])

plt.show()

#----------------------------------------------------------
# ---- earthquake rate ------------------------------------
#----------------------------------------------------------
ymax = 1.2*strend2/strend1
#fig = plt.figure(2, figsize=(12, 8))
fig = plt.figure(2, figsize=(12, 6))
ax1b = fig.add_subplot(131)
ax1b.set_xlabel(r'$(t-t_0)/T_{a1}$', fontsize=20)
ax1b.set_ylabel(r'$r / \chi_0 V \dot{\sigma}_{c1}$', fontsize=20)
ax1b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax1b.plot([(t[nstart]-t[nstep])/Ta1, (t[nt-1]-t[nstep])/Ta1], [1, 1], linewidth=1.0, ls='dotted', color='gray')

ax1b.plot((t[nstart:]-t[nstep])/Ta1 , r_tdsm[nstart:]/rinfty1, linewidth=2.0, ls='-', color='red', label=r'TDSM')
ax1b.plot((t[nstart:]-t[nstep])/Ta1 , r_rsm[nstart:]/rinfty1, linewidth=1.0, ls='--', color='blue', label=r'RSM')
ax1b.plot((t[nstart:]-t[nstep])/Ta1 , r_lcm[nstart:]/rinfty1, linewidth=1.0, ls='--', color='gray', label=r'LCM')
#plt.text((t[nstart+2]-t[nstep])/Ta1 , r_tdsm[nstart+2]/rinfty1 ,r'$\Delta\sigma/\dot{\sigma}_c=$'+'{:.0f}'.format(float(Ta1)), fontsize=20)

ax1b.plot((t[nstart:]-t[nstep])/Ta1 , (r_theo1[nstart:]+r_theo2[nstart:])/rinfty1, linewidth=2.5, ls='dotted', color='black')
# ---- theo1_ramp2 only for testing , to be deleted ---
ax1b.plot((t[nstart:]-t[nstep])/Ta1 , r_theo1_ramp2[nstart:]/rinfty1, linewidth=1.5, ls='dashed', color='orange')

plt.legend(loc='lower right',fontsize=20)
plt.figtext(0.09, 0.87, 'a)', fontsize=20)
plt.figtext(0.13, 0.82, r'$T_{a1}/T_{a2}=$'+'{:.1f}'.format(float(Ta1/Ta2)), fontsize=20)
if xrange_set:
    plt.xlim([tmin, tmax])
    plt.ylim([0, ymax])
    #ax1b.set_yticks([2,4,6,8,10])
    #ax1b.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax1b.set_xscale('log')
    #ax1b.set_yscale('log')
else:
    ax1b.plot((t[nstart:]-t[nstep])/Ta1 , r_theo1[nstart:]/rinfty1, linewidth=1.5, ls='--', color='orange')
    ax1b.plot((t[nstart:]-t[nstep])/Ta1 , r_theo2[nstart:]/rinfty1, linewidth=1.5, ls='--', color='green')

# ---- second case in right panel ------------------------------------
ax12b = fig.add_subplot(132)
ax12b.set_xlabel(r'$(t-t_0)/T_{a1}$', fontsize=20)
ax12b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax12b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax12b.plot([(t[nstart]-t[nstep])/Ta1, (t[nt-1]-t[nstep])/Ta12], [1, 1], linewidth=1.0, ls='dotted', color='gray')

ax12b.plot((t[nstart:]-t[nstep])/Ta1 , r_tdsm2[nstart:]/rinfty1, linewidth=2.0, ls='-', color='red', label=r'TDSM')
ax12b.plot((t[nstart:]-t[nstep])/Ta1 , r_rsm2[nstart:]/rinfty1, linewidth=1.0, ls='--', color='blue', label=r'RSM')
ax12b.plot((t[nstart:]-t[nstep])/Ta1 , r_lcm2[nstart:]/rinfty1, linewidth=1.0, ls='--', color='gray', label=r'LCM')

ax12b.plot((t[nstart:]-t[nstep])/Ta1 , (r_theo12[nstart:]+r_theo22[nstart:])/rinfty1, linewidth=2.5, ls='dotted', color='black')
# ---- theo2_ramp2 only for testing , to be deleted ---
ax12b.plot((t[nstart:]-t[nstep])/Ta1 , r_theo2_ramp2[nstart:]/rinfty1, linewidth=1.5, ls='dashed', color='orange')

#plt.legend(loc='upper right',fontsize=20)
#plt.figtext(0.485, 0.87, 'b)', fontsize=20)
#plt.figtext(0.58, 0.82, r'$T_{a1}/T_{a2}=$'+'{:.1f}'.format(float(Ta12/Ta22)), fontsize=20)
plt.figtext(0.365, 0.87, 'b)', fontsize=20)
plt.figtext(0.41, 0.82, r'$T_{a1}/T_{a2}=$'+'{:.1f}'.format(float(Ta12/Ta22)), fontsize=20)
if xrange_set:
    plt.xlim([tmin, tmax])
    plt.ylim([0, ymax])
    #ax12b.set_yticks([2,4,6,8,10])
    #ax12b.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
else:
    ax12b.plot((t[nstart:]-t[nstep])/Ta1 , r_theo12[nstart:]/rinfty1, linewidth=1.5, ls='--', color='orange')
    ax12b.plot((t[nstart:]-t[nstep])/Ta1 , r_theo22[nstart:]/rinfty1, linewidth=1.5, ls='--', color='green')

# ---- third case, ramp function  ------------------------------------
ax13b = fig.add_subplot(133)
ax13b.set_xlabel(r'$(t-t_0)/T_{a1}$', fontsize=20)
ax13b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax13b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax13b.plot([(t[nstart]-t[nstep])/Ta1, (t[nt-1]-t[nstep])/Ta12], [1, 1], linewidth=1.0, ls='dotted', color='gray')

ax13b.plot((t[nstart:]-t[nstep])/Ta1 , r_tdsm_ramp[nstart:]/rinfty1, linewidth=2.0, ls='-', color='red', label=r'TDSM')
ax13b.plot((t[nstart:]-t[nstep])/Ta1 , r_rsm_ramp[nstart:]/rinfty1, linewidth=1.0, ls='--', color='blue', label=r'RSM')
ax13b.plot((t[nstart:]-t[nstep])/Ta1 , r_lcm_ramp[nstart:]/rinfty1, linewidth=1.0, ls='--', color='gray', label=r'LCM')

ax13b.plot((t[nstart:]-t[nstep])/Ta1 , (r_theo1_ramp[nstart:]+r_theo2_ramp[nstart:]+r_theo3_ramp[nstart:])/rinfty1, linewidth=2.5, ls='dotted', color='black')
# ---- ramp2 model ---
ax13b.plot((t[nstart:]-t[nstep])/Ta1 , (r_theo1_ramp2[nstart:]+r_theo2_ramp2[nstart:]+r_theo3_ramp2[nstart:])/rinfty1, linewidth=1.5, ls='dashed', color='orange')

plt.figtext(0.635, 0.87, 'c)', fontsize=20)
if xrange_set:
    plt.xlim([tmin, tmax])
    plt.ylim([0, ymax])
    #ax13b.set_yticks([2,4,6,8,10])
    #ax13b.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
else:
    ax13b.plot((t[nstart:]-t[nstep])/Ta1 , r_theo1_ramp[nstart:]/rinfty1, linewidth=1.5, ls='--', color='orange')
    ax13b.plot((t[nstart:]-t[nstep])/Ta1 , r_theo2_ramp[nstart:]/rinfty1, linewidth=1.5, ls='--', color='green')

plt.show()
fig.savefig(str(pdf_file1)+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file1)+'.png', format='png',bbox_inches='tight')
