#--------------------------
# Plot Fig. 6a (see manuscript on TDSM ,Dahm, 2022) 
#--------------------------
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, BackgroundLoading
import matplotlib.pyplot as plt
import matplotlib.ticker
import math
import numpy as np

def scalsigma(stress):
    #c = 3.*np.log(stress)/2. + 1.9
    c = 2.*np.log(stress)/2. + 3.9
    return c

def gridrange(amin, amax, astep, rounding=math.ceil):
    amin = float(amin)
    amax = float(amax)
    astep = float(astep)
    na = int(rounding((amax-amin)/astep))
    amax = amin + (na-1) * astep
    a = np.linspace(amin, amax, na)
    return amin, amax, na, a

current_dir = Path(__file__).parent
config_file = current_dir / "config.toml"
pdf_file1 = current_dir / "plots/Dahm_fig6a"

print("set-up the tdsm, lcm and rsm class environments")
tdsm = TDSM(config=Config.open(config_file))
lcm  = LCM(config=Config.open(config_file))
trad = Traditional(config=Config.open(config_file))
rsm  = RSM(config=Config.open(config_file))

print("define model parameter for plotting, if different from config file")
hours = 3600.
days  = hours*24.
tstart = 0*days
#tend   = 200*days
tend   = 100*days
#deltat = 0.02*days  # = 72 hours
deltat = 0.01*days  # = 32 hours
nt = np.floor( (tend-tstart) / deltat).astype(int)

strend  = 2.5E-6
#strend  = 2.0E+6/(1500*deltat)  # = 5E-9 Pa/s
#strend  = 2.5E-9
strend  = 2.5E-5
#strend  = 2.5E-4
#sstep  =  [0.0001, 10., 20., 30., 40, 50]  # stress step in MPa
sstep  =  [0.0001, 3., 5., 10., 20., 30., 40, 50]  # stress step in MPa
#tstep  = 50*days         # time when stress step is acting
tstep  = 20*days
n1 = np.floor( (tstep-tstart) / deltat).astype(int)

chi0   = 5.E-3           # susceptibility to trigger earthquakes by unit stress increase
depthS = -1.E1           # skin depth in MPa (must be defined negativ)

deltaS    = -depthS/500.  # increment do discretize Coulomb stress axis
sigma_max = 10000.*deltaS # maximum depth on Coulomb axis (limit of integral)
precision = 18

print('deltat=',deltat/hours,' hours, from 0 to ',tend/hours,' hours, samples =',nt)
print('suscept. chi0=',chi0,' #/MPa')
print('stress step=[',sstep,'] MPa')
print('skin depth =',-depthS,' MPa')
print('tectonic trends =',strend*24.*hours,' MPa/days')
print('deltaS=',deltaS,' MPa, from 0 to ',sigma_max,' MPa')

#-------------------------------------------
#--------------- productivity
# ----  read real data ro comparison
fn_m55 = "../../text/lecture/seismologie/rate_state_seismicity/data_aftershocks/hainzl2008_fig2b_5.5-6.dat"
fn_m60 = "../../text/lecture/seismologie/rate_state_seismicity/data_aftershocks/hainzl2008_fig2b_6-6.5.dat"
fn_m65 = "../../text/lecture/seismologie/rate_state_seismicity/data_aftershocks/hainzl2008_fig2b_6.5-7.dat"
fn_m70 = "../../text/lecture/seismologie/rate_state_seismicity/data_aftershocks/hainzl2008_fig2b_7-7.5.dat"
fn_m75 = "../../text/lecture/seismologie/rate_state_seismicity/data_aftershocks/hainzl2008_fig2b_7.5-9.dat"
ex, ey = np.loadtxt(fn_m55, skiprows=1, usecols=(0,1), delimiter= ',', unpack=True)

t_m55 = 10**(ex)
r_m55 = 10**(ey)
print("m55",t_m55, r_m55)

ex, ey = np.loadtxt(fn_m60, skiprows=1, usecols=(0,1), delimiter= ',', unpack=True)
t_m60 = 10**(ex)
r_m60 = 10**(ey)
print("m60",t_m60, r_m60)

ex, ey = np.loadtxt(fn_m65, skiprows=1, usecols=(0,1), delimiter= ',', unpack=True)
t_m65 = 10**(ex)
r_m65 = 10**(ey)
print("m65",t_m65, r_m65)

ex, ey = np.loadtxt(fn_m70, skiprows=1, usecols=(0,1), delimiter= ',', unpack=True)
t_m70 = 10**(ex)
r_m70 = 10**(ey)
print("m70",t_m70, r_m70)


ex, ey = np.loadtxt(fn_m75, skiprows=1, usecols=(0,1), delimiter= ',', unpack=True)
t_m75 = 10**(ex)
r_m75 = 10**(ey)
print("m75",t_m75, r_m75)

#-----------------------------------------------------
#  simulate earthquake rate as a function of stress step
#-----------------------------------------------------
ns   = len(sstep)
val  = np.zeros(ns)
val1 = np.zeros(ns)
r_tdsm = np.zeros((ns,nt))
r_lcm  = np.zeros((ns,nt))
r_rsm  = np.zeros((ns,nt))
r_theo  = np.zeros((ns,nt))
n_tdsm = np.zeros((ns,nt-1))
n_lcm  = np.zeros((ns,nt-1))
n_rsm  = np.zeros((ns,nt-1))


nstep = np.floor( (tstep-tstart) / deltat).astype(int)+1

# ---- calculate equilibrium distribution of tectonic loading, to be used later to avoid initial oscillations
loading = BackgroundLoading(_config=tdsm.config, strend=strend, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz_background, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)

print('background')
print('len(cf)',len(cf))
print('len(chiz_background)',len(chiz_background))
for i in range(ns):

    loading = StepLoading(_config=tdsm.config, strend=strend, sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
    config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
    #plot(config, t, cf, r, xn)
    r_tdsm[i,:] = r[:]
    n_tdsm[i,:] = xn[:]

#    loading = StepLoading(_config=trad.config, strend=strend, sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
#    config, t, chiz, cf, r, xn = trad(loading=loading, chi0=chi0, deltat=deltat, tstart=tstart, tend=tend)
#    r_lcm[i,:] = r
#    n_lcm[i,:] = xn

#    loading = StepLoading(_config=rsm.config, strend=strend, sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
#    config, t, chiz, cf, r, xn = rsm(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
#    #plot(config, t, cf, r, xn)
#    r_rsm[i,:] = r
#    n_rsm[i,:] = xn

#-----------------------------------------------------
#  calculate calculate analytical approximations
#-----------------------------------------------------
rinfty = strend*chi0
Ta   = -depthS/strend
good = t >= (tstep-tstart)
#f1 = 0.15
f1 = 150
f2 = 0.2
for i in range(ns):
    r_theo[i,good] = (-depthS*np.exp(-(t[good]-tstep)*sstep[i]/(-deltat*depthS)) -depthS*np.exp(-(t[good]-tstep)/Ta) ) / (f1*np.exp(sstep[i]/depthS))
    r_theo[i,good] = +chi0*r_theo[i,good]/(deltat+(t[good]-tstep))
    r_theo[i,:] = r_theo[i,:] +rinfty

#-----------------------------------------------------
print("Plotting Fig. 6a")
#-----------------------------------------------------
Ta   = -depthS/strend
fa   = np.zeros(ns)
fa   = -np.asarray(sstep)/(depthS)
Ta0 = 0.0
scal = strend*chi0
tmin = -10.0
tmax = 120.0
xrange_set = False
xrange_set = True
lstyle = '-'

#if xrange_set:
#    nstart = nstep-0
#    nstart = nstep-1
#else:
#    nstep = 0
#    nstart = 0

# ---- Coulomb stress loading (input and interpolated input)
Ta = 1.0
scal = 1.0
#fig = plt.figure(1, figsize=(12, 8))
#fig , ax2a = plt.subplots(1, figsize=(8, 6))
fig = plt.figure(1, figsize=(8, 6))
ax2a = fig.add_subplot(111)
ax2a.set_xlabel(r'mainshock magnitude $M_W$', fontsize=20)
ax2a.set_ylabel(r'$r(t< 2.4 h)/r_{max}$', fontsize=20)
ax2a.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax2a.tick_params(axis = 'both', which = 'minor', labelsize = 18)

# ------------------ magnitude axis initialisieren und empirische Skalierung plotten  ----
magmin, magmax, nmag, magnitude = gridrange(5.25, 8.0, 0.25)
forecast1  = np.zeros(nmag)
alpha = 0.83
C = 4.E-7
forecast1 = C*10**(alpha*magnitude)
ax2a.plot(magnitude, forecast1, linewidth=1.0, ls='--', color='gray')

axes2 = ax2a.twiny()
#axes2.set_xlabel(r'$\sim \log \Delta \sigma_I^{3/2} + 1.9$', fontsize=20)
axes2.set_xlabel(r'$\sim \log \Delta \sigma_I + 3.9$', fontsize=20)


# ------------------ simulierte Productivity  ----
tlim = 0.1
tlim = 0.2  # intervall, indem aufsummiert wird (verschiebt Absolutwerte nach oben oder unten
tlim = 0.3  
#good = np.logical_and( t < (tlim-t_m75[0])*days +t[n1] , t > -t_m75[0] +t[n1] )
good = np.logical_and( t < (tlim-t_m75[0])*days +t[n1] , t > -t_m75[0]*days +t[n1] )
rmax1 = np.sum(r_tdsm[ns-1,good])
i=ns-1
rmax1 = np.sum(r_tdsm[i,good])
for i in range(1, ns):
    val[i] = np.sum(r_tdsm[i,good])
    val1[i] = np.sum(r_theo[i,good])
    ax2a.plot(scalsigma(sstep[i]), np.sum(r_tdsm[i,good])/rmax1, ls='None', ms=5, marker='o', mec= 'blue', mfc='gray', mew=2.0)

ax2a.plot(scalsigma(sstep[1:ns]), val[1:ns]/rmax1, ls='-', linewidth=1.5, color='blue')
#ax2a.plot(scalsigma(sstep), val/rmax1, ls='-', linewidth=1.5, color='blue')

# --- productivity from analytical solution of rsm model
ax2a.plot(scalsigma(sstep[1:ns]), val1[1:ns]/rmax1, ls='-', linewidth=1.5, color='green')
print("sstep[1:ns]",sstep[1:ns])
print("scal(sstep)",scalsigma(sstep[1:ns]))
print("val1       ",val1[1:ns]/rmax1)
print("val        ",val[1:ns]/rmax1)

# ---- plot real data observations, ISC catalog, after Hainzl and Marsan Fig. 2b
print("t_m55",t_m55.min(),t_m55.max(),r_m55.min(),r_m55.max())
print("t_m60",t_m60.min(),t_m60.max(),r_m60.min(),r_m60.max())
print("t_m65",t_m65.min(),t_m65.max(),r_m65.min(),r_m65.max())
print("t_m70",t_m70.min(),t_m70.max(),r_m70.min(),r_m70.max())
print("t_m75",t_m75.min(),t_m75.max(),r_m75.min(),r_m75.max())
nk = 3
rmax = np.sum(r_m75[0:nk])
ax2a.plot(5.75, np.sum(r_m55[0:nk])/rmax, ls='None', ms=15, marker='o', mec= 'black', mfc='None', mew=2.0)
ax2a.plot(6.25, np.sum(r_m60[0:nk])/rmax, ls='None', ms=15, marker='o', mec= 'black', mfc='None', mew=2.0)
ax2a.plot(6.75, np.sum(r_m65[0:nk])/rmax, ls='None', ms=15, marker='o', mec= 'black', mfc='None', mew=2.0)
ax2a.plot(7.25, np.sum(r_m70[0:nk])/rmax, ls='None', ms=15, marker='o', mec= 'black', mfc='None', mew=2.0)
ax2a.plot(7.75, np.sum(r_m75[0:nk])/rmax, ls='None', ms=15, marker='o', mec= 'black', mfc='None', mew=2.0)

#plt.legend(loc='upper right',fontsize=20)
plt.figtext(0.0, 0.93, 'a)', fontsize=20)
#plt.figtext(0.57, 0.70, r'$T_a=$'+'{:.1f}'.format(float(Ta/(24*hours)))+' days', fontsize=20)

if xrange_set:
    plt.yscale('log')
    ax2a.xaxis.label.set_color('black')
    ax2a.set_xlim(5.0, 8.0)
    axes2.set_xlim(5.0, 8.0)
    axes2.xaxis.label.set_color('blue')
    tkw = dict(size=4, width=1.5)
    axes2.tick_params(axis='x', colors='blue', **tkw)
    axes2.xaxis.set_tick_params(labelsize=18)
else:
    print('no range set')

plt.show()
fig.savefig(str(pdf_file1)+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file1)+'.png', format='png',bbox_inches='tight')
