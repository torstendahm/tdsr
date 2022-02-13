#--------------------------
# Plot Fig. 7 Morsleben (see manuscript on TDSM ,Dahm 2022) 
#--------------------------
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import ExternalFileLoading, BackgroundLoading
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mpl_colors
import matplotlib.cm as mpl_cm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import math
import numpy as np

km = 1000.
hours = 3600.
days = hours * 24.

current_dir = Path(__file__).parent
config_file = current_dir / "config.toml"
pdf_file1 = current_dir / "plots/Dahm_fig7"

infile     = 'data/stresschange_morsleben.dat'           # input stress loading
iascii = True
scal_cf = 1.E-6   # data provided in Pa, to be scaled to  MPa
scal_t = days     # time provilded in units of days, to be scaled to seconds
c_tstart = 0.0
fn_obsrate = 'data/observed_rate_southregion.txt'
V = 2250      # Volume of the study region in m^3

print("set-up the tdsm, lcm and rsm class environments")
tdsm = TDSM(config=Config.open(config_file))
#lcm  = LCM(config=Config.open(config_file))
trad = Traditional(config=Config.open(config_file))
rsm  = RSM(config=Config.open(config_file))

print("define model parameter for plotting, if different from config file")
hours = 3600.
tstart = -500*days
tend   = 440*days
deltat = 0.5*days
nt = int(math.ceil((tend - tstart) / deltat))

chi0  = [1.E6*0.00028*days/12, 1.E6*0.00052*days/12 , 1.E6*0.00065*days/12 ] # unit #/MPa m3
depthS  = [ -0.15 , -0.34 , -0.50 ] 

deltaS    = -np.asarray(depthS)/60.  # increment do discretize Coulomb stress axis
sigma_max = 3000.*deltaS # maximum depth on Coulomb axis (limit of integral)

Sshadow = -1.0      # stress shadow in MPa

print('deltat=',deltat/hours,' hours, from 0 to ',tend/hours,' hours')
print('suscept. chi0=',chi0,' #/MPa')
print('skin depth =',-np.asarray(depthS),' MPa')
print('deltaS=',deltaS,' MPa, from 0 to ',sigma_max,' MPa')

#-----------------------------------------------------
print("Calculate earthquake rates for read in file loading with tdsm, lcm and rsm")
#-----------------------------------------------------
r_tdsm  = np.zeros((3,nt))
r_lcm   = np.zeros((3,nt))
r_rsm   = np.zeros((3,nt))
n_tdsm  = np.zeros((3,nt-1))
n_lcm   = np.zeros((3,nt-1))
n_rsm   = np.zeros((3,nt-1))

# ---- External File Floading - calculate tdsm, rsm and traditional LCM
for i in range(3):
    loading = ExternalFileLoading(_config=tdsm.config, iascii=True, infile=infile, scal_t=scal_t, scal_cf=scal_cf, c_tstart=c_tstart, tstart=tstart, tend=tend, deltat=deltat)
    config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=False, chi0=chi0[i], depthS=depthS[i], deltaS=deltaS[i], sigma_max=sigma_max[i], deltat=deltat, tstart=tstart, tend=tend)
    print('tstart, tend, deltat=',tstart,tend,deltat)
    print('cf=',len(cf))
    print('t =',len(t),'  nt=',nt)
    #plot(config, t, cf, r, xn)
    r_tdsm[i,:] = r[:]
    n_tdsm[i,:] = xn[:]

    loading = ExternalFileLoading(_config=trad.config, iascii=True, infile=infile, scal_t=scal_t, scal_cf=scal_cf, c_tstart=c_tstart, tstart=tstart, tend=tend, deltat=deltat)
    config, t, cf_shad, cf, r, xn = trad(loading=loading, chi0=chi0[i], deltat=deltat, tstart=tstart, tend=tend)
    #plot(config, t, cf, r, xn)
    r_lcm[i,:] = r
    n_lcm[i,:] = xn

    #loading = ExternalFileLoading(_config=rsm.config, iascii=True, infile=infile, scal_t=scal_t, scal_cf=scal_cf, c_tstart=c_tstart, tstart=tstart, tend=tend, deltat=deltat)
    #config, t, chiz, cf, r, xn = rsm(loading=loading, chi0=chi0[i], depthS=depthS, deltat=deltat, tstart=tstart, tend=tend)
    #plot(config, t, cf, r, xn)
    #r_rsm[i,:] = r
    #n_rsm[i,:] = xn

#-----------------------------------------------------
print("Plotting Fig. 7")
#-----------------------------------------------------
n1 = 0
Ta = 1.0*days   # skalierung ist bereits auf Tage
scal = 1.*hours/deltat  # die beobachtete Rate wird in events / hour geplottet
# scal1 = 1.E6  # jetzt bereits auf MPa umgerechnet
scal1 = 1.0


xrange_set = True
#xrange_set = False
lstyle = ('-', '--', '--')
col    = ('red', 'orange', 'blue')
if xrange_set:
    n0 = 0
    t0 = t[n1]
    t0 = 0.0

    tmin = 0
    tmax = 100
    cmin = -10
    cmax = 200
else:
    n0 = 0
    t0 = 0.0

fig = plt.figure(1, figsize=(16, 16))

# ---- Coulomb stress loading (input and interpolated input
ax1a = fig.add_subplot(211)
mpl.rc('legend', fontsize=12)
mpl.rc('axes', labelsize=20)
ax1a.set_ylabel('$\sigma_c$ (MPa)', fontsize=20)
ax1a.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1a.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax1a.plot(t/Ta, cf_shad/scal1, linewidth=2.0, ls='--',color='darkgray')
ax1a.plot(t/Ta, cf/scal1, linewidth=2.0, ls='-',color='black', label="Coulomb stress")
print("time t=",t[0]/Ta,t[-1]/Ta,deltat/Ta)
ax1a.legend(loc='upper left')

if xrange_set:
    plt.xlim([tmin, tmax])
    plt.ylim([-0.2, 6.3])

# ---- observed and predicted earthquake rate ----
ax1b = fig.add_subplot(212)
ax1b.set_xlabel("$t$ (days)", fontsize=20)
ax1b.set_ylabel(r'$r(M_W>-5)\, (1/hours)$', fontsize=20)
ax1b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

# stress shadow and acquisition stops indicated ----
# ---------------- Anzahl der Beben als Funktion der Zeit
day1 = 990    # day when the stress shadow is exceeded at MH influx point
day1 = 1000   # day when the stress shadow is exceeded at MH influx point
gap1a = 1064.27
gap1b = 1072.03
gap2a = 1084.28
gap2b = 1098.98
gap3a = 1224.2
#ax1b.plot([day1, day1], [-50 , 150], linewidth=0.5, ls='--',color='gray')
#ax1b.plot([gap1a, gap1a], [-50 , 150], linewidth=0.5, ls='--',color='gray')
#ax1b.plot([gap1b, gap1b], [-50 , 150], linewidth=0.5, ls='--',color='gray')
#ax1b.plot([gap2a, gap2a], [-50 , 150], linewidth=0.5, ls='--',color='gray')
#ax1b.plot([gap2b, gap2b], [-50 , 150], linewidth=0.5, ls='--',color='gray')
#ax1b.plot([gap3a, gap3a], [-50 , 150], linewidth=0.5, ls='--',color='gray')

# --- observed events, compiled by Sebastian
rate_obs, t_obs = np.loadtxt(fn_obsrate, skiprows=2, usecols=(0,1), unpack=True)
good = t_obs <= gap3a
ax1b.plot(t_obs[good], rate_obs[good], linewidth=2.5, ls='-', color='black', label='observed')

# ----- new model with memory
for i in range(3):

    ax1b.plot((t[n0:nt-1]-t0)/Ta, r_tdsm[i,n0:nt-1]/scal, linewidth=1.7, ls=lstyle[i], color=col[i], label=r'$\delta\sigma=$%g $MPa,\chi_0=$%3.2g $kPa^{-1} \,m^{-3}$'% (-depthS[i],chi0[i]/(1.E3*V) ))

if xrange_set:
    plt.ylim([cmin, cmax])
    plt.xlim([tmin, tmax])

ax1b.legend(loc='upper left', fontsize=18)
plt.figtext(0.08, 0.88, 'c)', fontsize=20)

plt.show()
fig.savefig(str(pdf_file1)+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file1)+'.png', format='png',bbox_inches='tight')

# ------------- summary of source parameters --------------
for i in range(3):
    print(' ')
    print(' case ',i,' in ',col[i])
    chi0prime = chi0[i]*1000./ V
    Vchar = 1./(1000*chi0prime)  # in Vchar (in m3) entsteht 1 event bei Erhoehung um 1 MPa (=1000 kPa)
    Vchar = Vchar/1.E9           # in km3
    print('chi0 =',chi0[i]*1000.,' events / kPa,   chi0_prime=',chi0prime,' events / (kPa m^3)')
    print('chi0_prime=',chi0prime*1.E12,' events / (MPa km^3)')
    print('characterstic volume=',Vchar,' km^3')
    print('characterstic length=',np.power(Vchar,1./3.),' km')
    print('delta sigma=',-depthS[i]/(1.E3),' kPa')
    print('delta sigma=',-depthS[i]/(1.E6),' MPa')

# ----- save observed and predicted data ----
#np.savetxt(fn_out1,np.transpose([(t-t[n1])/Ta, cf/scal1]),fmt=['%g','%g'],delimiter=" ",header="t(days), sigma(MPa), V=%gm3 , chi0=%g , depthS=%gMPa" %(V,chi0,depthS/1.E6), comments="" )
