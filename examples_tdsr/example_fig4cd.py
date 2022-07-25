#--------------------------
# Plot Fig. 4c and d (see Dahm and Hainzl, 2022)
#--------------------------
import sys
sys.path.insert(0, "../")
from pathlib import Path

from tdsm import Config, TDSM, TDSR1, LCM, CFM, Traditional, RSM , RSD
from tdsm.plotting import plot
from tdsm.loading import StepLoading, FourPointLoading, BackgroundLoading, CyclicLoading
from tdsm.utils import Eq7

import numpy as np
import math
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent.parent
config_file = current_dir / "config.toml"
figname = current_dir / "../plots_tdsr/fig4cd"

print("set-up the tdsr, linear Coulomb failure and Rate and State models")
#tdsm = TDSM(config=Config.open(config_file))
tdsr = TDSR1(config=Config.open(config_file))
#lcm  = LCM(config=Config.open(config_file))
cfm  = CFM(config=Config.open(config_file))
#trad = Traditional(config=Config.open(config_file))
#rsm  = RSM(config=Config.open(config_file))
rsd  = RSD(config=Config.open(config_file))

#----------------------------------------
print("define model parameter for plotting, if different from config file")
#----------------------------------------
hours = 1.       #time unit
depthS = -1.0    # (=-dsig) skin depth in MPa (must be negativ, will be changed with next release)
t0 = 1.0         # mean failure time in sec for critically stressed source
strend = 1.0     # (=dotsigc) in MPa/timeunit: tectonic stressing rate before loading 
chi0 = 1.0       # susceptibility to trigger critically stressed sources if a unit step increas is applied (r0 = chi0*strend) 

ampsin = 5.0E0            # Amplitude of sinusoidal loading
Tsin = ampsin / strend
Tmax = 3.*Tsin

tstart = 0.0
tend   = Tmax 

deltat = 0.02*hours
t0 = deltat
nt = int(math.ceil((tend - tstart) / deltat))  # =NT = len(t)
taxis_log=0    # linear time axis discretization = default in config.toml
#ntlog=100
#tstartl = 0.01*hours

# ---- discretizing stress axis for integration with
deltaS    = -depthS/60.  # increment do discretize Coulomb stress axis
sigma_max = 3000.*deltaS # maximum depth on Coulomb axis (limit of integral)

#----------------------------------------
# Plot Properties
#----------------------------------------
# make a bigger global font size and sans-serif style
plt.rc('font', family='sans-serif')
plt.rc('font', size=12)
plt.rc('legend', fontsize=10)
fig, ax = plt.subplots(2, 1, figsize=(7,6) )
plt.subplots_adjust(hspace=0.0)

#----------------------------------------
# Calculate earthquake rates for cyclic loading
#----------------------------------------
loading = CyclicLoading(_config=tdsr.config, strend=strend, ampsin=ampsin, Tsin=Tsin, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = tdsr(loading=loading, chi0=chi0, t0=t0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, iX0switch=0, deltat=deltat, taxis_log=0, tstart=0, tend=tend)
cfs    = cf
r_tdsr = r

loading = CyclicLoading(_config=cfm.config, strend=strend, ampsin=ampsin, Tsin=Tsin, deltat=deltat, tstart=0, tend=tend)
config, t, chiz, cf, r, xn = cfm(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
cf_shad = chiz[0:nt]
r_cfm = r

loading = CyclicLoading(_config=rsd.config, strend=strend, ampsin=ampsin, Tsin=Tsin, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz, cf, r, xn = rsd(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
r_rsd = r

#----------------------------------------
# plot results 
#----------------------------------------
ax[0].plot(t, cf_shad, c='k', ls='dotted', lw=1, alpha=0.3)
ax[0].plot(t, cf, c='k')
ax[0].set_xlim(0, Tmax)
ax[0].set_ylim(1.1*np.min(cf), 1.1*np.max(cf))
ax[0].set_ylabel(r'$\sigma_c(t) - \sigma_c(0)$')
ax[0].set_xticks([])
ax[0].text(-0.18,.93, '(c)', fontsize=14, color='k', horizontalalignment='left', transform=ax[0].transAxes)

ax[1].plot(t, r_tdsr, c='b', lw=3, alpha=0.6, label='TDSR')
ax[1].plot(t, r_rsd, c='r', ls='dashed', lw=1, label='RS')
ax[1].plot(t, r_cfm, c='g', ls='dashed', lw=1, label='CF')

ax[1].set_xlim(0, Tmax)
ax[1].set_ylim(0, 1.07*np.max((np.max(r_tdsr), np.max(r_rsd), np.max(r_cfm))))
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Seismicity rate  $R$ / $r_0$')
ax[1].legend()
ax[1].text(-0.18,.93, '(d)', fontsize=14, color='k', horizontalalignment='left', transform=ax[1].transAxes)

plt.show()
fig.savefig(str(figname)+'.pdf', format='pdf',dpi=300,bbox_inches='tight')
fig.savefig(str(figname)+'.png', format='png',dpi=300,bbox_inches='tight')
print('\n\t OUTPUT: %s\n' % (figname))
