#--------------------------
# Plot Fig. 3a (see Dahm and Hainzl, 2022)
#--------------------------
import sys
sys.path.insert(0, "../")
from pathlib import Path

from tdsm import Config, TDSM, TDSR1, LCM, CFM, Traditional, RSM , RSD, RSD1
from tdsm.plotting import plot
from tdsm.loading import StepLoading, FourPointLoading, BackgroundLoading, CyclicLoading
from tdsm.utils import Eq5

import numpy as np
import math
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent.parent
config_file = current_dir / "config.toml"
figname = current_dir / "../plots/fig3a"

print("set-up the tdsr, linear Coulomb failure and Rate and State models")
tdsm = TDSM(config=Config.open(config_file))
tdsr = TDSR1(config=Config.open(config_file))
#lcm  = LCM(config=Config.open(config_file))
cfm  = CFM(config=Config.open(config_file))
#trad = Traditional(config=Config.open(config_file))
#rsm  = RSM(config=Config.open(config_file))
#rsd  = RSD(config=Config.open(config_file))
rsd  = RSD1(config=Config.open(config_file))

#----------------------------------------
print("define model parameter for plotting, if different from config file")
#----------------------------------------
hours = 1.       #time unit
depthS = -5.0    # (=-dsig) skin depth in MPa (must be negativ, will be changed with next release)
t0 = 1.0         # mean failure time in sec for critically stressed source
strend = 1.0     # (=dotsigc) in MPa/timeunit: tectonic stressing rate before loading 
chi0 = 1.0       # susceptibility to trigger critically stressed sources if a unit step increas is applied (r0 = chi0*strend) 

# ---- define Sshadow (Zmin) --------
Sshadow = [0.0, 30.0, 60.0]

tstart = 1.e-4*hours
tend   = 100.*hours

deltat = 0.1*hours     # obsolet if logarithmic time axis
nt = int(math.ceil((tend - tstart) / deltat))  # =NT = len(t) - overwritten if logarithmic time axis
taxis_log=0    # linear time axis discretization = default in config.toml
ntlog=1000
tstartl = 1.E-4 *hours

# ---- discretizing stress axis for integration with
deltaS    = -depthS/60.  # increment do discretize Coulomb stress axis
sigma_max = 3000.*deltaS # maximum depth on Coulomb axis (limit of integral)
#precision = 18   # wird bisher nicht beneotigt in tdsr

iX0switch = 1           # uniform distribution

#----------------------------------------
# Plot Properties
#----------------------------------------
Tmin = 0
Tmax = 70
cb = ['gray', 'b', 'g']
al = [0.3, 0.5, 0.7]
# make a bigger global font size and sans-serif style
plt.rc('font', family='sans-serif')
plt.rc('font', size=12)
plt.rc('legend', fontsize=10)
fig, ax = plt.subplots(1, 1, figsize=(7,5) )

#----------------------------------------
# Calculate earthquake rates for cyclic loading
#----------------------------------------
ns   = len(Sshadow)
#r_tdsm = np.zeros((ns,nt))
r_tdsr = np.zeros((ns,nt))
r_cfm  = np.zeros((ns,nt))
r_rsd  = np.zeros((ns,nt))

for k, Zmin in enumerate(Sshadow):
    loading = BackgroundLoading(_config=tdsr.config, strend=strend, taxis_log=taxis_log, ntlog=ntlog, deltat=deltat, tstart=tstartl, tend=tend)
    config, t, chiz, cf, r, xn = tdsr(loading=loading, chi0=chi0, t0=t0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, iX0switch=iX0switch, Sshadow=Sshadow[k], taxis_log=taxis_log, ntlog=ntlog, deltat=deltat, tstart=tstartl, tend=tend)
    r_tdsr[k,:] = r

#    loading = BackgroundLoading(_config=tdsm.config, strend=strend, taxis_log=taxis_log, ntlog=ntlog, deltat=deltat, tstart=tstartl, tend=tend)
#    config, t, chiz, cf, r, xn = tdsm(loading=loading, chi0=chi0, t0=t0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, iX0switch=iX0switch, Sshadow=Sshadow[k], taxis_log=taxis_log, ntlog=ntlog, deltat=deltat, tstart=tstartl, tend=tend)
#    r_tdsm[k,:] = r

    loading = BackgroundLoading(_config=cfm.config, strend=strend, taxis_log=taxis_log, ntlog=ntlog, deltat=deltat, tstart=tstartl, tend=tend)
    config, t, chiz, cf, r, xn = cfm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, Sshadow=Sshadow[k], taxis_log=taxis_log, ntlog=ntlog, deltat=deltat, tstart=tstartl, tend=tend)
    r_cfm[k,:] = r

    loading = BackgroundLoading(_config=rsd.config, strend=strend, taxis_log=taxis_log, ntlog=ntlog, deltat=deltat, tstart=tstartl, tend=tend)
    config, t, chiz, cf, r, xn = rsd(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, Sshadow=Sshadow[k], taxis_log=taxis_log, ntlog=ntlog, deltat=deltat, tstart=tstartl, tend=tend)
    r_rsd[k,:] = r

#----------------------------------------
# plot results 
#----------------------------------------
r0 = chi0*strend
for k, Zmin in enumerate(Sshadow):
    r_theo = Eq5(t, Sshadow[k], chi0, t0, -depthS, strend)
    ax.plot(t, r_cfm[k,:]/r0, c=cb[k], lw=1, alpha=0.3, ls='dashed', label=r'CF')
    ax.plot(t, r_rsd[k,:], c=cb[k], lw=1, alpha=0.9, ls='dotted', label=r'RS$_\mathrm{subcrit}$')
    #ax.plot(t, r_tdsm[k,:]/r0, c='black', ls='dotted', lw=5, alpha=0.3, zorder=2)
    ax.plot(t, r_tdsr[k,:]/r0, c=cb[k], lw=3, alpha=0.5, zorder=2, label=r'$\zeta_\mathrm{min}=%.0f$' % (Zmin))
    ax.plot(t, r_theo/r0, c=cb[k], lw=1, alpha=1.0, zorder=2, label='Eq.(5)')

ax.set_xlim(Tmin, Tmax)
ax.set_ylim(-0.2, 5.2)
ax.set_xlabel('Time')
ax.set_ylabel(r'Rate  $R \, / \, (\chi_0\, \dot\sigma_{c})$')
ax.legend()

plt.show()
fig.savefig(str(figname)+'.pdf', format='pdf',dpi=300,bbox_inches='tight')
fig.savefig(str(figname)+'.png', format='png',dpi=300,bbox_inches='tight')
print('\n\t OUTPUT: %s\n' % (figname))
