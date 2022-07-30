#--------------------------
# Example how to read in a Coulomb stress loading function
# seeSupplemenrtary figure in Dahm and Hainzl (2022) for complex CSF evolution
#--------------------------
import sys
sys.path.insert(0, "../")
from pathlib import Path
from tdsm import Config, TDSM, TDSR1, LCM, Traditional, CFM, RSM , RSD, RSD1
from tdsm.plotting import plot
from tdsm.loading import ExternalFileLoading, BackgroundLoading
import matplotlib.pyplot as plt
import matplotlib.ticker
import math
import numpy as np

current_dir = Path(__file__).parent.parent
config_file = current_dir / "config.toml"
pdf_file1 = current_dir / "../plots_tdsr/fig_CFS_readascii"
data_file = current_dir / "../data/CFSloading_synthetic.dat"

print("set-up the tdsm, lcm and rsm class environments")
tdsm = TDSM(config=Config.open(config_file))
tdsr = TDSR1(config=Config.open(config_file))
#lcm  = LCM(config=Config.open(config_file))
#trad = Traditional(config=Config.open(config_file))
#cfm  = CFM(config=Config.open(config_file))
# rsm  = RSM(config=Config.open(config_file))   # funktiooniert  nicht
# rsd  = RSD(config=Config.open(config_file))   # funktiooniert  nicht
rsd1  = RSD1(config=Config.open(config_file))


print("define model parameter for plotting, if different from config file")
hours = 1.0
tstart = 0.*hours
tend   = 10*hours
deltat = 0.001*hours
t0     = 1.0*deltat
nt = np.floor( (tend-tstart) / deltat).astype(int)

strend  = 1.0
chi0   = 1.0             # susceptibility to trigger earthquakes by unit stress increase
depthS = -1.0            # skin depth in MPa (must be defined negativ)

deltaS    = -depthS/500.  # increment do discretize Coulomb stress axis
sigma_max = 10000.*deltaS # maximum depth on Coulomb axis (limit of integral)
#precision = 18

iX0switch = 0           # steady state distribution

iascii = True
scal_cf = 1.0   # data provided in Pa, to be scaled to  MPa
scal_t  = 1.0    # time provilded in units of days, to be scaled to seconds
c_tstart = 0.0

#-----------------------------------------------------
print("Calculate stress loading function and save in directory ",data_file)
#-----------------------------------------------------
#T = 10.0
#NT = 10000
dsig = -depthS
t = np.linspace(tstart, tend, nt)
dS = np.zeros(nt)

N1 = int(0.1 * nt)
N2 = int(0.3 * nt)
N3 = int(0.5 * nt)
N4 = int(0.6 * nt)
N5 = int(0.8 * nt)

dS[:N1] = strend * deltat
dS[N2:N3] = 3 * strend * deltat
dS[N4:] = 0.5 * strend * deltat
dS[N1] += 5 * dsig
dS[N2] -= 5 * dsig
dS[N4] += 2 * dsig
S = np.cumsum(dS)

# External loading file is expected to have two headin lines (see tdsm/loading.py mthod ExternalFileLoading
# first column is Coulomb stress , secind column is time (note that scaling factors can be provided during call)
np.savetxt(data_file, np.transpose([S, t]), header='\n'.join(["Ascii loading file is expected to have two header lines", "1st column: Coulomb stress, 2nd column: time"]))

#-----------------------------------------------------
print("Calculate earthquake rates with tdsm, lcm and rsm")
#-----------------------------------------------------
cfs   = np.zeros(nt)

#r_tdsm = np.zeros(nt)
r_tdsr = np.zeros(nt)
#r_lcm  = np.zeros(nt)
#r_cfm  = np.zeros(nt)
r_rsd  = np.zeros(nt)
#r_rsd  = np.zeros(nt)

loading = ExternalFileLoading(_config=tdsr.config, iascii=True, infile=data_file, scal_t=scal_t, scal_cf=scal_cf, strend=strend, c_tstart=c_tstart, tstart=tstart, tend=tend, deltat=deltat)
config, t, chiz, cf, r_tdsr, xn = tdsr(loading=loading, chi0=chi0, t0=t0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, iX0switch=iX0switch, deltat=deltat, taxis_log=0, tstart=tstart, tend=tend)

loading = ExternalFileLoading(_config=rsd1.config, iascii=True, infile=data_file, scal_t=scal_t, scal_cf=scal_cf, strend=strend, c_tstart=c_tstart, tstart=tstart, tend=tend, deltat=deltat)
config, t, chiz, cf, r_rsd, xn = rsd1(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=tstart, tend=tend)

#-------------------------------------------------
# Plot results 
#-------------------------------------------------
# General Plot Properties
# make a bigger global font size and sans-serif style
plt.rc('font', family='sans-serif')
plt.rc('font', size=12)
plt.rc('legend', fontsize=10)
fig, ax = plt.subplots(2, 1, figsize=(7,6) )
plt.subplots_adjust(hspace=0.0)

ax[0].plot(t, S, c='k')
#ax[0].plot(t, cf, c='grey', lw=3, ls='dotted')
ax[0].set_xlim(tstart, tend)
ax[0].set_ylim(1.1*np.min(S), 1.1*np.max(S))
ax[0].set_ylabel(r'$\sigma_c(t) - \sigma_c(0)$')
ax[0].set_xticks([])
ax[1].plot(t, r_tdsr, c='b', lw=3, alpha=0.4, label='TDSR')
ax[1].plot(t, r_rsd, c='r', ls='dashed', lw=1, label='RS')
ax[1].set_xlim(tstart, tend)
ax[1].set_yscale('log')
#ax[1].set_ylim(0.5*np.min((np.min(r_tdsr), np.min(r_rsd))), 1.5*np.max((np.max(r_tdsr), np.max(r_rsd))))
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Seismicity rate  $R$ / $r_0$')
ax[1].legend()

plt.show()
fig.savefig(str(pdf_file1)+'.pdf', dpi=300, format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file1)+'.png', dpi=300, format='png',bbox_inches='tight')
print('\n\t OUTPUT: %s\n' % (pdf_file1))
