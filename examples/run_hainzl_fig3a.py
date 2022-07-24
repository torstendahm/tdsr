#--------------------------
# Fig. 3a  Dahm and Hainzl 2022
import sys
sys.path.insert(0, "../")
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, FourPointLoading, BackgroundLoading
import matplotlib.pyplot as plt
import math
import numpy as np

current_dir = Path(__file__).parent.parent
config_file = current_dir / "config.toml"
pdf_file1 = current_dir / "../plots/hainzl_fig3a"

t0 = 1.0
dsig = 5.0
X0 = 1.0
rotsigc = 1.0
r0 = X0 * dotsigc

Zminvalues = [0.0, 30.0, 60.0]

t = np.linspace(1e-4, 100, 1000)
S = dotsigc * t
Z = functions.Zvalues(S, 0, t0, dsig)

print("set-up the tdsr class environments")
tdsm = TDSM(config=Config.open(config_file))
tdsr = TDSR(config=Config.open(config_file))
lcm  = LCM(config=Config.open(config_file))
trad = Traditional(config=Config.open(config_file))
rsm  = RSM(config=Config.open(config_file))

print("define model parameter for plotting, if different from config file")
hours = 3600.

tstart = 0*hours
tend   = 200*hours
deltat = 0.1*hours

delval = 0.2*hours   # Parameter only to define other parameter, can be removed or replaced by numerical value
strend  = [0.25/(75*delval), 1.0/(75*delval), 1.75/(75*delval)]  # stressing rate in MPa/s

sstep  = 1.5E0           # stress step in MPa
tstep  = 5.0*hours       # time when stress step is acting
#tstep  = 50.0*hours     

chi0   = 1.E4            # susceptibility to trigger earthquakes by unit stress increase
depthS = -0.3            # skin depth in MPa (must be defined negativ)

deltaS    = -depthS/60.  # increment do discretize Coulomb stress axis
sigma_max = 3000.*deltaS # maximum depth on Coulomb axis (limit of integral)

print('deltat=',deltat/hours,' hours, from 0 to ',tend/hours,' hours')
print('suscept. chi0=',chi0,' #/MPa')
nstep = int(math.ceil((tstep - tstart) / deltat))-1
print('stress step=',sstep,' MPa, tstep=',tstep/hours, (nstep+0)*deltat/hours,nstep)
print('skin depth =',-depthS,' MPa')
print('tectonic trends = [',strend[0]*24.*hours,strend[1]*24.*hours, strend[2]*24.*hours,'] MPa/days')
print('deltaS=',deltaS,' MPa, from 0 to ',sigma_max,' MPa')

#-----------------------------------------------------
print("Calculate earthquake rates with tdsm, lcm and rsm")
#-----------------------------------------------------
cfs    = np.zeros((3,nt))
r_tdsm = np.zeros((3,nt))
r_theo = np.zeros((3,nt))
r_lcm  = np.zeros((3,nt))
r_rsm  = np.zeros((3,nt))
n_tdsm = np.zeros((3,nt-1))
n_lcm  = np.zeros((3,nt-1))
n_rsm  = np.zeros((3,nt-1))

for i in range(3):

    # ---- calculate equilibrium distribution of tectonic loading, to be used later to avoid initial oscillations
    loading = BackgroundLoading(_config=tdsm.config, strend=strend[i],deltat=deltat, tstart=tstart, tend=tend)
    config, t, chiz_background, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, deltat=deltat, tstart=0, tend=tend)

    loading = StepLoading(_config=tdsm.config, strend=strend[i], sstep=sstep, tstep=tstep, deltat=deltat)
    config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, deltat=deltat, tstart=0, tend=tend)
    #plot(config, t, cf, r, xn)
    cfs[i,:]    = cf[:]
    r_tdsm[i,:] = r[:]
    n_tdsm[i,:] = xn[:]
    #print('t[nstep]=',t[nstep]/hours,' h, tstep=',tstep/hours, (nstep+0)*deltat/hours,nstep)

    #loading = StepLoading(_config=lcm.config, strend=strend[1], sstep=sstep, tstep=tstep)
    #config, t, chiz, cf, r, xn = lcm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, deltat=deltat, tstart=0, tend=tend)

    loading = StepLoading(_config=trad.config, strend=strend[i], sstep=sstep, tstep=tstep, deltat=deltat)
    config, t, cf_shadow, cf, r, xn = trad(loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend)
    #plot(config, t, cf, r, xn)
    r_lcm[i,:] = r
    n_lcm[i,:] = xn

    loading = StepLoading(_config=rsm.config, strend=strend[i], sstep=sstep, tstep=tstep, deltat=deltat)
    config, t, cf_shadow, cf, r, xn = rsm(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
    #plot(config, t, cf, r, xn)
    r_rsm[i,:] = r
    n_rsm[i,:] = xn

    # calculate theoretical Omori curves
    good = t >= (tstep-tstart)
    dS = -depthS
    A = chi0*(dS)
    t0 = 0.1*delval  # t0 muss deutlich kleiner als deltat sein, damit die peak Amplituden richtig
    strend_background = strend[i]
    #t1 = deltat
    t1 = dS/(sstep/delval-strend_background)
    rback = chi0*dS*(1.-np.exp(-strend_background*t0/dS))/t0
    r_theo[i,good] = A/(t1+(t[good]-tstep))
    r_theo[i,good] = r_theo[i,good]*(1.-np.exp(-np.exp(-strend_background*(t[good]-tstep)/dS)*(t[good]-tstep+t1)/t0) )
    r_theo[i,:] = r_theo[i,:]+rback

#-----------------------------------------------------
print("Plotting Fig. 3b")
#-----------------------------------------------------
coltheo = 'black'
colnum  = 'orange'
thicktheo = 1.5
thicknum = 5.0

Ta   = -depthS/(np.asarray(strend))
Ta0 = 0.0
nstep = int(math.ceil((tstep - tstart) / deltat))-1
scal = chi0*strend[1]
tmin = -0.5
tmax = 5.5
xrange_set = False
xrange_set = True
lstyle = ('--', '-', '-.')
plot_theo = False
#plot_theo = True

fig = plt.figure(1, figsize=(12, 8))
# ---- Coulomb stress loading (input and interpolated input)
ax1a = fig.add_subplot(211)
ax1a.set_ylabel('$\sigma_c/\Delta\sigma$', fontsize=20)
ax1a.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1a.tick_params(axis = 'both', which = 'minor', labelsize = 18)

for i in range(3):
    ax1a.plot((t-t[nstep]-deltat)/Ta[1], cfs[i,:]/(sstep), linewidth=2.0, ls=lstyle[i], color='black', label=r'$T_a/T_{a1}=$'+'{:.2f}'.format(float(Ta[i]/Ta[1])))

plt.legend(loc='lower right',fontsize=20)
if xrange_set:
    plt.xlim([tmin, tmax])
    plt.ylim([-1, 3.3])

# ---- earthquake rate ------------------------------------
ax1b = fig.add_subplot(212)
ax1b.set_xlabel(r'$(t-t_s)/T_{a1}$', fontsize=20)
ax1b.set_ylabel(r'$r / \chi_0 V \dot{\sigma}_{c1}$', fontsize=20)
ax1b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

for i in [0,2]:
    ax1b.plot((t-t[nstep]-deltat)/Ta[1] , r_lcm[i,:]/scal, linewidth=1.5, ls=lstyle[i], color='lightgray')
    ax1b.plot((t-t[nstep]-deltat)/Ta[1] , r_tdsm[i,:]/scal, linewidth=thicknum, ls=lstyle[i], color=colnum)
    ax1b.plot((t-t[nstep]-deltat)/Ta[1] , r_rsm[i,:]/scal, linewidth=1.0, ls=lstyle[1], color='blue')
    if plot_theo:
        ax1b.plot((t-t[nstep])/Ta[1] , r_theo[i,:]/scal, linewidth=thicktheo, ls=lstyle[0], color=coltheo)

i=1
ax1b.plot((t-t[nstep]-deltat)/Ta[1] , r_lcm[i,:]/scal, linewidth=1.5, ls=lstyle[i], color='lightgray', label=r'LCM')
ax1b.plot((t-t[nstep]-deltat)/Ta[1] , r_tdsm[i,:]/scal, linewidth=thicknum, ls=lstyle[i], color=colnum, label=r'TDSM')
ax1b.plot((t-t[nstep]-deltat)/Ta[1] , r_rsm[i,:]/scal, linewidth=1.0, ls=lstyle[1], color='blue', label=r'RSM')
if plot_theo:
    ax1b.plot((t-t[nstep])/Ta[1] , r_theo[i,:]/scal, linewidth=thicktheo, ls=lstyle[0], color=coltheo)

plt.legend(loc='upper right',fontsize=20)
plt.figtext(0.06, 0.87, 'b)', fontsize=20)
plt.figtext(0.15, 0.82, r'$\Delta\sigma/\delta\sigma=$'+'{:d}'.format(int(-sstep/depthS)), fontsize=20)
if xrange_set:
    plt.xlim([tmin, tmax])
ax1b.set_yscale('log')

plt.show()
fig.savefig(str(pdf_file1)+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file1)+'.png', format='png',bbox_inches='tight')

#-----------------------------------------------------
print("Plotting Fig. B1a")
#-----------------------------------------------------
#fig = plt.figure(1, figsize=(12, 4))
fig = plt.figure(1, figsize=(12, 5.5))
ax1c = fig.add_subplot(111)
ax1c.set_xlabel(r'$(t-t_s)/T_{a1}$', fontsize=20)
ax1c.set_ylabel(r'$n-n_0$', fontsize=20)
ax1c.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1c.tick_params(axis = 'both', which = 'minor', labelsize = 18)

for i in [0,2]:
    ax1c.plot((t[0:-1]-t[nstep])/Ta[1] , n_tdsm[i,:]-n_tdsm[i,nstep], linewidth=thicknum, ls=lstyle[i], color=colnum)
    #ax1c.plot((t[0:-1]-t[nstep])/Ta[1] , n_rsm[i,:]-n_rsm[i,nstep], linewidth=1.0, ls=lstyle[1], color='blue')
    ax1c.plot((t[0:-1]-t[nstep])/Ta[1] , n_lcm[i,:]-n_lcm[i,nstep], linewidth=1.5, ls=lstyle[i], color='black')

i=1
ax1c.plot((t[0:-1]-t[nstep])/Ta[1] , n_rsm[i,:]-n_rsm[i,nstep], linewidth=1.0, ls=lstyle[1], color='blue', label=r'RSM ($T_a/T_{a1}=$'+'{:.0f})'.format(float(Ta[i]/Ta[1])))
ax1c.plot((t[0:-1]-t[nstep])/Ta[1] , n_tdsm[i,:]-n_tdsm[i,nstep], linewidth=thicknum, ls=lstyle[i], color=colnum, label=r'TDSM')
ax1c.plot((t[0:-1]-t[nstep])/Ta[1] , n_lcm[i,:]-n_lcm[i,nstep], linewidth=1.5, ls=lstyle[i], color='black', label=r'LCM')

plt.legend(loc='upper left',fontsize=20)
plt.figtext(0.02, 0.87, 'a)', fontsize=20)
if xrange_set:
    plt.xlim([tmin, tmax])
    plt.ylim([-10, 60])
    plt.ylim([-10, 120])
#ax1c.set_yscale('log')

plt.show()
fig.savefig(str(pdf_file2)+'.pdf', format='pdf', bbox_inches='tight')
fig.savefig(str(pdf_file2)+'.png', format='png', bbox_inches='tight')
