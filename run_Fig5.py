#--------------------------
# Plot Fig. 5 and B1b (see manuscript on TDSM ,Dahm 2022) 
#--------------------------
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, FourPointLoading, BackgroundLoading, CyclicLoading
import matplotlib.pyplot as plt
import math
import numpy as np

current_dir = Path(__file__).parent
config_file = current_dir / "config.toml"
pdf_file1 = current_dir / "plots/Dahm_fig5a"
pdf_file2 = current_dir / "plots/Dahm_figB1b"
pdf_file3 = current_dir / "plots/Dahm_fig5b"

print("set-up the tdsm, lcm and rsm class environments")
tdsm = TDSM(config=Config.open(config_file))
lcm  = LCM(config=Config.open(config_file))
trad = Traditional(config=Config.open(config_file))
rsm  = RSM(config=Config.open(config_file))


print("define model parameter for plotting, if different from config file")
hours = 3600.
tstart = 0*hours
tend   = 50*hours
deltat = 0.2*hours
nt = int(math.ceil((tend - tstart) / deltat))

#strend  = [0.25/(75*deltat), 1.0/(75*deltat), 1.75/(75*deltat)]  # stressing rate in MPa/s
strend = 1.75/(75*deltat) # stressing rate in MPa/s

ampsin = 1.0E0            # Amplitude of sinusoidal loading
Tsin   = 30_000           # Period ofsinusoidal loading

sstep  = [-0.5E0 , -1.0E0] # stress step in MPa
tstep  = 20.0*hours       # time when stress step is acting

chi0   = 1.E4            # susceptibility to trigger earthquakes by unit stress increase
depthS = -0.3            # skin depth in MPa (must be defined negativ)

deltaS    = -depthS/60.  # increment do discretize Coulomb stress axis
sigma_max = 3000.*deltaS # maximum depth on Coulomb axis (limit of integral)



print('deltat=',deltat/hours,' hours, from 0 to ',tend/hours,' hours')
print('suscept. chi0=',chi0,' #/MPa')
print('stress steps=[',sstep[0],sstep[1],'] MPa')
print('skin depth =',-depthS,' MPa')
print('tectonic trend =',strend,' MPa/days')
print('deltaS=',deltaS,' MPa, from 0 to ',sigma_max,' MPa')

#-----------------------------------------------------
print("Calculate earthquake rates  for step loading with tdsm, lcm and rsm")
#-----------------------------------------------------
cfs     = np.zeros((2,nt))
cf_shad = np.zeros((2,nt))
r_tdsm  = np.zeros((2,nt))
r_lcm   = np.zeros((2,nt))
r_rsm   = np.zeros((2,nt))
n_tdsm  = np.zeros((2,nt-1))
n_lcm   = np.zeros((2,nt-1))
n_rsm   = np.zeros((2,nt-1))

# ---- calculate equilibrium distribution of tectonic loading, to be used later to avoid initial oscillations
loading = BackgroundLoading(_config=tdsm.config, strend=strend)
config, t, chiz_background, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, deltat=deltat, tstart=0, tend=tend)

for i in range(2):

    loading = StepLoading(_config=tdsm.config, strend=strend, sstep=sstep[i], tstep=tstep)
    config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, deltat=deltat, tstart=0, tend=tend)
    #plot(config, t, cf, r, xn)
    cfs[i,:]    = cf[:]
    r_tdsm[i,:] = r[:]
    n_tdsm[i,:] = xn[:]

    #loading = StepLoading(_config=lcm.config, strend=strend, sstep=sstep[i], tstep=tstep)
    #config, t, chiz, cf, r, xn = lcm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, deltat=deltat, tstart=0, tend=tend)

    loading = StepLoading(_config=trad.config, strend=strend, sstep=sstep[i], tstep=tstep)
    config, t, chiz, cf, r, xn = trad(loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend)
    #plot(config, t, cf, r, xn)
    r_lcm[i,:] = r
    n_lcm[i,:] = xn
    cf_shad[i,:] = chiz[0:nt]

    loading = StepLoading(_config=rsm.config, strend=strend, sstep=sstep[i], tstep=tstep)
    config, t, chiz, cf, r, xn = rsm(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
    #plot(config, t, cf, r, xn)
    r_rsm[i,:] = r
    n_rsm[i,:] = xn

#-----------------------------------------------------
print("Calculate earthquake rates for cyclic loading with tdsm, lcm and rsm")
#-----------------------------------------------------

loading = CyclicLoading(_config=tdsm.config, strend=strend, ampsin=ampsin, Tsin=Tsin)
config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, deltat=deltat, tstart=0, tend=tend)
cfs_cycle    = cf
r_tdsm_cycle = r
n_tdsm_cycle = xn

loading = CyclicLoading(_config=trad.config, strend=strend, ampsin=ampsin, Tsin=Tsin)
config, t, chiz, cf, r, xn = trad(loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend)
cf_shad_cycle = chiz[0:nt]
r_lcm_cycle = r
n_lcm_cycle = xn

loading = CyclicLoading(_config=rsm.config, strend=strend, ampsin=ampsin, Tsin=Tsin)
config, t, chiz, cf, r, xn = rsm(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
r_rsm_cycle = r
n_rsm_cycle = xn

#-----------------------------------------------------
print("Plotting Fig. 5a")
#-----------------------------------------------------
fa   = np.zeros(len(sstep))
#fa   = (-np.asarray(sstep))/(strend*deltat)
fa   = np.asarray(sstep)/depthS
#Ta   = -depthS/(hours*strend)
#Ta   = -depthS/(deltat*strend)
Ta   = -depthS/strend
Ta0 = 0.0

nstep = int(math.ceil((tstep - tstart) / deltat))
scal = chi0*strend
tmin = -0.5
tmax = +7.5
smin = 3.0
smax = 15.0
tmin1 = -0.0
tmax1 = +17.0
smin1 = 3.0
smax1 = 15.0
xrange_set = True
#xrange_set = False
lstyle = ('--', '-')

fig = plt.figure(1, figsize=(12, 8))
# ---- Coulomb stress loading (input and interpolated input)
ax1a = fig.add_subplot(211)
ax1a.set_ylabel('$\sigma_c/\dot{\sigma}_c T_a$', fontsize=20)
ax1a.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1a.tick_params(axis = 'both', which = 'minor', labelsize = 18)

for i in range(2):
    ax1a.plot((t-t[nstep])/Ta, cf_shad[i,:]/(strend*Ta), linewidth=1.5, ls=lstyle[i], color='lightgray')
    ax1a.plot((t-t[nstep])/Ta, cfs[i,:]/(strend*Ta), linewidth=2.0, ls=lstyle[i], color='black', label=r'$\Delta\sigma/\delta\sigma=$'+'{:d}'.format(int(fa[i])))

plt.legend(loc='lower right',fontsize=20)
if xrange_set:
    plt.ylim([smin, smax])
    plt.xlim([tmin, tmax])

# ---- earthquake rate ------------------------------------
ax1b = fig.add_subplot(212)
ax1b.set_xlabel(r'$(t-t_0)/T_a$', fontsize=20)
ax1b.set_ylabel(r'$r \, / \chi_0 V \dot{\sigma}_c$', fontsize=20)
ax1b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

i=0
ax1b.plot((t-t[nstep])/Ta, r_lcm[i,:]/scal, linewidth=1.5, ls=lstyle[i], color='lightgray')
ax1b.plot((t-t[nstep])/Ta, r_tdsm[i,:]/scal, linewidth=2.0, ls=lstyle[i], color='red')
ax1b.plot((t-t[nstep])/Ta, r_rsm[i,:]/scal, linewidth=1.0, ls=lstyle[1], color='blue')

i=1
ax1b.plot((t-t[nstep])/Ta, r_lcm[i,:]/scal, linewidth=1.5, ls=lstyle[i], color='lightgray', label=r'LCM')
ax1b.plot((t-t[nstep])/Ta, r_tdsm[i,:]/scal, linewidth=2.0, ls=lstyle[i], color='red', label=r'TDSM')
ax1b.plot((t-t[nstep])/Ta, r_rsm[i,:]/scal, linewidth=1.0, ls=lstyle[1], color='blue', label=r'RSM')

plt.legend(loc='lower right',fontsize=20)
plt.figtext(0.02, 0.87, 'a)', fontsize=20)
#plt.figtext(0.15, 0.82, r'$\Delta \sigma/\dot{\sigma}_c T_a=$'+'{:d}'.format(int(depthS/(strend*Ta))), fontsize=20)
#plt.figtext(0.15, 0.82, r'$T_a=$'+'{:d}'.format(int(Ta/(hours)))+' hours', fontsize=20)
plt.figtext(0.15, 0.82, r'$T_a=\delta \sigma/\dot{\sigma}_c$', fontsize=20)
if xrange_set:
    plt.xlim([tmin, tmax])
#ax1b.set_yscale('log')

plt.show()
fig.savefig(str(pdf_file1)+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file1)+'.png', format='png',bbox_inches='tight')

#-----------------------------------------------------
print("Plotting Fig. B1b")
#-----------------------------------------------------
fig = plt.figure(1, figsize=(12, 4))
ax1c = fig.add_subplot(111)
ax1c.set_xlabel(r'$(t-t_0)/T_a$', fontsize=20)
ax1c.set_ylabel(r'$n-n_0$', fontsize=20)
ax1c.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1c.tick_params(axis = 'both', which = 'minor', labelsize = 18)

i=0
ax1c.plot((t[0:-1]-t[nstep])/Ta, n_lcm[i,:]-n_lcm[i,nstep], linewidth=1.5, ls=lstyle[i], color='black')
ax1c.plot((t[0:-1]-t[nstep])/Ta, n_tdsm[i,:]-n_tdsm[i,nstep], linewidth=2.0, ls=lstyle[i], color='red')
#ax1c.plot((t[0:-1]-t[nstep])/Ta, n_rsm[i,:]-n_rsm[i,nstep], linewidth=1.0, ls=lstyle[1], color='blue')

i=1
ax1c.plot((t[0:-1]-t[nstep])/Ta, n_lcm[i,:]-n_lcm[i,nstep], linewidth=1.5, ls=lstyle[i], color='black', label=r'LCM')
ax1c.plot((t[0:-1]-t[nstep])/Ta, n_tdsm[i,:]-n_tdsm[i,nstep], linewidth=2.0, ls=lstyle[i], color='red', label=r'TDSM')
#ax1c.plot((t[0:-1]-t[nstep])/Ta, n_rsm[i,:]-n_rsm[i,nstep], linewidth=1.0, ls=lstyle[1], color='blue', label=r'RSM')

plt.legend(loc='upper left',fontsize=20)
plt.figtext(0.02, 0.87, 'b)', fontsize=20)
if xrange_set:
    plt.xlim([tmin, tmax])
    plt.ylim([-2, 18])
ax1b.set_yscale('log')

plt.show()
fig.savefig(str(pdf_file2)+'.pdf', format='pdf', bbox_inches='tight')
fig.savefig(str(pdf_file2)+'.png', format='png', bbox_inches='tight')

#-----------------------------------------------------
print("Plotting Fig. 5b")
#-----------------------------------------------------
fa   = np.zeros(len(sstep))
fa   = (-np.asarray(sstep))/(strend*deltat)
Ta   = -depthS/strend
Ta0 = 0.0

nstep = int(math.ceil((Tsin/2. - tstart) / deltat))
scal = chi0*strend
lstyle = ('--', '-')

fig = plt.figure(1, figsize=(12, 8))
# ---- Coulomb stress loading (input and interpolated input)
ax1a = fig.add_subplot(211)
ax1a.set_ylabel('$\sigma_c/\dot{\sigma}_c T_a$', fontsize=20)
ax1a.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1a.tick_params(axis = 'both', which = 'minor', labelsize = 18)

i=1
ax1a.plot((t-t[nstep])/Ta, cf_shad_cycle/(strend*Ta), linewidth=1.5, ls=lstyle[i], color='lightgray')
ax1a.plot((t-t[nstep])/Ta, cfs_cycle/(strend*Ta), linewidth=2.0, ls=lstyle[i], color='black', label=r'$\Delta\sigma/\dot{\sigma}_c\Delta t=$'+'{:d}'.format(int(fa[i])))

#plt.legend(loc='lower right',fontsize=20)
if xrange_set:
    plt.xlim([tmin1, tmax1])
    plt.ylim([smin1, smax1])

# ---- earthquake rate ------------------------------------
ax1b = fig.add_subplot(212)
ax1b.set_xlabel(r'$(t-t_0)/T_a$', fontsize=20)
ax1b.set_ylabel(r'$r \, / \chi_0 V \dot{\sigma}_c$', fontsize=20)
ax1b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

i=1
ax1b.plot((t-t[nstep])/Ta, r_lcm_cycle/scal, linewidth=1.5, ls=lstyle[i], color='lightgray')
ax1b.plot((t-t[nstep])/Ta, r_tdsm_cycle/scal, linewidth=2.0, ls=lstyle[i], color='red')
ax1b.plot((t-t[nstep])/Ta, r_rsm_cycle/scal, linewidth=1.0, ls=lstyle[1], color='blue')

plt.legend(loc='lower right',fontsize=20)
plt.figtext(0.02, 0.87, 'b)', fontsize=20)
plt.figtext(0.15, 0.82, r'$T/T_a=$'+'{:d}'.format(int(Tsin/Ta)), fontsize=20)
if xrange_set:
    plt.xlim([tmin1, tmax1])

plt.show()
fig.savefig(str(pdf_file3)+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file3)+'.png', format='png',bbox_inches='tight')
