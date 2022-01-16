#--------------------------
# Plot Fig. 4 and B1a (see manuscript on TDSM ,Dahm, 2022) 
#--------------------------
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, BackgroundLoading
import matplotlib.pyplot as plt
import matplotlib.ticker
import math
import numpy as np

current_dir = Path(__file__).parent
config_file = current_dir / "config.toml"
pdf_file0 = current_dir / "plots/Dahm_fig4a"
pdf_file1 = current_dir / "plots/Dahm_fig4c"

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

#strend  = 3.0E-6
strend  = 2.5E-6
#sstep  =  [0.05, 0.55, 1.05, 1.55, 2.05, 2.55, 3.05, 3.55]  # stress step in MPa
sstep  =  [0.05, 0.45, 0.85, 1.25, 1.65, 2.05, 2.45, 2.85]  # stress step in MPa
sstep  =  [0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 3.0]  # stress step in MPa
#sstep  =  [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75]  # stress step in MPa
tstep  = 5.0*hours         # time when stress step is acting
#istep = np.floor( (tstep-tstart) / deltat).astype(int)
nstep = np.floor( (tstep-tstart) / deltat).astype(int)+1


chi0   = 1.E4            # susceptibility to trigger earthquakes by unit stress increase
depthS = -0.3            # skin depth in MPa (must be defined negativ)

#deltaS    = -depthS/60.  # increment do discretize Coulomb stress axis
#sigma_max = 3000.*deltaS # maximum depth on Coulomb axis (limit of integral)
deltaS    = -depthS/500.  # increment do discretize Coulomb stress axis
sigma_max = 10000.*deltaS # maximum depth on Coulomb axis (limit of integral)
precision = 18

print('deltat=',deltat,' sec = ',deltat/hours,' hours, from 0 to ',tend/hours,' hours, samples =',nt)
print('Ta   =',-depthS/strend,' sec = ',-depthS/(strend*hours),' hours')
print('suscept. chi0=',chi0,' #/MPa')
print('stress step=[',sstep,'] MPa')
print('skin depth =',-depthS,' MPa')
print('tectonic trends =',strend*24.*hours,' MPa/days,  dot(sigma)*deltat/deltaS=',-strend*deltat/depthS)
print('deltaS=',deltaS,' MPa, from 0 to ',sigma_max,' MPa')

#-----------------------------------------------------
print("Calculate earthquake rates with tdsm, lcm and rsm")
#-----------------------------------------------------
ns   = len(sstep)
cfs    = np.zeros((ns,nt))
r_tdsm_eq8 = np.zeros(nt)
K1         = np.zeros(nt)
romo       = np.zeros(nt)
r_tdsm = np.zeros((ns,nt))
r_lcm  = np.zeros((ns,nt))
r_rsm  = np.zeros((ns,nt))
n_tdsm = np.zeros((ns,nt-1))
n_lcm  = np.zeros((ns,nt-1))
n_rsm  = np.zeros((ns,nt-1))

r_tdsm_theo = np.zeros((ns,nt))
r_rsm_theo  = np.zeros((ns,nt))

# ---- initial dfecay according to equation 8 for zero tectonic stress rate
loading = BackgroundLoading(_config=tdsm.config, strend=0.0, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz_background, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
r_tdsm_eq8 = r

# ---- calculate equilibrium distribution of tectonic loading, to be used later to avoid initial oscillations
# loading = BackgroundLoading(_config=tdsm.config, strend=strend)
loading = BackgroundLoading(_config=tdsm.config, strend=strend, deltat=deltat, tstart=tstart, tend=tend)
config, t, chiz_background, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)

print('background')
print('len(cf)',len(cf))
print('len(chiz_background)',len(chiz_background))
rinfty = strend*chi0
Ta   = -depthS/strend
good = t >= (tstep-tstart)
for i in range(ns):

    loading = StepLoading(_config=tdsm.config, strend=strend, sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
    config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
    #plot(config, t, cf, r, xn)
    cfs[i,:]    = cf[:]
    r_tdsm[i,:] = r[:]
    n_tdsm[i,:] = xn[:]

    loading = StepLoading(_config=trad.config, strend=strend, sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
    config, t, chiz, cf, r, xn = trad(loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend)
    #plot(config, t, cf, r, xn)
    r_lcm[i,:] = r
    n_lcm[i,:] = xn

    loading = StepLoading(_config=rsm.config, strend=strend, sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
    config, t, chiz, cf, r, xn = rsm(loading=loading, chi0=chi0, depthS=depthS, deltat=deltat, tstart=0, tend=tend)
    #plot(config, t, cf, r, xn)
    r_rsm[i,:] = r
    n_rsm[i,:] = xn

    # calculate analytical approximations
    f1 = 0.15
    f1 = 150
    #f1 = 0.0
    f2 = 0.2
    K = -f2*depthS+f1*sstep[i]
    K = -depthS*np.exp(-(t[good]-tstep)/Ta)/(1.-np.exp(sstep[i]/depthS))
    K = -depthS*deltat/(Ta*(1.+(np.exp(sstep[i]/depthS)-1.)*1.))
    K1[good] = -depthS*np.exp(-(t[good]-tstep)*sstep[i]/(-deltat*depthS)) -depthS*np.exp(-(t[good]-tstep)/Ta)
    #K1[good] = -sstep[i]*np.exp(-(t[good]-tstep)*sstep[i]/(-deltat*depthS))/depthS -depthS*np.exp(-(t[good]-tstep)/Ta)/sstep[i]
    #K1[good] = (-depthS*np.exp(-(t[good]-tstep)*sstep[i]/(-deltat*depthS)))/(-depthS*np.exp(-(t[good]-tstep)/Ta))
    K1[good] = (-depthS*np.exp(-(t[good]-tstep)/Ta)) / (f1*np.exp(sstep[i]/depthS))
    K1[good] = (-depthS*np.exp(-(t[good]-tstep)*sstep[i]/(-deltat*depthS)) -depthS*np.exp(-(t[good]-tstep)/Ta) ) / (f1*np.exp(sstep[i]/depthS))
    #K1[good] = (-depthS*np.exp(-(t[good]-tstep)*sstep[i]/(-deltat*depthS)))
    #K1[nstep-1:nstep+2] = f1*sstep[i] +K1[nstep-1:nstep+2]
    #K = -depthS/( 1.-np.exp(+sstep[i]/depthS) )
    #K = -depthS*np.exp(+sstep[i]/depthS)
    #r_tdsm_theo[i,good] = +chi0*K*np.exp(-(t[good]-tstep)/Ta)/(deltat+(t[good]-tstep))
    r_tdsm_theo[i,good] = +chi0*K1[good]/(deltat+(t[good]-tstep))
    #r_tdsm_theo[i,good] = chi0*(-depthS*np.exp(-(t[good]-tstep)/Ta)/(deltat+(t[good]-tstep)) +chi0*fact_adhoc*sstep[i])*np.exp(-(t[good]-tstep)/Ta)/(deltat+(t[good]-tstep))
    #r_tdsm_theo[i,good] = +chi0*(-depthS+sstep[i])*np.exp(-(t[good]-tstep)/Ta)
    f3 = 0.01
    #r_tdsm_theo[i,nstep+1:nstep+2] = r_tdsm_theo[i,nstep+1:nstep+2] +chi0*f3*sstep[i]/deltat
    r_tdsm_theo[i,:] = r_tdsm_theo[i,:] +rinfty
    print('i=',i,' sstep=',sstep[i],'K=',K1[nstep-2:nstep+4],' theo minmax=',np.amin(r_tdsm_theo[i,nstep-2:nstep+4])/rinfty, np.amax(r_tdsm_theo[i,nstep-2:nstep+4])/rinfty)

    r_rsm_theo[i,good] = 1./(1.+(np.exp(sstep[i]/depthS)-1.)*np.exp(-(t[good]-tstep)/Ta) )
    r_rsm_theo[i,:] = rinfty*r_rsm_theo[i,:]
    #print('i=',i,' sstep=',sstep[i],' rsm  minmax=',np.amin(r_rsm_theo[i,good])/rinfty, np.amax(r_rsm_theo[i,good])/rinfty)
    print('i=',i,' sstep=',sstep[i],' rsm  minmax=',np.amin(r_rsm_theo[i,nstep-2:nstep+4])/rinfty, np.amax(r_rsm_theo[i,nstep-2:nstep+4])/rinfty)

#----------------------------
# theoretical spontaneous decay of fully filled house
#----------------------------
A = chi0*(-depthS)
t1 = deltat
romo[good] = A/(t1+(t[good]-t[nstep]))
#romo = A/(t1+t)
#for i in range(nstep,nt):
#    #romo[i] = A/(t1+((t[i]-t[nstep])-deltat))
#    romo[i] = A/(t1+(t[i]-t[nstep]))

#-----------------------------------------------------
print("Plotting Fig. 4b")
#-----------------------------------------------------
Ta   = -depthS/strend
scal1 = -chi0*depthS/deltat
tmin = -10.0
tmax = 120.0
xrange_set = False
xrange_set = True
lstyle = '-'
if xrange_set:
    nstart = nstep-1
else:
    nstep = 0
    nstart = 0

fig = plt.figure(2, figsize=(6, 8))
ax = fig.add_subplot(111)
#ax.set_xlabel(r'$(t-t_0)/T_a$', fontsize=20)
#ax.set_xlabel(r'$(t-t_0)/T_a$', fontsize=20)
ax.set_xlabel(r'$t/\Delta t$', fontsize=20)
ax.set_ylabel(r'$r \Delta t / \chi_0 V$', fontsize=20)
ax.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax.plot([(t[0])/deltat,(t[-1])/deltat], [0, 0], linewidth=1.0, ls='-', color='black')
ax.plot((t[nstart:]-t[nstep])/deltat, r_tdsm_eq8[0:-nstart]/scal1, linewidth=2.5, ls='-', color='red', label=r'numerical')
ax.plot((t[nstart:]-t[nstep])/deltat, romo[nstart:]/scal1, linewidth=3.5, color='black', ls='--', label=r'equation 8')

plt.figtext(0.0, 0.87, 'a)', fontsize=20)
#plt.figtext(0.57, 0.83, r'$T_a=$'+'{:.1f}'.format(float(Ta/(24*hours)))+' days', fontsize=20)
#plt.figtext(0.57, 0.83, r'$\delta\sigma=$'+'{:.1f}'.format(float(-depthS))+' MPa', fontsize=20)
plt.figtext(0.57, 0.68, r'$\delta\sigma=$'+'{:.1f}'.format(float(-depthS))+' MPa', fontsize=20)
plt.legend(loc='upper right',fontsize=20)
if xrange_set:
    #plt.xlim([tmin, tmax])
    #plt.ylim([0.9, 200])
    #ax1b.set_yticks([1,10,50,100,200])
    plt.xlim([-5, 50])

plt.show()
fig.savefig(str(pdf_file0)+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file0)+'.png', format='png',bbox_inches='tight')

#-----------------------------------------------------
print("Plotting Fig. 4c")
#-----------------------------------------------------
Ta   = -depthS/strend
fa   = np.zeros(ns)
fa   = -np.asarray(sstep)/(depthS)
Ta0 = 0.0
rinfty = strend*chi0
tmin = -10.0
tmax = 120.0
xrange_set = False
xrange_set = True
lstyle = '-'

if xrange_set:
    nstart = nstep-0
    nstart = nstep-1
else:
    nstep = 0
    nstart = 0

# ---- Coulomb stress loading (input and interpolated input)
#fig = plt.figure(1, figsize=(12, 8))
fig = plt.figure(1, figsize=(6, 8))
ax1a = fig.add_subplot(111)
ax1a.set_ylabel('$\sigma_c/\dot{\sigma}_cT_a$', fontsize=20)
ax1a.set_xlabel(r'$(t-t_0)/T_a$', fontsize=20)
ax1a.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1a.tick_params(axis = 'both', which = 'minor', labelsize = 18)

for i in range(ns):
    ax1a.plot((t[nstart:]-t[nstep])/Ta, cfs[i,nstart:]/(Ta*strend), linewidth=2.0, ls=lstyle, color='black', label=r'$\Delta\sigma/\delta\sigma=$'+'{:.1f}'.format(float(fa[i])))

plt.legend(loc='lower right',fontsize=20)
if xrange_set:
    #plt.xlim([tmin, tmax])
    print(' to be defined')
plt.show()

# ---- earthquake rate ------------------------------------
fig = plt.figure(2, figsize=(6, 8))
ax1b = fig.add_subplot(111)
ax1b.set_xlabel(r'$(t-t_0)/T_a$', fontsize=20)
ax1b.set_ylabel(r'$r / \chi_0 V \dot{\sigma}_c$', fontsize=20)
ax1b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax1b.plot([(t[nstart]-t[nstep])/Ta, (t[nt-1]-t[nstep])/Ta], [1, 1], linewidth=1.0, ls='-', color='black')
for i in range(ns-1):
    #ax1b.plot((t[nstart:]-t[nstep])/Ta , r_lcm[i,nstart:]/rinfty, linewidth=1.5, ls=lstyle, color='lightgray')
    ax1b.plot((t[nstart:]-t[nstep])/Ta , r_tdsm[i,nstart:]/rinfty, linewidth=2.0, ls='-', color='red')
    ax1b.plot((t[nstart:]-t[nstep])/Ta , r_rsm[i,nstart:]/rinfty, linewidth=1.0, ls='--', color='blue')


    ax1b.plot((t[nstart:]-t[nstep])/Ta , r_tdsm_theo[i,nstart:]/rinfty, linewidth=3.5, ls='--', color='yellow')
    ax1b.plot((t[nstart:]-t[nstep])/Ta , r_rsm_theo[i,nstart:]/rinfty, linewidth=1.5, ls='dotted', color='green')

i=ns-1
#ax1b.plot((t[nstart:]-t[nstep])/Ta , r_lcm[i,nstart:]/rinfty, linewidth=1.5, ls=lstyle, color='lightgray', label=r'LCM')
ax1b.plot((t[nstart:]-t[nstep])/Ta , r_tdsm[i,nstart:]/rinfty, linewidth=2.0, ls='-', color='red', label=r'TDSM')
ax1b.plot((t[nstart:]-t[nstep])/Ta , r_rsm[i,nstart:]/rinfty, linewidth=1.0, ls='--', color='blue', label=r'RSM')
plt.text((t[nstart+2]-t[nstep])/Ta , r_tdsm[i,nstart+2]/rinfty ,r'$\Delta\sigma/\delta\sigma=$'+'{:.0f}'.format(float(fa[i])), fontsize=20)

ax1b.plot((t[nstart:]-t[nstep])/Ta , r_tdsm_theo[i,nstart:]/rinfty, linewidth=3.5, ls='--', color='yellow')
ax1b.plot((t[nstart:]-t[nstep])/Ta , r_rsm_theo[i,nstart:]/rinfty, linewidth=3.5, ls='dotted', color='green')

for i in range(ns-1):
    #ax1b.plot((t[nstart+2]-t[nstep])/Ta , r_tdsm[i,nstart+2]/rinfty , marker='o' , color='black', markersize=2)
    plt.text((t[nstart+2]-t[nstep])/Ta , r_tdsm[i,nstart+2]/rinfty ,r'{:.1f}'.format(float(fa[i])), fontsize=20)

plt.legend(loc='upper right',fontsize=20)
plt.figtext(0.0, 0.87, 'c)', fontsize=20)
plt.figtext(0.57, 0.70, r'$T_a=$'+'{:.1f}'.format(float(Ta/(24*hours)))+' days', fontsize=20)
#plt.text((t[nstart+2]-t[nstep])/Ta , r_tdsm[i,nstart+2]/rinfty , r'$\Delta\sigma/\delta\sigma=$', fontsize=20)
if xrange_set:
    #plt.xlim([tmin, tmax])
    plt.ylim([0.9, 200])
    ax1b.set_xscale('log')
    ax1b.set_yscale('log')
    #ax1b.set_yticks([1,10,50,100,200])
    ax1b.set_yticks([1,10,50,100])
    ax1b.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
else:
    plt.ylim([0, 2])

plt.show()
fig.savefig(str(pdf_file1)+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file1)+'.png', format='png',bbox_inches='tight')
