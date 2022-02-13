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
tend   = 600*hours
tend   = 20000*hours
#deltat = 0.2*hours
deltat = 600.0
deltat = 18000.0
#nt = int(math.ceil((tend - tstart) / deltat))
nt = np.floor( (tend-tstart) / deltat).astype(int)

#strend  = 3.0E-6
#strend  = 2.5E-6
#strend  = 0.5E-6
strend  = 1.0E-8
##sstep  =  [0.05, 0.55, 1.05, 1.55, 2.05, 2.55, 3.05, 3.55]  # stress step in MPa
##sstep  =  [0.05, 0.45, 0.85, 1.25, 1.65, 2.05, 2.45, 2.85]  # stress step in MPa
sstep  =  [0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 3.0]  # stress step in MPa
sstep  =  [0.30, 0.60, 0.90, 1.20, 1.50, 1.80, 3.0]  # stress step in MPa
#sstep  =  [0.60, 0.90, 1.20, 1.50, 1.80, 2.10, 2.40, 2.70, 3.0]  # stress step in MPa
##sstep  =  [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75]  # stress step in MPa
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
print('tectonic trends =',strend,' MPa/s =',strend*24.*hours,' MPa/days,  dot(sigma)*deltat/deltaS=',-strend*deltat/depthS)
print('deltaS=',deltaS,' MPa, from 0 to ',sigma_max,' MPa')

#-----------------------------------------------------
print("Calculate earthquake rates with tdsm, lcm and rsm")
#-----------------------------------------------------
ns   = len(sstep)
cfs    = np.zeros((ns,nt))
r_tdsm_eq8 = np.zeros(nt)
K1         = np.zeros(nt)
romo       = np.zeros(nt)
romo1      = np.zeros(nt)
r_tdsm = np.zeros((ns,nt))
r_lcm  = np.zeros((ns,nt))
r_rsm  = np.zeros((ns,nt))
n_tdsm = np.zeros((ns,nt-1))
n_lcm  = np.zeros((ns,nt-1))
n_rsm  = np.zeros((ns,nt-1))

r_tdsm_theo = np.zeros((ns,nt))
prod        = np.zeros(ns)      # estimated from integration of theoretical r(t)
ctheo       = np.zeros(ns)      # argument of integration kernel 
r_rsm_theo  = np.zeros((ns,nt))
prod_rsm_theo  = np.zeros(ns)

# ---- initial dfecay according to equation 8 for zero tectonic stress rate
#strend_background = 0.0
strend_background = 0.00001
#strend_background = 0.0001
#strend_background = 0.0002
#loading = BackgroundLoading(_config=tdsm.config, strend=0.0, deltat=deltat, tstart=tstart, tend=tend)
loading = BackgroundLoading(_config=tdsm.config, strend=strend_background, deltat=deltat, tstart=tstart, tend=tend)
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
good1 = t >= (tstep-tstart)-deltat
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
    #f1 = 0.15
    #f1 = 150
    #f2 = 0.2
    #K = -f2*depthS+f1*sstep[i]
    #K = -depthS*np.exp(-(t[good]-tstep)/Ta)/(1.-np.exp(sstep[i]/depthS))
    #K = -depthS*deltat/(Ta*(1.+(np.exp(sstep[i]/depthS)-1.)*1.))
    #K1[good] = -depthS*np.exp(-(t[good]-tstep)*sstep[i]/(-deltat*depthS)) -depthS*np.exp(-(t[good]-tstep)/Ta)
    #K1[good] = (-depthS*np.exp(-(t[good]-tstep)/Ta)) / (f1*np.exp(sstep[i]/depthS))
    #K1[good] = (-depthS*np.exp(-(t[good]-tstep)*sstep[i]/(-deltat*depthS)) -depthS*np.exp(-(t[good]-tstep)/Ta) ) / (f1*np.exp(sstep[i]/depthS))
    #r_tdsm_theo[i,good] = +chi0*K1[good]/(deltat+(t[good]-tstep))
    #f3 = 0.01
    #r_tdsm_theo[i,:] = r_tdsm_theo[i,:] +rinfty
    #print('i=',i,' sstep=',sstep[i],'K=',K1[nstep-2:nstep+4],' theo minmax=',np.amin(r_tdsm_theo[i,nstep-2:nstep+4])/rinfty, np.amax(r_tdsm_theo[i,nstep-2:nstep+4])/rinfty)

    r_rsm_theo[i,:] = rinfty
    r_rsm_theo[i,good] = rinfty/(1.+(np.exp(sstep[i]/depthS)-1.)*np.exp(-(t[good]-tstep)/Ta) )
    prod_rsm_theo[i] = rinfty/(1.+(np.exp(sstep[i]/depthS)-1.)*np.exp(-(t[nstep+1]-tstep)/Ta))
    print('i=',i,' sstep=',sstep[i],' rsm  minmax=',np.amin(r_rsm_theo[i,nstep-2:nstep+4])/rinfty, np.amax(r_rsm_theo[i,nstep-2:nstep+4])/rinfty,' Productivity=',prod_rsm_theo[i])

    #----------------------------------------------
    # ----- calculate own (new) theoretical Omori curves ----
    # t1 bestimmt die Hoehe des peaks
    # t0 bestimmt die Laenge der Abfallkurve, aehnlich wie Ta - mit trial and error auf 0.14*deltat 
    #----------------------------------------------
    dS = -depthS
    strend_tmp = strend
    #strend_tmp = 5*strend
    t0 = 0.1*deltat  # t0 muss deutlich kleiner als deltat sein, damit die peak Amplituden richtig
    #t1 = delta

    factor = 380.+(10.-380.)*(sstep[i]-sstep[0])/(sstep[-1]-sstep[0])
    #print(i,'sstep=',sstep[0],sstep[-1])
    #print('factor=',factor)
    t1 = dS/(sstep[i]/deltat-strend_tmp)
    #t1 = dS/(sstep[i]/(0.1*deltat)-strend_tmp)
    #t1 = dS/(sstep[i]/(10*deltat)-strend_tmp)
    #t1 = dS/(sstep[i]/(380*deltat)-strend_tmp)
    t1 = dS/(sstep[i]/(factor*deltat)-strend_tmp)

    print('i=',i,' sstep=',sstep[i],' deltat=',deltat,' dS=',dS,' strend=',strend)
    a = dS/(sstep[i]/deltat-strend_tmp)
    print('t1 own=', a)
    print('t1 rsm=', 1./prod_rsm_theo[i],a*prod_rsm_theo[i]) 
    b= np.exp(-sstep[i]/dS)*dS/strend/(1.-np.exp(-sstep[i]/dS))
    print('t1 c  =', b, b*prod_rsm_theo[i],a*prod_rsm_theo[i] )

    #p2 = np.exp(-sigma_max/dS)/t0
    #p2 = np.exp(-0.35*sigma_max/dS)/t0   # 
    #p2 = np.exp(-0.40*sigma_max/dS)/t0   # 
    #p2 = np.exp(-0.50*sigma_max/dS)/t0   # mit 0.5 besserer fit wenn strend>0

    rback = rinfty
    #rback = chi0*dS*(1.-np.exp(-strend*t0/dS))/t0

    print('# estimate factor prod  und   c  by integration') 
    nev_lcm = 1.0*chi0*(sstep[i])-rback*deltat
    v = dS*deltat/sstep[i]
    #values = np.linspace(v/500., v, 5, endpoint=True)
    values = np.linspace(1.0*v, 800*v, 5000, endpoint=True)
    for k, t1 in enumerate(values):
        t0 = 0.14*deltat
        p2 = np.exp(-0.50*sigma_max/dS)  # kleine Korrektur durch endliche Integration
        r_tdsm_theo[i,good1] = chi0*dS/(t1+(t[good1]-tstep))
        r_tdsm_theo[i,good1] = r_tdsm_theo[i,good1]*(np.exp(-p2*(t[good1]-t[nstep]+t1)/t0)-np.exp(-np.exp(-strend_tmp*(t[good1]-tstep)/dS)*(t[good1]-tstep+t1)/t0) )
        nev = np.trapz(r_tdsm_theo[i,good1], dx=deltat)
        if nev <= nev_lcm:
            print('k=',k,' c=',t1,' nev=',nev,' n_lcm=',nev_lcm)
            prod[i] = t1
            break

    #print(' t1=',t1,' nev=',nev)
    r_tdsm_theo[i,:] = r_tdsm_theo[i,:]+rback

#----------------------------
# theoretical spontaneous decay of fully filled house
# Fazit der Tests: 
# 1. Die Loesung erklaert den Peak und background und Ta. Der detailliete Verlauf fuer kleine Zeiten weicht etwas ab, evtl. durch Approx bei Herleitung
# 2. Variieren von Ta veraendert die Abklingzeit wie gewuenscht. 
# 3. Die Peakhoehe wird fuer ein volles Haus reproduziert, wenn t1 = deltat. Wir t1 > deltat gewaehlt, dann wird die Peak HOehe kleiner
#    Die Abklingform (Ta) bleiben aber gleich, ebenso der BAckground. Man kann durch Steurung von t1 die Response fuer unterschiedliche Depletion Stadien
#    simulieren. t1 = depthS / ( Delta Sigma / Delta t - strend )  
# 4. t0 muss etwa 0.1 Deltat sein, damit Peak Amplituden richtig bestimmt
# 5. Die Integration in the Theorie geht bis zeta = infty. Im Programm aber nur bis etwa 2 * depthS
#----------------------------
A = chi0*(-depthS)
t1 = deltat   # t1 muss auf deltat gesetzt werden, damit uebereinstimmung mit numerischer Simulation
#t1 = 5*deltat   # t1 muss auf deltat gesetzt werden, damit uebereinstimmung mit numerischer Simulation
t0 = 0.1*deltat  # t0 muss deutlich kleiner als deltat sein, damit die peak Amplituden richtig 
dS = -depthS
strend_tmp = strend_background
#strend_tmp = 1E+0*strend_background
p2 = np.exp(-sigma_max/dS)/t0
#p2 = np.exp(-0.35*sigma_max/dS)/t0   # mit 0.35 erhalte ich perfekten fit wenn strend>0

# 1. time dependent component of spontaneous depletion under steady rate (see appendix)
romo[good] = A/(t1+(t[good]-t[nstep]))
#romo1[good] = romo[good]*(1.-np.exp(-np.exp(-strend_background*(t[good]-t[nstep])/(-depthS))*(t[good]-t[nstep]+t1)/(t1)) )
romo1[good] = romo[good]*(np.exp(-p2*(t[good]-t[nstep]+t1))-np.exp(-np.exp(-strend_tmp*(t[good]-t[nstep])/dS)*(t[good]-t[nstep]+t1)/t0) )
# 2. time independent component of steady rate loading (see appendix)
# rback = chi0*strend_background
rback = chi0*(-depthS)*(1.-np.exp(-strend_background*t0/dS))/t0
romo1[good] = romo1[good]+rback

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
    nstart = 1

fig = plt.figure(2, figsize=(6, 8))
ax = fig.add_subplot(111)
#ax.set_xlabel(r'$(t-t_0)/T_a$', fontsize=20)
#ax.set_xlabel(r'$(t-t_0)/T_a$', fontsize=20)
ax.set_xlabel(r'$t/\Delta t$', fontsize=20)
ax.set_ylabel(r'$r \Delta t / \chi_0 V \delta \sigma$', fontsize=20)
ax.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax.plot([(t[0])/deltat,(t[-1])/deltat], [0, 0], linewidth=1.0, ls='-', color='black')
ax.plot((t[nstart:]-t[nstep])/deltat, r_tdsm_eq8[0:-nstart]/scal1, linewidth=2.5, ls='-', color='red', label=r'numerical')
ax.plot((t[nstart:]-t[nstep])/deltat, romo1[nstart:]/scal1, linewidth=3.5, color='black', ls='--', label=r'$equation (9)$')
#ax.plot((t[nstart:]-t[nstep])/deltat, romo[nstart:]/scal1, linewidth=3.5, color='green', ls='--')

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
ax1b.set_xlabel(r'$(t-t_s)/T_a$', fontsize=20)
ax1b.set_ylabel(r'$r / \chi_0 V \dot{\sigma}_c$', fontsize=20)
ax1b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax1b.plot([(t[nstart]-t[nstep])/Ta, (t[nt-1]-t[nstep])/Ta], [1, 1], linewidth=1.0, ls='-', color='black')
for i in range(ns-1):
    #ax1b.plot((t[nstart:]-t[nstep])/Ta , r_lcm[i,nstart:]/rinfty, linewidth=1.5, ls=lstyle, color='lightgray')
    ax1b.plot((t[nstart:]-t[nstep])/Ta , r_tdsm[i,nstart:]/rinfty, linewidth=2.0, ls='-', color='red')
    #ax1b.plot((t[nstart:]-t[nstep])/Ta , r_rsm[i,nstart:]/rinfty, linewidth=1.0, ls='--', color='blue')

    ax1b.plot((t[nstart:]-t[nstep]+deltat)/Ta , r_tdsm_theo[i,nstart:]/rinfty, linewidth=1.5, ls='--', color='blue')
    #ax1b.plot((t[nstart:]-t[nstep])/Ta , r_rsm_theo[i,nstart:]/rinfty, linewidth=1.5, ls='dashed', color='yellow')
    #ax1b.plot(t[nstep-5]/Ta , prod_rsm_theo[i]/rinfty, ls='None', ms=7, marker='o', mec= 'darkgreen', mfc='None', mew=2.0)
    #ax1b.plot(deltat/Ta, prod_rsm_theo[i]/rinfty, ls='None', ms=7, marker='o', mec= 'darkgreen', mfc='None', mew=2.0)

i=ns-1
#ax1b.plot((t[nstart:]-t[nstep])/Ta , r_lcm[i,nstart:]/rinfty, linewidth=1.5, ls=lstyle, color='lightgray', label=r'LCM')
ax1b.plot((t[nstart:]-t[nstep])/Ta , r_tdsm[i,nstart:]/rinfty, linewidth=2.0, ls='-', color='red', label=r'TDSM')
#ax1b.plot((t[nstart:]-t[nstep])/Ta , r_rsm[i,nstart:]/rinfty, linewidth=1.5, ls='--', color='blue', label=r'RSM')
ax1b.plot((t[nstart:]-t[nstep])/Ta , r_rsm[2,nstart:]/rinfty, linewidth=1.5, ls='dotted', color='green', label=r'RSM')
plt.text((t[nstart+2]-t[nstep])/Ta , 1.3*r_tdsm[i,nstart+2]/rinfty ,r'$\Delta\sigma/\delta\sigma=$'+'{:.0f}'.format(float(fa[i])), fontsize=20)
plt.text((t[nstart+2]-t[nstep])/Ta , 0.97*r_tdsm[i,nstart+2]/rinfty ,r'$t_1/T_a=$'+'{:.4f}'.format(float(prod[i]/Ta)), fontsize=20)

ax1b.plot((t[nstart:]-t[nstep]+deltat)/Ta , r_tdsm_theo[i,nstart:]/rinfty, linewidth=1.5, ls='--', color='blue', label='eq. (12)')
#ax1b.plot((t[nstart:]-t[nstep])/Ta , r_rsm_theo[i,nstart:]/rinfty, linewidth=3.5, ls='--', color='yellow')
#ax1b.plot(deltat/Ta , prod_rsm_theo[i]/rinfty, ls='None', ms=7, marker='o', mec= 'darkgreen', mfc='None', mew=2.0)

for i in range(ns-1):
    #ax1b.plot((t[nstart+2]-t[nstep])/Ta , r_tdsm[i,nstart+2]/rinfty , marker='o' , color='black', markersize=2)
    #plt.text((t[nstart+2]-t[nstep])/Ta , r_tdsm[i,nstart+2]/rinfty ,r'{:.1f}'.format(float(fa[i])), fontsize=20)
    plt.text((t[nstart+2]-t[nstep])/Ta , r_tdsm[i,nstart+2]/rinfty ,r'{:.1f};{:2.3f}'.format(float(fa[i]),float(prod[i]/Ta)), fontsize=20)

plt.legend(loc='upper right',fontsize=20)
plt.figtext(0.0, 0.87, 'c)', fontsize=20)
plt.figtext(0.55, 0.63, r'$T_a=$'+'{:.0f}'.format(float(Ta/(24*hours)))+' days', fontsize=20)
#plt.text((t[nstart+2]-t[nstep])/Ta , r_tdsm[i,nstart+2]/rinfty , r'$\Delta\sigma/\delta\sigma=$', fontsize=20)
if xrange_set:
    #plt.xlim([tmin, tmax])
    #plt.ylim([0.9, 200])
    plt.ylim([0.9, 200])
    plt.ylim([0.9, 1500])
    ax1b.set_xscale('log')
    ax1b.set_yscale('log')
    #ax1b.set_yticks([1,10,50,100,200])
    ax1b.set_yticks([1,10,50,100])
    ax1b.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
else:
    #plt.ylim([0, 2])
    print(' to be defined')

plt.show()
fig.savefig(str(pdf_file1)+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file1)+'.png', format='png',bbox_inches='tight')
