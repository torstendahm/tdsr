#--------------------------
# Plot Fig. 5a (see manuscript on TDSM ,Dahm, 2022) 
#--------------------------
import sys
sys.path.insert(0, "../")
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, BackgroundLoading
import matplotlib.pyplot as plt
import matplotlib.ticker
import math
import numpy as np

current_dir = Path(__file__).parent.parent
config_file = current_dir / "config.toml"
pdf_file1 = current_dir / "../plots/Dahm_fig5a"

# ----  real data for comparison
fn_m55 = "../data/hainzl2008_fig2b_5.5-6.dat"
fn_m60 = "../data/hainzl2008_fig2b_6-6.5.dat"
fn_m65 = "../data/hainzl2008_fig2b_6.5-7.dat"
fn_m70 = "../data/hainzl2008_fig2b_7-7.5.dat"
fn_m75 = "../data/hainzl2008_fig2b_7.5-9.dat"

print("set-up the tdsm, lcm and rsm class environments")
tdsm = TDSM(config=Config.open(config_file))
lcm  = LCM(config=Config.open(config_file))
trad = Traditional(config=Config.open(config_file))
rsm  = RSM(config=Config.open(config_file))


print("define model parameter for plotting, if different from config file")
for_test_only = True
for_test_only = False

hours = 3600.
days = hours * 24.
tstart = 0*days
tend   = 1000*days
deltat = 0.04*days
nt = np.floor( (tend-tstart) / deltat).astype(int)

if for_test_only:
    strend  = [ 1.0E-8 , 1.0E-8 ]   # test mit brut force auf kleinere stress rate zu gehen
    sstep  =  [ 0.57E0, 5.9E0 ]  # stress step in MPa
else:
    strend  = [ 5.0E-10 , 5.0E-10 ]
    sstep  =  [ 0.0270E0, 0.3E0 ]

tstep  = 5.0*days       # time when stress step is acting
nstep = np.floor( (tstep-tstart) / deltat).astype(int)+1


if for_test_only:
    deltaS    = 1./300.      # increment do discretize Coulomb stress axis
    sigma_max = 20000.*deltaS # maximum depth on Coulomb axis (limit of integral)
    precision = 25
    chi0   = 5.00E4            # susceptibility to trigger earthquakes by unit stress increase
    depthS = [-0.072 , -0.55]    
else:
    deltaS    = 1./3000.     # increment do discretize Coulomb stress axis
    sigma_max = 50000.*deltaS # maximum depth on Coulomb axis (limit of integral)
    precision = 25
    chi0   = 1.00E6
    depthS = [-0.0036 , -0.0275]    

print('deltat=',deltat,' sec = ',deltat/hours,' hours, from 0 to ',tend/hours,' hours, samples =',nt)
print('Ta   =',-np.asarray(depthS)/np.asarray(strend),' sec = ',-np.asarray(depthS)/(np.asarray(strend)*hours),' hours')
print('suscept. chi0=',chi0,' #/MPa')
print('stress step=[',sstep,'] MPa')
print('skin depth =',-np.asarray(depthS),' MPa')
print('tectonic trends =',strend,' MPa/s =',np.asarray(strend)*24.*hours,' MPa/days =',np.asarray(strend)*24.*hours*365.,' MPa/yr')
print('dot(sigma)*deltat/deltaS=',-np.asarray(strend)*deltat/np.asarray(depthS))
print('deltaS=',deltaS,' MPa, from 0 to ',sigma_max,' MPa')

#-------------------------------------------
# read real data for comparison 
#-------------------------------------------
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
print("Calculate earthquake rates with tdsm, lcm and rsm")
#-----------------------------------------------------
ns   = len(sstep)
cfs    = np.zeros((ns,nt))
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

rinfty      = np.zeros(ns)
Ta          = np.zeros(ns)

# ---- calculate equilibrium distribution of tectonic loading, to be used later to avoid initial oscillations
for i in range(ns):
    loading = BackgroundLoading(_config=tdsm.config, strend=strend[i], deltat=deltat, tstart=tstart, tend=tend)
    config, t, chiz_background, cf, r, xn = tdsm(loading=loading, chi0=chi0, depthS=depthS[i], deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)

    print('background')
    print('len(cf)',len(cf))
    print('len(chiz_background)',len(chiz_background))
    rinfty[i] = strend[i]*chi0
    Ta[i]   = -depthS[i]/strend[i]


    good = t >= (tstep-tstart)
    good1 = t >= (tstep-tstart)-deltat

    loading = StepLoading(_config=tdsm.config, strend=strend[i], sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
    config, t, chiz, cf, r, xn = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS[i], deltaS=deltaS, sigma_max=sigma_max, precision=precision, deltat=deltat, tstart=tstart, tend=tend)
    cfs[i,:]    = cf[:]
    r_tdsm[i,:] = r[:]
    n_tdsm[i,:] = xn[:]

    loading = StepLoading(_config=trad.config, strend=strend[i], sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
    config, t, chiz, cf, r, xn = trad(loading=loading, chi0=chi0, deltat=deltat, tstart=0, tend=tend)
    r_lcm[i,:] = r
    n_lcm[i,:] = xn

    loading = StepLoading(_config=rsm.config, strend=strend[i], sstep=sstep[i], tstep=tstep, deltat=deltat, tstart=tstart, tend=tend)
    config, t, chiz, cf, r, xn = rsm(loading=loading, chi0=chi0, depthS=depthS[i], deltat=deltat, tstart=0, tend=tend)
    r_rsm[i,:] = r
    n_rsm[i,:] = xn

    #----------------------------------------------
    # ----- calculate own (new) theoretical Omori curves ----
    # t1 bestimmt die Hoehe des peaks
    # t0 bestimmt die Laenge der Abfallkurve, aehnlich wie Ta - mit trial and error auf 0.14*deltat 
    #----------------------------------------------
    dS = -depthS[i]
    strend_tmp = strend[i]
    t0 = 0.1*deltat  # t0 muss deutlich kleiner als deltat sein, damit die peak Amplituden richtig

    factor = 380.+(10.-380.)*(sstep[i]-sstep[0])/(sstep[-1]-sstep[0])
    t1 = dS/(sstep[i]/(factor*deltat)-strend_tmp)

    print('i=',i,' sstep=',sstep[i],' deltat=',deltat,' dS=',dS,' strend=',strend)
    a = dS/(sstep[i]/deltat-strend_tmp)
    print('t1 own=', a)
    print('t1 rsm=', 1./prod_rsm_theo[i],a*prod_rsm_theo[i]) 
    b= np.exp(-sstep[i]/dS)*dS/strend/(1.-np.exp(-sstep[i]/dS))
    print('t1 c  =', b, b*prod_rsm_theo[i],a*prod_rsm_theo[i] )
    rback = rinfty[i]

    print('# estimate factor prod  und   c  by integration') 
    nev_lcm = 1.0*chi0*(sstep[i])-rback*deltat
    v = dS*deltat/sstep[i]
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

    r_tdsm_theo[i,:] = r_tdsm_theo[i,:]+rback

#-----------------------------------------------------
print("Plotting Fig. 6c")
#-----------------------------------------------------
Ta   = -np.asarray(depthS)/np.asarray(strend)
#fa   = np.zeros(ns)
fa   = -np.asarray(sstep)/np.asarray(depthS)
Ta0 = 1.0
tmin = 0.01
if for_test_only:
    tmax = 1000.0
else:
    tmax = 1000.0
xrange_set = False
xrange_set = True
lstyle = '-'

if xrange_set:
    nstart = nstep-0
    nstart = nstep-1
else:
    nstart = nstep-0
    nstart = nstep-1

fig = plt.figure(2, figsize=(6, 8))
ax1b = fig.add_subplot(111)
ax1b.set_xlabel(r'$(t-t_s)$ days', fontsize=20)
ax1b.set_ylabel(r'$r$ (1/days)', fontsize=20)
ax1b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

print('tstart=',nstart,nstep,t[nstart],t[nstep])
print('infty =',rinfty[0],rinfty[1])
print('Ta0 =',Ta0,'  Ta=',Ta)
print('t     =',len(t),t[0],t[-1:])
print('r_tdsm=',len(r_tdsm[0,:]),r_tdsm[0,0],r_tdsm[0,-1:])
col = ['blue' , 'green' ]
line = ['-' , '-' ]
tshift = [t_m55[0], t_m75[0] ]
for i in range(ns):

    ax1b.plot((t[nstart:]-t[nstep])/(Ta0*days)+tshift[i] , r_tdsm[i,nstart:], linewidth=2.0, ls='-', color=col[i])

#---- sum of masurements between m=5.5 and M=6.5 (M6 +-0.5) and m=7.0 and M=8.0 (M7.5 +-0.5)
if for_test_only:
    plt.plot(t_m55[0:18], (r_m55[0:18]+r_m60[0:18])/2., ls='None', ms=15, marker='o', mec= 'blue', mfc='None', mew=2.0, label=r'$M=6\pm 0.5$')
    plt.plot(t_m75, (r_m70+r_m75)/2., ls='None', ms=15, marker='o', mec= 'green', mfc='None', mew=2.0, label=r'$M\geq 7$')
else:
    plt.plot(t_m55[0:18], (r_m55[0:18]+r_m60[0:18])/2., ls='None', ms=15, marker='o', mec= 'blue', mfc='None', mew=2.0, label=r'$M=6\pm 0.5$')
    plt.plot(t_m75, (r_m70+r_m75)/2., ls='None', ms=15, marker='o', mec= 'green', mfc='None', mew=2.0, label=r'$M\geq 7$')

plt.legend(loc='upper right',fontsize=20)
plt.figtext(0.0, 0.87, 'a)', fontsize=20)
plt.figtext(0.55, 0.70, r'$\delta\sigma=$'+'{:.0f}'.format(float(-depthS[1]*1000))+r'kPa', fontsize=20, color='green')
plt.figtext(0.55, 0.65, r'$T_a=$'+'{:.0f}'.format(float((-depthS[1]/strend[1])/days))+'days', fontsize=20, color='green')
plt.figtext(0.135, 0.18, r'$\delta\sigma=$'+'{:.0f}'.format(float(-depthS[0]*1000))+r'kPa', fontsize=20, color='blue')
plt.figtext(0.135, 0.13, r'$T_a=$'+'{:.0f}'.format(float((-depthS[0]/strend[0])/days))+'days', fontsize=20, color='blue')
if xrange_set:
    plt.xlim([tmin, tmax])
    if for_test_only:
        plt.ylim([0.0003, 50])
    else:
        plt.ylim([0.0003, 50])
    ax1b.set_xscale('log')
    ax1b.set_yscale('log')
else:
    print(' to be defined')

plt.show()
fig.savefig(str(pdf_file1)+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(str(pdf_file1)+'.png', format='png',bbox_inches='tight')
