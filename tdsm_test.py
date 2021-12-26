################################
# Coulomb Failure Model 
# New approach to integrate the probability function along Coulomb axis 
# rate und ratez (chi, chiz, p, pz) to compare CFM versus RSM)
# T. Dahm, 15.06.2020
################################
import math
import numpy as num
from scipy import special as spl

import matplotlib
##matplotlib.use('PS')
##matplotlib.use('TkAgg')
matplotlib.use('Qt4Agg')
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mpl_colors
import matplotlib.cm as mpl_cm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

km = 1000.
hours = 3600.
days = hours * 24.

# ----------------- stress loading parameter and output file -------------
# faengt bei S0 an, d.h. S0 ist im Regelfall 0
#deltat = 1*hours  # sampleinterval in t in seconds
#n1 = 10          # fist phase must be larger at least 2
#n2 = 11          # second phase at least n1+1, better n1+2 as central differences taken
deltat = 0.2*hours  # sampleinterval in t in seconds
n1 = 50          # fist phase must be larger at least 2
n2 = 51          # second phase at least n1+1, better n1+2 as central differences taken
sc0 = 0.0        # Coulomb stress values
sc1 = 0.5
#sc2 = 1.5
#sc3 = 2.5
sc2 = 1.0
sc3 = 2.0
sc3 = 2.0
sc2b = 2.5
sc3b = 3.5
sc2c = 0.8
sc3c = 1.8
# In dem Beispiel gibt es 150 Samples (30 /0.2 ). Tectonic rate wird aus ersten 50 sample berechnet, 
tectstressrate = (sc1-sc0)/(n1*deltat)
chi0 = 10.        # susceptibility
Sshadow = -0.5    # initial stress shadow
Sshadow = 0.0    # initial stress shadow
file = 'CFMprob_ramp2_paper_chi%2.0f' % (chi0)
file2 = 'CFMprob_ramp2_paper_totalnumber_chi%2.0f' % (chi0)

# ----------------- trigger probability and chi distribution -------------
#deltaS = 0.02     # sampling of Coulomb stress axis
#sigma_max = 2.
deltaS = 0.01     # sampling of Coulomb stress axis
#deltaS = 0.1     # sampling of Coulomb stress axis
sigma_max = 10.0

#---------------- boundary layer nucleation effect -----------------------
depthS = -0.2     # Tiefe der Einflussfunktion 
#depthS = -0.5     # Tiefe der Einflussfunktion 
#depthS = -0.7     # Tiefe der Einflussfunktion 
#depthS = -1.0     # Tiefe der Einflussfunktion 
#depthS = -2.0     # Tiefe der Einflussfunktion 

# mit arr() wird definiert, dass np.array() immer floating point variablen
def arr(l):
    return num.array(l, dtype=num.float)

def cfun(kappa, eta, lamb_u, lamb, alpha, shear, scal):
    term1 = (lamb_u-lamb)*(lamb+2.*shear)
    term2 = alpha*alpha*(lamb_u+2.*shear)
    c = term1/(scal*term2)
    return c

def gridrange(amin, amax, astep, rounding=math.ceil):
    amin = float(amin)
    amax = float(amax)
    astep = float(astep)
    na = int(rounding((amax-amin)/astep))
    amax = amin + (na-1) * astep
    a = num.linspace(amin, amax, na)
    return amin, amax, na, a

def CMrate(r, c, no_k, dS, CFS):
    R = np.zeros(1+len(dS))
    R[0] = r
    print('c:%f range: %d %d' % (c,len(dS),len(CFS)))
    for i in range(len(dS)):
#        print('dS(%d): %e of %d' % (i, dS[i], len(dS)))
        if (no_k == 0):
            if dS[i] > 0 and i>0:
                R[1+i] = (r + c * dS[i])
            else:
                R[1+i] = r
        else:
            if dS[i] > 0 and i>0 and (CFS[i] >= np.amax(CFS[:i])):
                R[1+i] = (r + c * dS[i])
            else:
                R[1+i] = r
    return R

def scc(tau_max, sign_max, pf, fric):
    '''
    Coulomb stress change = Coulomb failure stress 
    Definition: compressive stress is negative

    Input: 
    tau_rake: shear stress in direction in given plane in direction slip vector [Pa]
    tau_max : shear stress on 45 deg plane [Pa]
    sign: normal stress on plane [Pa] (compressive stress is negative)
    sign_max: normal stress for 45 deg plane [Pa] 
    pf  : fluid pore pressure [MPa] (here positive press ure is compression)
    fric: friction coefficient (dimensionless)

    return:    scc [Pa]
    '''

    scc = num.sqrt(1.+fric*fric)*tau_max +fric*(sign_max+pf)

    return scc

def scf(tau, sign, pf, fric):
    '''
    Coulomb stress change = Coulomb failure stress 
    Definition: compressive stress is negative

    Input: 
    tau: shear stress in plane [Pa]
    sign: normal stress on plane [Pa] (compressive stress is negative)
    pf  : fluid pore pressure [MPa] (here positive press ure is compression)
    fric: friction coefficient (dimensionless)

    return:    scf [Pa]
    '''

    scf       = tau + fric*(sign +pf)

    return scf

def scalc(sig, strike, dip, rake):
    '''
    Version according to my script seismology3
    using Aki's coordinate system: 0=x --> north
                                   1=y --> east
                                   2=z --> depth (i.e. positive values)

    Input:
    sig:       stress tensor [Pa], 3x3 numpy array
    strike:    strike of the earthquake mechanism [degree], float
    dip:       dip of the earthquake mechanism [degree], float
    rake:      rake of the earthquake mechanism [degree], float

    return:  tau_rake, tau_plane, tau_max, sign, sign_max 
             tau_rake [Pa]: shear stress in given plane in rake direction 
             tau_plane[Pa]: max shear stress in given plane
             taum_max [Pa]: max shear stress on arbitrary plane
             sign [Pa]:     normal stress on given plane
             sign_max [Pa]: normal stress on plane with max shear stress
    '''
    if sig.shape != (3, 3):
        print('input shape of sig not correct')

    deg2rad = num.pi / 180.0
    st0 = deg2rad * strike
    di0 = deg2rad * dip
    ra0 = deg2rad * rake

    nvec = num.zeros((3))
    slip = num.zeros((3))
    trac = num.zeros((3))
    taus = num.zeros((3))

    nvec[0] = -num.sin(di0)*num.sin(st0)
    nvec[1] = +num.sin(di0)*num.cos(st0)
    nvec[2] = -num.cos(di0)

    slip[0] = +num.cos(ra0)*num.cos(st0)+num.cos(di0)*num.sin(ra0)*num.sin(st0)
    slip[1] = +num.cos(ra0)*num.sin(st0)-num.cos(di0)*num.sin(ra0)*num.cos(st0)
    slip[2] = -num.sin(di0)*num.sin(ra0)

    trac = num.dot(sig, nvec.T)
    sign = num.dot(trac,nvec.T)
    taus = num.cross( nvec, num.cross(trac,nvec) )

    evalues, evecs = num.linalg.eig(stress)
    s1 = num.real( evalues.max() )
    s3 = num.real( evalues.min() )
    #shear_max = (num.real(evalues[0]) - num.real(evalues[1])) / 2.0    
    tau_max  = (s1 - s3) / 2.0    
    sign_max = (s1 + s3) / 2.0    
    #print(iz,iy,ix,it, "eigenwerte=", evalues, "vec=", evecs[:,0], evecs[:,1], evecs[:,2])

    #sign1 = 0.0
    #tau_rake1 = 0.0
    #for j in range(3):
    #    for i in range(3):
    #        sign1 += nvec[i] * sig[i][j] * nvec[j]
    #        tau_rake1 += slip[i] * sig[i][j] * nvec[j]

    tau_rake  = num.dot(taus, slip.T)
    tau_plane = num.linalg.norm(taus)

    return  tau_rake, tau_plane, tau_max, sign, sign_max 

############# main code ###########################
# ------------------ Coulomb stress axis , p und chi initialisieren  ----
smin, smax, nsigma, sigma = gridrange(-sigma_max, +sigma_max, deltaS)
chi = num.zeros(nsigma)
chiz = num.zeros(nsigma)
chib = num.zeros(nsigma)
chic = num.zeros(nsigma)
chizb = num.zeros(nsigma)
chizc = num.zeros(nsigma)
p   = num.zeros(nsigma)
pz  = num.zeros(nsigma)

# ------------------ stress loading function cf definieren ----
tmin, tmax, nt, t = gridrange(0*hours, 30*hours, deltat)
cf    = num.zeros(nt)
cfb   = num.zeros(nt)    # Zum Vergleich - stress step is doppelt so hoch
cfc   = num.zeros(nt)    # stress step ist kleiner
rate  = num.zeros(nt)
ratez  = num.zeros(nt)
rateb  = num.zeros(nt)
ratec  = num.zeros(nt)
ratezb  = num.zeros(nt)
ratezc  = num.zeros(nt)
rate1  = num.zeros(nt)
chi1  = num.zeros(nt)
neq   = num.zeros(nt-1)
neqb  = num.zeros(nt-1)
neqc  = num.zeros(nt-1)
neqz  = num.zeros(nt-1)
neqzb = num.zeros(nt-1)
neqzc = num.zeros(nt-1)
neq1  = num.zeros(nt-1)

for i in range(0,n1):
    cf[i] = sc0 + i*(sc1-sc0)/float(n1-1)
    #cfb[i] = sc0 + i*(sc1-sc0)/float(n1-1)
    #cfc[i] = sc0 + i*(sc1-sc0)/float(n1-1)
    cfb[i] = cf[i]
    cfc[i] = cf[i]
    #print('1',i,cf[i],i,sc1-sc0,n1)
for i in range(n1,n2):
    cf[i] = cf[n1-1] + (i-n1+1)*(sc2-sc1)/float(n2-n1)
    cfb[i] = cfb[n1-1] + (i-n1+1)*(sc2b-sc1)/float(n2-n1)
    cfc[i] = cfc[n1-1] + (i-n1+1)*(sc2c-sc1)/float(n2-n1)
    #print('2',i,cf[i],(i-n1+1),sc2-sc1,n2-n1)
for i in range(n2,nt):
    cf[i] = cf[n2-1] + (i-n2+1)*(sc3-sc2)/float(nt-n2)
    cfb[i] = cfb[n2-1] + (i-n2+1)*(sc3b-sc2b)/float(nt-n2)
    cfc[i] = cfc[n2-1] + (i-n2+1)*(sc3c-sc2c)/float(nt-n2)
    #print('3',i,cf[i],(i-n2+1),sc3-sc2,nt-n2)

# ---------------- Neue Ansatz fuer Bebenrate, deltan=rate, durch Integration ----
p   = num.heaviside( sigma , 1 )
pz   = num.heaviside( sigma , 1 )
# ---- probability function wird nach unten modifiziert durch boundary layer influence zone
ndepth  = num.around(depthS/deltaS,0).astype(num.int)
nzero = int(nsigma/2)
#for k in range(0,-2*ndepth):
#    #print('k=',k,' nzero+ndepth+k=',nzero+ndepth+k,pz[nzero+ndepth+k])
#    #pz[nzero+ndepth+k] = 0.0 +float(-k)/ndepth                     # linear
#    #pz[nzero+ndepth+k] = 0.0 +(1.-num.cos(-0.5*num.pi*k/ndepth ) )  # 1-cos() first cycle
#    # x = 1.-(1.-2.*num.pi*(-ndepth-k)/ndepth)*num.exp(+2.*num.pi*(-ndepth-k)/ndepth) # Brunes Moment
#    # pz[nzero+ndepth+k] =1.-x
#    pz[nzero+2*ndepth+k] = num.exp(+2.*num.pi*(-2*ndepth-k)/ndepth) # einfache e Abfall
#    #pz[nzero+ndepth+k] = num.exp(+2.*num.pi*(-ndepth-k)/ndepth) # einfache e Abfall
#    #print('k=',k,' nzero+ndepth+k=',nzero+ndepth+k,pz[nzero+ndepth+k])

for k in range(0,-6*ndepth):
    pz[nzero+6*ndepth+k] = num.exp((-6*ndepth-k)/ndepth) # einfache e Abfall

print('nzero=',nzero,' ndepth=',ndepth)

chi = num.heaviside( -sigma , 0 )
chib = num.heaviside( -sigma , 0 )
chic = num.heaviside( -sigma , 0 )
chiz = num.heaviside( -sigma , 0 )
chizb = num.heaviside( -sigma , 0 )
chizc = num.heaviside( -sigma , 0 )
print(sigma)
print(p)
print(pz)
print(chi)

# ...... initiale stress shadow falls noetig
nshift  = num.around(Sshadow/deltaS,0).astype(num.int)
chi = num.roll(chi,nshift)
chiz = num.roll(chi,nshift)
if (nshift > 0):
   chi[0:nshift] = 1.0
   chiz[0:nshift] = 1.0
elif (nshift < 0):
   chi[nshift:] = 0.0
   chiz[nshift:] = 0.0

#print('shadow',nshift,Sshadow,deltaS)
#print(chi)

# ...... incremental approach  (Achtung: Rundungsfehler bei shifts testweise durch resid beruecksichtigt)
rate[0] = 0.0 
ratez[0] = 0.0 
resid = 0.
#for i in range(1,3):
for i in range(1,nt):
    deltacf = cf[i]-(cf[i-1]-resid)
    nshift  = num.around(deltacf/deltaS,0).astype(num.int)
    resid = deltacf - nshift*deltaS
    chi = num.roll(chi,nshift)
    chiz = num.roll(chiz,nshift)
    if (nshift > 0):
       chi[0:nshift] = 1.0
       chiz[0:nshift] = 1.0
    elif (nshift < 0):
       chi[nshift:] = 0.0
       chiz[nshift:] = 0.0
    rate[i] = num.trapz(chi*p)*deltaS  
    ratez[i] = num.trapz(chiz*pz)*deltaS  

    print(i,nshift,deltacf,deltaS,' rate= ',rate[i])
    #print(chi)

    # ---- Anpassen von chi (memory effect)
    chi = chi*(1.-p)
    chiz = chiz*(1.-pz)

rate = rate*chi0/deltat
ratez = ratez*chi0/deltat

# ----- das selbe nochmal fuer stress step mit doppelter magnitude
rateb[0] = 0.0 
ratezb[0] = 0.0 
residb = 0.
for i in range(1,nt):
    deltacfb = cfb[i]-(cfb[i-1]-residb)
    nshiftb  = num.around(deltacfb/deltaS,0).astype(num.int)
    residb = deltacfb - nshiftb*deltaS
    chib = num.roll(chib,nshiftb)
    chizb = num.roll(chizb,nshiftb)
    if (nshiftb > 0):
       chib[0:nshiftb] = 1.0
       chizb[0:nshiftb] = 1.0
    elif (nshiftb < 0):
       chib[nshiftb:] = 0.0
       chizb[nshiftb:] = 0.0
    rateb[i] = num.trapz(chib*p)*deltaS  
    ratezb[i] = num.trapz(chizb*pz)*deltaS  

    # ---- Anpassen von chi (memory effect)
    chib = chib*(1.-p)
    chizb = chizb*(1.-pz)

rateb = rateb*chi0/deltat
ratezb = ratezb*chi0/deltat

# ----- das selbe nochmal fuer stress step mit kleinerer magnitude
ratec[0] = 0.0 
ratezc[0] = 0.0 
residc = 0.
for i in range(1,nt):
    #deltacfc = cfc[i]-(cfc[i-1]-residc)
    deltacfc = cfc[i]-(cfc[i-1]-residc)
    nshiftc  = num.around(deltacfc/deltaS,0).astype(num.int)
    residc = deltacfc - nshiftc*deltaS
    chic =  num.roll(chic, nshiftc)
    chizc = num.roll(chizc,nshiftc)
    if (nshiftc > 0):
       chic[0:nshiftc] = 1.0
       chizc[0:nshiftc] = 1.0
    elif (nshiftc < 0):
       chic[nshiftc:] = 0.0
       chizc[nshiftc:] = 0.0
    ratec[i]  = num.trapz(chic*p)*deltaS  
    ratezc[i] = num.trapz(chizc*pz)*deltaS  

    # ---- Anpassen von chi (memory effect)
    chic = chic*(1.-p)
    chizc = chizc*(1.-pz)

ratec = ratec*chi0/deltat
ratezc = ratezc*chi0/deltat

# ---------------- Traditionelle Berechnung mit if Funktion ----
S0 = -Sshadow
chi1[0] = cf[0]-Sshadow
for i in range(1,nt-1):
    if (cf[i] >= cf[i-1] and cf[i] >= S0):
        S0 = cf[i]
        rate1[i] = chi0*(cf[i]-cf[i-1])/deltat
    else:
        rate1[i] = 0.0/deltat
    chi1[i] = S0

# ---------------- Anzahl der Beben als Funktion der Zeit
c0 = 0.
neq[0]  = c0
neqb[0]  = c0
neqc[0]  = c0
neqz[0] = c0
neqzb[0] = c0
neqzc[0] = c0
for i in range(1,nt-2):
    neq[i] = num.trapz(rate[0:i+1])
    neqb[i] = num.trapz(rateb[0:i+1])
    neqc[i] = num.trapz(ratec[0:i+1])
    neqz[i] = num.trapz(ratez[0:i+1])
    neqzb[i] = num.trapz(ratezb[0:i+1])
    neqzc[i] = num.trapz(ratezc[0:i+1])

c0 = 0.
neq1[0] = c0
for i in range(1,nt-2):
    neq1[i] = num.trapz(rate1[0:i+1])


print(t)
print(len(cf))
print(cf)
print(rate1)

#----------------------------
# Ergebnisse plotten
#-----------------------------
Ta = deltat*(num.exp(1.) -1.)
scal = -depthS*chi0/deltat
#fig = plt.figure(1, figsize=(16, 12))
fig = plt.figure(1, figsize=(12, 8))
tmin = -5
tmax = 25

# ---- Coulomb stress loading (input and interpolated input
ax1a = fig.add_subplot(211)
#labelpos(ax1a, 2, 1.5)
#ax1a.set_title('n1/n2='+str(n1)+'/'+str(n2)+', S0/S1/S2/S3='+str(sc0)+'/'+str(sc1)+'/'+str(sc2)+'/'+str(sc3)+r', $\chi=$'+str(chi0)+r'$Pa^{-1}$', fontsize=22)
#ax1a.set_xlabel("$t$ (hours)", fontsize=20)
#ax1a.set_ylabel('$\sigma_c$ (Pa)', fontsize=20)
ax1a.set_ylabel('$\sigma_c/\dot{\sigma}_c\Delta t$', fontsize=20)
ax1a.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1a.tick_params(axis = 'both', which = 'minor', labelsize = 18)

fa = (sc2-sc1)/(tectstressrate*deltat)
fb = (sc2b-sc1)/(tectstressrate*deltat)
fc = (sc2c-sc1)/(tectstressrate*deltat)
#ax1a.plot((t-t[n1])/Ta, chi1, linewidth=2.0, ls='-',color='darkgray',marker='o',markersize=2.5)
ax1a.plot((t-t[n1])/Ta, cfc/(tectstressrate*deltat), linewidth=2.0, ls='-.', color='black', label=r'$\Delta\sigma/\dot{\sigma}_c\Delta t=$'+'{:d}'.format(num.int(fc)))
ax1a.plot((t-t[n1])/Ta, cf/(tectstressrate*deltat), linewidth=2.0, ls='-', color='black', label=r'$\Delta\sigma/\dot{\sigma}_c\Delta t=$'+'{:d}'.format(num.int(fa)))
ax1a.plot((t-t[n1])/Ta, cfb/(tectstressrate*deltat), linewidth=2.0, ls='--', color='black', label=r'$\Delta\sigma/\dot{\sigma}_c\Delta t=$'+'{:d}'.format(num.int(fb)))

plt.legend(loc='lower right',fontsize=20)
plt.xlim([tmin, tmax])

# ---- Coulomb stress change (at mid sample positions)
ax1b = fig.add_subplot(212)
fa = -depthS/(tectstressrate*deltat)
fb = 0.0
#ax1b.set_title("earthquake rate", fontsize=22)
#ax1b.set_xlabel("$t$ (hours)", fontsize=20)
#ax1b.set_ylabel(r'$r \, (1/hours)$', fontsize=20)
ax1b.set_xlabel(r'$(t-t_0)/\Delta t(e-1)$', fontsize=20)
#ax1b.set_ylabel(r'$r \cdot \Delta t / \chi_0 \delta \sigma$', fontsize=20)
ax1b.set_ylabel(r'$r / \chi_0 V \delta \dot{\sigma}$', fontsize=20)
ax1b.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

#ax1b.plot((t[0:nt-1]-t[n1])/Ta, rate1[0:nt-1]/scal, linewidth=1.0, color='gray', ls='--',)
ax1b.plot((t[0:nt-1]-t[n1])/Ta, rate[0:nt-1]/scal, linewidth=1.5, color='lightgray',  label=r'$\delta \dot{\sigma}/\dot{\sigma}_c=$'+'{:.1f}'.format(num.int(fb)))
ax1b.plot((t[40:nt-1]-t[n1])/Ta, ratez[40:nt-1]/scal, linewidth=3.5, color='red', label=r'$\delta \dot{\sigma}/\dot{\sigma}_c=$'+'{:d}'.format(num.int(fa)))

ax1b.plot((t[40:nt-1]-t[n1])/Ta, ratezb[40:nt-1]/scal, linewidth=2.0, color='red', ls='--')
ax1b.plot((t[40:nt-1]-t[n1])/Ta, ratezc[40:nt-1]/scal, linewidth=2.0, color='red', ls='-.')

plt.legend(loc='upper right',fontsize=20)
plt.figtext(0.02, 0.87, 'a)', fontsize=20)
plt.xlim([tmin, tmax])
ax1b.set_yscale('log')

# ----- plot and save figure ----
fig.savefig(file+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(file+'.png', format='png',resolution=100.0, oversample=2.0, size=None, width=None, height=None, psconvert=True, bbox_inches='tight')

plt.show()

# ------- plot absolute number of earthquakes for own and Coulomb failure models
fig = plt.figure(1, figsize=(12, 4))
ax1c = fig.add_subplot(111)
#ax1c.set_title("earthquake rate", fontsize=22)
#ax1c.set_xlabel("$t$ (hours)", fontsize=20)
#ax1c.set_ylabel(r'$r \, (1/hours)$', fontsize=20)
ax1c.set_xlabel(r'$(t-t_0)/\Delta t(e-1)$', fontsize=20)
ax1c.set_ylabel(r'$n-n[-1]$', fontsize=20)
ax1c.tick_params(axis = 'both', which = 'major', labelsize = 20)
ax1c.tick_params(axis = 'both', which = 'minor', labelsize = 18)

ax1c.plot((t[40:nt-1]-t[n1])/Ta, neqz[40:nt-1]-neqz[n1-1], linewidth=2.5, color='red', label=r'$\delta \dot{\sigma}/\dot{\sigma}_c=$'+'{:d}'.format(num.int(fa    )))

ax1c.plot((t[40:nt-1]-t[n1])/Ta, neqzb[40:nt-1]-neqzb[n1-1], linewidth=2.5, color='red', ls='--')
ax1c.plot((t[40:nt-1]-t[n1])/Ta, neqzc[40:nt-1]-neqzc[n1-1], linewidth=2.5, color='red', ls='-.')

#ax1c.plot((t[0:nt-1]-t[n1])/Ta, rate1[0:nt-1]/scal, linewidth=1.0, color='gray', ls='--',)
ax1c.plot((t[40:nt-1]-t[n1])/Ta, neq[40:nt-1]-neq[n1-1], linewidth=1.5, color='gray',  label=r'$\delta \dot{\sigma}/\dot{\sigma}_c=$'+'{:.1f}'.format(num.    int(fb)))
ax1c.plot((t[40:nt-1]-t[n1])/Ta, neqb[40:nt-1]-neqb[n1-1], linewidth=1.5, color='gray', ls='--')
ax1c.plot((t[40:nt-1]-t[n1])/Ta, neqc[40:nt-1]-neqc[n1-1], linewidth=1.5, color='gray', ls='-.')
#ax1c.plot((t[40:nt-1]-t[n1])/Ta, ratec[40:nt-1]-ratec[n1-1], linewidth=1.5, color='gray', ls='-.')

plt.legend(loc='lower right',fontsize=20)
plt.figtext(0.02, 0.87, 'a)', fontsize=20)
#tmin = -5
#tmax = 25
plt.xlim([tmin, tmax])
#ax1c.set_yscale('log')
plt.show()

# ----- plot and save figure ----
fig.savefig(file2+'.pdf', format='pdf',bbox_inches='tight')
fig.savefig(file2+'.png', format='png',resolution=250.0, oversample=2.0, size=None, width=None, height=None, psconvert=True, bbox_inches='tight')

#cf = cfun(kappa, eta, lamb_u, lamb, alpha, shear, scal)
