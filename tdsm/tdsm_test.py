################################
# Time Dependent Seismicity Model - Plotting results
# T. Dahm, R. Dahm 26.12.2021
################################
'''
    Add description here

    Input:
    return:  tau_rake, tau_plane, tau_max, sign, sign_max 
'''

import math
import numpy as num
from scipy import special as spl
import pickle

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
deltat = 0.2*hours  # sampleinterval in t in seconds
n1 = 50          # fist phase must be larger at least 2
n2 = 51          # second phase at least n1+1, better n1+2 as central differences taken
sc0 = 0.0        # Coulomb stress values
sc1 = 0.5
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

def gridrange(amin, amax, astep, rounding=math.ceil):
    amin = float(amin)
    amax = float(amax)
    astep = float(astep)
    na = int(rounding((amax-amin)/astep))
    amax = amin + (na-1) * astep
    a = num.linspace(amin, amax, na)
    return amin, amax, na, a

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
    cfb[i] = cf[i]
    cfc[i] = cf[i]
for i in range(n1,n2):
    cf[i] = cf[n1-1] + (i-n1+1)*(sc2-sc1)/float(n2-n1)
    cfb[i] = cfb[n1-1] + (i-n1+1)*(sc2b-sc1)/float(n2-n1)
    cfc[i] = cfc[n1-1] + (i-n1+1)*(sc2c-sc1)/float(n2-n1)
for i in range(n2,nt):
    cf[i] = cf[n2-1] + (i-n2+1)*(sc3-sc2)/float(nt-n2)
    cfb[i] = cfb[n2-1] + (i-n2+1)*(sc3b-sc2b)/float(nt-n2)
    cfc[i] = cfc[n2-1] + (i-n2+1)*(sc3c-sc2c)/float(nt-n2)

# ---------------- Neue Ansatz fuer Bebenrate, deltan=rate, durch Integration ----
p   = num.heaviside( sigma , 1 )
pz   = num.heaviside( sigma , 1 )
# ---- probability function wird nach unten modifiziert durch boundary layer influence zone
ndepth  = num.around(depthS/deltaS,0).astype(num.int)
nzero = int(nsigma/2)
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

# ----- write meta data 
# dump all variables and time series
out = dict(
     meta = dict(
                 chi0 = chi0,
                 Sshadow = Sshadow,
                 deltaS = deltaS,
                 sigma_max = sigma_max,
                 depthS = depthS,
                 nt = nt, 
                 n1 = n1, 
                 deltat = deltat,
                 sc1 = sc1,
                 sc2 = sc2,
                 sc2b = sc2b, 
                 sc2c = sc2c, 
                 tectstressrate = tectstressrate
                 ),
     data = [t, cf, cfb, cfc, rate, ratez, ratezb, ratezc, neq, neqb, neqc, neqz, neqzb, neqzc]
        )

fn_out = 'results/tdsm.out'
pickle.dump(out, open(fn_out, 'wb'))
