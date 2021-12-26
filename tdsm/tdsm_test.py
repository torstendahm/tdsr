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

# ----------------- physical model parameter ----------------------------
chi0 = 10.        # susceptibility
depthS = -0.2     # Tiefe der Einflussfunktion 

# ----------------- parameter to modify starting state
Sshadow = 0.0    # initial stress shadow , usually set to zero

# ----------------- discretisation parameter -----------------------------
deltat = 0.2*hours  # sampleinterval in t in seconds
deltaS = 0.01       # sampling of Coulomb stress axis
sigma_max = 10.0

# ----------------- stress loading models -------------
# 1. step function overlaying steady stressing rate
model1 = True
if model1: 
    n1 = 50          # fist phase must be larger at least 2
    n2 = 51          # second phase at least n1+1, better n1+2 as central differences taken
    sc0 = 0.0        # Coulomb stress values
    sc1 = 0.5
    sc2 = 1.0
    sc3 = 2.0
    # In dem Beispiel gibt es 150 Samples (30 /0.2 ). Tectonic rate wird aus ersten 50 sample berechnet, 
    tectstressrate = (sc1-sc0)/(n1*deltat)

# ----------------- model to be used for forecast -- only one of three can be true  -------------
tdsm = True
#tdsm = False

#lcm  = True
lcm  = False

#traditional = True
traditional = False

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
#chi = num.zeros(nsigma)
chiz = num.zeros(nsigma)
#p   = num.zeros(nsigma)
pz  = num.zeros(nsigma)

# ------------------ stress loading function cf definieren ----
tmin, tmax, nt, t = gridrange(0*hours, 30*hours, deltat)
cf    = num.zeros(nt)
ratez  = num.zeros(nt)
chi1  = num.zeros(nt)
neqz  = num.zeros(nt-1)

for i in range(0,n1):
    cf[i] = sc0 + i*(sc1-sc0)/float(n1-1)
for i in range(n1,n2):
    cf[i] = cf[n1-1] + (i-n1+1)*(sc2-sc1)/float(n2-n1)
for i in range(n2,nt):
    cf[i] = cf[n2-1] + (i-n2+1)*(sc3-sc2)/float(nt-n2)

pz   = num.heaviside( sigma , 1 )

#------------------------------------------------
# ---------------- Neue Ansatz fuer Bebenrate, deltan=rate, durch Integration ----
#------------------------------------------------
if tdsm:
    # ---- probability function wird nach unten modifiziert durch boundary layer influence zone
    ndepth  = num.around(depthS/deltaS,0).astype(num.int)
    nzero = int(nsigma/2)
    for k in range(0,-6*ndepth):
        pz[nzero+6*ndepth+k] = num.exp((-6*ndepth-k)/ndepth) # einfache e Abfall

    print('nzero=',nzero,' ndepth=',ndepth)

    chiz = num.heaviside( -sigma , 0 )
    print(sigma)
    print(pz)

    # ...... initiale stress shadow falls noetig
    nshift  = num.around(Sshadow/deltaS,0).astype(num.int)
    chiz = num.roll(chiz,nshift)
    if (nshift > 0):
       chiz[0:nshift] = 1.0
    elif (nshift < 0):
       chiz[nshift:] = 0.0

    # ...... incremental approach  (Achtung: Rundungsfehler bei shifts testweise durch resid beruecksichtigt)
    ratez[0] = 0.0 
    resid = 0.
    for i in range(1,nt):
        deltacf = cf[i]-(cf[i-1]-resid)
        nshift  = num.around(deltacf/deltaS,0).astype(num.int)
        resid = deltacf - nshift*deltaS
        chiz = num.roll(chiz,nshift)
        if (nshift > 0):
           chiz[0:nshift] = 1.0
        elif (nshift < 0):
           chiz[nshift:] = 0.0
        ratez[i] = num.trapz(chiz*pz)*deltaS  
    
        # ---- Anpassen von chi (memory effect)
        chiz = chiz*(1.-pz)
    
    ratez = ratez*chi0/deltat

#------------------------------------------------
# ------------- linear Coulomb failure model ---
#------------------------------------------------
if lcm:
    chiz = num.heaviside( -sigma , 0 )
    print(sigma)
    print(pz)

    # ...... initiale stress shadow falls noetig
    nshift  = num.around(Sshadow/deltaS,0).astype(num.int)
    chiz = num.roll(chiz,nshift)
    if (nshift > 0):
       chiz[0:nshift] = 1.0
    elif (nshift < 0):
       chiz[nshift:] = 0.0

    # ...... incremental approach  (Achtung: Rundungsfehler bei shifts testweise durch resid beruecksichtigt)
    ratez[0] = 0.0 
    resid = 0.
    for i in range(1,nt):
        deltacf = cf[i]-(cf[i-1]-resid)
        nshift  = num.around(deltacf/deltaS,0).astype(num.int)
        resid = deltacf - nshift*deltaS
        chiz = num.roll(chiz,nshift)
        if (nshift > 0):
           chiz[0:nshift] = 1.0
        elif (nshift < 0):
           chiz[nshift:] = 0.0
        ratez[i] = num.trapz(chiz*pz)*deltaS  
    
        # ---- Anpassen von chi (memory effect)
        chiz = chiz*(1.-pz)
    
    ratez = ratez*chi0/deltat


#------------------------------------------------
# ------------- linear Coulomb failure model auf traditionelle Art berechnet ---
#------------------------------------------------
if traditional:
    # ---------------- Traditionelle Berechnung mit if Funktion ----
    S0 = -Sshadow
    chi1[0] = cf[0]-Sshadow
    for i in range(1,nt-1):
        if (cf[i] >= cf[i-1] and cf[i] >= S0):
            S0 = cf[i]
            ratez[i] = chi0*(cf[i]-cf[i-1])/deltat
        else:
            ratez[i] = 0.0/deltat
        chiz[i] = S0

# ---------------- Anzahl der Beben als Funktion der Zeit
c0 = 0.
neqz[0] = c0
for i in range(1,nt-2):
    neqz[i] = num.trapz(ratez[0:i+1])

print(t)
print(len(cf))
print(cf)

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
                 sc0 = sc0,
                 sc1 = sc1,
                 sc2 = sc2,
                 sc3 = sc3,
                 tectstressrate = tectstressrate
                 ),
     data = [t, cf, ratez, neqz]
        )

fn_out = 'results/tdsm.out'
pickle.dump(out, open(fn_out, 'wb'))
