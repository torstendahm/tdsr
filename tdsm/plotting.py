################################
# Time Dependent Seismicity Model - Plotting results
# T. Dahm, R. Dahm 26.12.2021
################################

def plot(config, t, cf, ratez, neqz):
    from matplotlib import pyplot as plt
    plt.plot(cf)
    plt.show()
    plt.close()
    plt.plot(ratez)
    plt.show()
    plt.close()
    plt.plot(neqz)
    plt.show()
    plt.close()

def plot_paper(config, t, cf, ratez, neqz):
    import math
    import numpy as num
    import pickle
    import matplotlib
    from matplotlib import pyplot as plt
    import matplotlib as mpl
    import matplotlib.colors as mpl_colors
    import matplotlib.cm as mpl_cm
    from matplotlib.colors import BoundaryNorm
    from matplotlib.ticker import MaxNLocator

    nt = len(t)
    Ta = config.deltat*(num.exp(1.) -1.)
    scal = -config.depthS*config.chi0/config.deltat
    #fig = plt.figure(1, figsize=(16, 12))
    fig = plt.figure(1, figsize=(12, 8))
    tmin = -5
    tmax = 25

    # ---- Coulomb stress loading (input and interpolated input
    ax1a = fig.add_subplot(211)
    ax1a.set_ylabel('$\sigma_c/\dot{\sigma}_c\Delta t$', fontsize=20)
    ax1a.tick_params(axis = 'both', which = 'major', labelsize = 20)
    ax1a.tick_params(axis = 'both', which = 'minor', labelsize = 18)
    
    tectstressrate = config.loading.stress_rate
    fa = (config.loading.sc2-config.loading.sc1)/(tectstressrate * config.deltat)
    ax1a.plot((t-t[config.loading.n1])/Ta, cf/(tectstressrate*config.deltat), linewidth=2.0, ls='-', color='black', label=r'$\Delta\sigma/\dot{\sigma}_c\Delta t=$'+'{:d}'.format(num.int(fa)))

    plt.legend(loc='lower right',fontsize=20)
    plt.xlim([tmin, tmax])

    # ---- Coulomb stress change (at mid sample positions)
    ax1b = fig.add_subplot(212)
    fa = -config.depthS/(tectstressrate*config.deltat)
    fb = 0.0
    ax1b.set_xlabel(r'$(t-t_0)/\Delta t(e-1)$', fontsize=20)
    ax1b.set_ylabel(r'$r / \chi_0 V \delta \dot{\sigma}$', fontsize=20)
    ax1b.tick_params(axis = 'both', which = 'major', labelsize = 20)
    ax1b.tick_params(axis = 'both', which = 'minor', labelsize = 18)

    ax1b.plot((t[40:nt-1]-t[config.loading.n1])/Ta, ratez[40:nt-1]/scal, linewidth=3.5, color='red', label=r'$\delta \dot{\sigma}/\dot{\sigma}_c=$'+'{:d}'.format(num.int(fa)))

    plt.legend(loc='upper right',fontsize=20)
    plt.figtext(0.02, 0.87, 'a)', fontsize=20)
    plt.xlim([tmin, tmax])
    ax1b.set_yscale('log')

    # ----- plot and save figure ----
    # fig.savefig('plots/'+file+'.pdf', format='pdf',bbox_inches='tight')
    # fig.savefig('plots/'+file+'.png', format='png',resolution=100.0, oversample=2.0, size=None, width=None, height=None, psconvert=True, bbox_inches='tight')

    plt.show()

    # ------- plot absolute number of earthquakes for own and Coulomb failure models
    fig = plt.figure(1, figsize=(12, 4))
    ax1c = fig.add_subplot(111)
    ax1c.set_xlabel(r'$(t-t_0)/\Delta t(e-1)$', fontsize=20)
    ax1c.set_ylabel(r'$n-n[-1]$', fontsize=20)
    ax1c.tick_params(axis = 'both', which = 'major', labelsize = 20)
    ax1c.tick_params(axis = 'both', which = 'minor', labelsize = 18)

    ax1c.plot((t[40:nt-1]-t[config.loading.n1])/Ta, neqz[40:nt-1]-neqz[config.loading.n1-1], linewidth=2.5, color='red', label=r'$\delta \dot{\sigma}/\dot{\sigma}_c=$'+'{:d}'.format(num.int(fa    )))

    plt.legend(loc='lower right',fontsize=20)
    plt.figtext(0.02, 0.87, 'a)', fontsize=20)
    plt.xlim([tmin, tmax])
    plt.show()

    # ----- plot and save figure ----
    # fig.savefig('plots/'+file2+'.pdf', format='pdf',bbox_inches='tight')
    # fig.savefig('plots/'+file2+'.png', format='png',resolution=250.0, oversample=2.0, size=None, width=None, height=None, psconvert=True, bbox_inches='tight')
