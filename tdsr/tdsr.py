"""
TDSR (Time Dependent Stress Response) seismicity models
====================================
TDSR is a python tool to simulate time dependent earthquake rates for
different stress loading scenarios.
Earthquake rate is definded as number of events per time unit with
magnitudes above a completeness magnitude.
Stress loading is understood as Coulomb stress as a function of time.
The stress loading is assumed homogeneous within the rock volume for which
the simulation is performed. The details of the TDSR seismicity model
and examles are described in Dahm and Hainzl (2022), JGR, http://doi.org/10.1029/2022JB024443  . 

Different loading scenarios are supported and defined by loading
classes which are imported from tdsr.loading.
Loading classes include

    * step in Coulomb stress :class:`tdsr.loading.StepLoading`
    * a constant background stress rate :class:`tdsr.loading.BackgroundLoading`
    * a change in stress rate :class:`tdsr.loading.TrendchangeLoading`
    * a cyclic stress rate superposed to a constant background trend
      :class:`tdsr.loading.CyclicLoading`
    * a stress curve defined by 4 points :class:`tdsr.loading.FourPointLoading`
    * or a ramp like loading scenario :class:`tdsr.loading.RampLoading`

Alternatively, the loading time function can be read from an
external file with :class:`tdsr.loading.CustomLoading`.
Stress loading functions can also be defined in the calling python
script and directly passed to the seismicity model classes
(example not yet provided).

Additional to simulations of the new TDSR model, the toolbox can
compare results with earthquake rates calculated from other seismicity
models, using the same loading functions and parameter settings as TDSR.
All models are defined as classes and can be loaded from tdsr
(e.g. TDSR1 for the TDSR model, provided by tdsr.tdsr.TDSR1).
Other classes included so far are rate and state
seismicity models (RSD, RSD1) and linear Coulomb failure
models (CFM, Traditional).

An elementary plotting class is supported, which is imported from
tdsr.plotting. However, the examples provided in the toolbox use
their own plotting solutions.

Default settings of input parameter can be defined in config.toml.
Default values can be overwritten when calling the seismicity model
classes or  loading scenarios (see examples).

The examples collected in the subdirectory ''examples'' reproduce figures
published in Dahm and Hainzl (2022).
All theory of the TDSR model and the examples are described in
Dahm and Hainzl (2022).
The software toolbox in python can be cited as
Dahm, T.; Hainzl, S.; Dahm, R.A. (2022): Time-dependent stress response seismicity models (TDSR). GFZ Data Services. https://doi.org/10.5880/GFZ.2.1.2022.002 .
The TDSR  github project here, however, is intended as an open source
project that may change and improve over the years.
Please send us your suggestions or addings and contribute to the
git project if interested.
"""

from copy import deepcopy
from typing import Optional, Tuple

import numpy as np
import numpy.typing as npt

from tdsr.config import Config
from tdsr.loading import Loading
from tdsr.exceptions import MissingParameter
from tdsr.utils import (
    DEBUG,
    X0gaussian,
    X0steady,
    X0uniform,
    Zvalues,
    gridrange,
    gridrange_log,
    pf,
    shifted,
)

Result = Tuple[
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
]


class LCM(object):
    """
    Class of the Linear Coulomb Failure Model (LCM) used as base
    for calculating time dependent seismicity with class TDSR.
    Note that TDSR is outdated and was replaced by TDSR1.
    For simulations with the linear Coulomb Failure model use
    classes "CFM" or "Traditional".
    """

    def __init__(self, config: Optional[Config] = None) -> None:
        """
        Add description here

        Parameters
        ---------
        config
            Optional config to use
        """
        self.config = config or Config()

    def __call__(
        self,
        chi0: Optional[float] = None,
        t0: Optional[float] = None,
        chiz: Optional[float] = None,
        depthS: Optional[float] = None,
        Sshadow: Optional[float] = None,
        iX0: Optional[str] = None,
        Zmean: Optional[float] = None,
        Zstd: Optional[float] = None,
        equilibrium: Optional[bool] = None,
        deltat: Optional[float] = None,
        tstart: Optional[float] = None,
        tend: Optional[float] = None,
        taxis_log: Optional[bool] = None,
        ntlog: Optional[int] = None,
        deltaS: Optional[float] = None,
        sigma_max: Optional[int] = None,
        precision: Optional[int] = None,
        loading: Optional[Loading] = None,
    ) -> Result:
        config = deepcopy(self.config)
        config.merge(
            dict(
                chi0=chi0,
                t0=t0,
                depthS=depthS,
                Sshadow=Sshadow,
                iX0=iX0,
                Zmean=Zmean,
                Zstd=Zstd,
                equilibrium=equilibrium,
                deltat=deltat,
                tstart=tstart,
                tend=tend,
                taxis_log=taxis_log,
                ntlog=ntlog,
                deltaS=deltaS,
                sigma_max=sigma_max,
                precision=precision,
            )
        )
        if loading is not None:
            config.loading = loading
        # if config.equilibrium:
        if chiz is not None:
            self._prepare(config)
            chiz_background = chiz
            # print('chiz_background=',chiz_background)
            self.chiz = chiz_background
        else:
            # print('chiz is none')
            self._prepare(config)
        return self._compute(config)

    def _prepare(self, config: Config) -> None:
        (self.smin, self.smax, self.nsigma, self.sigma, self.dZ) = gridrange(
            -config.sigma_max, +config.sigma_max, config.deltaS
        )
        if config.taxis_log:
            (self.tmin, self.tmax, self.nt, self.t, self.dt) = gridrange_log(
                config.tstart, config.tend, config.ntlog
            )
        else:
            (self.tmin, self.tmax, self.nt, self.t, self.dt) = gridrange(
                config.tstart, config.tend, config.deltat
            )

        self.chiz = np.zeros(self.nsigma)
        #  self.pz will be overridden by TDSR subclass
        self.pz = np.heaviside(self.sigma, 1)

        loading = config.loading
        if loading is None:
            raise MissingParameter("missing loading function")
        self.cf = loading.values(length=self.nt)

    def _compute(self, config: Config) -> Result:
        if config.equilibrium:
            # chiz = np.heaviside(-self.sigma, 0)
            # raise ValueError("the option to calculate or read the equilibrium function chiz is not yet implemented")
            print("chiz bereits uebergeben, daher nicht ueberschrieben mit Heaviside")
        else:
            self.chiz = np.heaviside(-self.sigma, 0)

        # nshift = np.around(config.Sshadow / config.deltaS, 0).astype(int)
        nshift = np.around(-1.0 * config.Sshadow / config.deltaS, 0).astype(int)
        self.chiz = shifted(self.chiz, nshift)

        ratez = np.zeros(self.nt)
        ratez[0] = 0.0
        resid = 0.0
        for i in range(1, self.nt):
            # optional to be changed: config.deltaS may be replaced by self.dZ[i] if sigma axis not discretized with equal sampling (see gridrange output)
            deltacf = self.cf[i] - (self.cf[i - 1] - resid)
            nshift = np.around(deltacf / config.deltaS, 0).astype(int)
            resid = deltacf - nshift * config.deltaS
            # shift chiz (memory effect)
            self.chiz = shifted(self.chiz, nshift)

            # ratez[i] = np.trapz(self.chiz * self.pz) * config.deltaS
            ratez[i] = (
                np.trapz(config.chi0 * self.chiz * self.pz * config.deltaS) / self.dt[i]
            )
            # cut off chiz
            self.chiz = self.chiz * (1.0 - self.pz)

        # ratez = ratez * config.chi0 / config.deltat
        # ratez = ratez
        neqz = np.zeros(self.nt - 1)
        # neqz[0] = 0.0
        # replace loop by numpy.cumsum to directly calculate  the cumulative sum
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
            # neqz[i] = np.trapz(ratez[0 : i + 1] , self.t[0:i+1])  # 24 july 2022: beruecksichtigt, dass zeitreihe nicht gelichabstaendig sein kann
        return self.t, self.chiz, self.cf, ratez, neqz


class TDSR(LCM):
    """
    TDSR class realisation builts on LCM and is deprecated.
    Use TDSR1 instead.
    """

    def _compute(self, config: Config) -> Result:
        # use exponential decay for pz
        ndepth = np.around(config.depthS / config.deltaS, 0).astype(int)
        nzero = int(self.nsigma / 2)
        window = int(config.precision * ndepth)
        self.pz[nzero + window : nzero] = np.exp(-np.arange(window, 0) / ndepth)
        return super()._compute(config)


class TDSR1(object):
    """
    TDSR1 simulates the time dependent stress rate seismicity model TDSR.
    Key parameters are

        * ``chi0``: susceptibility to trigger earthquakes after unit
          step increase
        * ``t0``: mean time to failure for critically stressed sources
        * ``depthS``: skin depth parameter in exponential trigger functions

    Different source distributions ``iX0`` can be considered:

        * If ``iX0=="uniform"``,  a uniform disribution of stress stages
          at seismic sources is assumed.
        * If ``iX0=="gaussian"`` a gaussian distribution of stress stages
          is assumed. The mean stress level of the Gaussian peak is defined
          by ``Zmean`` and the standard deviation by ``Zstd``.
        * If ``iX0switch=="equilibrium"`` a steady state distribution
          of sources is assumed, using ``config.strend``.

    All three cases of stress distributions can be modified by a shift of
    ``Sshadow`` on the stress axis to simulate a subcritical stress state.
    If ``taxis_log==False`` an equally space time sampling is assumed
    with interval ``deltat``.
    """

    def __init__(self, config: Optional[Config] = None) -> None:
        """
        Add description here

        Parameters
        ---------
        config
            Optional config to use
        """
        self.config = config or Config()

    def __call__(
        self,
        chi0: Optional[float] = None,
        t0: Optional[float] = None,
        chiz: Optional[float] = None,
        depthS: Optional[float] = None,
        Sshadow: Optional[float] = None,
        iX0: Optional[str] = None,
        Zmean: Optional[float] = None,
        Zstd: Optional[float] = None,
        equilibrium: Optional[bool] = None,
        deltat: Optional[float] = None,
        tstart: Optional[float] = None,
        tend: Optional[float] = None,
        taxis_log: Optional[bool] = None,
        ntlog: Optional[int] = None,
        deltaS: Optional[float] = None,
        sigma_max: Optional[int] = None,
        precision: Optional[int] = None,
        loading: Optional[Loading] = None,
    ) -> Result:
        config = deepcopy(self.config)
        config.merge(
            dict(
                chi0=chi0,
                t0=t0,
                depthS=depthS,
                Sshadow=Sshadow,
                iX0=iX0,
                Zmean=Zmean,
                Zstd=Zstd,
                equilibrium=equilibrium,
                deltat=deltat,
                tstart=tstart,
                tend=tend,
                taxis_log=taxis_log,
                ntlog=ntlog,
                deltaS=deltaS,
                sigma_max=sigma_max,
                precision=precision,
            )
        )
        if loading is not None:
            config.loading = loading
        # if config.equilibrium:
        if chiz is not None:
            self._prepare(config)
            chiz_background = chiz
            # print('chiz_background=',chiz_background)
            self.chiz = chiz_background
        else:
            # print('chiz is none')
            self._prepare(config)
        return self._compute(config)

    def _prepare(self, config: Config) -> None:
        (self.smin, self.smax, self.nsigma, self.sigma, self.dZ) = gridrange(
            -config.sigma_max, +config.sigma_max, config.deltaS
        )
        #print("taxis_log=", config.taxis_log, " ntlog=", config.ntlog)
        if config.taxis_log:
            (self.tmin, self.tmax, self.nt, self.t, self.dt) = gridrange_log(
                config.tstart, config.tend, config.ntlog
            )
        else:
            (self.tmin, self.tmax, self.nt, self.t, self.dt) = gridrange(
                config.tstart, config.tend, config.deltat
            )
        # self.chiz = np.zeros(self.nsigma)
        # self.pz = np.heaviside(self.sigma, 1)
        loading = config.loading
        if loading is None:
            raise MissingParameter("missing loading function")
        self.cf = loading.values(length=self.nt)

    def _compute(self, config: Config) -> Result:
        ratez = np.zeros(self.nt)
        # ratez[0] = 0.0
        dsig = -config.depthS
        t0 = config.t0
        X0 = config.chi0
        # nicht gut geloest, macht nur Sinn wenn Background
        # mit step verwendet wird (tstep=0)
        if config.taxis_log:
            # Zmin = config.Sshadow
            # nicht gut geloest, sstep bei tstep=0 sollte mit sstep
            # definiert werden - noch zu aendern
            Zmin = config.loading.sstep
        else:
            Zmin = 0.0
        # print('Zmin ',Zmin,' config.loading.strend=',config.loading.strend,' chi0=',config.chi0,' t0=',t0,' X0=',X0,' dsig=',dsig)
        # dt = np.ediff1d(self.t, to_end=self.t[-1]-self.t[-2])  # wird bereits in  gridrange berechnet
        # Z = functions.Zvalues(self.cf, 0, t0, dsig) # t0 kann  raus, da nicht benutzt
        # Z = Zvalues(self.cf, 0.0, 0.0, dsig)
        Z = Zvalues(self.cf, Zmin, 0.0, dsig)
        # Z = 0.04*Z  # falls es mit dem zu grossen Range und der groben Diskretisierung der zeta Achse zu Problemen kommt
        dZ = np.ediff1d(Z, to_end=Z[-1] - Z[-2])
        # print('smin=',np.amin(self.cf),' smax=',np.amax(self.cf),' ns=',len(self.cf))
        # print('zmin=',np.amin(Z),' zmax=',np.amax(Z),' nz=',len(Z))
        # print('zvalues ',Z)
        dS = np.ediff1d(self.cf, to_end=self.cf[-1] - self.cf[-2])

        if config.iX0.lower() == "equilibrium":
            # steady state equilibrium before loading starts
            r0 = config.chi0 * config.loading.strend
            X = X0steady(Z + Zmin, r0, config.t0, -config.depthS, config.loading.strend)

        elif config.iX0.lower() == "uniform":
            # uniform Distribution of stress states (e.g. Fig. 3a)
            X = X0uniform(Z, config.Sshadow, config.chi0)

        elif config.iX0.lower() == "gaussian":
            # gaussian distribution before loading starts
            if DEBUG:
                print(
                    "Gaussian distribution with " "Zmean=",
                    config.Zmean,
                    " Zstd=",
                    config.Zstd,
                )
            X = X0gaussian(Z + Zmin, config.Zmean, config.Zstd, config.chi0)

        else:
            raise InvalidParameter(
                'iX0 must be of of "equilibrium", "uniform", "gaussian", ',
                "but got ",
                config.iX0.lower(),
            )

        self.chiz = X
        # print('xmin=',np.amin(X),' xmax=',np.amax(X),' nx=',len(X))
        # print('i=0 1/pf=',1./pf(Z, config.t0, -config.depthS)[0:3],' ... ',1./pf(Z, config.t0, -config.depthS)[-3:-1])
        # for i in range(self.nt):
        # for i in range(1, self.nt):
        for i in range(self.nt):
            # dX = self.chiz/ tf(Z, t0, config.deltaS) * self.dt[i]
            # dX = X / tf(Z, t0, -config.depthS) * self.dt[i]
            # dX = X / tf(Z, config.t0, -config.depthS) * self.dt[i]
            dX = X * pf(Z, config.t0, -config.depthS) * self.dt[i]
            dX[(dX > X)] = X[
                (dX > X)
            ]  # wenn diese Zeile entfaellt, dann muss nicht mit dt multipliziert werden
            ratez[i] = (
                np.sum(dX * dZ) / self.dt[i]
            )  # oben wird dX mit dt multipliziert, hier dividiert.
            # ratez[i] = np.trapz(dX * dZ) / self.dt[i]  # oben wird dX mit dt multipliziert, hier dividiert.
            # print('len(Z)=',len(Z),' len(dS)=',len(dS),' nt=',self.nt)
            Z -= dS[i]
            X -= dX

        # print('i=nt-1 1/pf=',1./pf(Z, config.t0, -config.depthS)[0:3],' ... ',1./pf(Z, config.t0, -config.depthS)[-3:-1])
        # print('rmin=',np.amin(ratez),' rmax=',np.amax(ratez),' nx=',len(ratez))
        neqz = np.zeros(self.nt - 1)
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return self.t, self.chiz, self.cf, ratez, neqz


class Traditional(LCM):
    """Linear Coulomb Failure Model (LCM) realisation using a simple (traditional) approach."""

    def _compute(self, config: Config) -> Result:
        ratez = np.zeros(self.nt)
        cf_shad = np.zeros(self.nt)
        S0 = +config.Sshadow
        cf_shad[0] = self.cf[0] - config.Sshadow
        for i in range(1, self.nt - 1):
            if self.cf[i] >= self.cf[i - 1] and self.cf[i] >= S0:
                S0 = self.cf[i]
                ratez[i] = self.cf[i] - self.cf[i - 1]
            else:
                ratez[i] = 0.0
            cf_shad[i] = S0
        ratez = ratez * config.chi0 / config.deltat
        neqz = np.zeros(self.nt - 1)
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return self.t, cf_shad, self.cf, ratez, neqz


class CFM(LCM):
    """The class Coulomb Failure Model (CFM ) is ismilar to Traditional but coded different."""

    def _compute(self, config: Config) -> Result:
        ratez = np.zeros(self.nt)
        cf_shad = np.zeros(self.nt)
        S0 = +config.Sshadow
        cf_shad[0] = self.cf[0] - config.Sshadow
        for i in range(1, self.nt - 1):
            if self.cf[i] >= self.cf[i - 1] and self.cf[i] >= S0:
                S0 = self.cf[i]
                ratez[i] = self.cf[i] - self.cf[i - 1]
            else:
                ratez[i] = 0.0
            cf_shad[i] = S0
        ratez = ratez * config.chi0 / config.deltat
        neqz = np.zeros(self.nt - 1)
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return self.t, cf_shad, self.cf, ratez, neqz


class RSM(LCM):
    """Rate and State Model (RCM) class documentation.
    RSM estimates the time dependent seimicity response for a given stress loading scenario.
    Theory is described in Dietrich (1994), JGR.
    """

    def _compute(self, config: Config) -> Result:
        cf_shad = np.zeros(self.nt)
        S0 = -config.Sshadow
        # self.chiz[0] = self.cf[0] - config.Sshadow
        cf_shad[0] = self.cf[0] - config.Sshadow
        ratez = np.zeros(self.nt)
        dS = np.ediff1d(self.cf, to_end=config.loading.strend)
        rinfty = config.chi0 * config.loading.strend
        Asig = -config.depthS
        #print("Asig", Asig)
        gamma = 1.0
        ratez[0] = 1.0
        for i in range(1, self.nt):
            dum = gamma / config.loading.strend
            # Cattania, PhD Eq.(6.2), Dieterich JGR 1994, Eq.(17):
            gamma = (dum - config.deltat / dS[i - 1]) * np.exp(
                -dS[i - 1] / Asig
            ) + config.deltat / dS[i - 1]
            gamma *= config.loading.strend
            ratez[i] = 1.0 / gamma
            cf_shad[i] = S0
            # self.chiz[i] = S0
        ratez = rinfty * ratez
        neqz = np.zeros(self.nt - 1)
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return self.t, cf_shad, self.cf, ratez, neqz


class RSD(LCM):
    """Rate and State Model a la Dietrich (RSrate_Dietrich) class documentation.
    RSrate_Dietrich estimates the time dependent seimicity response for a given stress loading scenario.
    Theory is described in Dietrich (1994), JGR.
    """

    def _compute(self, config: Config) -> Result:
        cf_shad = np.zeros(self.nt)
        ratez = np.zeros(self.nt)  # nt = len(S)
        gamma = 1.0
        ratez[0] = 1.0

        # dt = np.ediff1d(t, to_end=t[-1]-t[-2])
        # dS = np.ediff1d(S, to_end=S[-1]-S[-2])
        dt = np.ediff1d(self.t, to_end=config.loading.strend)
        dS = np.ediff1d(self.cf, to_end=config.loading.strend)

        rinfty = config.chi0 * config.loading.strend
        Asig = -config.depthS
        #print("Asig", Asig)
        gamma = 1.0
        ratez[0] = 1.0
        for i in range(1, self.nt):
            dum = gamma / config.loading.strend
            # Cattania, PhD Eq.(6.2), Dieterich JGR 1994, Eq.(17):
            gamma = (dum - config.deltat / dS[i - 1]) * np.exp(
                (-dS[i - 1] + Asig * np.log(config.loading.strend)) / Asig
            ) + config.loading.strend * config.deltat / dS[i - 1]
            ratez[i] = 1.0 / gamma
            # cf_shad[i] = S0
        ratez = rinfty * ratez
        neqz = np.zeros(self.nt - 1)
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return self.t, cf_shad, self.cf, ratez, neqz


class RSD1(LCM):
    """
    Seismicity response of RS for a continuous stress evolution S(t) for times t
    according to Heimisson & Segall, JGR 2018, Eq.(20) & Eq.(29) & Eq.(34)
    """

    def _compute(self, config: Config) -> Result:
        cf_shad = np.zeros(self.nt)
        S0 = +config.Sshadow
        cf_shad[0] = self.cf[0] - config.Sshadow
        t1 = np.min(self.t[(self.cf > config.Sshadow)])
        if t1 > self.t[0]:
            i0 = np.argmax(self.t[(self.t < t1)])
            tb = np.interp(config.Sshadow, self.cf[i0 : i0 + 2], self.t[i0 : i0 + 2])
        else:
            tb = self.t[0]
        ti = self.t[(self.t >= tb)] - tb
        Si = self.cf[(self.t >= tb)] - config.Sshadow
        dt = np.ediff1d(ti, to_end=ti[-1] - ti[-2])
        Asig = -config.depthS
        ta = Asig / config.loading.strend
        r0 = config.chi0 * config.loading.strend
        K = np.exp(Si / Asig)
        integK = np.cumsum(K * dt)
        ratez = np.zeros(len(self.t))
        ratez[(self.t >= tb)] = r0 * K / (1.0 + integK / ta)

        neqz = np.zeros(self.nt - 1)
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return self.t, cf_shad, self.cf, ratez, neqz
