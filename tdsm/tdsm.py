"""
TDSM
====================================
Add some documentation about TDSM here
"""

from copy import deepcopy
from typing import Optional, Tuple

import numpy as np
import numpy.typing as npt

from tdsm.config import Config
from tdsm.loading import Loading
from tdsm.utils import gridrange, shifted

Result = Tuple[
    Config,
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
]


class LCM(object):
    """LCM class documentation."""

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
        chiz: Optional[float] = None,
        depthS: Optional[float] = None,
        Sshadow: Optional[float] = None,
        equilibrium: Optional[bool] = None,
        deltat: Optional[float] = None,
        tstart: Optional[float] = None,
        tend: Optional[float] = None,
        deltaS: Optional[float] = None,
        sigma_max: Optional[int] = None,
        precision: Optional[int] = None,
        loading: Optional[Loading] = None,
    ) -> Result:
        config = deepcopy(self.config)
        config.merge(
            dict(
                chi0=chi0,
                depthS=depthS,
                Sshadow=Sshadow,
                equilibrium=equilibrium,
                deltat=deltat,
                tstart=tstart,
                tend=tend,
                deltaS=deltaS,
                sigma_max=sigma_max,
                precision=precision,
            )
        )
        if loading is not None:
            config.loading = loading
        #if config.equilibrium:
        if chiz is not None:
            self._prepare(config)
            chiz_background = chiz
            #print('chiz_background=',chiz_background)
            self.chiz = chiz_background
        else:
            #print('chiz is none')
            self._prepare(config)
        return self._compute(config)

    def _prepare(self, config: Config) -> None:
        (self.smin, self.smax, self.nsigma, self.sigma) = gridrange(
            -config.sigma_max, +config.sigma_max, config.deltaS
        )
        (self.tmin, self.tmax, self.nt, self.t) = gridrange(
            config.tstart, config.tend, config.deltat
        )
        self.chiz = np.zeros(self.nsigma)
        self.pz = np.heaviside(self.sigma, 1)
        loading = config.loading
        if loading is None:
            raise ValueError("missing loading function")
        self.cf = loading.values(length=self.nt)

    def _compute(self, config: Config) -> Result:
        if config.equilibrium:
            #chiz = np.heaviside(-self.sigma, 0)
            #raise ValueError("the option to calculate or read the equilibrium function chiz is not yet implemented")
            print("chiz bereits uebergeben, daher nicht ueberschrieben mit Heaviside")
        else: 
            self.chiz = np.heaviside(-self.sigma, 0)

        nshift = np.around(config.Sshadow / config.deltaS, 0).astype(int)
        self.chiz = shifted(self.chiz, nshift)

        ratez = np.zeros(self.nt)
        ratez[0] = 0.0
        resid = 0.0
        for i in range(1, self.nt):
            deltacf = self.cf[i] - (self.cf[i - 1] - resid)
            nshift = np.around(deltacf / config.deltaS, 0).astype(int)
            resid = deltacf - nshift * config.deltaS
            # shift chiz (memory effect)
            self.chiz = shifted(self.chiz, nshift)
            ratez[i] = np.trapz(self.chiz * self.pz) * config.deltaS  # type: ignore
            # cut off chiz
            self.chiz = self.chiz * (1.0 - self.pz)

        ratez = ratez * config.chi0 / config.deltat
        neqz = np.zeros(self.nt - 1)
        # neqz[0] = 0.0
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return config, self.t, self.chiz, self.cf, ratez, neqz


class TDSM(LCM):
    def _compute(self, config: Config) -> Result:
        # use exponential decay for pz
        ndepth = np.around(config.depthS / config.deltaS, 0).astype(int)
        nzero = int(self.nsigma / 2)
        window = int(config.precision * ndepth)
        self.pz[nzero + window : nzero] = np.exp(-np.arange(window, 0) / ndepth)
        return super()._compute(config)


class Traditional(LCM):
    def _compute(self, config: Config) -> Result:
        ratez = np.zeros(self.nt)
        S0 = -config.Sshadow
        self.chiz[0] = self.cf[0] - config.Sshadow
        for i in range(1, self.nt - 1):
            if self.cf[i] >= self.cf[i - 1] and self.cf[i] >= S0:
                S0 = self.cf[i]
                #ratez[i] = config.chi0 * (self.cf[i] - self.cf[i - 1]) / config.deltat
                ratez[i] = (self.cf[i] - self.cf[i - 1])
            else:
                #ratez[i] = 0.0 / config.deltat
                ratez[i] = 0.0
            self.chiz[i] = S0
        ratez = ratez * config.chi0 / config.deltat
        neqz = np.zeros(self.nt - 1)
        # neqz[0] = 0.0
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return config, self.t, self.chiz, self.cf, ratez, neqz

class RSM(LCM):
    def _compute(self, config: Config) -> Result:
        S0 = -config.Sshadow
        self.chiz[0] = self.cf[0] - config.Sshadow
        ratez = np.zeros(self.nt)
        dS = np.ediff1d(self.cf, to_end=config.loading.strend)
        rinfty = config.chi0*config.loading.strend
        Asig = -config.depthS
        print('Asig',Asig)
        gamma = 1.0
        ratez[0] = 1.0
        for i in range(1, self.nt):
            dum = gamma / config.loading.strend
            # Cattania, PhD Eq.(6.2), Dieterich JGR 1994, Eq.(17):
            gamma = (dum - config.deltat/dS[i-1]) * np.exp(-dS[i-1]/Asig) + config.deltat/dS[i-1]
            gamma *= config.loading.strend
            ratez[i] = 1.0 / gamma
            self.chiz[i] = S0
        ratez = rinfty*ratez
        neqz = np.zeros(self.nt - 1)
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return config, self.t, self.chiz, self.cf, ratez, neqz

