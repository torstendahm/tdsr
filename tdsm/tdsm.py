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
from tdsm.constants import HOURS
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
        depthS: Optional[float] = None,
        Sshadow: Optional[float] = None,
        deltat: Optional[float] = None,
        tstart: Optional[float] = None,
        tend: Optional[float] = None,
        deltaS: Optional[float] = None,
        sigma_max: Optional[int] = None,
        precision: Optional[int] = None,
    ) -> Result:
        config = deepcopy(self.config)
        config.merge(
            dict(
                chi0=chi0,
                depthS=depthS,
                Sshadow=Sshadow,
                deltat=deltat,
                tstart=tstart,
                tend=tend,
                deltaS=deltaS,
                sigma_max=sigma_max,
                precision=precision,
            )
        )
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
        chiz = np.heaviside(-self.sigma, 0)
        nshift = np.around(config.Sshadow / config.deltaS, 0).astype(int)
        chiz = shifted(chiz, nshift)

        ratez = np.zeros(self.nt)
        ratez[0] = 0.0
        resid = 0.0
        for i in range(1, self.nt):
            deltacf = self.cf[i] - (self.cf[i - 1] - resid)
            nshift = np.around(deltacf / config.deltaS, 0).astype(int)
            resid = deltacf - nshift * config.deltaS
            # shift chiz (memory effect)
            chiz = shifted(chiz, nshift)
            ratez[i] = np.trapz(chiz * self.pz) * config.deltaS  # type: ignore
            # cut off chiz
            chiz = chiz * (1.0 - self.pz)

        ratez = ratez * config.chi0 / config.deltat
        neqz = np.zeros(self.nt - 1)
        # neqz[0] = 0.0
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return config, self.t, self.cf, ratez, neqz


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
                ratez[i] = config.chi0 * (self.cf[i] - self.cf[i - 1]) / config.deltat
            else:
                ratez[i] = 0.0 / config.deltat
            self.chiz[i] = S0
        ratez = ratez * config.chi0 / config.deltat
        neqz = np.zeros(self.nt - 1)
        # neqz[0] = 0.0
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])  # type: ignore
        return config, self.t, self.cf, ratez, neqz
