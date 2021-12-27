import numpy as np
from typing import Optional
from tdsm.utils import shifted, gridrange
from tdsm.config import Config
from tdsm.constants import HOURS
from copy import deepcopy


class LCM(object):
    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        assert self.config.loading is not None

    def __call__(
        self,
        chi0=None,
        depthS=None,
        Sshadow=None,
        deltat=None,
        tstart=None,
        tend=None,
        deltaS=None,
        sigma_max=None,
        precision=None,
    ):
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

    def _prepare(self, config: Config):
        (self.smin, self.smax, self.nsigma, self.sigma) = gridrange(
            -config.sigma_max, +config.sigma_max, config.deltaS
        )
        (self.tmin, self.tmax, self.nt, self.t) = gridrange(
            config.tstart, config.tend, config.deltat
        )
        self.chiz = np.zeros(self.nsigma)
        self.pz = np.heaviside(self.sigma, 1)
        self.cf = config.loading.values(length=self.nt)

    def _compute(self, config: Config):
        chiz = np.heaviside(-self.sigma, 0)
        nshift = np.around(config.Sshadow / config.deltaS, 0).astype(int)
        chiz = shifted(chiz, nshift)

        ratez = np.zeros(self.nt)
        ratez[0] = 0.0
        resid = 0.0
        for i in range(1, self.nt):
            deltacf = self.cf[i] - (self.cf[i - 1] - resid)
            nshift = np.around(deltacf / config.deltaS, 0).astype(np.int)
            resid = deltacf - nshift * config.deltaS
            # shift chiz (memory effect)
            chiz = shifted(chiz, nshift)
            ratez[i] = np.trapz(chiz * self.pz) * config.deltaS
            # cut off chiz
            chiz = chiz * (1.0 - self.pz)

        ratez = ratez * config.chi0 / config.deltat
        neqz = np.zeros(self.nt - 1)
        # neqz[0] = 0.0
        for i in range(1, self.nt - 2):
            neqz[i] = np.trapz(ratez[0 : i + 1])
        return config, self.t, self.cf, ratez, neqz


class TDSM(LCM):
    def _compute(self, config: Config):
        # use exponential decay for pz
        ndepth = np.around(config.depthS / config.deltaS, 0).astype(int)
        nzero = int(self.nsigma / 2)
        window = int(config.precision * ndepth)
        print(window)
        self.pz[nzero + window : nzero] = np.exp(-np.arange(window, 0) / ndepth)
        return super()._compute(config)


class Traditional(LCM):
    def _compute(self, config: Config):
        ratez = np.zeros(self.nt)
        S0 = -Sshadow
        self.chi1[0] = cf[0] - Sshadow
        for i in range(1, nt - 1):
            if cf[i] >= cf[i - 1] and cf[i] >= S0:
                S0 = cf[i]
                ratez[i] = chi0 * (cf[i] - cf[i - 1]) / deltat
            else:
                ratez[i] = 0.0 / deltat
            chiz[i] = S0
