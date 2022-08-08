from tdsr.loading.loading import Loading

from typing import TYPE_CHECKING
import numpy as np
import numpy.typing as npt
from tdsr.types import Number
from tdsr.utils import DEBUG


if TYPE_CHECKING:
    from tdsr.config import Config


class BackgroundLoading(Loading):
    __name__: str = "Background"

    def __init__(
        self,
        _config: "Config",
        strend: Number = 7.0e-5,
        tstart: Number = 0.0,
        tend: Number = 86400.0,
        deltat: Number = 720.0,
        taxis_log: Number = 0,
        ntlog: Number = 1000,
        sstep: Number = 0.0,
    ):
        self.config = _config
        self.strend = strend
        self.tstart = tstart
        self.tend = tend
        self.deltat = deltat
        self.taxis_log = taxis_log
        self.ntlog = ntlog
        self.sstep = sstep

    @property
    def stress_rate(self) -> float:
        # return (self.sc1 - self.sc0) / (self.n1 * self.deltat)
        return (self.sc3 - self.sc0) / (self.tend - self.tstart)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        # sc0 = self.sstep
        # spaeter evtl wieder auskommentieren:
        # so ist es gleich wie TDSR  fuer step
        sc0 = 0.0
        if self.taxis_log == 1:
            nt = self.ntlog
            tvalues = np.logspace(np.log10(self.tstart), np.log10(self.tend), nt)
            cf = sc0 + tvalues * self.strend
        else:
            nt = length
            # sc3 = self.tend*self.strend
            sc3 = sc0 + (self.tend - self.tstart) * self.strend
            cf = np.linspace(sc0, sc3, num=nt)

        print(
            "in background: tstart, tend, dottau=", self.tstart, self.tend, self.strend
        )
        print("taxis_log=", self.taxis_log)
        print("ntlog=", self.ntlog)
        print("nt=", length, nt, " len(cf)=", len(cf))
        print("min and max cf=", np.amin(cf), np.amax(cf))
        from pprint import pprint

        if DEBUG:
            pprint(
                dict(
                    sc0=sc0,
                    sc3=sc3,
                    nt=nt,
                )
            )
        return np.hstack(
            [
                cf,
            ]
        )
