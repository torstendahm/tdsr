from tdsr.loading.loading import Loading

from typing import TYPE_CHECKING, Optional
import numpy as np
import numpy.typing as npt
from tdsr.types import Number
from tdsr.utils import DEBUG


if TYPE_CHECKING:
    from tdsr.config import Config


class StepLoading(Loading):
    __name__: str = "Step"

    def __init__(
        self,
        _config: "Config",
        strend: Number = 7.0e-5,
        tstep: Optional[Number] = None,
        sstep: Number = 1.0,
        tstart: Number = 0.0,
        tend: Number = 86400.0,
        taxis_log: Number = 0,
        ntlog: Number = 1000,
        deltat: Number = 720.0,
    ):
        self.config = _config
        self.strend = strend
        self.tstep = tstep or self.config.tend / 2
        self.sstep = sstep
        self.tstart = tstart
        self.tend = tend
        self.taxis_log = taxis_log
        self.ntlog = ntlog
        self.deltat = deltat

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        if self.taxis_log == 1:
            print("taxis_log=1 (logarithmic scale), step is at tstep=0 and start>0")
            print(
                "Steploading not configured for logarithmic t-axis. "
                "Use BackgroundLoading with sstep=sstep and tdsr() with "
                "Sshadow=sstep"
            )
            exit()
            nt = self.ntlog
            tvalues = np.logspace(np.log10(self.tstart), np.log10(self.tend), nt)
            seg2 = sc0 + tvalues * self.strend
        else:
            n1 = np.floor((self.tstep - self.tstart) / self.deltat).astype(int) + 1
            nt = length
            sinterval = self.deltat * self.strend
            ninter1 = n1
            ninter2 = nt - n1 - 1
            sc0 = 0.0
            sc1 = float(ninter1) * sinterval
            sc2 = sc1 + self.sstep
            sc3 = sc2 + float(ninter2) * sinterval
            seg1 = np.linspace(sc0, sc1, num=n1 + 1)[0:-1]
            seg3 = np.linspace(sc2, sc3, num=nt - n1)
            # print('n1=',n1,' nt=',nt,' ninter1=',ninter1,' ninter2=',ninter2)
            # print(' test segment 1')
            # print(seg1[-3:])
            # print('slope=',seg1[-1]-seg1[-2])
            # print(' ')
            # print(' test segment 3')
            # print(seg3[0:3])
            # print('slope=',seg3[1]-seg3[0])
            # print(' ')
            # print('len1=',len(seg1),' len3=',len(seg3))
            # exit()

        from pprint import pprint

        if DEBUG:
            pprint(
                dict(
                    sc0=sc0,
                    sc1=sc1,
                    sc2=sc2,
                    sc3=sc3,
                    n1=n1,
                    nt=nt,
                )
            )
        if not (n1 >= 0 and n1 + 1 <= nt):
            raise ValueError("tstep must be greater than zero and smaller than tend")
        return np.hstack([seg1, seg3])