from tdsr.loading.loading import Loading

from typing import TYPE_CHECKING, Optional
import numpy as np
import numpy.typing as npt
from tdsr.types import Number
from tdsr.utils import DEBUG
from tdsr.exceptions import InvalidParameter


if TYPE_CHECKING:
    from tdsr.config import Config


class StepLoading(Loading):
    """
    Class StepLoading (short name "Step") can define a stress step over a linear trend of Coulomb stress loading. An equal distance sampling of the stress function with sampling interval of ``deltat`` is generated if ``taxis_log==False``. Controlling parameter are the background stress trend (``strend``), the time when the stress step is applied (``tstep``), and the magnitude of the stress step (``sstep``). If ``taxis_log==True`` a logarithic time sampling is considered. Then, the time of the step is set to ``t=tstart``, and the  number of sampling points (``ntlog``) is specified instead of the sampling interval.
    """

    __name__: str = "Step"

    def __init__(
        self,
        strend: Number = 7.0e-5,
        tstep: Optional[Number] = None,
        sstep: Number = 1.0,
        tstart: Number = 0.0,
        # tstep: Optional[Number] = None,
        tend: Optional[Number] = None,
        # tend: Number = 86400.0,
        taxis_log: bool = False,
        ntlog: Number = 1000,
        deltat: Number = 720.0,
        config: Optional["Config"] = None,
    ):
        self.config = config
        self.strend = strend
        self.tend = tend or self.config.tend
        # todo
        assert self.tend is not None
        self.tstep = tstep or self.tend / 2
        self.sstep = sstep
        self.tstart = tstart
        self.taxis_log = taxis_log
        self.ntlog = ntlog
        self.deltat = deltat

    @property
    def stress_rate(self) -> float:
        """the loading stress rate"""
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        """the loading values"""
        if self.taxis_log:
            print("logarithmic scale used, step is at tstep=0 and start>0")
            raise InvalidParameter(
                "Steploading not configured for logarithmic t-axis. "
                "Use BackgroundLoading with sstep=sstep and tdsr() with "
                "Sshadow=sstep"
            )
            # exit()
            # nt = self.ntlog
            # tvalues = np.logspace(np.log10(self.tstart), np.log10(self.tend), nt)
            # seg2 = sc0 + tvalues * self.strend
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
            raise InvalidParameter(
                "tstep must be greater than zero and smaller than tend"
            )
        return np.hstack([seg1, seg3])
