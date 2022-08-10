from tdsr.loading.loading import Loading

from typing import TYPE_CHECKING, Optional
import numpy as np
import numpy.typing as npt
from tdsr.types import Number
from tdsr.utils import DEBUG
from tdsr.exceptions import InvalidParameter


if TYPE_CHECKING:
    from tdsr.config import Config


class TrendchangeLoading(Loading):
    """
    Class TrendchangeLoading (short  name "Trendchange") defines a loading function where stress rate changes abruptly from strend to strend2 at time tstep. A constant sampling interval is defined by deltat between tstart and tend.
    """

    __name__: str = "Trendchange"

    def __init__(
        self,
        strend: Number = 7.0e-5,
        strend2: Number = 7.0e-4,
        tstep: Optional[Number] = None,
        tstart: Number = 0.0,
        tend: Number = 86400.0,
        deltat: Number = 720.0,
        config: Optional["Config"] = None,
    ):
        self.config = config
        self.strend = strend
        self.strend2 = strend2
        self.tstep = tstep or self.config.tend / 2
        self.tstart = tstart
        self.tend = tend
        self.deltat = deltat

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        n1 = np.floor((self.tstep - self.tstart) / self.deltat).astype(int) + 1
        nt = length
        sinterval1 = self.deltat * self.strend
        sinterval2 = self.deltat * self.strend2
        ninter1 = n1
        ninter2 = nt - n1 - 1
        sc0 = 0.0
        sc1 = float(ninter1) * sinterval1
        sc2 = sc1 + float(ninter2) * sinterval2
        seg1 = np.linspace(sc0, sc1, num=n1 + 1)[0:-1]
        seg2 = np.linspace(sc1, sc2, num=nt - n1)

        from pprint import pprint

        if DEBUG:
            pprint(
                dict(
                    sc0=sc0,
                    sc1=sc1,
                    sc2=sc2,
                    n1=n1,
                    nt=nt,
                )
            )
        if not (n1 >= 0 and n1 + 1 <= nt):
            raise InvalidParameter(
                "tstep must be greater than zero and smaller than tend"
            )
        return np.hstack([seg1, seg2])
