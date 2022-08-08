from tdsr.loading.loading import Loading

from typing import TYPE_CHECKING, Optional
import numpy as np
import numpy.typing as npt
from tdsr.types import Number
from tdsr.utils import DEBUG


if TYPE_CHECKING:
    from tdsr.config import Config


class CyclicLoading(Loading):
    """
    Class CyclicLoading (short  name "Cycle") defines a sinusoidal stress loading function with period Tsin and amplitude ampsin. A linear trend is superposed (strend). An equal distance sampling of the stress function with sampling interval of "deltat" is generated.
    """

    __name__: str = "Cycle"

    def __init__(
        self,
        strend: Number = 7.0e-5,
        ampsin: Number = 0.2e0,
        Tsin: Number = 43200,
        tstart: Number = 0.0,
        # tend: Number = 86400.0,
        tend: Optional[Number] = None,
        deltat: Number = 720.0,
        config: Optional["Config"] = None,
    ):
        self.config = config
        self.strend = strend
        self.ampsin = ampsin
        self.Tsin = Tsin
        self.tstart = tstart
        # todo
        self.tend = tend or self.config.tend
        self.deltat = deltat

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        sc0 = 0.0
        sc3 = (self.tend - self.tstart) * self.strend
        nt = length
        from pprint import pprint

        seg1 = np.linspace(sc0, sc3, num=nt)
        temp = np.linspace(self.tstart, self.tend, num=nt)
        seg2 = self.ampsin * (1.0 - np.cos(2.0 * np.pi * temp / self.Tsin))
        seg = seg1 + seg2

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
                seg,
            ]
        )
