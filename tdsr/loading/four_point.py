from tdsr.loading.loading import Loading
from typing import TYPE_CHECKING, Optional
import numpy as np
import numpy.typing as npt
from tdsr.types import Number


if TYPE_CHECKING:
    from tdsr.config import Config


class FourPointLoading(Loading):
    """
    FourPointLoading (short name "4points") can define a Coulomb stress
    loading by 4 stress values at subsequent times.
    Equally-spaced stress function with a sampling interval of "deltat" is
    generated. The first Coulomb stress value (sc0) is assumed at time tstart.
    The last one (sc3) at time tend.
    The sample  times of the two stress values in between (sc1 and sc2)
    are defined by integer number of the sampling interval,
    n1 and n2, respectively.
    """

    __name__: str = "4points"

    def __init__(
        self,
        n1: int = 50,
        n2: int = 51,
        deltat: Number = 0.0,
        sc0: Number = 0.0,
        sc1: Number = 0.5,
        sc2: Number = 1.0,
        sc3: Number = 2.0,
        config: Optional["Config"] = None,
    ):
        self.n1 = n1
        self.n2 = n2
        self.deltat = float(deltat)
        self.sc0 = float(sc0)
        self.sc1 = float(sc1)
        self.sc2 = float(sc2)
        self.sc3 = float(sc3)

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        return np.hstack(
            [
                np.linspace(self.sc0, self.sc1, num=self.n1),
                np.linspace(self.sc1, self.sc2, num=self.n2 - self.n1 + 1)[1:],
                np.linspace(self.sc2, self.sc3, num=length - self.n2 + 1)[1:],
            ]
        )
