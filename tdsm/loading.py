from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Dict, Optional, Type

import numpy as np
import numpy.typing as npt
#from typing import Type as npt

if TYPE_CHECKING:
    from tdsm.config import Config
    from tdsm.tdsm import Result

from tdsm.utils import DEBUG, Number


class Loading(ABC):
    __name__: str = ""

    def __init__(self) -> None:
        pass

    @property
    def name(self) -> str:
        return self.__name__

    @abstractmethod
    def values(self, length: int) -> npt.NDArray[np.float64]:
        pass

    @property
    @abstractmethod
    def stress_rate(self) -> float:
        pass


class FourPointLoading(Loading):
    __name__: str = "4points"

    def __init__(
        self,
        _config: "Config",
        n1: int = 50,
        n2: int = 51,
        deltat: Number = 0.0,
        sc0: Number = 0.0,
        sc1: Number = 0.5,
        sc2: Number = 1.0,
        sc3: Number = 2.0,
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


class StepLoading(Loading):
    __name__: str = "Step"

    def __init__(
        self,
        _config: "Config",
        strend: Number = 7.0E-5,
        tstep: Optional[Number] = None,
        sstep: Number = 1.0,
    ):
        self.config = _config
        self.strend = strend
        self.tstep = tstep or self.config.tend / 2
        self.sstep = sstep

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        sc0 = 0.0
        sc1 = (self.tstep - self.config.tstart) * self.strend
        sc2 = sc1 + self.sstep
<<<<<<< HEAD
        sc2plus = sc2 + self.config.deltat * self.strend
        sc3 = sc2plus + (self.config.tend - self.tstep - self.config.deltat) * self.strend
        #sc3 = sc2 + (self.config.tend - self.tstep - self.config.deltat) * self.strend
=======
        sc3 = sc2 + (self.config.tend - self.tstep - self.config.deltat) * self.strend
>>>>>>> refs/remotes/origin/master
        n1 = np.floor(self.tstep / self.config.deltat).astype(int)
        n2 = n1 + 1
        nt = length
        from pprint import pprint

        if DEBUG:
            pprint(
                dict(
                    sc0=sc0,
                    sc1=sc1,
                    sc2=sc2,
                    sc3=sc3,
                    n1=n1,
                    n2=n2,
                    nt=nt,
                )
            )
        if not (n1 >= 0 and n2 <= nt):
            raise ValueError("tstep must be greater than zero and smaller than tend")
        return np.hstack(
            [
                np.linspace(sc0, sc1, num=n1),
                np.linspace(sc1, sc2, num=n2 - n1 + 1)[1:],
                #np.linspace(sc2, sc3, num=nt - n2 + 1)[1:],
                np.linspace(sc2plus, sc3, num=nt - n2 + 1)[1:],
            ]
        )

class BackgroundLoading(Loading):
    __name__: str = "Background"

    def __init__(
        self,
        _config: "Config",
<<<<<<< HEAD
        strend: Number = 7.0E-5,
=======
        strend: Number = 0.00001,
>>>>>>> refs/remotes/origin/master
    ):
        self.config = _config
        self.strend = strend

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        sc0 = 0.0
        sc3 = (self.config.tend - self.config.tstart) * self.strend
        nt = length
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
                np.linspace(sc0, sc3, num=nt),
            ]
        )



LOADING: Dict[str, Type[Loading]] = {
    "step": StepLoading,
    "4points": FourPointLoading,
    "background": BackgroundLoading,
}
