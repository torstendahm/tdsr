from abc import ABC, abstractmethod
from typing import Dict, Type

import numpy as np
import numpy.typing as npt

from tdsm.utils import Number


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


LOADING: Dict[str, Type[Loading]] = {
    "step": StepLoading,
    "4points": FourPointLoading,
}
