import numpy as np
from tdsm.utils import Number, PathLike
from abc import ABC, abstractmethod


class Loading(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def values(self, length: int):
        pass

    @property
    @abstractmethod
    def stress_rate(self):
        pass


class StepLoading(Loading):
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
        self.deltat = deltat
        self.sc0 = sc0
        self.sc1 = sc1
        self.sc2 = sc2
        self.sc3 = sc3

    @property
    def stress_rate(self):
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int):
        return np.hstack(
            [
                np.linspace(self.sc0, self.sc1, num=self.n1),
                np.linspace(self.sc1, self.sc2, num=self.n2 - self.n1 + 1)[1:],
                np.linspace(self.sc2, self.sc3, num=length - self.n2 + 1)[1:],
            ]
        )

    # @classmethod
    # def default(cls):
    #     cls(
    #         n1=50,
    #         n2=51,
    #         sc0=0.0,
    #         sc1=0.5,
    #         sc2=1.0,
    #         sc3=2.0,
    #     )


LOADING = dict(
    step=StepLoading,
)
