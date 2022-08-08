from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt


class Loading(ABC):
    """
    Loading abstract base class to be implemented by concrete
    loadings.
    """

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
