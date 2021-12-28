import os
from typing import Callable, Tuple, Union

import numpy as np
import numpy.typing as npt

Number = Union[int, float]
PathLike = Union[str, os.PathLike]


def shifted(x: npt.NDArray[np.float64], nshift: int) -> npt.NDArray[np.float64]:
    # other option: use padding and cut
    x = np.roll(x, nshift)
    if nshift > 0:
        x[0:nshift] = 1.0
    elif nshift < 0:
        x[nshift:] = 0.0
    return x


def gridrange(
    amin: Number,
    amax: Number,
    astep: Number,
    rounding: Callable[[float], float] = np.ceil,
) -> Tuple[float, float, int, npt.NDArray[np.float64]]:
    amin = float(amin)
    amax = float(amax)
    astep = float(astep)
    na = int(rounding((amax - amin) / astep))
    amax = amin + (na - 1) * astep
    a = np.linspace(amin, amax, na)
    return amin, amax, na, a
