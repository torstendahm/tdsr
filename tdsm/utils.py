import os
from typing import Union
import numpy as np

Number = Union[int, float]
PathLike = Union[str, os.PathLike]


def shifted(x: np.array, nshift: int) -> np.array:
    # other option: use padding and cut
    x = np.roll(x, nshift)
    if nshift > 0:
        x[0:nshift] = 1.0
    elif nshift < 0:
        x[nshift:] = 0.0
    return x


def gridrange(amin: Number, amax: Number, astep: Number, rounding=np.ceil):
    amin = float(amin)
    amax = float(amax)
    astep = float(astep)
    na = int(rounding((amax - amin) / astep))
    amax = amin + (na - 1) * astep
    a = np.linspace(amin, amax, na)
    return amin, amax, na, a