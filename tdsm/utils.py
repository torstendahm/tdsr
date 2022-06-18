import os
#import sys
import pickle as pkl
from typing import TYPE_CHECKING, Callable, Tuple, Union

import numpy as np
import numpy.typing as npt
#from typing import Type as npt

from scipy.stats import norm

if TYPE_CHECKING:
    from tdsm.tdsm import Result

DEBUG = os.environ.get("DEBUG") is not None
Number = Union[int, float]
PathLike = Union[str, os.PathLike]


def save(result: "Result", filename: PathLike) -> None:
    with open(filename, "wb") as f:
        pkl.dump(result, f)


def load(filename: PathLike) -> "Result":
    with open(filename, "rb") as f:
        result: Result = pkl.load(f)
        return result


##### Discretization ##############################

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

def Zvalues(S, Sstep, t0, dsig):
    NZ = 1000
    dS = np.ediff1d(S, to_end=S[-1]-S[-2])
    Smax = np.maximum(0, Sstep + np.max(np.cumsum(dS)))
    Z1 = -Smax - 10 * dsig
    Z2 = Smax + 10 * dsig
    Z = np.linspace(Z1, Z2, NZ)
    return Z

def shifted(x: npt.NDArray[np.float64], nshift: int) -> npt.NDArray[np.float64]:
    # other option: use padding and cut
    x = np.roll(x, nshift)
    if nshift > 0:
        x[0:nshift] = 1.0
    elif nshift < 0:
        x[nshift:] = 0.0
    return x

##### INITIAL CONDITIONS #####################

def X0steady(Z, r0, t0, dsig, dotsigc):
    '''
    Steady state initial stress distribution 
    --> R(t)=r0 for stressing rate dotS = dotsigc
    Eq.(6) of Dahm & Hainzl (2022)
    '''
    X = (r0 / dotsigc) * np.exp(-np.exp(-Z/dsig) * dsig /(t0 * dotsigc))   
    return X

def X0uniform(Z, Zmin, X0):
    '''
    Uniform initial stress disribution for Z>=Zmin
    '''
    X = X0 * np.heaviside(Z-Zmin, 1)
    return X

def X0gaussian(Z, Zmean, Zstd, X0):
    '''
    Normal initial stress disribution with mean Zmean and standard deviation Zstd
    '''
    X = X0 * norm.pdf(Z, loc=Z0mean, scale=Z0std)
    return X

##### THEORETICAL FUNCTIONS #######################

def Eq5(t, Zmin, X0, t0, dsig, dotsigc):
    '''
    Constant stressing (with rate dotsigc) given an initially uniform stress distribution
    Eq.(5) in Dahm & Hainzl (2022)
    '''
    ta = dsig / dotsigc 
    fac = X0 * dotsigc
    nominator = 1.0 - np.exp(-(ta/t0) * (1 - np.exp(-t/ta)) * np.exp(-(Zmin - dotsigc*t)/dsig))
    denominator = 1.0 - np.exp(-t/ta)
    R = fac * nominator / denominator
    return R
 
def Eq7(t, dS, r0, dsig, dotsigc, dotsigca):
    '''
    Stress step dS at time t=0 given an initial steady state state related constant rate r0 for dotsigc
    with constant stressing rate dotsigca for t>0
    Eq.(8) in Dahm & Hainzl (2022)
    '''
    f = dotsigca / dotsigc
    R = f * r0 / ((f * np.exp(-dS/dsig) - 1) * np.exp(-dotsigca*t/dsig) + 1)
    return R

def Eq8(t, dS, r0, dsig, dotsigc):
    '''
    Stress step dS at time t=0 given an initial steady state state related constant rate r0 for dotsigc
    without stress changes for t>0
    Eq.(7) in Dahm & Hainzl (2022)
    '''
    R = r0 / (np.exp(-dS/dsig) + (dotsigc/dsig) * t)
    return R
