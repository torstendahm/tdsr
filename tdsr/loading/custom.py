from tdsr.loading.loading import Loading

from pathlib import Path
from typing import TYPE_CHECKING, Optional
import numpy as np
import numpy.typing as npt
from tdsr.utils import gridrange, DEBUG
from tdsr.types import Number, PathLike
from tdsr.exceptions import InvalidParameter, MissingParameter


if TYPE_CHECKING:
    from tdsr.config import Config


class CustomLoading(Loading):
    """
    CustomLoading (short  name "custom") is used to read an arbitrary stress loading file from disk or via an argument. Both ascii and binary files can be loaded if the formatting is correct. See examples for further explanations
    """

    __name__: str = "Custom"

    def __init__(
        self,
        file: Optional[PathLike] = None,
        data: Optional[PathLike] = None,
        strend: Number = 7.0e-5,
        tstart: Number = 0.0,
        tend: Number = 86400.0,
        taxis_log: bool = False,
        deltat: Number = 720.0,
        scal_t: Number = 3600 * 24,
        scal_cf: Number = 1.0e-6,
        c_tstart: Number = 0.0,
        config: Optional["Config"] = None,
    ):
        self.config = config
        self.strend = strend
        self.tstart = tstart
        self.tend = tend
        self.taxis_log = taxis_log
        self.deltat = deltat
        self.c_tstart = c_tstart
        self.scal_cf = scal_cf
        self.scal_t = scal_t

        if file is not None:
            if not Path(file).is_file():
                raise FileNotFoundError  # ("input file does not exist")
            with open(file, "rb") as f:
                self.data = np.loadtxt(
                    f,
                    skiprows=2,
                    usecols=(0, 1),
                    unpack=False,
                )
        elif data is not None:
            if data.shape[1] != 2:
                raise InvalidParameter("data shape must be length x 2 (time, stress")
            self.data = data
        else:
            raise MissingParameter("custom loading missing file or data argument")

    @property
    def stress_rate(self) -> float:
        # return (self.sc1 - self.sc0) / (self.n1 * self.deltat)
        raise NotImplementedError

    def values(self, length: int) -> npt.NDArray[np.float64]:

        if self.taxis_log:
            raise InvalidParameter(
                "logarithmic time samples not possible when "
                "reading in stress loading function. Set taxis_log=False."
            )

        # scaled time
        t_obs = self.data[:, 0] * self.scal_t
        # scaled stress
        cf_obs = self.data[:, 1] * self.scal_cf

        tmin_obs = t_obs[0]
        tmax_obs = t_obs[-1]

        if self.tstart > tmin_obs:
            print("tmin_obs=", tmin_obs, " tstart=", self.tstart)
            raise InvalidParameter(
                "tstart must be smaller than the begin time of series read in"
            )
        tmin, tmax, nt, t, dt = gridrange(self.tstart, self.tend, self.deltat)

        if self.tstart < tmin_obs:
            t_temp = np.insert(t_obs, 0, self.tstart, axis=None)
            c_temp = np.insert(cf_obs, 0, self.c_tstart, axis=None)
            cf = np.interp(t, t_temp, c_temp)
            if DEBUG:
                print("t_temp =", len(t_temp), " len(cf_temp)=", len(c_temp))
        else:
            cf = np.interp(t, t_obs, cf_obs)

        if DEBUG:
            print("tstart, tend, deltat=", self.tstart, self.tend, self.deltat)
            print("t =", len(t), " len(cf)=", len(cf))
            print("t_ob   =", len(t_obs), " len(cf_ob  )=", len(cf_obs))

        return np.hstack(
            [
                cf,
            ]
        )
