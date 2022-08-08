from tdsr.loading.loading import Loading

from pathlib import Path
from typing import TYPE_CHECKING, Optional
import numpy as np
import numpy.typing as npt
from tdsr.utils import gridrange
from tdsr.types import Number


if TYPE_CHECKING:
    from tdsr.config import Config


class ExternalFileLoading(Loading):
    __name__: str = "InFile"

    def __init__(
        self,
        _config: "Config",
        strend: Number = 7.0e-5,
        tstart: Number = 0.0,
        tend: Number = 86400.0,
        taxis_log: int = 0,
        deltat: Number = 720.0,
        scal_t: Number = 3600 * 24,
        scal_cf: Number = 1.0e-6,
        c_tstart: Number = 0.0,
        iascii: Optional[bool] = None,
        infile: Optional = None,
    ):
        self.config = _config
        self.strend = strend
        self.iascii = iascii
        self.infile = infile
        self.tstart = tstart
        self.tend = tend
        self.taxis_log = taxis_log
        self.deltat = deltat
        self.c_tstart = c_tstart
        self.scal_cf = scal_cf
        self.scal_t = scal_t

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:

        if self.taxis_log != 0:
            print(
                "logarithmic time samples not possible when "
                "reading in stress loading function. Set taxis_log = 0"
            )
            exit()
        if self.iascii:
            if not Path(self.infile).is_file():
                raise ValueError("input file does not exist")
            cf_obs, t_obs = np.loadtxt(
                self.infile, skiprows=2, usecols=(0, 1), unpack=True
            )
            t_obs = t_obs * self.scal_t  # scale time
            cf_obs = cf_obs * self.scal_cf  # scale stress
            tmin_obs = t_obs[0]
            # tmax_obs = t_obs[-1]
            # time series stops with end time of series read in
            # self.tend = tmax_obs
            # nt = len(t_obs)
            # dt = (tmax-tmin)/nt
        else:
            raise ValueError("so far only ascii format supported for input file")

        # insert one value before start sample of stress change read in
        if self.tstart > tmin_obs:
            print("tmin_obs=", tmin_obs, " tstart=", self.tstart)
            raise ValueError(
                "tstart must be smaller than the begin time of series read in"
            )
        tmin, tmax, nt, t, dt = gridrange(self.tstart, self.tend, self.deltat)

        if self.tstart < tmin_obs:
            t_temp = np.insert(t_obs, 0, self.tstart, axis=None)
            c_temp = np.insert(cf_obs, 0, self.c_tstart, axis=None)
            # cf = np.interp(t[0:-1], t_temp, c_temp)
            cf = np.interp(t, t_temp, c_temp)
            # print('t_temp =',len(t_temp),' len(cf_temp)=',len(c_temp))
        else:
            cf = np.interp(t, t_obs, cf_obs)
        # print('tstart, tend, deltat=',self.tstart,self.tend,self.deltat)
        # print('t =',len(t),' len(cf)=',len(cf))
        # print('t_ob   =',len(t_obs),' len(cf_ob  )=',len(cf_obs))

        return np.hstack(
            [
                cf,
            ]
        )
