from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Dict, Optional, Type

import numpy as np
import numpy.typing as npt
#from typing import Type as npt
#from scipy import interpolate
import os.path

if TYPE_CHECKING:
    from tdsm.config import Config
    from tdsm.tdsm import Result

from tdsm.utils import DEBUG, Number, gridrange


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
        tstart: Number = 0.,
        tend: Number = 86400.,
        deltat: Number = 720. 
    ):
        self.config = _config
        self.strend = strend
        self.tstep = tstep or self.config.tend / 2
        self.sstep = sstep
        self.tstart = tstart
        self.tend   = tend
        self.deltat = deltat

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        #n1 = np.floor( (self.tstep - self.config.tstart) / self.config.deltat).astype(int) +1
        n1 = np.floor( (self.tstep - self.tstart) / self.deltat).astype(int) +1
        nt = length
        sinterval = self.deltat * self.strend
        ninter1 = n1
        ninter2 = nt-n1-1
        sc0 = 0.0
        sc1 = float(ninter1) * sinterval
        sc2 = sc1 + self.sstep
        sc3 = sc2 + float(ninter2) * sinterval
        seg1 = np.linspace(sc0,sc1, num=n1+1)[0:-1]
        seg3 = np.linspace(sc2,sc3, num=nt-n1)
        #print('n1=',n1,' nt=',nt,' ninter1=',ninter1,' ninter2=',ninter2)
        #print(' test segment 1')
        #print(seg1[-3:])
        #print('slope=',seg1[-1]-seg1[-2])
        #print(' ')
        #print(' test segment 3')
        #print(seg3[0:3])
        #print('slope=',seg3[1]-seg3[0])
        #print(' ')
        #print('len1=',len(seg1),' len3=',len(seg3))
        #exit()

        from pprint import pprint

        if DEBUG:
            pprint(
                dict(
                    sc0=sc0,
                    sc1=sc1,
                    sc2=sc2,
                    sc3=sc3,
                    n1=n1,
                    nt=nt,
                )
            )
        if not (n1 >= 0 and n1+1 <= nt):
            raise ValueError("tstep must be greater than zero and smaller than tend")
        return np.hstack(
            [
                seg1 , seg3 
                # np.linspace(sc0, sc1, num=n1),
                #np.linspace(sc2, sc3, num=nt - n1 + 1)[1:],
            ]
        )

class CyclicLoading(Loading):
    __name__: str = "Cycle"

    def __init__(
        self,
        _config: "Config",
        strend: Number = 7.0E-5,
        ampsin: Number = 0.2E0,
        Tsin: Number = 43200,
        tstart: Number = 0.,
        tend: Number = 86400.,
        deltat: Number = 720. 
    ):
        self.config = _config
        self.strend = strend
        self.ampsin = ampsin
        self.Tsin = Tsin
        self.tstart = tstart
        self.tend   = tend
        self.deltat = deltat

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        sc0 = 0.0
        sc3 = (self.tend - self.tstart) * self.strend
        nt = length
        from pprint import pprint

        seg1 = np.linspace(sc0, sc3, num=nt)
        temp = np.linspace(self.tstart, self.tend, num=nt)
        seg2 = self.ampsin*(1.-np.cos(2.*np.pi*temp/self.Tsin))
        seg = seg1 + seg2

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
                seg, 
            ]
        )

class BackgroundLoading(Loading):
    __name__: str = "Background"

    def __init__(
        self,
        _config: "Config",
        strend: Number = 7.0E-5,
        tstart: Number = 0.,
        tend: Number = 86400.,
        deltat: Number = 720. 
    ):
        self.config = _config
        self.strend = strend
        self.tstart = tstart
        self.tend   = tend
        self.deltat = deltat

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        sc0 = 0.0
        sc3 = (self.tend - self.tstart) * self.strend
        print('in background: tstart, tend, dottau=', self.tstart, self.tend, self.strend)
        print('nt=', length)
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

class TrendchangeLoading(Loading):
    __name__: str = "Trendchange"

    def __init__(
        self,
        _config: "Config",
        strend: Number = 7.0E-5,
        strend2: Number = 7.0E-4,
        tstep: Optional[Number] = None,
        tstart: Number = 0.,
        tend: Number = 86400.,
        deltat: Number = 720. 
    ):
        self.config = _config
        self.strend = strend
        self.strend2 = strend2
        self.tstep = tstep or self.config.tend / 2
        self.tstart = tstart
        self.tend   = tend
        self.deltat = deltat

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        n1 = np.floor( (self.tstep - self.tstart) / self.deltat).astype(int) +1
        nt = length
        sinterval1 = self.deltat * self.strend
        sinterval2 = self.deltat * self.strend2
        ninter1 = n1
        ninter2 = nt-n1-1
        sc0 = 0.0
        sc1 = float(ninter1) * sinterval1
        sc2 = sc1 + float(ninter2) * sinterval2
        seg1 = np.linspace(sc0,sc1, num=n1+1)[0:-1]
        seg2 = np.linspace(sc1,sc2, num=nt-n1)

        from pprint import pprint

        if DEBUG:
            pprint(
                dict(
                    sc0=sc0,
                    sc1=sc1,
                    sc2=sc2,
                    n1=n1,
                    nt=nt,
                )
            )
        if not (n1 >= 0 and n1+1 <= nt):
            raise ValueError("tstep must be greater than zero and smaller than tend")
        return np.hstack(
            [
                seg1 , seg2 
            ]
        )

class RampLoading(Loading):
    __name__: str = "Ramp"

    def __init__(
        self,
        _config: "Config",
        strend: Number = 7.0E-5,
        strend2: Number = 7.0E-4,
        strend3: Number = 7.0E-5,
        nsample2: Number = 20,
        tstep: Optional[Number] = None,
        tstart: Number = 0.,
        tend: Number = 86400.,
        deltat: Number = 720. 
    ):
        self.config = _config
        self.strend = strend
        self.strend2 = strend2
        self.strend3 = strend3
        self.nsample2 = nsample2
        self.tstep = tstep or self.config.tend / 2
        self.tstart = tstart
        self.tend   = tend
        self.deltat = deltat

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:
        nt = length
        n1 = np.floor( (self.tstep - self.tstart) / self.deltat).astype(int) +1
        n2 = self.nsample2
        n3 = nt-(n1+n2)
        sinterval1 = self.deltat * self.strend
        sinterval2 = self.deltat * self.strend2
        sinterval3 = self.deltat * self.strend3
        ninter1 = n1
        ninter2 = n2
        ninter3 = n3+1
        sc0 = 0.0
        sc1 = float(ninter1) * sinterval1
        sc2 = sc1 + float(ninter2-1) * sinterval2
        sc3 = sc2 + float(ninter3-1) * sinterval3
        seg1 = np.linspace(sc0,sc1, num=n1+1)[0:-1]
        seg2 = np.linspace(sc1,sc2, num=ninter2)
        seg3 = np.linspace(sc2,sc3, num=ninter3)[1:]
        #print(' nt=',nt,' n1=',n1,' n2=',n2,' n3=',n3,' sum=',n1+n2+n3)
        #print(' test segment 1')
        #print(seg1[-3:])
        #print('slope=',seg1[-1]-seg1[-2])
        #print(' ')
        #print(' test segment 2')
        #print(seg2[0:3],seg2[-2:])
        #print('slope=',seg2[1]-seg2[0])
        #print(' ')
        #print(' test segment 3')
        #print(seg3[0:3])
        #print('slope=',seg3[1]-seg3[0])
        #print(' ')
        #print('len1=',len(seg1),' len2=',len(seg2),' len3=',len(seg3))

        from pprint import pprint

        if DEBUG:
            pprint(
                dict(
                    sc0=sc0,
                    sc1=sc1,
                    sc2=sc2,
                    sc3=sc3,
                    n1=n1,
                    nt=nt,
                )
            )
        if not (n1 >= 0 and n1+1 <= nt):
            raise ValueError("tstep must be greater than zero and smaller than tend")
        return np.hstack(
            [
                seg1 , seg2 , seg3 
            ]
        )

class ExternalFileLoading(Loading):
    __name__: str = "InFile"

    def __init__(
        self,
        _config: "Config",
        tstart: Number = 0.,
        tend: Number = 86400.,
        deltat: Number = 720., 
        scal_t: Number = 3600*24,
        scal_cf: Number = 1.E-6,
        c_tstart: Number = 0.0,
        iascii: Optional[bool] = None,
        infile: Optional = None
    ):
        self.config = _config
        self.iascii  = iascii
        self.infile = infile
        self.tstart = tstart
        self.tend   = tend
        self.deltat = deltat
        self.c_tstart = c_tstart
        self.scal_cf = scal_cf
        self.scal_t  = scal_t

    @property
    def stress_rate(self) -> float:
        return (self.sc1 - self.sc0) / (self.n1 * self.deltat)

    def values(self, length: int) -> npt.NDArray[np.float64]:

        if self.iascii:
            if not os.path.isfile(self.infile):
                raise ValueError("input file does not exist")
            cf_obs, t_obs = np.loadtxt(self.infile, skiprows=2, usecols=(0,1), unpack=True)
            t_obs = t_obs*self.scal_t     # scale time 
            cf_obs = cf_obs*self.scal_cf  # scale stress
            tmin_obs = t_obs[0]
            tmax_obs = t_obs[-1]
            self.tend     = tmax_obs # time series stops with end time of series read in
            #nt = len(t_obs)
            #dt = (tmax-tmin)/nt
        else:
            raise ValueError("so far only ascii format supported for input file")

        # insert one value before start sample of stress change read in
        if self.tstart >= tmin_obs:
            raise ValueError("tstart must be smaller than the begin time of series read in")
        tmin, tmax, nt, t = gridrange(self.tstart, self.tend, self.deltat)
        #cf    = np.zeros(nt)

        t_temp = np.insert(t_obs,0,self.tstart,axis=None)
        c_temp = np.insert(cf_obs,0,self.c_tstart,axis=None)
        #f = interpolate.interp1d(t_temp,c_temp, kind='linear')
        #cf = f(t)
        cf = np.interp(t[0:-1], t_temp, c_temp)
        print('tstart, tend, deltat=',self.tstart,self.tend,self.deltat)
        print('t =',len(t))
        print('cf=',len(cf))

        return np.hstack(
            [
                cf, 
            ]
        )

LOADING: Dict[str, Type[Loading]] = {
    "step": StepLoading,
    "4points": FourPointLoading,
    "background": BackgroundLoading,
    "cycle": CyclicLoading,
    "trendchange": TrendchangeLoading,
    "ramp": RampLoading,
    "infile": ExternalFileLoading,
}
