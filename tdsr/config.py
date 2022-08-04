################################
# Time Dependent Seismicity Model - Config
# T. Dahm, R. Dahm 26.12.2021
################################

from typing import Any, Dict, Optional

import toml

from tdsr.constants import HOURS
from tdsr.loading import LOADING, Loading, StepLoading
from tdsr.utils import Number, PathLike


class Config(object):
    def __init__(
        self,
        chi0: Number = 10000.0,
        depthS: Number = -0.5,
        Sshadow: Number = 0.0,
        t0: Number = 0.2 * HOURS,
        deltat: Number = 0.2 * HOURS,
        tstart: Number = 0 * HOURS,
        tend: Number = 30 * HOURS,
        taxis_log: int = 0,
        ntlog: int = 1000,
        deltaS: Number = 0.0125,
        iX0switch: int = 0,
        Zmean: float = 0.0,
        Zstd: float = 1.0,
        equilibrium: bool = False,
        sigma_max: int = 25,
        precision: int = 18,
        loading: Optional[Loading] = None,
    ) -> None:
        self.chi0 = float(chi0)
        self.depthS = float(depthS)
        self.Sshadow = float(Sshadow)
        self.Zmean = float(Zmean)
        self.Zstd = float(Zstd)
        self.equilibrium = bool(equilibrium)
        self.t0 = float(t0)
        self.deltat = float(deltat)
        self.tstart = float(tstart)
        self.tend = float(tend)
        self.taxis_log = int(taxis_log)
        self.ntlog = int(ntlog)
        self.deltaS = float(deltaS)
        self.sigma_max = float(sigma_max)
        self.precision = int(precision)
        self.iX0switch = int(iX0switch)
        self.loading = loading
        if loading is None:
            self.loading = StepLoading(_config=self)

    def merge(self, config: Dict[str, Any]) -> None:
        config = {k: v for k, v in config.items() if v is not None}
        if "chi0" in config:
            self.chi0 = config["chi0"]
        if "t0" in config:
            self.t0 = config["t0"]
        if "depthS" in config:
            self.depthS = config["depthS"]
        if "Sshadow" in config:
            self.Sshadow = config["Sshadow"]
        if "equilibrium" in config:
            self.equilibrium = config["equilibrium"]
        if "deltat" in config:
            self.deltat = config["deltat"]
        if "tstart" in config:
            self.tstart = config["tstart"]
        if "tend" in config:
            self.tend = config["tend"]
        if "taxis_log" in config:
            self.taxis_log = config["taxis_log"]
        if "ntlog" in config:
            self.ntlog = config["ntlog"]
        if "deltaS" in config:
            self.deltaS = config["deltaS"]
        if "sigma_max" in config:
            self.sigma_max = config["sigma_max"]
        if "iX0switch" in config:
            self.iX0switch = config["iX0switch"]
        if "Zmean" in config:
            self.Zmean = config["Zmean"]
        if "Zstd" in config:
            self.Zstd = config["Zstd"]
        if "precision" in config:
            self.precision = config["precision"]

    @classmethod
    def open(cls, config_file: PathLike) -> "Config":
        config = cls()
        parsed_config = dict(toml.load(config_file))
        config.merge(parsed_config)
        try:
            loading_typ = parsed_config["use_loading"]
            loading_params = parsed_config["loading"][loading_typ]
            loading_cls = LOADING[loading_typ]
            config.loading = loading_cls(_config=config, **loading_params)
        except KeyError:
            # todo: warn about missing loading
            pass
        return config


DEFAULT_CONFIG = Config()
