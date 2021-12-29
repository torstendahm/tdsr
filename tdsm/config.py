from typing import Any, Dict, Optional

import toml

from tdsm.constants import HOURS
from tdsm.loading import LOADING, Loading, StepLoading
from tdsm.utils import Number, PathLike


class Config(object):
    def __init__(
        self,
        chi0: Number = 10.0,
        depthS: Number = -0.2,
        Sshadow: Number = 0.0,
        deltat: Number = 0.2 * HOURS,
        tstart: Number = 0 * HOURS,
        tend: Number = 30 * HOURS,
        deltaS: Number = 0.01,
        equilibrium: bool = False,
        sigma_max: int = 10,
        precision: int = 12,
        loading: Optional[Loading] = None,
    ) -> None:
        self.chi0 = float(chi0)
        self.depthS = float(depthS)
        self.Sshadow = float(Sshadow)
        self.Equilibrium = bool(equilibrium)
        self.deltat = float(deltat)
        self.tstart = float(tstart)
        self.tend = float(tend)
        self.deltaS = float(deltaS)
        self.sigma_max = sigma_max
        self.precision = precision
        self.loading = loading
        if loading is None:
            self.loading = StepLoading(_config=self)

    def merge(self, config: Dict[str, Any]) -> None:
        config = {k: v for k, v in config.items() if v is not None}
        if "chi0" in config:
            self.chi0 = config["chi0"]
        if "depthS" in config:
            self.depthS = config["depthS"]
        if "Sshadow" in config:
            self.depthS = config["Sshadow"]
        if "equilibrium" in config:
            self.equilibrium = config["equilibrium"]
        if "deltat" in config:
            self.deltat = config["deltat"]
        if "tstart" in config:
            self.tstart = config["tstart"]
        if "tend" in config:
            self.tend = config["tend"]
        if "deltaS" in config:
            self.deltaS = config["deltaS"]
        if "sigma_max" in config:
            self.sigma_max = config["sigma_max"]
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
