import toml
from pprint import pprint
from typing import Optional, Dict
from tdsm.constants import HOURS
from tdsm.utils import Number, PathLike
from tdsm.loading import Loading, StepLoading, LOADING


class Config(object):
    def __init__(
        self,
        chi0=10.0,
        depthS=-0.2,
        Sshadow=0.0,
        deltat=0.2 * HOURS,
        tstart=0 * HOURS,
        tend=30 * HOURS,
        deltaS=0.01,
        sigma_max=10.0,
        precision=12,
        loading: Optional[Loading] = None,
    ) -> None:
        self.chi0 = chi0
        self.depthS = depthS
        self.Sshadow = Sshadow
        self.deltat = deltat
        self.tstart = tstart
        self.tend = tend
        self.deltaS = deltaS
        self.sigma_max = sigma_max
        self.precision = precision
        if loading is None:
            self.loading = StepLoading(deltat=deltat)

    def merge(self, config: Dict[str, str]) -> None:
        config = {k: v for k, v in config.items() if v is not None}
        if "chi0" in config:
            self.chi0 = config["chi0"]
        if "depthS" in config:
            self.depthS = config["depthS"]
        if "Sshadow" in config:
            self.depthS = config["Sshadow"]
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
    def open(cls, config_file: PathLike):
        config = cls()
        parsed_config = toml.load(config_file)

        loading: Optional[Loading] = None
        try:
            loading_typ = parsed_config["use_loading"]
            loading_params = parsed_config["loading"][loading_typ]
            loading_cls = LOADING[loading_typ]
            config.loading = loading_cls(**loading_params)
        except KeyError:
            # todo: warn about missing loading
            pass

        config.merge(parsed_config)
        return config
