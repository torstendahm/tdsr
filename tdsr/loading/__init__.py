################################
# Time Dependent Seismicity Model - Loading implementations
# T. Dahm, R. Dahm 26.12.2021
################################

from typing import Dict, Type
from tdsr.loading.loading import Loading
from tdsr.loading.four_point import FourPointLoading
from tdsr.loading.step import StepLoading
from tdsr.loading.cyclic import CyclicLoading
from tdsr.loading.background import BackgroundLoading
from tdsr.loading.trend_change import TrendchangeLoading
from tdsr.loading.ramp import RampLoading
from tdsr.loading.custom import CustomLoading

LOADING: Dict[str, Type[Loading]] = {
    "step": StepLoading,
    "4points": FourPointLoading,
    "background": BackgroundLoading,
    "cycle": CyclicLoading,
    "trendchange": TrendchangeLoading,
    "ramp": RampLoading,
    "custom": CustomLoading,
}
