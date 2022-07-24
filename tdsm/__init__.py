# -*- coding: utf-8 -*-

"""Top-level package for tdsm."""

__author__ = """torstendahm"""
__email__ = "torsten.dahm@gfz-potsdam.de"
__version__ = "0.0.1"

from tdsm.config import Config
from tdsm.tdsm import LCM, TDSM, TDSR1, Result, Traditional, CFM, RSM, RSD
from tdsm.utils import load, save

__all__ = [
    "Result",
    "Config",
    "save",
    "load",
    "TDSM",
    "TDSR1",
    "LCM",
    "Traditional",
    "CFM",
    "RSD",
    "RSM",
]
