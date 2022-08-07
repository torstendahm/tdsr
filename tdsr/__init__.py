# -*- coding: utf-8 -*-

"""Top-level package for tdsr."""

__author__ = """torstendahm"""
__email__ = "torsten.dahm@gfz-potsdam.de"

from tdsr.version import version as __version__
from tdsr.config import Config
from tdsr.tdsr import CFM, LCM, RSD, RSD1, RSM, TDSR, TDSR1, Result, Traditional
from tdsr.utils import load, save

__all__ = [
    "Result",
    "Config",
    "save",
    "load",
    "TDSR",
    "TDSR1",
    "LCM",
    "Traditional",
    "CFM",
    "RSD",
    "RSD1",
    "RSM",
]
