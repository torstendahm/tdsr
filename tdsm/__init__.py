# -*- coding: utf-8 -*-

"""Top-level package for tdsm."""

__author__ = """torstendahm"""
__email__ = "torsten.dahm@gfz-potsdam.de"
__version__ = "0.0.1"

from tdsm.tdsm import LCM, TDSM, Result, Traditional
from tdsm.utils import load, save

__all__ = [
    "Result",
    "save",
    "load",
    "TDSM",
    "LCM",
    "Traditional",
]
