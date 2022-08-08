#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""TDSR exceptions"""


class TDSRException(Exception):
    pass


class InvalidParameter(TDSRException):
    def __init__(self, message: str) -> None:
        super().__init__(message)


class MissingParameter(TDSRException):
    def __init__(self, message: str) -> None:
        super().__init__(message)
