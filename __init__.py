# -*- coding: utf-8 -*-

"""Top-level package for tdsm."""

__author__ = """torstendahm"""
__email__ = "torsten.dahm@gfz-potsdam.de"
__version__ = "0.0.1"

from proto_compile.options import BaseCompilerOptions, CompilerOptions, CompileTarget
from proto_compile.proto_compile import compile, compile_grpc_web, compile_python_grpc
from proto_compile.versions import Target

__all__ = [
    "compile",
    "compile_grpc_web",
    "compile_python_grpc",
    "Target",
    "BaseCompilerOptions",
    "CompileTarget",
    "CompilerOptions",
]
