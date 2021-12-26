# -*- coding: utf-8 -*-

"""Console script for tdsm."""
import os
import sys
import typing

import click

import proto_compile.proto_compile as compiler
import proto_compile.versions as versions
from proto_compile.options import BaseCompilerOptions
from proto_compile.utils import PathLike
from proto_compile.versions import Target


def assert_valid_dir(
    ctx: click.core.Context, param: click.core.Parameter, value: str
) -> str:
    if os.path.exists(value) and os.path.isdir(value):
        return value
    raise click.BadParameter(
        'Directory "' + value + '" does not exist.', ctx=ctx, param=param
    )


def optional_assert_valid_dir(
    ctx: click.core.Context, param: click.core.Parameter, value: str
) -> str:
    if value == "":
        return value
    return assert_valid_dir(ctx, param, value)


base_proto_parent_dir_help = (
    "base proto parent dir used for protoc -I=<base_proto_parent_dir>. ",
    "Must be a valid directory that contains the proto files in <proto_source_dir>",
)


@click.group()
@click.argument("proto-source-dir", callback=assert_valid_dir, type=click.Path())
@click.argument("output-dir", type=click.Path())
@click.option(
    "--minimal-include-dir",
    default=False,
    type=bool,
    help=str(
        "determine the shortest proto source include path that contains all protos in the given source dir"
    ),
)
@click.option(
    "--clear-output-dirs",
    "-clear",
    is_flag=True,
    default=False,
    help=str("whether to clear the output directories before compilation"),
)
@click.option(
    "--verbosity",
    default=0,
    help=str("level of verbosity when printing to stdout (the higher the more output)"),
)
@click.option(
    "--protoc-version",
    default=versions.DEFAULT_PROTOC_VERSION,
    help="protoc version to use (default is %s)" % versions.DEFAULT_PROTOC_VERSION,
)
@click.pass_context
def proto_compile(
    ctx: click.Context,
    proto_source_dir: str,
    output_dir: str,
    minimal_include_dir: bool,
    clear_output_dirs: bool,
    verbosity: int,
    protoc_version: str,
) -> None:
    ctx.ensure_object(dict)
    ctx.obj["COMPILER_OPTIONS"] = BaseCompilerOptions(
        proto_source_dir=proto_source_dir,
        output_dir=output_dir,
        minimal_include_dir=minimal_include_dir,
        clear_output_dirs=clear_output_dirs,
        verbosity=verbosity,
        protoc_version=protoc_version,
    )


@proto_compile.command()
@click.option(
    "--js_out_options",
    default="import_style=commonjs,binary",
    help=str("options for the javascript proto compiler"),
)
@click.option(
    "--grpc_web_out_options",
    default="import_style=typescript,mode=grpcwebtext",
    help=str("options for the grpc web proto compiler"),
)
@click.option(
    "--grpc_web_plugin_version",
    default=versions.DEFAULT_PLUGIN_VERSIONS[Target.GRPC_WEB],
    help="grpc web plugin version to use (default is %s)"
    % versions.DEFAULT_PLUGIN_VERSIONS[Target.GRPC_WEB],
)
@click.pass_context
def grpc_web(
    ctx: click.Context,
    js_out_options: str,
    grpc_web_out_options: str,
    grpc_web_plugin_version: str,
) -> int:
    """ compile using the grpc-web preset """
    try:
        compiler.compile_grpc_web(
            options=ctx.obj["COMPILER_OPTIONS"],
            js_out_options=js_out_options,
            grpc_web_out_options=grpc_web_out_options,
            grpc_web_plugin_version=grpc_web_plugin_version,
        )
    except Exception as e:  # pragma: no cover
        raise click.ClickException(str(e))
    return 0


@proto_compile.command()
@click.option(
    "--py_out_options", default=None, help=str("options for the python proto compiler"),
)
@click.option(
    "--py_output_dir",
    default=None,
    help=str("separate output dir for the python generated files"),
)
@click.option(
    "--py_grpc_out_options",
    default=None,
    help=str("options for the python grpc proto compiler"),
)
@click.option(
    "--py_grpc_output_dir",
    default=None,
    help=str("separate output dir for the python grpc generated files"),
)
@click.pass_context
def python_grpc(
    ctx: click.Context,
    py_out_options: typing.Optional[str],
    py_output_dir: typing.Optional[PathLike],
    py_grpc_out_options: typing.Optional[str],
    py_grpc_output_dir: typing.Optional[PathLike],
) -> int:
    """ compile using the python grpc preset """
    try:
        compiler.compile_python_grpc(
            options=ctx.obj["COMPILER_OPTIONS"],
            py_out_options=py_out_options,
            py_output_dir=py_output_dir,
            py_grpc_out_options=py_grpc_out_options,
            py_grpc_output_dir=py_grpc_output_dir,
        )
    except Exception as e:  # pragma: no cover
        raise click.ClickException(str(e))
    return 0


if __name__ == "__main__":
    sys.exit(proto_compile(obj=dict()))  # pragma: no cover
