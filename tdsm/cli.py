# -*- coding: utf-8 -*-

"""Console script for tdsm."""
import os
import sys
from pathlib import Path
from typing import Optional
from tdsm.utils import PathLike
from tdsm.config import Config
from tdsm import TDSM, LCM, Traditional
import click


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

DEFAULT_CONFIG = Config()


@click.group(invoke_without_command=True)
@click.argument("input", callback=assert_valid_dir, type=click.Path())
@click.argument("output", default=None, type=click.Path())
@click.option(
    "-c",
    "--config",
    default=None,
    type=click.Path(),
    help=str("configuration file path"),
)
@click.option(
    "--chi0",
    default=None,
    type=float,
    help="chi0 (default %s)" % (DEFAULT_CONFIG.chi0),
)
@click.option(
    "--depthS",
    default=None,
    type=float,
    help="depthS",
)
@click.option(
    "--Sshadow",
    default=None,
    type=float,
    help="Sshadow",
)
@click.option(
    "--deltat",
    default=None,
    type=float,
    help="deltat",
)
@click.option(
    "--tstart",
    default=None,
    type=float,
    help="tstart",
)
@click.option(
    "--tend",
    default=None,
    type=float,
    help="tend",
)
@click.option(
    "--deltaS",
    default=None,
    type=float,
    help="deltaS",
)
@click.option(
    "--sigma_max",
    default=None,
    type=int,
    help="sigma_max",
)
@click.option(
    "--precision",
    default=None,
    type=int,
    help="precision",
)
@click.option(
    "--precision",
    default=None,
    type=float,
    help="precision",
)
@click.option(
    "--loading",
    default=None,
    type=str,
    help="loading",
)
@click.pass_context
def tdsm(
    ctx: click.Context,
    input: PathLike,
    output: Optional[PathLike],
    config: Optional[PathLike],
    chi0: Optional[float],
    depths: Optional[float],
    sshadow: Optional[float],
    deltat: Optional[float],
    tstart: Optional[float],
    tend: Optional[float],
    deltas: Optional[float],
    sigma_max: Optional[int],
    precision: Optional[int],
    loading: Optional[str],
) -> None:
    """ TDSM method """
    # check if output dir is available
    print(output)

    if ctx.invoked_subcommand is None:
        if config is not None:
            conf = Config.open(config)
        else:
            conf = Config()
        conf.merge(
            dict(
                hi0=chi0,
                depthS=depths,
                Sshadow=sshadow,
                deltat=deltat,
                tstart=tstart,
                tend=tend,
                deltaS=deltas,
                sigma_max=sigma_max,
                precision=precision,
            )
        )
        print(conf)
        tdsm = TDSM(config=conf)
        result = tdsm()
        print(result)
    else:
        ctx.ensure_object(dict)
        ctx.obj["CONFIG"] = BaseParams(
            input_dir=input_dir,
            output_dir=output_dir,
        )


@tdsm.command()
@click.pass_context
def lcm(
    ctx: click.Context,
) -> int:
    """ LCM method """
    return 0


@tdsm.command()
@click.pass_context
def traditional(
    ctx: click.Context,
) -> int:
    """ traditional method """
    return 0


if __name__ == "__main__":
    sys.exit(tdsm(obj=dict()))  # pragma: no cover
