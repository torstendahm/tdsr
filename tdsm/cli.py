# -*- coding: utf-8 -*-

"""Console script for tdsm."""
import os
import sys
from pathlib import Path
from typing import Optional

import click

from tdsm import LCM, TDSM, Traditional
from tdsm.config import Config
from tdsm.utils import PathLike


def optional_valid_dir_or_file(
    ctx: click.core.Context, param: click.core.Parameter, value: str
) -> str:
    return value
    # if value == "":
    #     return value
    # if os.path.exists(value) and os.path.isdir(value) and not os.path.isfile(value):
    # return value
    # raise click.BadParameter(
    #     "%s is not a file or directory." % (value), ctx=ctx, param=param
    # )


DEFAULT_CONFIG = Config()


def get_output_file(out: PathLike, name: str) -> PathLike:
    _out = Path(out)
    is_dir = _out.suffix == ""
    if is_dir:
        return _out / f"{name}.pkl"
    return _out


def check_output_file(
    file: PathLike, force: bool = False, ask: bool = True
) -> PathLike:
    _file = Path(file)
    if _file.exists() and _file.is_file():
        if not force:
            if not ask or (
                ask
                and not click.confirm(
                    "%s already exists. Do you want to overwrite it?" % _file,
                    default=False,
                )
            ):
                raise FileExistsError("%s already exists" % _file)
    # create parent directories
    _file.parent.mkdir(parents=True, exist_ok=True)
    return _file


@click.group(invoke_without_command=True)
@click.option(
    "-i",
    "--input",
    callback=optional_valid_dir_or_file,
    default=None,
    help="input loading stress time series",
    type=click.Path(),
)
@click.option(
    "-o",
    "--output",
    callback=optional_valid_dir_or_file,
    default=None,
    help="output file",
    type=click.Path(),
)
@click.option(
    "-c",
    "--config",
    default=None,
    type=click.Path(),
    help="configuration file path",
)
@click.option(
    "-f",
    "--force",
    default=False,
    type=bool,
    help="force overwrite output results",
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
    help="depthS (default %s)" % (DEFAULT_CONFIG.depthS),
)
@click.option(
    "--Sshadow",
    default=None,
    type=float,
    help="Sshadow (default %s)" % (DEFAULT_CONFIG.Sshadow),
)
@click.option(
    "--deltat",
    default=None,
    type=float,
    help="deltat (default %s)" % (DEFAULT_CONFIG.deltat),
)
@click.option(
    "--tstart",
    default=None,
    type=float,
    help="tstart (default %s)" % (DEFAULT_CONFIG.tstart),
)
@click.option(
    "--tend",
    default=None,
    type=float,
    help="tend (default %s)" % (DEFAULT_CONFIG.tend),
)
@click.option(
    "--deltaS",
    default=None,
    type=float,
    help="deltaS (default %s)" % (DEFAULT_CONFIG.deltaS),
)
@click.option(
    "--sigma_max",
    default=None,
    type=int,
    help="sigma_max (default %s)" % (DEFAULT_CONFIG.sigma_max),
)
@click.option(
    "--precision",
    default=None,
    type=int,
    help="precision (default %s)" % (DEFAULT_CONFIG.precision),
)
@click.option(
    "--loading",
    default=None,
    type=str,
    help="loading (default %s)"
    % ("None" if DEFAULT_CONFIG.loading is None else DEFAULT_CONFIG.loading.name),
)
@click.pass_context
def tdsm(
    ctx: click.Context,
    input: PathLike,
    output: Optional[PathLike],
    config: Optional[PathLike],
    force: bool,
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
    """TDSM method"""
    output_path: Optional[PathLike] = None
    input_path: Optional[PathLike] = None
    if output:
        output_path = Path(output).resolve()
        output_file = get_output_file(output_path, name="tdsm")
        output_file = check_output_file(output_file, force=force)
        print(output_file)

    if config is not None:
        conf = Config.open(config)
    else:
        conf = Config()
    conf.merge(
        dict(
            chi0=chi0,
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

    if ctx.invoked_subcommand is None:
        tdsm = TDSM(config=conf)
        result = tdsm()
        print(result)
    else:
        ctx.ensure_object(dict)
        ctx.obj["PARAMS"] = dict(
            config=conf,
            input_path=input_path,
            output_path=output_path,
            force=force,
        )


@tdsm.command()
@click.pass_context
def lcm(
    ctx: click.Context,
) -> int:
    """LCM method"""
    params = ctx.obj["PARAMS"]
    config = params["config"]
    output_file = get_output_file(params["output_path"], name="lcm")
    output_file = check_output_file(output_file, force=params["force"])
    return 0


@tdsm.command()
@click.pass_context
def traditional(
    ctx: click.Context,
) -> int:
    """traditional method"""
    params = ctx.obj["PARAMS"]
    config = params["config"]
    output_file = get_output_file(params["output_path"], name="lcm")
    output_file = check_output_file(output_file, force=params["force"])
    return 0


if __name__ == "__main__":
    sys.exit(tdsm(obj=dict()))  # pragma: no cover
