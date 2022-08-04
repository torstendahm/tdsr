# -*- coding: utf-8 -*-

################################
# Time Dependent Seismicity Model - CLI for tdsr
# T. Dahm, R. Dahm 26.12.2021
################################


"""Console script for tdsr."""
import sys
from pathlib import Path
from typing import Optional

import click

from tdsr import LCM, TDSR, Traditional, save
from tdsr.config import Config, DEFAULT_CONFIG
from tdsr.utils import PathLike


def optional_valid_dir_or_file(
    ctx: click.core.Context, param: click.core.Parameter, value: str
) -> str:
    if value is None or value == "":
        return value
    if Path(value).exists() and Path(value).is_file():
        return value
    raise click.BadParameter("%s is not a file" % (value), ctx=ctx, param=param)


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
                    "%s already exists. Use -f to force overwriting. Do you want to overwrite it?"
                    % _file,
                    default=False,
                )
            ):
                raise click.Abort("%s already exists" % _file)
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
    "--force",
    is_flag=True,
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
def tdsr(
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
    """tdsr method"""
    output_file: Optional[PathLike] = None
    input_file: Optional[PathLike] = None
    if output:
        # try:
        output_path = Path(output).resolve()
        output_file = get_output_file(output_path, name="tdsr")
        output_file = check_output_file(output_file, force=force)
        # except FileExistsError:
        #     raise Click

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
        tdsr = TDSR(config=conf)
        result = tdsr()
        if output_file:
            save(result, output_file)
            print("saved to ", output_file)
    else:
        ctx.ensure_object(dict)
        ctx.obj["PARAMS"] = dict(
            config=conf,
            input_file=input_file,
            output_file=output_file,
            force=force,
        )


@tdsr.command()
@click.pass_context
def lcm(
    ctx: click.Context,
) -> int:
    """LCM method"""
    params = ctx.obj["PARAMS"]
    config = params["config"]
    output_file = get_output_file(params["output_path"], name="lcm")
    output_file = check_output_file(output_file, force=params["force"])
    lcm = LCM(config=config)
    result = lcm()
    print(result)
    return 0


@tdsr.command()
@click.pass_context
def traditional(
    ctx: click.Context,
) -> int:
    """Traditional method"""
    params = ctx.obj["PARAMS"]
    config = params["config"]
    output_file = get_output_file(params["output_path"], name="lcm")
    output_file = check_output_file(output_file, force=params["force"])
    traditional = Traditional(config=config)
    result = traditional()
    print(result)
    return 0


if __name__ == "__main__":
    sys.exit(tdsr(obj=dict()))  # pragma: no cover
