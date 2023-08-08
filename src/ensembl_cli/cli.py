import os
import pathlib
import shutil

from pprint import pprint

import click
import wakepy.keep

from ensembl_cli import __version__
from ensembl_cli.download import _cfg, download_aligns, download_species
from ensembl_cli.util import read_config


def listpaths(dirname, glob_pattern):
    """return path to all files matching glob_pattern"""
    fns = [str(p) for p in pathlib.Path(dirname).glob(glob_pattern)]
    return fns if fns else None


def sorted_by_size(local_path, dbnames, debug=False):
    """returns dbnames ordered by directory size"""
    join = os.path.join
    getsize = os.path.getsize
    size_dbnames = []
    for dbname in dbnames:
        path = join(local_path, dbname.name)
        size = sum(getsize(join(path, fname)) for fname in os.listdir(path))
        size_dbnames.append([size, dbname])
    size_dbnames.sort()

    if debug:
        pprint(size_dbnames)

    sizes, dbnames = zip(*size_dbnames)
    return dbnames


# defining some of the options
_cfgpath = click.option(
    "-c",
    "--configpath",
    default=_cfg,
    type=pathlib.Path,
    help="path to config file specifying databases, only "
    "species or compara at present",
)
_verbose = click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="causes stdout/stderr from rsync download to be " "written to screen",
)
_numprocs = click.option(
    "-n",
    "--numprocs",
    type=int,
    default=1,
    help="number of processes to use for download",
)
_force = click.option(
    "-f",
    "--force_overwrite",
    is_flag=True,
    help="drop existing database if it exists prior to " "installing",
)
_debug = click.option(
    "-d",
    "--debug",
    is_flag=True,
    help="maximum verbosity, and reduces number of downloads",
)
_dbrc_out = click.option(
    "-o",
    "--outpath",
    type=pathlib.Path,
    help="path to directory to export all rc contents",
)
_release = click.option("-r", "--release", type=int, help="Ensembl release number")


@click.group()
@click.version_option(__version__)
def main():
    """admin tools for an Ensembl MySQL installation"""
    pass


@main.command()
@_cfgpath
@_debug
@_verbose
def download(configpath, debug, verbose):
    """download databases from Ensembl using rsync, can be done in parallel"""
    if configpath.name == _cfg:
        click.secho(
            "WARN: using the built in demo cfg, will write to /tmp", fg="yellow"
        )

    config = read_config(configpath)

    with wakepy.keep.running():
        download_species(config, debug, verbose)
        download_aligns(config, debug, verbose)

    click.secho(f"Downloaded to {config.staging_path}", fg="green")


@main.command()
@_cfgpath
@_force
@_verbose
def install(configpath, force_overwrite, verbose):
    """create the local db's"""
    from ensembl_cli.install import (
        local_install_compara,
        local_install_genomes,
    )

    if configpath.name == _cfg:
        click.secho(
            "WARN: using the built in demo cfg, will write to /tmp", fg="yellow"
        )

    config = read_config(configpath)
    if force_overwrite:
        shutil.rmtree(config.install_path, ignore_errors=True)

    with wakepy.keep.running():
        local_install_genomes(config)
        local_install_compara(config)

    click.secho(f"Contents installed to {str(config.install_path)!r}", fg="green")


@main.command()
@_dbrc_out
def exportrc(outpath):
    """exports the rc directory to the nominated path

    setting an environment variable ENSEMBLDBRC with this path
    will force its contents to override the default ensembl_cli settings"""
    from ensembl_cli.util import ENSEMBLDBRC

    shutil.copytree(ENSEMBLDBRC, outpath)
    # remove the python module file
    for fn in pathlib.Path(outpath).glob("__init__.py*"):
        fn.unlink()
    click.secho(f"Contents written to {outpath}", fg="green")


if __name__ == "__main__":
    main()
