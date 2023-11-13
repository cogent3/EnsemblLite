import pathlib
import shutil

import click
import wakepy.keep

from rich.progress import track
from trogon import tui

from ensembl_lite import __version__
from ensembl_lite._config import (
    DOWNLOADED_CONFIG_NAME,
    INSTALLED_CONFIG_NAME,
    read_config,
    read_installed_cfg,
    write_installed_cfg,
)
from ensembl_lite.download import (
    _cfg,
    download_aligns,
    download_homology,
    download_species,
    get_species_for_alignments,
)


# defining some of the options
_cfgpath = click.option(
    "-c",
    "--configpath",
    default=_cfg,
    type=pathlib.Path,
    help="path to config file specifying databases, only "
    "species or compara at present",
)
_download = click.option(
    "-d",
    "--download",
    type=pathlib.Path,
    help="path to local download directory, contains a cfg file",
)
_installation = click.option(
    "--installation",
    type=pathlib.Path,
    help="path to local installation directory",
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


@tui()
@click.group()
@click.version_option(__version__)
def main():
    """tools for obtaining and interrogating subsets of https://ensembl.org genomic data"""
    pass


@main.command(no_args_is_help=True)
@_cfgpath
@_debug
@_verbose
def download(configpath, debug, verbose):
    """download data from Ensembl's ftp site"""
    if configpath.name == _cfg:
        click.secho(
            "WARN: using the built in demo cfg, will write to /tmp", fg="yellow"
        )

    config = read_config(configpath)
    if not any((config.species_dbs, config.align_names)):
        click.secho("No genomes, no alignments specified")
        exit(1)

    if not config.species_dbs:
        species = get_species_for_alignments(
            host=config.host,
            remote_path=config.remote_path,
            release=config.release,
            align_names=config.align_names,
        )
        config.update_species(species)

    if verbose:
        print(config.species_dbs)

    config.write()
    with wakepy.keep.running():
        download_species(config, debug, verbose)
        download_homology(config, debug, verbose)
        download_aligns(config, debug, verbose)

    click.secho(f"Downloaded to {config.staging_path}", fg="green")


@main.command(no_args_is_help=True)
@_download
@_force
@_verbose
def install(download, force_overwrite, verbose):
    """create the local representations of the data"""
    from ensembl_lite.install import (
        local_install_compara,
        local_install_genomes,
        local_install_homology,
    )

    configpath = download / DOWNLOADED_CONFIG_NAME
    config = read_config(configpath)
    if verbose:
        print(f"{config.install_path=}")

    if force_overwrite:
        shutil.rmtree(config.install_path, ignore_errors=True)

    config.install_path.mkdir(parents=True, exist_ok=True)
    write_installed_cfg(config)

    with wakepy.keep.running():
        local_install_genomes(config, force_overwrite=force_overwrite)
        local_install_compara(config, force_overwrite=force_overwrite)
        local_install_homology(config, force_overwrite=force_overwrite)

    click.secho(f"Contents installed to {str(config.install_path)!r}", fg="green")


@main.command(no_args_is_help=True)
@_dbrc_out
def exportrc(outpath):
    """exports sample config and species table to the nominated path

    setting an environment variable ENSEMBLDBRC with this path
    will force its contents to override the default ensembl_lite settings"""
    from ensembl_lite.util import ENSEMBLDBRC

    shutil.copytree(ENSEMBLDBRC, outpath)
    # we assume all files starting with alphabetical characters are valid
    for fn in pathlib.Path(outpath).glob("*"):
        if not fn.stem.isalpha():
            if fn.is_file():
                fn.unlink()
            else:
                # __pycache__ directory
                shutil.rmtree(fn)
    click.secho(f"Contents written to {outpath}", fg="green")


def _get_installed_config_path(ctx, param, path):
    """path to installed.cfg"""
    path = pathlib.Path(path)
    if path.name == INSTALLED_CONFIG_NAME:
        return path

    path = path / INSTALLED_CONFIG_NAME
    if not path.exists():
        click.secho(f"{str(path)} missing", fg="red")
        exit(1)
    return path


_installed = click.option(
    "-i",
    "--installed",
    required=True,
    callback=_get_installed_config_path,
    help="string pointing to installation",
)


@main.command(no_args_is_help=True)
@_installed
@click.option(
    "-o", "--outpath", required=True, type=pathlib.Path, help="path to write json file"
)
@click.option(
    "-r",
    "--relationship",
    type=click.Choice(["ortholog_one2one", "ortholog_one2many"]),
    default="ortholog_one2one",
    help="type of homology",
)
def homologs(installed, outpath, relationship):
    """exports all homolog groups of type relationship in json format"""
    import json

    from cogent3 import open_

    from ensembl_lite._homologydb import HomologyDb

    config = read_installed_cfg(installed)
    db_path = config.install_homologies / "homologies.sqlitedb"
    db = HomologyDb(source=db_path)
    related = list(db.get_related_groups(relationship))
    with open_(outpath, mode="wt") as out:
        json.dump(related, out)


@main.command(no_args_is_help=True)
@_installed
def installed(installed):
    """show what is installed"""
    from cogent3 import make_table

    from ensembl_lite.species import Species
    from ensembl_lite.util import rich_display

    config = read_installed_cfg(installed)

    genome_dir = config.genomes_path
    if genome_dir.exists():
        species = [fn.name for fn in genome_dir.glob("*")]
        data = {"species": [], "common name": []}
        for name in species:
            cn = Species.get_common_name(name, level="ignore")
            if not cn:
                continue
            data["species"].append(name)
            data["common name"].append(cn)

        table = make_table(data=data, title="Installed genomes")
        rich_display(table)

    # TODO as above
    compara_aligns = config.aligns_path
    if compara_aligns.exists():
        align_names = [
            fn.stem for fn in compara_aligns.glob("*") if not fn.name.startswith(".")
        ]
        table = make_table(
            data={"align name": align_names}, title="Installed whole genome alignments"
        )
        rich_display(table)


def _species_names_from_csv(ctx, param, species):
    """returns species names"""
    if species is not None:
        species = [s.strip().lower() for s in species.split(",")]
    return species


_species = click.option(
    "--species",
    required=True,
    callback=_species_names_from_csv,
    help="Single species name, or multiple (comma separated).",
)
_outdir = click.option(
    "--outdir",
    type=pathlib.Path,
    required=True,
    default="gene_metadata.tsv",
    help="Output file name.",
)
_limit = click.option(
    "--limit",
    type=int,
    default=None,
    help="Limit to this number of genes.",
    show_default=True,
)


@main.command(no_args_is_help=True)
@_installed
@_species
@_outdir
@_limit
def dump_genes(installed, species, outdir, limit):
    """Dump meta data table for genes from one species to <species>-<release>.gene_metadata.tsv"""
    from cogent3 import make_table
    from cogent3.core.annotation_db import GffAnnotationDb

    from ensembl_lite.species import Species

    config = read_installed_cfg(installed)
    species = species[0]
    path = config.installed_genome(species=species)
    if not path.exists():
        click.secho(f"{species!r} not in {str(installed.parent)!r}", fg="red")
        exit(1)

    # TODO: this filename should be defined in one place
    path = path / "features.gff3db"
    if not path.exists():
        click.secho(f"{path.name!r} is missing", fg="red")
        exit(1)

    annot_db = GffAnnotationDb(source=path)
    rows = []
    columns = [
        "name",
        "seqid",
        "source",
        "biotype",
        "start",
        "end",
        "score",
        "strand",
        "phase",
    ]
    for i, record in track(enumerate(annot_db.get_records_matching(biotype="gene"))):
        rows.append([record[c] for c in columns])
        if i == limit:
            break

    table = make_table(header=columns, data=rows)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"{path.parent.stem}-{config.release}-gene_metadata.tsv"
    table.write(outpath)
    click.secho(f"Finished wrote {str(outpath)!r}!", fg="green")


if __name__ == "__main__":
    main()
