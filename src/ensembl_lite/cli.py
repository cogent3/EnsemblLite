import pathlib
import shutil

import click


try:
    from wakepy.keep import running as keep_running
except ImportError:
    from ensembl_lite._util import fake_wake as keep_running

from trogon import tui

from ensembl_lite import __version__
from ensembl_lite import _config as elt_config
from ensembl_lite import _download as elt_download
from ensembl_lite._species import Species
from ensembl_lite._util import PathType


try:
    # trap flaky behaviour on linux
    with keep_running():
        ...

except NotImplementedError:
    from ensembl_lite._util import fake_wake as keep_running


def _get_installed_config_path(ctx, param, path) -> PathType:
    """path to installed.cfg"""
    path = pathlib.Path(path)
    if path.name == elt_config.INSTALLED_CONFIG_NAME:
        return path

    path = path / elt_config.INSTALLED_CONFIG_NAME
    if not path.exists():
        click.secho(f"{str(path)} missing", fg="red")
        exit(1)
    return path


def _values_from_csv(ctx, param, value) -> list[str] | None:
    if value is None:
        return

    return [f.strip() for f in value.split(",")]


def _species_names_from_csv(ctx, param, species) -> list[str] | None:
    """returns species names"""
    species = _values_from_csv(ctx, param, species)
    if species is None:
        return

    db_names = []
    for name in species:
        try:
            db_name = Species.get_ensembl_db_prefix(name)
        except ValueError:
            click.secho(f"ERROR: unknown species {name!r}", fg="red")
            exit(1)

        db_names.append(db_name)

    return db_names


# defining some of the options
_cfgpath = click.option(
    "-c",
    "--configpath",
    default=elt_download._cfg,
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
_installed = click.option(
    "-i",
    "--installed",
    required=True,
    callback=_get_installed_config_path,
    help="string pointing to installation",
)
_outpath = click.option(
    "-o", "--outpath", required=True, type=pathlib.Path, help="path to write json file"
)
_outdir = click.option(
    "-od", "--outdir", required=True, type=pathlib.Path, help="path to write files"
)
_align_name = click.option(
    "--align_name",
    default=None,
    help="Ensembl name of the alignment or a glob pattern, e.g. '*primates*'",
)
_ref = click.option("--ref", default=None, help="Reference species.")
_ref_genes_file = click.option(
    "--ref_genes_file",
    default=None,
    type=click.Path(resolve_path=True, exists=True),
    help=".csv or .tsv file with a header containing a stableid column",
)
_limit = click.option(
    "--limit",
    type=int,
    default=None,
    help="Limit to this number of genes.",
    show_default=True,
)

_verbose = click.option(
    "-v",
    "--verbose",
    is_flag=True,
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
_nprocs = click.option(
    "-np",
    "--num_procs",
    type=int,
    default=1,
    help="number of procs to use, defaults to 1",
)


_outdir = click.option(
    "--outdir",
    type=pathlib.Path,
    default=".",
    help="Output directory name.",
    show_default=True,
)

_species = click.option(
    "--species",
    required=True,
    callback=_species_names_from_csv,
    help="Single species name, or multiple (comma separated).",
)

_mask_features = click.option(
    "--mask_features",
    callback=_values_from_csv,
    help="biotypes to mask (comma separated).",
)


@tui()
@click.group()
@click.version_option(__version__)
def main():
    """tools for obtaining and interrogating subsets of https://ensembl.org genomic data"""
    pass


@main.command(no_args_is_help=True)
@_dbrc_out
def exportrc(outpath):
    """exports sample config and species table to the nominated path

    setting an environment variable ENSEMBLDBRC with this path
    will force its contents to override the default ensembl_lite settings"""
    from ensembl_lite._util import ENSEMBLDBRC

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


@main.command(no_args_is_help=True)
@_cfgpath
@_debug
@_verbose
def download(configpath, debug, verbose):
    """download data from Ensembl's ftp site"""
    if configpath.name == elt_download._cfg:
        # todo is this statement correct if we're seting a root dir now?
        click.secho(
            "WARN: using the built in demo cfg, will write to /tmp", fg="yellow"
        )
    config = elt_config.read_config(configpath, root_dir=pathlib.Path(".").resolve())

    if verbose:
        print(config)

    if not any((config.species_dbs, config.align_names)):
        click.secho("No genomes, no alignments specified", fg="red")
        exit(1)

    if not config.species_dbs:
        species = elt_download.get_species_for_alignments(
            host=config.host,
            remote_path=config.remote_path,
            release=config.release,
            align_names=config.align_names,
        )
        config.update_species(species)

    if verbose:
        print(config.species_dbs)

    config.write()
    with keep_running():
        elt_download.download_species(config, debug, verbose)
        elt_download.download_homology(config, debug, verbose)
        elt_download.download_aligns(config, debug, verbose)

    click.secho(f"Downloaded to {config.staging_path}", fg="green")


@main.command(no_args_is_help=True)
@_download
@_nprocs
@_force
@_verbose
def install(download, num_procs, force_overwrite, verbose):
    """create the local representations of the data"""
    from rich import progress

    from ensembl_lite._install import (
        local_install_alignments,
        local_install_genomes,
        local_install_homology,
    )

    configpath = download / elt_config.DOWNLOADED_CONFIG_NAME
    config = elt_config.read_config(configpath)
    if verbose:
        print(f"{config.install_path=}")

    if force_overwrite:
        shutil.rmtree(config.install_path, ignore_errors=True)

    config.install_path.mkdir(parents=True, exist_ok=True)
    elt_config.write_installed_cfg(config)
    with keep_running():
        with progress.Progress(
            progress.TextColumn("[progress.description]{task.description}"),
            progress.BarColumn(),
            progress.TaskProgressColumn(),
            progress.TimeRemainingColumn(),
            progress.TimeElapsedColumn(),
        ) as progress:
            local_install_genomes(
                config,
                force_overwrite=force_overwrite,
                max_workers=num_procs,
                verbose=verbose,
                progress=progress,
            )
            # On test cases, only 30% speedup from running install homology data
            # in parallel due to overhead of pickling the data, but considerable
            # increase in memory. So, run in serial to avoid memory issues since
            # it's reasonably fast anyway. (At least until we have
            # a more robust solution.)
            local_install_homology(
                config,
                force_overwrite=force_overwrite,
                max_workers=num_procs,
                verbose=verbose,
                progress=progress,
            )
            local_install_alignments(
                config,
                force_overwrite=force_overwrite,
                max_workers=num_procs,
                verbose=verbose,
                progress=progress,
            )

    click.secho(f"Contents installed to {str(config.install_path)!r}", fg="green")


@main.command(no_args_is_help=True)
@_installed
def installed(installed):
    """show what is installed"""
    from cogent3 import make_table

    from ensembl_lite._species import Species
    from ensembl_lite._util import rich_display

    config = elt_config.read_installed_cfg(installed)

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


@main.command(no_args_is_help=True)
@_installed
@_species
def species_summary(installed, species):
    """genome summary data for a species"""
    from ._genomedb import (
        _ANNOTDB_NAME,
        get_species_summary,
        load_annotations_for_species,
    )
    from ._util import rich_display

    config = elt_config.read_installed_cfg(installed)
    if species is None:
        click.secho("ERROR: a species name is required", fg="red")
        exit(1)

    if len(species) > 1:
        click.secho(f"ERROR: one species at a time, not {species!r}", fg="red")
        exit(1)

    species = species[0]
    path = config.installed_genome(species=species) / _ANNOTDB_NAME
    if not path.exists():
        click.secho(f"{species!r} not in {str(config.install_path.parent)!r}", fg="red")
        exit(1)

    annot_db = load_annotations_for_species(path=path)
    summary = get_species_summary(annot_db=annot_db, species=species)
    rich_display(summary)


@main.command(no_args_is_help=True)
@_installed
@_outdir
@_align_name
@_ref
@_ref_genes_file
@_mask_features
@_limit
@_force
@_verbose
def alignments(
    installed,
    outdir,
    align_name,
    ref,
    ref_genes_file,
    mask_features,
    limit,
    force_overwrite,
    verbose,
):
    """dump alignments for named genes"""
    from cogent3 import load_table

    from ensembl_lite._aligndb import AlignDb, write_alignments
    from ensembl_lite._genomedb import load_genome, update_stableid_prefixes
    from ensembl_lite._species import Species

    # todo support genomic coordinates, e.g. coord_name:start-stop:strand, for
    #  a reference species

    if not ref:
        click.secho(
            "ERROR: must specify a reference genome",
            fg="red",
        )
        exit(1)

    if force_overwrite:
        shutil.rmtree(outdir, ignore_errors=True)

    outdir.mkdir(parents=True, exist_ok=True)

    config = elt_config.read_installed_cfg(installed)
    # update the prefixes
    update_stableid_prefixes(config)
    align_path = config.path_to_alignment(align_name)
    if align_path is None:
        click.secho(
            f"{align_name!r} does not match any alignments under {str(config.aligns_path)!r}",
            fg="red",
        )
        exit(1)

    # load the gene stable ID's
    table = load_table(ref_genes_file)
    if "stableid" not in table.columns:
        click.secho(
            f"'stableid' column missing from {str(ref_genes_file)!r}",
            fg="red",
        )
        exit(1)

    align_db = AlignDb(source=align_path)
    ref_species = Species.get_ensembl_db_prefix(ref)
    if ref_species not in align_db.get_species_names():
        click.secho(
            f"species {ref!r} not in the alignment",
            fg="red",
        )
        exit(1)

    # get all the genomes
    if verbose:
        print(f"working on species {align_db.get_species_names()}")

    genomes = {
        sp: load_genome(config=config, species=sp)
        for sp in align_db.get_species_names()
    }

    write_alignments(
        align_db=align_db,
        genomes=genomes,
        limit=limit,
        mask_features=mask_features,
        outdir=outdir,
        ref_species=ref_species,
        stableids=table.columns["stableid"],
    )

    click.secho("Done!", fg="green")


@main.command(no_args_is_help=True)
@_installed
@_outpath
@click.option(
    "-r",
    "--relationship",
    type=click.Choice(["ortholog_one2one"]),
    default="ortholog_one2one",
    help="type of homology",
)
@_ref
@_nprocs
@_limit
@_force
@_verbose
def homologs(
    installed, outpath, relationship, ref, num_procs, limit, force_overwrite, verbose
):
    """exports all homolog groups of type relationship in fasta format"""
    from rich.progress import Progress, track

    from ensembl_lite._genomedb import load_genome
    from ensembl_lite._homologydb import (
        _HOMOLOGYDB_NAME,
        collect_seqs,
        load_homology_db,
    )

    if ref is None:
        click.secho("ERROR: a reference species name is required, use --ref", fg="red")
        exit(1)

    if force_overwrite:
        shutil.rmtree(outpath, ignore_errors=True)

    outpath.mkdir(parents=True, exist_ok=True)

    config = elt_config.read_installed_cfg(installed)
    Species.update_from_file(config.genomes_path / "species.tsv")
    # we all the protein coding gene IDs from the reference species
    genome = load_genome(config=config, species=ref)
    if verbose:
        print(f"loaded genome for {ref}")
    gene_ids = list(genome.get_ids_for_biotype(biotype="gene"))
    if verbose:
        print(f"found {len(gene_ids)} gene IDs for {ref}")
    db = load_homology_db(path=config.homologies_path / _HOMOLOGYDB_NAME)
    related = []
    for gid in track(gene_ids, description="Homolog search"):
        if rel := db.get_related_to(gene_id=gid, relationship_type=relationship):
            related.append(rel)

        if limit and len(related) >= limit:
            break

    if verbose:
        print(f"Found {len(related)} homolog groups")
    # todo create a directory data store writer and write all output to
    #  that. This requires homolog_group has a .source attribute
    get_seqs = collect_seqs(config=config)
    with Progress(transient=False) as progress:
        reading = progress.add_task(total=len(related), description="Extracting  🧬")
        for i, seqs in enumerate(
            get_seqs.as_completed(
                related,
                parallel=True,
                show_progress=False,
                par_kw=dict(max_workers=num_procs),
            )
        ):
            progress.update(reading, advance=1)
            if not seqs:
                if verbose:
                    print(f"{seqs=}")
                continue
            if not seqs.obj.seqs:
                if verbose:
                    print(f"{seqs.obj.seqs=}")
                continue

            # todo also need to be writing out a logfile, plus a meta data table of
            #  gene IDs and location info
            txt = [seq.to_fasta() for seq in seqs.obj.seqs]
            outname = outpath / f"seqcoll-{i}.fasta"
            with outname.open(mode="w") as outfile:
                outfile.write("".join(txt))


@main.command(no_args_is_help=True)
@_installed
@_species
@_outdir
@_limit
def dump_genes(installed, species, outdir, limit):
    """Dump meta data table for genes from one species to <species>-<release>.gene_metadata.tsv"""
    from ensembl_lite._genomedb import (
        _ANNOTDB_NAME,
        get_gene_table_for_species,
        load_annotations_for_species,
    )

    config = elt_config.read_installed_cfg(installed)
    if species is None:
        click.secho("ERROR: a species name is required", fg="red")
        exit(1)

    if len(species) > 1:
        click.secho(f"ERROR: one species at a time, not {species!r}", fg="red")
        exit(1)

    path = config.installed_genome(species=species[0]) / _ANNOTDB_NAME
    if not path.exists():
        click.secho(f"{species!r} not in {str(config.install_path.parent)!r}", fg="red")
        exit(1)

    annot_db = load_annotations_for_species(path=path)
    path = annot_db.source
    table = get_gene_table_for_species(annot_db=annot_db, limit=limit)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"{path.parent.stem}-{config.release}-gene_metadata.tsv"
    table.write(outpath)
    click.secho(f"Finished: wrote {str(outpath)!r}!", fg="green")


if __name__ == "__main__":
    main()
