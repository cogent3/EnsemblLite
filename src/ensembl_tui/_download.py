import pathlib
import re
import shutil
import typing

import cogent3
from rich.progress import Progress

from ensembl_tui import _config as eti_config
from ensembl_tui import _ftp_download as eti_ftp
from ensembl_tui import _mysql_core_attr as eti_db_attr
from ensembl_tui import _mysql_ingest as eti_db_ingest
from ensembl_tui import _name as eti_name
from ensembl_tui import _site_map as eti_site_map
from ensembl_tui import _species as eti_species
from ensembl_tui import _util as eti_util

DEFAULT_CFG = eti_util.get_resource_path("sample.cfg")

_invalid_seq = re.compile("(dna_(sm|rm)|(toplevel|primary_assembly).fa.gz)")


def valid_seq_file(name: str) -> bool:
    """unmasked genomic DNA sequences"""
    return _invalid_seq.search(name) is None


class valid_gff3_file:  # noqa: N801
    """whole genome gff3"""

    def __init__(self, release: str) -> None:
        self._valid = re.compile(f"([.]{release}[.]gff3[.]gz|README|CHECKSUMS)")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


def _remove_tmpdirs(path: eti_util.PathType) -> None:
    """delete any tmp dirs left over from unsuccessful runs"""
    tmpdirs = [p for p in path.glob("tmp*") if p.is_dir()]
    for tmpdir in tmpdirs:
        shutil.rmtree(tmpdir)


def get_core_db_dirnames(config: eti_config.Config) -> dict[str, str]:
    """maps species name to ftp path to mysql core dbs"""
    all_db_names = list(
        eti_ftp.listdir(config.host, f"{config.remote_release_path}/mysql"),
    )
    selected_species = {}
    for db_name in all_db_names:
        if "_core_" not in db_name:
            continue
        db = eti_name.EnsemblDbName(db_name.rsplit("/", maxsplit=1)[1])
        if db.species in config.species_dbs and db.db_type == "core":
            selected_species[db.prefix] = db_name
    return selected_species


def get_remote_mysql_paths(db_name: str) -> list[str]:
    return [f"{db_name}/{name}" for name in eti_db_attr.make_mysqldump_names()]


def make_core_db_templates(
    config: eti_config.Config,
    sp_db_map: dict[str, str],
    progress: Progress | None = None,
) -> None:
    """creates duckdb db files for importing Ensembl mysql data

    Parameters
    ----------
    config
        eti configuration
    sp_db_map
        mapping of species to mysql db names
    progress
        rich.progress context manager for tracking progress

    Notes
    -----
    Communicates with the Ensembl MySQL server to infer the table schema's.
    """
    table_names = eti_db_attr.get_all_tables()
    if progress is not None:
        msg = "Making db templates"
        make_templates = progress.add_task(
            total=len(table_names),
            description=msg,
        )

    template_dest = config.staging_template_path
    # get one species db name which wqe use to infer the db schema
    db_name = next(iter(sp_db_map.values())).split("/")[-1]
    template_dest.mkdir(parents=True, exist_ok=True)
    for table_name in table_names:
        eti_db_ingest.make_table_template(
            dest_dir=template_dest,
            db_name=db_name,
            table_name=table_name,
        )
        if progress is not None:
            progress.update(make_templates, description=msg, advance=1)


def download_species(
    config: eti_config.Config,
    debug: bool,
    verbose: bool,
    progress: Progress | None = None,
) -> None:
    """download seq and annotation data"""
    remote_template = f"{config.remote_release_path}/" + "{}"
    site_map = eti_site_map.get_site_map(config.host)
    if verbose:
        eti_util.print_colour(
            text=f"DOWNLOADING\n  ensembl release={config.release}",
            colour="green",
        )
        eti_util.print_colour(
            text="\n".join(f"  {d}" for d in config.species_dbs),
            colour="green",
        )
        eti_util.print_colour(
            text=f"\nWRITING to output path={config.staging_genomes}\n",
            colour="green",
        )

    sp_db_map = get_core_db_dirnames(config)

    # create the duckdb templates for the tables, if they don't exist
    make_core_db_templates(config, sp_db_map, progress=progress)

    msg = "Downloading genomes"
    if progress is not None:
        species_download = progress.add_task(
            total=len(config.species_dbs),
            description=msg,
        )

    patterns = {"fasta": valid_seq_file, "gff3": valid_gff3_file(config.release)}
    for key in config.species_dbs:
        db_prefix = eti_species.Species.get_ensembl_db_prefix(key)
        local_root = config.staging_genomes / db_prefix
        local_root.mkdir(parents=True, exist_ok=True)
        # getting genome sequences
        remote = site_map.get_seqs_path(db_prefix)
        remote_dir = remote_template.format(remote)
        remote_paths = list(
            eti_ftp.listdir(config.host, path=remote_dir, pattern=patterns["fasta"]),
        )
        if verbose:
            eti_util.print_colour(text=f"{remote_paths=}", colour="yellow")

        if debug:
            # we need the checksum files
            paths = [p for p in remote_paths if eti_util.is_signature(p)]
            # but fewer data files, to reduce time for debugging
            remote_paths = [p for p in remote_paths if not eti_util.dont_checksum(p)]
            remote_paths = remote_paths[:4] + paths

        dest_path = config.staging_genomes / db_prefix / "fasta"
        dest_path.mkdir(parents=True, exist_ok=True)
        # cleanup previous download attempts
        _remove_tmpdirs(dest_path)
        icon = "🧬🧬"
        eti_ftp.download_data(
            host=config.host,
            local_dest=dest_path,
            remote_paths=remote_paths,
            description=f"{db_prefix[:10]}... {icon}",
            do_checksum=True,
            progress=progress,
        )

        # getting the annotations from mysql tables
        remote_dir = sp_db_map[db_prefix]
        remote_paths = get_remote_mysql_paths(remote_dir)
        dest_path = config.staging_genomes / db_prefix / "mysql"
        dest_path.mkdir(parents=True, exist_ok=True)
        # cleanup previous download attempts
        _remove_tmpdirs(dest_path)
        icon = "📚"
        eti_ftp.download_data(
            host=config.host,
            local_dest=dest_path,
            remote_paths=remote_paths,
            description=f"{db_prefix[:10]}... {icon}",
            do_checksum=True,
            progress=progress,
        )

        if progress is not None:
            progress.update(species_download, description=msg, advance=1)


class valid_compara_align:  # noqa: N801
    """whole genome alignment data"""

    def __init__(self) -> None:
        self._valid = re.compile("([.](emf|maf)[.]gz|README|MD5SUM)")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


def download_aligns(
    config: eti_config.Config,
    debug: bool,
    verbose: bool,
    progress: Progress | None = None,
) -> None:
    """download whole genome alignments"""
    if not config.align_names:
        return

    site_map = eti_site_map.get_site_map(config.host)
    remote_template = (
        f"{config.remote_path}/release-{config.release}/{site_map.alignments_path}/"
        "{}"
    )

    msg = "Downloading alignments"
    if progress is not None:
        align_download = progress.add_task(
            total=len(config.species_dbs),
            description=msg,
        )

    valid_compara = valid_compara_align()
    for align_name in config.align_names:
        remote_path = remote_template.format(align_name)
        remote_paths = list(eti_ftp.listdir(config.host, remote_path, valid_compara))
        if verbose:
            print(remote_paths)

        if debug:
            # we need the checksum files
            paths = [p for p in remote_paths if eti_util.is_signature(p)]
            remote_paths = [p for p in remote_paths if not eti_util.is_signature(p)]
            remote_paths = remote_paths[:4] + paths

        local_dir = config.staging_aligns / align_name
        local_dir.mkdir(parents=True, exist_ok=True)
        _remove_tmpdirs(local_dir)
        eti_ftp.download_data(
            host=config.host,
            local_dest=local_dir,
            remote_paths=remote_paths,
            description=f"{align_name[:10]}...",
            do_checksum=True,
            progress=progress,
        )

        if progress is not None:
            progress.update(align_download, description=msg, advance=1)

    return


class valid_compara_homology:  # noqa: N801
    """homology tsv files"""

    def __init__(self) -> None:
        self._valid = re.compile("([.]tsv|[.]tsv[.]gz|README|MD5SUM)$")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


def download_homology(
    config: eti_config.Config,
    debug: bool,
    verbose: bool,
    progress: Progress | None = None,
) -> None:
    """downloads tsv homology files for each genome"""
    if not config.homologies:
        return

    site_map = eti_site_map.get_site_map(config.host)
    remote_template = (
        f"{config.remote_path}/release-{config.release}/{site_map.homologies_path}/"
        "{}"
    )

    local = config.staging_homologies

    msg = "Downloading homology"
    if progress is not None:
        species_download = progress.add_task(
            total=len(config.species_dbs),
            description=msg,
        )

    for db_name in config.db_names:
        remote_path = remote_template.format(db_name)
        remote_paths = list(
            eti_ftp.listdir(config.host, remote_path, valid_compara_homology()),
        )
        if verbose:
            print(f"{remote_path=}", f"{remote_paths=}", sep="\n")

        if debug:
            # we need the checksum files
            remote_paths = [p for p in remote_paths if not eti_util.is_signature(p)]
            remote_paths = remote_paths[:4]

        local_dir = local / db_name
        local_dir.mkdir(parents=True, exist_ok=True)
        _remove_tmpdirs(local_dir)
        eti_ftp.download_data(
            host=config.host,
            local_dest=local_dir,
            remote_paths=remote_paths,
            description=f"{db_name[:10]}...",
            do_checksum=False,  # no checksums for species homology files
            progress=progress,
        )

        if progress is not None:
            progress.update(species_download, description=msg, advance=1)

    return


def download_ensembl_tree(
    host: str,
    remote_path: str,
    release: str,
    tree_fname: str,
) -> cogent3.core.tree.PhyloNode:
    """loads a tree from Ensembl"""
    site_map = eti_site_map.get_site_map(host)
    url = f"https://{host}/{remote_path}/release-{release}/{site_map.trees_path}/{tree_fname}"
    return cogent3.load_tree(url)


def get_ensembl_trees(host: str, remote_path: str, release: str) -> list[str]:
    """returns trees from ensembl compara"""
    site_map = eti_site_map.get_site_map(host)
    path = f"{remote_path}/release-{release}/{site_map.trees_path}"
    return list(
        eti_ftp.listdir(host=host, path=path, pattern=lambda x: x.endswith(".nh")),
    )


def get_species_for_alignments(
    host: str,
    remote_path: str,
    release: str,
    align_names: typing.Iterable[str],
) -> dict[str, list[str]]:
    """return the species for the indicated alignments"""
    ensembl_trees = get_ensembl_trees(
        host=host,
        remote_path=remote_path,
        release=release,
    )
    aligns_trees = eti_util.trees_for_aligns(align_names, ensembl_trees)
    species = {}
    for tree_path in aligns_trees.values():
        tree = download_ensembl_tree(
            host=host,
            remote_path=remote_path,
            release=release,
            tree_fname=pathlib.Path(tree_path).name,
        )
        # dict structure is {common name: db prefix}, just use common name
        species |= {n: ["core"] for n in eti_species.species_from_ensembl_tree(tree)}
    return species
