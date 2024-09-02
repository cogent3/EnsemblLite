import pathlib
import re
import shutil

import click
from cogent3 import load_tree
from rich.progress import Progress

from ensembl_tui import _config as elt_config
from ensembl_tui import _ftp_download as elt_ftp
from ensembl_tui import _site_map as elt_site_map
from ensembl_tui import _species as elt_species
from ensembl_tui import _util as elt_util

_cfg = elt_util.get_resource_path("sample.cfg")

_invalid_seq = re.compile("(dna_(sm|rm)|(toplevel|primary_assembly).fa.gz)")


def valid_seq_file(name: str) -> bool:
    """unmasked genomic DNA sequences"""
    return _invalid_seq.search(name) is None


class valid_gff3_file:
    """whole genome gff3"""

    def __init__(self, release: str) -> None:
        self._valid = re.compile(f"([.]{release}[.]gff3[.]gz|README|CHECKSUMS)")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


def _remove_tmpdirs(path: elt_util.PathType):
    """delete any tmp dirs left over from unsuccessful runs"""
    tmpdirs = [p for p in path.glob("tmp*") if p.is_dir()]
    for tmpdir in tmpdirs:
        shutil.rmtree(tmpdir)


def download_species(
    config: elt_config.Config,
    debug: bool,
    verbose: bool,
    progress: Progress | None = None,
):
    """download seq and gff data"""
    remote_template = f"{config.remote_path}/release-{config.release}/" + "{}"
    site_map = elt_site_map.get_site_map(config.host)
    if verbose:
        click.secho(f"DOWNLOADING\n  ensembl release={config.release}", fg="green")
        click.secho("\n".join(f"  {d}" for d in config.species_dbs), fg="green")
        click.secho(f"\nWRITING to output path={config.staging_genomes}\n", fg="green")

    msg = "Downloading genomes"
    if progress is not None:
        species_download = progress.add_task(
            total=len(config.species_dbs),
            description=msg,
        )

    patterns = dict(fasta=valid_seq_file, gff3=valid_gff3_file(config.release))
    for key in config.species_dbs:
        db_prefix = elt_species.Species.get_ensembl_db_prefix(key)
        local_root = config.staging_genomes / db_prefix
        local_root.mkdir(parents=True, exist_ok=True)
        for subdir in ("fasta", "gff3"):
            if subdir == "fasta":
                remote = site_map.get_seqs_path(db_prefix)
            else:
                remote = site_map.get_annotations_path(db_prefix)

            remote_dir = remote_template.format(remote)
            remote_paths = list(
                elt_ftp.listdir(config.host, path=remote_dir, pattern=patterns[subdir]),
            )
            if verbose:
                print(f"{remote_paths=}")
            if debug:
                # we need the checksum files
                paths = [p for p in remote_paths if elt_util.is_signature(p)]
                # but fewer data files, to reduce time for debugging
                remote_paths = [
                    p for p in remote_paths if not elt_util.dont_checksum(p)
                ]
                remote_paths = remote_paths[:4] + paths

            dest_path = config.staging_genomes / db_prefix / subdir
            dest_path.mkdir(parents=True, exist_ok=True)
            _remove_tmpdirs(dest_path)
            icon = "🧬🧬" if subdir == "fasta" else "📚"
            elt_ftp.download_data(
                host=config.host,
                local_dest=dest_path,
                remote_paths=remote_paths,
                description=f"{db_prefix[:10]}... {icon}",
                do_checksum=True,
                progress=progress,
            )

        if progress is not None:
            progress.update(species_download, description=msg, advance=1)


class valid_compara_align:
    """whole genome alignment data"""

    def __init__(self) -> None:
        self._valid = re.compile("([.](emf|maf)[.]gz|README|MD5SUM)")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


def download_aligns(
    config: elt_config.Config,
    debug: bool,
    verbose: bool,
    progress: Progress | None = None,
):
    """download whole genome alignments"""
    if not config.align_names:
        return

    site_map = elt_site_map.get_site_map(config.host)
    remote_template = (
        f"{config.remote_path}/release-{config.release}/{site_map.alignments_path}/"
        + "{}"
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
        remote_paths = list(elt_ftp.listdir(config.host, remote_path, valid_compara))
        if verbose:
            print(remote_paths)

        if debug:
            # we need the checksum files
            paths = [p for p in remote_paths if elt_util.is_signature(p)]
            remote_paths = [p for p in remote_paths if not elt_util.is_signature(p)]
            remote_paths = remote_paths[:4] + paths

        local_dir = config.staging_aligns / align_name
        local_dir.mkdir(parents=True, exist_ok=True)
        _remove_tmpdirs(local_dir)
        elt_ftp.download_data(
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


class valid_compara_homology:
    """homology tsv files"""

    def __init__(self) -> None:
        self._valid = re.compile("([.]tsv[.]gz|README|MD5SUM)")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


def download_homology(
    config: elt_config.Config,
    debug: bool,
    verbose: bool,
    progress: Progress | None = None,
):
    """downloads tsv homology files for each genome"""
    if not config.homologies:
        return

    site_map = elt_site_map.get_site_map(config.host)
    remote_template = (
        f"{config.remote_path}/release-{config.release}/{site_map.homologies_path}/"
        + "{}"
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
            elt_ftp.listdir(config.host, remote_path, valid_compara_homology()),
        )
        if verbose:
            print(remote_paths)

        if debug:
            # we need the checksum files
            remote_paths = [p for p in remote_paths if not elt_util.is_signature(p)]
            remote_paths = remote_paths[:4]

        local_dir = local / db_name
        local_dir.mkdir(parents=True, exist_ok=True)
        _remove_tmpdirs(local_dir)
        elt_ftp.download_data(
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


def download_ensembl_tree(host: str, remote_path: str, release: str, tree_fname: str):
    """loads a tree from Ensembl"""
    site_map = elt_site_map.get_site_map(host)
    url = f"https://{host}/{remote_path}/release-{release}/{site_map.trees_path}/{tree_fname}"
    return load_tree(url)


def get_ensembl_trees(host: str, remote_path: str, release: str) -> list[str]:
    """returns trees from ensembl compara"""
    site_map = elt_site_map.get_site_map(host)
    path = f"{remote_path}/release-{release}/{site_map.trees_path}"
    return list(
        elt_ftp.listdir(host=host, path=path, pattern=lambda x: x.endswith(".nh")),
    )


def get_species_for_alignments(
    host: str,
    remote_path: str,
    release: str,
    align_names: list[str],
) -> dict[str, list[str]]:
    """return the species for the indicated alignments"""
    ensembl_trees = get_ensembl_trees(
        host=host,
        remote_path=remote_path,
        release=release,
    )
    aligns_trees = elt_util.trees_for_aligns(align_names, ensembl_trees)
    species = {}
    for tree_path in aligns_trees.values():
        tree_path = pathlib.Path(tree_path)
        tree = download_ensembl_tree(
            host=host,
            remote_path=remote_path,
            release=release,
            tree_fname=tree_path.name,
        )
        # dict structure is {common name: db prefix}, just use common name
        species |= {
            n: ["core"] for n in elt_species.species_from_ensembl_tree(tree).keys()
        }
    return species
