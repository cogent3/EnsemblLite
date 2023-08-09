from __future__ import annotations

import os
import re
import shutil

import click

from ensembl_cli.ftp_download import download_data, listdir
from ensembl_cli.species import Species
from ensembl_cli.util import (
    Config,
    dont_checksum,
    get_resource_path,
    is_signature,
)


_cfg = get_resource_path("ensembldb_download.cfg")

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


def _remove_tmpdirs(path: os.PathLike):
    """delete any tmp dirs left over from unsuccessful runs"""
    tmpdirs = [p for p in path.glob("tmp*") if p.is_dir()]
    for tmpdir in tmpdirs:
        shutil.rmtree(tmpdir)


def download_species(config: Config, debug: bool, verbose: bool):
    """download seq and gff data"""
    remote_template = f"{config.remote_path}/release-{config.release}/" + "{}/{}"

    if verbose:
        click.secho(f"DOWNLOADING\n  ensembl release={config.release}", fg="green")
        click.secho("\n".join(f"  {d}" for d in config.species_dbs), fg="green")
        click.secho(f"\nWRITING to output path={config.staging_path}\n", fg="green")

    patterns = dict(fasta=valid_seq_file, gff3=valid_gff3_file(config.release))
    for key in config.species_dbs:
        db_prefix = Species.get_ensembl_db_prefix(key)
        local_root = config.staging_path / db_prefix
        local_root.mkdir(parents=True, exist_ok=True)
        for subdir in ("fasta", "gff3"):
            path = remote_template.format(subdir, db_prefix)
            path = f"{path}/dna" if subdir == "fasta" else path
            dest_path = config.staging_path / db_prefix / subdir
            dest_path.mkdir(parents=True, exist_ok=True)
            remote_paths = list(
                listdir(config.host, path=path, pattern=patterns[subdir])
            )
            if verbose:
                print(f"{remote_paths=}")
            if debug:
                # we need the checksum files
                paths = [p for p in remote_paths if is_signature(p)]
                # but fewer data files, to reduce time for debugging
                remote_paths = [p for p in remote_paths if not dont_checksum(p)]
                remote_paths = remote_paths[:4] + paths

            _remove_tmpdirs(dest_path)
            download_data(
                host=config.host,
                local_dest=dest_path,
                remote_paths=remote_paths,
                description=f"{db_prefix[:5]}.../{subdir}",
                do_checksum=True,
            )

    return


class valid_compara_align:
    """whole genome alignment data"""

    def __init__(self) -> None:
        self._valid = re.compile("([.]emf[.]gz|README|MD5SUM)")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


def download_aligns(config: Config, debug: bool, verbose: bool):
    """download whole genome alignments"""
    remote_template = (
        f"{config.remote_path}/release-{config.release}/emf/ensembl-compara/multiple_alignments/"
        + "{}"
    )
    valid_compara = valid_compara_align()
    local = config.staging_path / "compara"
    for align_name in config.align_names:
        remote_path = remote_template.format(align_name)
        remote_paths = list(listdir(config.host, remote_path, valid_compara))
        if verbose:
            print(remote_paths)

        if debug:
            # we need the checksum files
            paths = [p for p in remote_paths if is_signature(p)]
            remote_paths = [p for p in remote_paths if not is_signature(p)]
            remote_paths = remote_paths[:4] + paths

        local_dir = local / align_name
        local_dir.mkdir(parents=True, exist_ok=True)
        _remove_tmpdirs(local_dir)
        download_data(
            host=config.host,
            local_dest=local_dir,
            remote_paths=remote_paths,
            description=f"compara/{align_name[:5]}...",
            do_checksum=True,
        )

    return
