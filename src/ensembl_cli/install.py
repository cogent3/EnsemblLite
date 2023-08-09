import os
import shutil

from cogent3 import load_annotations, load_seq, open_
from rich.progress import track
from unsync import unsync

from ensembl_cli.util import Config


@unsync(cpu_bound=True)
def _install_one_seq(src, dest_dir):
    seq = load_seq(src, moltype="dna", label_to_name=lambda x: x.split()[0])
    with open_(dest_dir / f"{seq.name}.fa.gz", mode="wt") as outfile:
        outfile.write(seq.to_fasta(block_size=int(1e9)))
    return True


@unsync(cpu_bound=True)
def _install_one_annotations(src, dest):
    if dest.exists():
        return True

    _ = load_annotations(path=src, write_path=dest)
    return True


def _install_gffdb(src_dir: os.PathLike, dest_dir: os.PathLike):
    src_dir = src_dir / "gff3"
    dest = dest_dir / "features.gff3db"
    paths = list(src_dir.glob("*.gff3.gz"))
    return [_install_one_annotations(path, dest) for path in paths]


def _install_seqs(src_dir: os.PathLike, dest_dir: os.PathLike):
    src_dir = src_dir / "fasta"
    paths = list(src_dir.glob("*.fa.gz"))
    return [_install_one_seq(path, dest_dir) for path in paths]


def local_install_genomes(config: Config, force_overwrite: bool):
    if force_overwrite:
        shutil.rmtree(config.install_path, ignore_errors=True)

    # we create the local installation
    config.install_path.mkdir(parents=True, exist_ok=True)
    # we create subdirectories for each species
    for db_name in config.db_names:
        sp_dir = config.install_path / db_name
        sp_dir.mkdir(parents=True, exist_ok=True)

    # for each species, we identify the download and dest paths for annotations
    tasks = []
    for db_name in config.db_names:
        src_dir = config.staging_path / db_name
        dest_dir = config.install_path / db_name
        tasks.extend(_install_seqs(src_dir, dest_dir))

    # we now load the individual gff3 files and write to annotation db's
    for db_name in config.db_names:
        src_dir = config.staging_path / db_name
        dest_dir = config.install_path / db_name
        tasks.extend(_install_gffdb(src_dir, dest_dir))
    # we do all tasks in one go
    _ = [
        t.result()
        for t in track(tasks, description="Installing genomes...", transient=True)
    ]


def local_install_compara(config: Config, force_overwrite: bool) -> Config:
    if force_overwrite:
        shutil.rmtree(config.install_path, ignore_errors=True)

    # we create the local installation
    config.install_path.mkdir(parents=True, exist_ok=True)

    # create the tasks for each alignment name
    return config
