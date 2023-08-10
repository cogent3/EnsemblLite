import os
import shutil
import typing

from cogent3 import load_annotations, load_seq, make_seq, open_
from cogent3.util import parallel as PAR
from rich.progress import track
from unsync import unsync

from ensembl_cli import maf
from ensembl_cli.aligndb import AlignDb
from ensembl_cli.convert import seq_to_gap_coords
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

    return


def seq2gaps(record: dict):
    seq = make_seq(record.pop("seq"))
    record["gap_spans"] = seq_to_gap_coords(seq)
    return record


def _load_one_align(path: os.PathLike) -> typing.Iterable[dict]:
    records = []
    for block_id, align in enumerate(maf.parse(path)):
        converted = []
        for maf_name, seq in align.items():
            record = maf_name.to_dict()
            record["block_id"] = block_id
            record["source"] = path.name
            record["seq"] = seq
            converted.append(seq2gaps(record))
        records.extend(converted)
    return records


def local_install_compara(config: Config, force_overwrite: bool):
    if force_overwrite:
        shutil.rmtree(config.install_path, ignore_errors=True)

    for align_name in config.align_names:
        src_dir = config.staging_path / "compara" / align_name
        dest_dir = config.install_path / "compara"
        dest_dir.mkdir(parents=True, exist_ok=True)
        # write out to a db with align_name
        db = AlignDb(source=(dest_dir / f"{align_name}.sqlitedb"))
        records = []
        paths = list(src_dir.glob(f"{align_name}*maf*"))
        max_workers = min(len(paths), 10)
        for result in track(
            PAR.as_completed(_load_one_align, paths, max_workers=max_workers),
            transient=True,
            description="Installing compara...",
            total=len(paths),
        ):
            records.extend(result)

        db.add_records(records=records)
        db.close()

    return
