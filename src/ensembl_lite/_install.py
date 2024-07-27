from __future__ import annotations

import shutil
import typing

from rich.progress import Progress

from ensembl_lite._align import AlignDb
from ensembl_lite._config import Config
from ensembl_lite._genome import (
    _ANNOTDB_NAME,
    fasta_to_hdf5,
    make_annotation_db,
)
from ensembl_lite._homology import (
    _HOMOLOGYDB_NAME,
    HomologyDb,
    compressor,
    inflate,
    load_homologies,
    pickler,
)
from ensembl_lite._maf import load_align_records
from ensembl_lite._species import SPECIES_NAME, Species
from ensembl_lite._util import PathType, get_iterable_tasks


def _make_src_dest_annotation_paths(
    src_dir: PathType, dest_dir: PathType
) -> list[tuple[PathType, PathType]]:
    src_dir = src_dir / "gff3"
    dest = dest_dir / _ANNOTDB_NAME
    paths = list(src_dir.glob("*.gff3.gz"))
    return [(path, dest) for path in paths]


def local_install_genomes(
    config: Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
    progress: typing.Optional[Progress] = None,
):
    if force_overwrite:
        shutil.rmtree(config.install_genomes, ignore_errors=True)
    # we create the local installation
    config.install_genomes.mkdir(parents=True, exist_ok=True)
    # we create subdirectories for each species
    for db_name in list(config.db_names):
        sp_dir = config.install_genomes / db_name
        sp_dir.mkdir(parents=True, exist_ok=True)

    # for each species, we identify the download and dest paths for annotations
    db_names = list(config.db_names)
    if max_workers:
        max_workers = min(len(db_names) + 1, max_workers)

    if verbose:
        print(f"genomes {max_workers=}")

    # we load the individual gff3 files and write to annotation db's
    src_dest_paths = []
    for db_name in config.db_names:
        src_dir = config.staging_genomes / db_name
        dest_dir = config.install_genomes / db_name
        src_dest_paths.extend(_make_src_dest_annotation_paths(src_dir, dest_dir))

    msg = "Installing features 📚"
    if progress is not None:
        writing = progress.add_task(total=len(src_dest_paths), description=msg)

    tasks = get_iterable_tasks(
        func=make_annotation_db, series=src_dest_paths, max_workers=max_workers
    )
    for db_name, prefixes in tasks:
        if verbose:
            print(f"{db_name=} {prefixes=}")

        if prefixes:
            for prefix in prefixes:
                Species.add_stableid_prefix(db_name, prefix)

        if progress is not None:
            progress.update(writing, description=msg, advance=1)

    species_table = Species.to_table()
    species_table.write(config.install_genomes / SPECIES_NAME)
    if verbose:
        print("Finished installing features ")

    msg = "Installing  🧬🧬"
    if progress is not None:
        writing = progress.add_task(total=len(db_names), description=msg, advance=0)
    # we parallelise across databases
    writer = fasta_to_hdf5(config=config)
    tasks = get_iterable_tasks(func=writer, series=db_names, max_workers=max_workers)
    for result in tasks:
        if not result:
            print(result)
            raise RuntimeError

        if progress is not None:
            progress.update(writing, description=msg, advance=1)

    if verbose:
        print("Finished installing sequences ")
    return


def local_install_alignments(
    config: Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
    progress: typing.Optional[Progress] = None,
):
    if force_overwrite:
        shutil.rmtree(config.install_aligns, ignore_errors=True)

    aln_loader = load_align_records(set(config.db_names))

    for align_name in config.align_names:
        src_dir = config.staging_aligns / align_name
        dest_dir = config.install_aligns
        dest_dir.mkdir(parents=True, exist_ok=True)
        # write out to a db with align_name
        db = AlignDb(source=(dest_dir / f"{align_name}.sqlitedb"))
        records = []
        paths = list(src_dir.glob(f"{align_name}*maf*"))

        if max_workers and max_workers > 1:
            # we adjust the maximum workers to the number of paths
            max_workers = min(len(paths) + 1, max_workers or 0)

        if verbose:
            print(f"{max_workers=}")

        series = get_iterable_tasks(
            func=aln_loader, series=paths, max_workers=max_workers
        )

        msg = "Installing alignments"
        if progress is not None:
            writing = progress.add_task(total=len(paths), description=msg, advance=0)

        for result in series:
            if not result:
                print(result)
                raise RuntimeError

            records.extend(result)

            if progress is not None:
                progress.update(writing, description=msg, advance=1)

        db.add_records(records=records)
        db.close()

    if verbose:
        print("Finished installing alignments")

    return


def local_install_homology(
    config: Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
    progress: typing.Optional[Progress] = None,
):
    if force_overwrite:
        shutil.rmtree(config.install_homologies, ignore_errors=True)

    config.install_homologies.mkdir(parents=True, exist_ok=True)

    outpath = config.install_homologies / _HOMOLOGYDB_NAME
    db = HomologyDb(source=outpath)

    dirnames = []
    for sp in config.db_names:
        path = config.staging_homologies / sp
        dirnames.extend(list(path.glob("*.tsv.gz")))

    if max_workers:
        max_workers = min(len(dirnames) + 1, max_workers)
    else:
        max_workers = 1

    if verbose:
        print(f"homologies {max_workers=}")

    loader = load_homologies(allowed_species=set(config.db_names))
    if max_workers > 1:
        loader = loader + pickler + compressor

    msg = "Installing homologies"
    if progress is not None:
        writing = progress.add_task(total=len(dirnames), description=msg, advance=0)

    tasks = get_iterable_tasks(func=loader, series=dirnames, max_workers=max_workers)
    for result in tasks:
        if max_workers > 1:
            # reconstitute the blosc compressed data
            result = inflate(result)

        for rel_type, records in result.items():
            db.add_records(records=records, relationship_type=rel_type)

        if progress is not None:
            progress.update(writing, description=msg, advance=1)

    no_records = len(db) == 0
    db.close()
    if no_records:
        outpath.unlink()

    if verbose:
        print("Finished installing homologies")
