from __future__ import annotations

import shutil

from rich.progress import Progress, track

from ensembl_lite._aligndb import AlignDb
from ensembl_lite._config import Config
from ensembl_lite._genomedb import (
    _ANNOTDB_NAME,
    _load_one_annotations,
    fasta_to_hdf5,
)
from ensembl_lite._homologydb import (
    HomologyDb,
    compressor,
    inflate,
    load_homologies,
    pickler,
)
from ensembl_lite._maf import load_align_records
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

    with Progress(transient=False) as progress:
        msg = "Installing  🧬 features"
        writing = progress.add_task(total=len(src_dest_paths), description=msg)
        tasks = get_iterable_tasks(
            func=_load_one_annotations, series=src_dest_paths, max_workers=max_workers
        )
        for _ in tasks:
            progress.update(writing, description=msg, advance=1)

    if verbose:
        print("Finished installing features ")

    with Progress(transient=True) as progress:
        writing = progress.add_task(
            total=len(db_names), description="Installing  🧬", advance=0
        )
        # we parallelise across databases
        writer = fasta_to_hdf5(config=config)
        tasks = get_iterable_tasks(
            func=writer, series=db_names, max_workers=max_workers
        )
        for result in tasks:
            if not result:
                print(result)
                raise RuntimeError
            progress.update(writing, description="Installing  🧬", advance=1)

    if verbose:
        print("Finished installing sequences ")
    return


def local_install_alignments(
    config: Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
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

        for result in track(
            series,
            transient=False,
            description="Installing alignments",
            total=len(paths),
        ):
            if not result:
                print(result)
                raise RuntimeError
            records.extend(result)

        db.add_records(records=records)
        db.close()

    if verbose:
        print("Finished installing homologies")

    return


def local_install_homology(
    config: Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
):
    if force_overwrite:
        shutil.rmtree(config.install_homologies, ignore_errors=True)

    config.install_homologies.mkdir(parents=True, exist_ok=True)

    outpath = config.install_homologies / "homologies.sqlitedb"
    db = HomologyDb(source=outpath)

    dirnames = []
    for sp in config.db_names:
        path = config.staging_homologies / sp
        dirnames.extend(list(path.glob("*.tsv.gz")))

    if max_workers:
        max_workers = min(len(dirnames) + 1, max_workers)

    if verbose:
        print(f"homologies {max_workers=}")

    loader = load_homologies(allowed_species=set(config.db_names))
    if max_workers > 1:
        loader = loader + pickler + compressor

    with Progress(transient=False) as progress:
        msg = "Installing homologies"
        writing = progress.add_task(total=len(dirnames), description=msg, advance=0)
        tasks = get_iterable_tasks(
            func=loader, series=dirnames, max_workers=max_workers
        )
        for result in tasks:
            if max_workers > 1:
                # reconstitute the blosc compressed data
                result = inflate(result)

            for rel_type, records in result.items():
                db.add_records(records=records, relationship_type=rel_type)
            progress.update(writing, description=msg, advance=1)

    no_records = len(db) == 0
    db.close()
    if no_records:
        outpath.unlink()

    if verbose:
        print("Finished installing homologies")
