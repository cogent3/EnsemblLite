from __future__ import annotations

import shutil
import typing

import numpy

from cogent3 import load_annotations, make_seq, open_
from cogent3.app.composable import LOADER, define_app
from cogent3.app.typing import IdentifierType
from cogent3.parse.fasta import MinimalFastaParser
from cogent3.parse.table import FilteringParser
from cogent3.util.io import iter_splitlines
from rich.progress import Progress, track

from ensembl_lite import _maf
from ensembl_lite._aligndb import AlignDb, AlignRecord
from ensembl_lite._config import _COMPARA_NAME, Config
from ensembl_lite._genomedb import _ANNOTDB_NAME, _SEQDB_NAME, SeqsDataHdf5
from ensembl_lite._homologydb import HomologyDb
from ensembl_lite._species import Species
from ensembl_lite._util import PathType, get_iterable_tasks


def _rename(label: str) -> str:
    return label.split()[0]


@define_app
class fasta_to_hdf5:
    def __init__(self, config: Config):
        self.config = config

    def main(self, db_name: str) -> bool:
        src_dir = self.config.staging_genomes / db_name
        dest_dir = self.config.install_genomes / db_name

        seq_store = SeqsDataHdf5(
            source=dest_dir / _SEQDB_NAME,
            species=Species.get_species_name(db_name),
            mode="w",
        )

        src_dir = src_dir / "fasta"
        for path in src_dir.glob("*.fa.gz"):
            for label, seq in MinimalFastaParser(iter_splitlines(path)):
                seqid = _rename(label)
                seq_store.add_record(seqid=seqid, seq=seq)
                del seq

        seq_store.close()

        return True


def _get_seqs(src: PathType) -> list[tuple[str, str]]:
    with open_(src) as infile:
        data = infile.read().splitlines()
    name_seqs = list(MinimalFastaParser(data))
    return [(_rename(name), seq) for name, seq in name_seqs]


def _load_one_annotations(src_dest: tuple[PathType, PathType]) -> bool:
    src, dest = src_dest
    if dest.exists():
        return True

    db = load_annotations(path=src, write_path=dest)
    db.db.close()
    del db
    return True


def _make_src_dest_annotation_paths(
    src_dir: PathType, dest_dir: PathType
) -> list[tuple[PathType, PathType]]:
    src_dir = src_dir / "gff3"
    dest = dest_dir / _ANNOTDB_NAME
    paths = list(src_dir.glob("*.gff3.gz"))
    return [(path, dest) for path in paths]


T = tuple[PathType, list[tuple[str, str]]]


def _prepped_seqs(
    src_dir: PathType, dest_dir: PathType, progress: Progress, max_workers: int
) -> T:
    src_dir = src_dir / "fasta"
    paths = list(src_dir.glob("*.fa.gz"))
    dest = dest_dir / _SEQDB_NAME
    all_seqs = []

    common_name = Species.get_common_name(src_dir.parent.name)
    msg = f"📚🗜️ {common_name} seqs"
    load = progress.add_task(msg, total=len(paths))
    tasks = get_iterable_tasks(func=_get_seqs, series=paths, max_workers=max_workers)
    for result in tasks:
        all_seqs.extend(result)
        progress.update(load, advance=1, description=msg)

    progress.update(load, visible=False)
    return dest, all_seqs


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

    with Progress(transient=True) as progress:
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


def seq2gaps(record: dict) -> AlignRecord:
    seq = make_seq(record.pop("seq"))
    indel_map, _ = seq.parse_out_gaps()
    if indel_map.num_gaps:
        record["gap_spans"] = numpy.array(
            [indel_map.gap_pos, indel_map.get_gap_lengths()], dtype=numpy.int32
        ).T
    else:
        record["gap_spans"] = numpy.array([], dtype=numpy.int32)
    return AlignRecord(**record)


@define_app(app_type=LOADER)
class _load_one_align:
    def __init__(self, species: set[str] | None = None):
        self.species = species or {}

    def main(self, path: IdentifierType) -> typing.Iterable[dict]:
        records = []
        for block_id, align in enumerate(_maf.parse(path)):
            converted = []
            for maf_name, seq in align.items():
                if maf_name.species not in self.species:
                    continue
                record = maf_name.to_dict()
                record["block_id"] = f"{path.name}-{block_id}"
                record["source"] = path.name
                record["seq"] = seq
                converted.append(seq2gaps(record))
            records.extend(converted)
        return records


def local_install_compara(
    config: Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
):
    if force_overwrite:
        shutil.rmtree(config.install_path / _COMPARA_NAME, ignore_errors=True)

    aln_loader = _load_one_align(set(config.db_names))

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
            transient=True,
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


class LoadHomologies:
    def __init__(self, allowed_species: set):
        self._allowed_species = allowed_species
        # map the Ensembl columns to HomologyDb columns

        self.src_cols = (
            "homology_type",
            "species",
            "gene_stable_id",
            "protein_stable_id",
            "homology_species",
            "homology_gene_stable_id",
            "homology_protein_stable_id",
        )
        self.dest_col = (
            "relationship",
            "species_1",
            "gene_id_1",
            "prot_id_1",
            "species_2",
            "gene_id_2",
            "prot_id_2",
            "source",
        )
        self._reader = FilteringParser(
            row_condition=self._matching_species, columns=self.src_cols, sep="\t"
        )

    def _matching_species(self, row):
        return {row[1], row[4]} <= self._allowed_species

    def __call__(self, paths: typing.Iterable[PathType]) -> list:
        final = []
        for path in paths:
            rows = list(self._reader(iter_splitlines(path)))
            header = rows.pop(0)
            assert list(header) == list(self.src_cols), (header, self.src_cols)
            rows = [r + [path.name] for r in rows]
            final.extend(rows)

        return final


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
        dirnames.append(list(path.glob("*.tsv.gz")))

    loader = LoadHomologies(allowed_species=set(config.db_names))
    if max_workers:
        max_workers = min(len(dirnames) + 1, max_workers)

    if verbose:
        print(f"homologies {max_workers=}")

    with Progress(transient=True) as progress:
        msg = "Installing homologies"
        writing = progress.add_task(total=len(dirnames), description=msg, advance=0)
        tasks = get_iterable_tasks(
            func=loader, series=dirnames, max_workers=max_workers
        )
        for rows in tasks:
            db.add_records(records=rows, col_order=loader.dest_col)
            del rows
            progress.update(writing, description=msg, advance=1)

    no_records = len(db) == 0
    db.close()
    if no_records:
        outpath.unlink()

    if verbose:
        print("Finished installing homologies")
