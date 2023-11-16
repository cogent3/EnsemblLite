import os
import shutil
import typing

from cogent3 import load_annotations, load_seq, make_seq, open_
from cogent3.parse.table import FilteringParser
from cogent3.util import parallel as PAR
from rich.progress import track
from unsync import unsync

from ensembl_lite import maf
from ensembl_lite._aligndb import AlignDb
from ensembl_lite._config import _COMPARA_NAME, Config
from ensembl_lite._genomedb import (
    _ANNOTDB_NAME,
    _SEQDB_NAME,
    CompressedGenomeSeqsDb,
)
from ensembl_lite._homologydb import HomologyDb
from ensembl_lite.convert import seq_to_gap_coords
from ensembl_lite.util import elt_compress_it


@unsync(cpu_bound=True)
def _install_one_seq(src: os.PathLike) -> typing.Tuple[str, bytes]:
    seq = load_seq(src, moltype="dna", label_to_name=lambda x: x.split()[0])
    return seq.name, elt_compress_it(str(seq))


@unsync(cpu_bound=True)
def _install_one_annotations(src: os.PathLike, dest: os.PathLike) -> bool:
    if dest.exists():
        return True

    _ = load_annotations(path=src, write_path=dest)
    return True


def _install_gffdb(src_dir: os.PathLike, dest_dir: os.PathLike) -> list[bool]:
    src_dir = src_dir / "gff3"
    dest = dest_dir / _ANNOTDB_NAME
    paths = list(src_dir.glob("*.gff3.gz"))
    return [_install_one_annotations(path, dest) for path in paths]


T = typing.Tuple[os.PathLike, typing.List[typing.Tuple[str, bytes]]]


def _install_seqs(src_dir: os.PathLike, dest_dir: os.PathLike) -> T:
    src_dir = src_dir / "fasta"
    paths = list(src_dir.glob("*.fa.gz"))
    dest = dest_dir / _SEQDB_NAME
    return dest, [_install_one_seq(path) for path in paths]


def local_install_genomes(config: Config, force_overwrite: bool):
    if force_overwrite:
        shutil.rmtree(config.install_genomes, ignore_errors=True)

    # we create the local installation
    config.install_genomes.mkdir(parents=True, exist_ok=True)
    # we create subdirectories for each species
    for db_name in config.db_names:
        sp_dir = config.install_genomes / db_name
        sp_dir.mkdir(parents=True, exist_ok=True)

    # for each species, we identify the download and dest paths for annotations
    # our tasks here are the load/compress steps
    tasks = {}
    for db_name in config.db_names:
        src_dir = config.staging_genomes / db_name
        dest_dir = config.install_genomes / db_name
        dest, tsks = _install_seqs(src_dir, dest_dir)
        tasks[dest] = tsks

    for dest, tsks in tasks.items():
        db = CompressedGenomeSeqsDb(source=dest, species=dest.parent.name)
        records = [
            tsk.result()
            for tsk in track(
                tsks,
                description=f"Installing {dest.parent.name} seqs...",
                transient=True,
            )
        ]
        db.add_compressed_records(records=records)

    # we now load the individual gff3 files and write to annotation db's
    tasks = []
    for db_name in config.db_names:
        src_dir = config.staging_genomes / db_name
        dest_dir = config.install_genomes / db_name
        tasks.extend(_install_gffdb(src_dir, dest_dir))

    # we do all tasks in one go
    _ = [
        t.result()
        for t in track(tasks, description="Installing annotations...", transient=True)
    ]

    db.close()

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
        shutil.rmtree(config.install_path / _COMPARA_NAME, ignore_errors=True)

    for align_name in config.align_names:
        src_dir = config.staging_aligns / align_name
        dest_dir = config.install_aligns
        dest_dir.mkdir(parents=True, exist_ok=True)
        # write out to a db with align_name
        db = AlignDb(source=(dest_dir / f"{align_name}.sqlitedb"))
        records = []
        paths = list(src_dir.glob(f"{align_name}*maf*"))
        max_workers = min(len(paths), 10)
        for result in track(
            PAR.as_completed(_load_one_align, paths, max_workers=max_workers),
            transient=True,
            description="Installing aligns...",
            total=len(paths),
        ):
            records.extend(result)

        db.add_records(records=records)
        db.close()

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

    def __call__(self, paths: typing.Iterable[os.PathLike]) -> list:
        final = []
        for path in paths:
            with open_(path) as infile:
                # we bulk load because it's faster than the default line-by-line
                # iteration on a file
                data = infile.read().splitlines()

            rows = list(self._reader(data))
            header = rows.pop(0)
            assert list(header) == list(self.src_cols), (header, self.src_cols)
            rows = [r + [path.name] for r in rows]
            final.extend(rows)

        return final


def local_install_homology(config: Config, force_overwrite: bool):
    if force_overwrite:
        shutil.rmtree(config.install_homologies, ignore_errors=True)

    config.install_homologies.mkdir(parents=True, exist_ok=True)

    outpath = config.install_homologies / "homologies.sqlitedb"
    db = HomologyDb(source=outpath)

    dirnames = [config.staging_homologies / sp for sp in config.db_names]
    loader = LoadHomologies(allowed_species=set(config.db_names))
    # On test cases, only 30% speedup from running in parallel due to overhead
    # of pickling the data, but considerable increase in memory. So, run
    # in serial to avoid memory issues since it's reasonably fast anyway.
    for dirname in track(
        dirnames,
        transient=True,
        description="Installing homologies...",
    ):
        rows = loader(dirname.glob("*.tsv.gz"))
        db.add_records(records=rows, col_order=loader.dest_col)
        del rows

    no_records = len(db) == 0
    db.close()
    if no_records:
        outpath.unlink()
