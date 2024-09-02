import pathlib
import typing
from collections import defaultdict
from dataclasses import dataclass

import h5py
import numpy
from cogent3.app.composable import define_app
from cogent3.core.alignment import Aligned, Alignment
from cogent3.core.location import _DEFAULT_GAP_DTYPE, IndelMap

from ensembl_tui import _genome as elt_genome
from ensembl_tui import _storage_mixin as elt_mixin
from ensembl_tui import _util as elt_util

_no_gaps = numpy.array([], dtype=_DEFAULT_GAP_DTYPE)

GAP_STORE_SUFFIX = "indels-hdf5_blosc2"
ALIGN_STORE_SUFFIX = "align_coords-sqlitedb"


@dataclass(slots=True)
class AlignRecord:
    """a record from an AlignDb

    Notes
    -----
    Can return fields as attributes or like a dict using the field name as
    a string.
    """

    source: str
    block_id: int
    species: str
    seqid: str
    start: int
    stop: int
    strand: str
    gap_spans: numpy.ndarray

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, item, value):
        setattr(self, item, value)

    def __eq__(self, other):
        attrs = "block_id", "species", "seqid", "start", "stop", "strand"
        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return (self.gap_spans == other.gap_spans).all()

    def __hash__(self):
        return hash(
            (
                self.block_id,
                self.species,
                self.seqid,
                self.start,
                self.stop,
                self.strand,
            ),
        )

    @property
    def gap_data(self):
        if len(self.gap_spans):
            gap_pos, gap_lengths = self.gap_spans.T
        else:
            gap_pos, gap_lengths = _no_gaps.copy(), _no_gaps.copy()

        return gap_pos, gap_lengths


ReturnType = tuple[str, tuple]  # the sql statement and corresponding values


class GapStore(elt_mixin.Hdf5Mixin):
    # store gap data from aligned sequences
    def __init__(
        self,
        source: elt_util.PathType,
        align_name: str | None = None,
        mode: str = "r",
        in_memory: bool = False,
    ):
        self.source = pathlib.Path(source)
        self.mode = "w-" if mode == "w" else mode
        h5_kwargs = (
            dict(
                driver="core",
                backing_store=False,
            )
            if in_memory
            else {}
        )
        try:
            self._file = h5py.File(source, mode=self.mode, **h5_kwargs)
        except OSError:
            print(f"{source=}")
            raise

        if "r" not in self.mode and "align_name" not in self._file.attrs:
            assert align_name
            self._file.attrs["align_name"] = align_name
        if (
            align_name
            and (file_species := self._file.attrs.get("align_name", None)) != align_name
        ):
            raise ValueError(f"{self.source.name!r} {file_species!r} != {align_name}")
        self.align_name = self._file.attrs["align_name"]

    def add_record(self, *, index: int, gaps: numpy.ndarray):
        # dataset names must be strings
        index = str(index)
        if index in self._file:
            stored = self._file[index]
            if (gaps == stored).all():
                # already seen this index
                return
            # but it's different, which is a problem
            raise ValueError(f"{index!r} already present but with different gaps")
        self._file.create_dataset(
            name=index,
            data=gaps,
            chunks=True,
            **elt_util._HDF5_BLOSC2_KWARGS,
        )
        self._file.flush()

    def get_record(self, *, index: int) -> numpy.ndarray:
        return self._file[str(index)][:]


# TODO add a table and methods to support storing the species tree used
#  for the alignment and for getting the species tree
class AlignDb(elt_mixin.SqliteDbMixin):
    table_name = "align"
    _align_schema = {
        "id": "INTEGER PRIMARY KEY",  # used to uniquely identify gap_spans in bound GapStore
        "source": "TEXT",  # the file path
        "block_id": "INTEGER",  # the tree id from MAF
        "species": "TEXT",
        "seqid": "TEXT",
        "start": "INTEGER",
        "stop": "INTEGER",
        "strand": "TEXT",
    }

    _index_columns = {"align": ("id", "block_id", "seqid", "start", "stop")}

    def __init__(self, *, source=":memory:", mode="a"):
        """
        Parameters
        ----------
        source
            location to store the db, defaults to in memory only
        """
        # note that data is destroyed
        source = pathlib.Path(source)
        self.source = source
        if source.name == ":memory:":
            gap_path = "memory"
            kwargs = dict(in_memory=True)
        else:
            gap_path = source.parent / f"{source.stem}.{GAP_STORE_SUFFIX}"
            kwargs = dict(in_memory=False)

        self.gap_store = GapStore(
            source=gap_path,
            align_name=source.stem,
            mode=mode,
            **kwargs,
        )
        self._db = None
        self._init_tables()

    def add_records(self, records: typing.Sequence[AlignRecord]):
        # bulk insert
        col_order = [
            row[1]
            for row in self.db.execute(
                f"PRAGMA table_info({self.table_name})",
            ).fetchall()
            if row[1] != "id"
        ]

        # we need to identify block_id's that have already been used
        block_ids = tuple({r.block_id for r in records})
        val_placeholder = ", ".join("?" * len(block_ids))
        sql = f"SELECT DISTINCT(block_id) from {self.table_name} WHERE block_id IN ({val_placeholder})"
        used = {r[0] for r in self.db.execute(sql, block_ids).fetchall()}

        val_placeholder = ", ".join("?" * len(col_order))
        sql = f"INSERT INTO {self.table_name} ({', '.join(col_order)}) VALUES ({val_placeholder}) RETURNING id"

        for i in range(len(records)):
            if records[i].block_id in used:
                continue

            index = self.db.execute(sql, [records[i][c] for c in col_order]).fetchone()
            index = index["id"]
            self.gap_store.add_record(index=index, gaps=records[i].gap_spans)

    def _get_block_id(
        self,
        *,
        species,
        seqid: str,
        start: int | None,
        stop: int | None,
    ) -> list[str]:
        sql = f"SELECT block_id from {self.table_name} WHERE species = ? AND seqid = ?"
        values = species, seqid
        if start is not None and stop is not None:
            # as long as start or stop are within the record start/stop, it's a match
            sql = f"{sql} AND ((start <= ? AND ? < stop) OR (start <= ? AND ? < stop))"
            values += (start, start, stop, stop)
        elif start is not None:
            # the aligned segment overlaps start
            sql = f"{sql} AND start <= ? AND ? < stop"
            values += (start, start)
        elif stop is not None:
            # the aligned segment overlaps stop
            sql = f"{sql} AND start <= ? AND ? < stop"
            values += (stop, stop)

        return self.db.execute(sql, values).fetchall()

    def get_records_matching(
        self,
        *,
        species,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
    ) -> typing.Iterable[AlignRecord]:
        # make sure python, not numpy, integers
        start = None if start is None else int(start)
        stop = None if stop is None else int(stop)

        # We need the block IDs for all records for a species whose coordinates
        # lie in the range (start, stop). We then search for all records with
        # each block id. We return full records.
        # Client code is responsible for creating Aligned sequence instances
        # and the Alignment.

        # TODO: there's an issue here with records being duplicated, solved
        #   for now by making AlignRecord hashable and using a set for block_ids
        block_ids = {
            r["block_id"]
            for r in self._get_block_id(
                species=species,
                seqid=seqid,
                start=start,
                stop=stop,
            )
        }
        values = ", ".join("?" * len(block_ids))
        sql = f"SELECT * from {self.table_name} WHERE block_id IN ({values})"
        results = defaultdict(set)
        for record in self.db.execute(sql, tuple(block_ids)).fetchall():
            record = {k: record[k] for k in record.keys()}
            index = record.pop("id")
            record["gap_spans"] = self.gap_store.get_record(index=index)
            results[record["block_id"]].add(AlignRecord(**record))

        return results.values()

    def get_species_names(self) -> list[str]:
        """return the list of species names"""
        return list(self.get_distinct("species"))


def get_alignment(
    align_db: AlignDb,
    genomes: dict,
    ref_species: str,
    seqid: str,
    ref_start: int | None = None,
    ref_end: int | None = None,
    namer: typing.Callable | None = None,
    mask_features: list[str] | None = None,
) -> typing.Iterable[Alignment]:
    """yields cogent3 Alignments"""

    if ref_species not in genomes:
        raise ValueError(f"unknown species {ref_species!r}")

    align_records = align_db.get_records_matching(
        species=ref_species,
        seqid=seqid,
        start=ref_start,
        stop=ref_end,
    )
    # sample the sequences
    for block in align_records:
        # we get the gaps corresponding to the reference sequence
        # and convert them to a IndelMap instance. We then convert
        # the ref_start, ref_end into align_start, align_end. Those values are
        # used for all other species -- they are converted into sequence
        # coordinates for each species -- selecting their sequence,
        # building the Aligned instance, and selecting the annotation subset.
        for align_record in block:
            if align_record.species == ref_species and align_record.seqid == seqid:
                # ref_start, ref_end are genomic positions and the align_record
                # start / stop are also genomic positions
                genome_start = align_record.start
                genome_end = align_record.stop
                gap_pos, gap_lengths = align_record.gap_data
                gaps = IndelMap(
                    gap_pos=gap_pos,
                    gap_lengths=gap_lengths,
                    parent_length=genome_end - genome_start,
                )

                # We use the IndelMap object to identify the alignment
                # positions the ref_start / ref_end correspond to. The alignment
                # positions are used below for slicing each sequence in the
                # alignment.

                # make sure the sequence start and stop are within this
                # aligned block
                seq_start = max(ref_start or genome_start, genome_start)
                seq_end = min(ref_end or genome_end, genome_end)
                # make these coordinates relative to the aligned segment
                if align_record.strand == "-":
                    # if record is on minus strand, then genome stop is
                    # the alignment start
                    seq_start, seq_end = genome_end - seq_end, genome_end - seq_start
                else:
                    seq_start = seq_start - genome_start
                    seq_end = seq_end - genome_start

                align_start = gaps.get_align_index(seq_start)
                align_end = gaps.get_align_index(seq_end)
                break
        else:
            raise ValueError(f"no matching alignment record for {ref_species!r}")

        seqs = {}
        for align_record in block:
            record_species = align_record.species
            genome = genomes[record_species]
            # We need to convert the alignment coordinates into sequence
            # coordinates for this species.
            genome_start = align_record.start
            genome_end = align_record.stop
            gap_pos, gap_lengths = align_record.gap_data
            gaps = IndelMap(
                gap_pos=gap_pos,
                gap_lengths=gap_lengths,
                parent_length=genome_end - genome_start,
            )

            # We use the alignment indices derived for the reference sequence
            # above
            seq_start = gaps.get_seq_index(align_start)
            seq_end = gaps.get_seq_index(align_end)
            seq_length = seq_end - seq_start
            if align_record.strand == "-":
                # if it's neg strand, the alignment start is the genome stop
                seq_start = gaps.parent_length - seq_end

            s = genome.get_seq(
                seqid=align_record.seqid,
                start=genome_start + seq_start,
                stop=genome_start + seq_start + seq_length,
                namer=namer,
                with_annotations=False,
            )
            # we now trim the gaps for this sequence to the sub-alignment
            gaps = gaps[align_start:align_end]

            if align_record.strand == "-":
                s = s.rc()

            if not namer:
                strand_symbol = -1 if align_record.strand == "-" else 1
                s.name = f"{s.name}:{strand_symbol}"

            aligned = Aligned(gaps, s)
            if aligned.name not in seqs:
                seqs[aligned.name] = aligned
            elif str(aligned) == str(seqs[aligned.name]):
                print(f"duplicated {s.name}")

        aln = Alignment(list(seqs.values()))
        aln.annotation_db = genome.annotation_db
        if mask_features:
            aln = aln.with_masked_annotations(biotypes=mask_features)

        yield aln


def _add_alignments(*alns, sep="?") -> Alignment:
    """concatenates alignments using sep as spacer"""
    defaults = ["?" * len(aln) for aln in alns]
    all_names = set()
    for aln in alns:
        all_names.update(set(aln.names))

    result = {n: [] for n in all_names}
    for aln, default in zip(alns, defaults, strict=False):
        data = aln.to_dict()
        for name in all_names:
            result[name].append(data.get(name, default))

    result = {n: sep.join(data) for n, data in result.items()}
    return Alignment(data=result, moltype=aln.moltype)


@define_app
class construct_alignment:
    """reassemble an alignment that maps to a given genomic segment

    If the segment spans multiple alignments these are joinded using
    the sep character.
    """

    def __init__(
        self,
        align_db: AlignDb,
        genomes: dict[str, elt_genome.Genome],
        mask_features: list[str] | None = None,
        sep: str = "?",
    ) -> None:
        self._align_db = align_db
        self._genomes = genomes
        self._mask_features = mask_features
        self._sep = sep

    def main(self, segment: elt_genome.genome_segment) -> list[Alignment]:
        results = []
        for aln in get_alignment(
            self._align_db,
            self._genomes,
            segment.species,
            segment.seqid,
            segment.start,
            segment.stop,
            mask_features=self._mask_features,
        ):
            aln.info.source = segment.source
            results.append(aln)

        return results
