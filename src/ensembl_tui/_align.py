import pathlib
import typing
import uuid
from collections import defaultdict
from dataclasses import dataclass

import h5py
import numpy
import typing_extensions
from cogent3.app.composable import define_app
from cogent3.core.alignment import Aligned, Alignment
from cogent3.core.location import _DEFAULT_GAP_DTYPE, IndelMap

from ensembl_tui import _genome as eti_genome
from ensembl_tui import _storage_mixin as eti_storage
from ensembl_tui import _util as eti_util

_no_gaps = numpy.array([], dtype=_DEFAULT_GAP_DTYPE)

GAP_STORE_SUFFIX = "indels-hdf5_blosc2"
ALIGN_STORE_SUFFIX = "align_coords-sqlitedb"

ALIGN_ATTR_SCHEMA = (
    "align_id INTEGER PRIMARY KEY DEFAULT nextval('align_id_seq')",
    "source TEXT",
    "block_id BIGINT",
    "species TEXT",
    "seqid TEXT",
    "start INTEGER",
    "stop INTEGER",
    "strand TEXT",
)
ALIGN_ATTR_COLS = eti_util.make_column_constant(ALIGN_ATTR_SCHEMA)

VT = str | int | numpy.ndarray

# We store alignment data in two parts: the coordinates and the gaps
# this split was motivated by the limitation of sqlite3 - storing everying
# in a row meant querying would require loading a row into memory. With the
# move to parquet it may be possible to just add gap data as a column in the
# parquet file with BLOB type.


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

    def __getitem__(self, item: str) -> VT:
        return getattr(self, item)

    def __setitem__(self, item: str, value: VT) -> None:
        setattr(self, item, value)

    def __eq__(self, other: typing_extensions.Self) -> bool:
        attrs = "block_id", "species", "seqid", "start", "stop", "strand"
        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return (self.gap_spans == other.gap_spans).all()

    def __hash__(self) -> int:
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
    def gap_data(self) -> tuple[numpy.ndarray, numpy.ndarray]:
        if len(self.gap_spans):
            gap_pos, gap_lengths = self.gap_spans.T
        else:
            gap_pos, gap_lengths = _no_gaps.copy(), _no_gaps.copy()

        return gap_pos, gap_lengths


ReturnType = tuple[str, tuple]  # the sql statement and corresponding values


class GapStore(eti_storage.Hdf5Mixin):
    """store gap data from aligned sequences"""

    def __init__(
        self,
        source: eti_util.PathType,
        align_name: str | None = None,
        mode: str = "r",
        in_memory: bool = False,
    ) -> None:
        """
        Parameters
        ----------
        source
            path to the file
        align_name
            the Ensembl alignment name, e.g. '10_primates.epo'
        mode
            file open mode
        in_memory
            an in-memory HDF5 file
        """
        in_memory = in_memory or "memory" in str(source)
        # h5py requires a unique file path even for in-memory stores
        source = uuid.uuid4().hex if in_memory else source

        self.source = pathlib.Path(source)
        self.mode = "w-" if mode == "w" else mode
        h5_kwargs = (
            {
                "driver": "core",
                "backing_store": False,
            }
            if in_memory
            else {}
        )
        try:
            self._file = h5py.File(source, mode=self.mode, **h5_kwargs)
        except OSError as err:
            msg = f"h5py could not open {source=}"
            raise OSError(msg) from err

        if "r" not in self.mode and "align_name" not in self._file.attrs:
            if not align_name:
                msg = "align_name not present in GapStore attrs"
                raise ValueError(msg)
            self._file.attrs["align_name"] = align_name
        if (
            align_name
            and (file_species := self._file.attrs.get("align_name", None)) != align_name
        ):
            msg = f"{self.source.name!r} {file_species!r} != {align_name}"
            raise ValueError(msg)
        self.align_name = self._file.attrs["align_name"]

    def add_record(self, *, index: int | str, gaps: numpy.ndarray) -> None:
        # dataset names must be strings
        index = str(index)
        if index in self._file:
            stored = self._file[index]
            if (gaps == stored).all():
                # already seen this index
                return
            # but it's different, which is a problem
            msg = f"{index!r} already present but with different gaps"
            raise ValueError(msg)
        self._file.create_dataset(
            name=index,
            data=gaps,
            chunks=True,
            **eti_util.HDF5_BLOSC2_KWARGS,
        )
        self._file.flush()

    def get_record(self, *, index: int | str) -> numpy.ndarray:
        return self._file[str(index)][:]  # type: ignore


# TODO add a table and methods to support storing the species tree used
#  for the alignment and for getting the species tree
@dataclass(slots=True)
class AlignDb(eti_storage.DuckdbParquetBase):
    _tables: tuple[str, ...] = ("align_blocks",)
    gap_store: GapStore | None = None

    def _post_init(self) -> None:
        source = self.source
        if source.name == ":memory:":
            gap_path = "memory"
            kwargs = {"in_memory": True}
        else:
            gap_path = source / f"{source.stem}.{GAP_STORE_SUFFIX}"
            kwargs = {"in_memory": False}

        if self.gap_store is not None:
            return

        self.gap_store = GapStore(
            source=gap_path,
            align_name=source.stem,
            mode="r",
            **kwargs,
        )

    def _get_block_id(
        self,
        *,
        species: str,
        seqid: str,
        start: int | None,
        stop: int | None,
    ) -> list[str]:
        sql = f"SELECT block_id from {self._tables[0]} WHERE species = ? AND seqid = ?"
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

        return self.conn.sql(sql, params=values).fetchall()

    def get_records_matching(
        self,
        *,
        species: str,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
    ) -> typing.Iterable[set[AlignRecord]]:
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
            r[0]
            for r in self._get_block_id(
                species=species,
                seqid=seqid,
                start=start,
                stop=stop,
            )
        }
        if not block_ids:
            return []
        col_order = ", ".join(ALIGN_ATTR_COLS)
        values = ", ".join("?" * len(block_ids))
        sql = f"SELECT {col_order} from {self._tables[0]} WHERE block_id IN ({values})"
        results = defaultdict(set)
        for record in self.conn.sql(sql, params=tuple(block_ids)).fetchall():
            data = dict(zip(ALIGN_ATTR_COLS, record, strict=True))
            index = data.pop("align_id")
            data["gap_spans"] = self.gap_store.get_record(index=index)
            results[data["block_id"]].add(AlignRecord(**data))

        return results.values()

    def get_species_names(self) -> list[str]:
        """return the list of species names"""
        return list(self.get_distinct("species"))

    def get_distinct(self, field: str) -> list[str]:
        return [
            r[0]
            for r in self.conn.sql(
                f"SELECT DISTINCT {field} from {self._tables[0]}",
            ).fetchall()
        ]

    def num_records(self) -> int:
        return self.conn.sql(f"SELECT COUNT(*) from {self._tables[0]}").fetchone()[0]


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
            msg = f"no matching alignment record for {ref_species!r}"
            raise ValueError(msg)

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


@define_app
class construct_alignment:  # noqa: N801
    """reassemble an alignment that maps to a given genomic segment

    If the segment spans multiple alignments these are joinded using
    the sep character.
    """

    def __init__(
        self,
        align_db: AlignDb,
        genomes: dict[str, eti_genome.Genome],
        mask_features: list[str] | None = None,
        sep: str = "?",
    ) -> None:
        self._align_db = align_db
        self._genomes = genomes
        self._mask_features = mask_features
        self._sep = sep

    def main(self, segment: eti_genome.genome_segment) -> list[Alignment]:
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
