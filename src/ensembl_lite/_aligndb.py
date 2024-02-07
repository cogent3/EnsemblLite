from __future__ import annotations

import typing

from collections import defaultdict
from dataclasses import dataclass

import numpy

from cogent3.core.alignment import Alignment

from ensembl_lite._db_base import SqliteDbMixin, _compressed_array_proxy


@dataclass(slots=True)
class AlignRecord:
    """a record from an AlignDb

    Notes
    -----
    Can return fields as attributes or like a dict using the field name as
    a string.
    """

    source: str
    block_id: str
    species: str
    seqid: str
    start: int
    end: int
    strand: str
    gap_spans: numpy.ndarray

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, item, value):
        setattr(self, item, value)

    def __eq__(self, other):
        attrs = "source", "block_id", "species", "seqid", "start", "end", "strand"
        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return (self.gap_spans == other.gap_spans).all()


ReturnType = tuple[str, tuple]  # the sql statement and corresponding values


# todo add a table and methods to support storing the species tree used
#  for the alignment and for getting the species tree
class AlignDb(SqliteDbMixin):
    table_name = "align"
    _align_schema = {
        "source": "TEXT",  # the file path
        "block_id": "TEXT",  # <source file path>-<alignment number>
        "species": "TEXT",
        "seqid": "TEXT",
        "start": "INTEGER",
        "end": "INTEGER",
        "strand": "TEXT",
        "gap_spans": "compressed_array",
    }

    def __init__(self, *, source=":memory:"):
        """
        Parameters
        ----------
        source
            location to store the db, defaults to in memory only
        """
        # note that data is destroyed
        self.source = source
        self._db = None
        self._init_tables()

    def add_records(self, records: typing.Sequence[AlignRecord]):
        # bulk insert
        col_order = [
            row[1]
            for row in self.db.execute(
                f"PRAGMA table_info({self.table_name})"
            ).fetchall()
        ]
        for i in range(len(records)):
            records[i].gap_spans = _compressed_array_proxy(records[i].gap_spans)
            records[i] = [records[i][c] for c in col_order]

        val_placeholder = ", ".join("?" * len(col_order))
        sql = f"INSERT INTO {self.table_name} ({', '.join(col_order)}) VALUES ({val_placeholder})"
        self.db.executemany(sql, records)

    def _get_block_id(
        self,
        *,
        species,
        seqid: str,
        start: int | None,
        end: int | None,
    ) -> list[str]:
        sql = f"SELECT block_id from {self.table_name} WHERE species = ? AND seqid = ?"
        values = species, seqid
        if start is not None and end is not None:
            # as long as start or end are within the record start/end, it's a match
            sql = f"{sql} AND ((start <= ? AND ? < end) OR (start <= ? AND ? < end))"
            values += (start, start, end, end)
        elif start is not None:
            # the aligned segment overlaps start
            sql = f"{sql} AND start <= ? AND ? < end"
            values += (start, start)
        elif end is not None:
            # the aligned segment overlaps end
            sql = f"{sql} AND start <= ? AND ? < end"
            values += (end, end)

        return self.db.execute(sql, values).fetchall()

    def get_records_matching(
        self,
        *,
        species,
        seqid: str,
        start: int | None = None,
        end: int | None = None,
    ) -> typing.Iterable[AlignRecord]:
        # make sure python, not numpy, integers
        start = None if start is None else int(start)
        end = None if end is None else int(end)

        # We need the block IDs for all records for a species whose coordinates
        # lie in the range (start, end). We then search for all records with
        # each block id. We return full records.
        # Client code is responsible for creating Aligned sequence instances
        # and the Alignment.

        block_ids = [
            r["block_id"]
            for r in self._get_block_id(
                species=species, seqid=seqid, start=start, end=end
            )
        ]

        values = ", ".join("?" * len(block_ids))
        sql = f"SELECT * from {self.table_name} WHERE block_id IN ({values})"
        results = defaultdict(list)
        for record in self.db.execute(sql, block_ids).fetchall():
            record = {k: record[k] for k in record.keys()}
            results[record["block_id"]].append(AlignRecord(**record))

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
) -> typing.Generator[Alignment]:
    """yields cogent3 Alignments"""
    from ensembl_lite.convert import gap_coords_to_seq

    if ref_species not in genomes:
        raise ValueError(f"unknown species {ref_species!r}")

    align_records = align_db.get_records_matching(
        species=ref_species, seqid=seqid, start=ref_start, end=ref_end
    )

    # sample the sequences
    for block in align_records:
        # we get the gaps corresponding to the reference sequence
        # and convert them to a GapPosition instance. We then convert
        # the ref_start, ref_end into align_start, align_end. Those values are
        # used for all other species -- they are converted into sequence
        # coordinates for each species -- selecting their sequence,
        # building the Aligned instance, and selecting the annotation subset.
        for align_record in block:
            if align_record.species == ref_species and align_record.seqid == seqid:
                # ref_start, ref_end are genomic positions and the align_record
                # start / end are also genomic positions
                genome_start = align_record.start
                genome_end = align_record.end
                gaps = GapPositions(
                    align_record.gap_spans, seq_length=genome_end - genome_start
                )

                # We use the GapPosition object to identify the alignment
                # positions the ref_start / ref_end correspond to. The alignment
                # positions are used below for slicing each sequence in the
                # alignment.

                # make sure the sequence start and end are within this
                # aligned block
                seq_start = max(ref_start or genome_start, genome_start)
                seq_end = min(ref_end or genome_end, genome_end)
                # make these coordinates relative to the aligned segment
                if align_record.strand == "-":
                    # if record is on minus strand, then genome end is
                    # the alignment start
                    seq_start, seq_end = genome_end - seq_end, genome_end - seq_start
                else:
                    seq_start = seq_start - genome_start
                    seq_end = seq_end - genome_start

                align_start = gaps.from_seq_to_align_index(seq_start)
                align_end = gaps.from_seq_to_align_index(seq_end)
                break
        else:
            raise ValueError(f"no matching alignment record for {ref_species!r}")

        seqs = []
        for align_record in block:
            record_species = align_record.species
            genome = genomes[record_species]
            # We need to convert the alignment coordinates into sequence
            # coordinates for this species.
            genome_start = align_record.start
            genome_end = align_record.end
            gaps = GapPositions(
                align_record.gap_spans, seq_length=genome_end - genome_start
            )

            # We use the alignment indices derived for the reference sequence
            # above
            seq_start = gaps.from_align_to_seq_index(align_start)
            seq_end = gaps.from_align_to_seq_index(align_end)
            seq_length = seq_end - seq_start
            if align_record.strand == "-":
                # if it's neg strand, the alignment start is the genome end
                seq_start = gaps.seq_length - seq_end

            s = genome.get_seq(
                seqid=align_record.seqid,
                start=genome_start + seq_start,
                end=genome_start + seq_start + seq_length,
                namer=namer,
            )
            # we now trim the gaps for this sequence to the sub-alignment
            gaps = gaps[align_start:align_end]

            if align_record.strand == "-":
                s = s.rc()

            aligned = gap_coords_to_seq(gaps.gaps, s)
            seqs.append(aligned)

        yield Alignment(seqs)


def _gap_spans_cum_lengths(
    gaps: numpy.ndarray,
) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """returns 1D arrays in alignment coordinates of
    gap start, gap end, cumulative gap length"""
    if not len(gaps):
        r = numpy.array([], dtype=gaps.dtype)
        return r, r, r

    cumsum = gaps[:, 1].cumsum()
    sum_to_prev = 0
    gap_starts = numpy.empty(gaps.shape[0], dtype=gaps.dtype)
    gap_ends = numpy.empty(gaps.shape[0], dtype=gaps.dtype)
    for i, (p, l) in enumerate(gaps):
        gap_start = sum_to_prev + p
        gap_end = gap_start + l
        sum_to_prev = cumsum[i]
        gap_starts[i] = gap_start
        gap_ends[i] = gap_end

    return numpy.array(gap_starts), numpy.array(gap_ends), cumsum


@dataclass(slots=True)
class GapPositions:
    """records gap insertion index and length

    Notes
    -----
    This very closely parallels the cogent3.core.location.Map class,
    but is more memory efficient. When that class has been updated,
    this can be removed.
    """

    # 2D numpy int array,
    # each row is a gap
    # column 0 is sequence index of gap **relative to the alignment**
    # column 1 is gap length
    gaps: numpy.ndarray
    # length of the underlying sequence
    seq_length: int

    def __post_init__(self):
        if not len(self.gaps):
            # can get a zero length array with shape != (0, 0)
            # e.g. by slicing gaps[:0], but since there's no data
            # we force it to have zero elements on both dimensions
            self.gaps = self.gaps.reshape((0, 0))

        # make gap array immutable
        self.gaps.flags.writeable = False

    def __getitem__(self, item: slice) -> typing.Self:
        # we're assuming that this gap object is associated with a sequence
        # that will also be sliced. Hence, we need to shift the gap insertion
        # positions relative to this newly sliced sequence.
        if item.step:
            raise NotImplementedError(
                f"{type(self).__name__!r} does not support strides"
            )
        start = item.start or 0
        stop = item.stop or len(self)
        gaps = self.gaps.copy()
        if start < 0 or stop < 0:
            raise NotImplementedError(
                f"{type(self).__name__!r} does not support negative indexes"
            )

        # spans is in alignment indices
        # [(gap start, gap end), ...]
        gap_starts, gap_ends, cum_lengths = _gap_spans_cum_lengths(gaps)

        if not len(gaps) or stop < gap_starts[0] or start >= gap_ends[-1]:
            return type(self)(
                gaps=numpy.array([], dtype=gaps.dtype), seq_length=stop - start
            )

        # second column of spans is gap ends
        # which gaps does it fall between
        l = numpy.searchsorted(gap_ends, start, side="left")
        if gap_starts[l] <= start < gap_ends[l]:
            # start is within a gap
            begin = l
            begin_diff = start - gap_starts[l]
            gaps[l, 1] -= begin_diff
            shift = start - cum_lengths[l - 1] - begin_diff if l else gaps[l, 0]
        elif start == gap_ends[l]:
            # at gap boundary
            begin = l + 1
            shift = start - cum_lengths[l]
        else:
            # not within a gap
            begin = l
            shift = start - cum_lengths[l - 1] if l else start

        # start search for end from l index
        r = numpy.searchsorted(gap_ends[l:], stop, side="right") + l
        if r == len(gaps):
            # stop is after last gap
            end = r
        elif gap_starts[r] < stop <= gap_ends[r]:
            # within gap
            end = r + 1
            end_diff = gap_ends[r] - stop
            gaps[r, 1] -= end_diff
        else:
            end = r

        result = gaps[begin:end]
        result[:, 0] -= shift
        if not len(result):
            # no gaps
            seq_length = stop - start
        else:
            seq_length = self.from_align_to_seq_index(
                stop
            ) - self.from_align_to_seq_index(start)

        return type(self)(gaps=result, seq_length=seq_length)

    def __len__(self):
        total_gaps = self.gaps[:, 1].sum() if len(self.gaps) else 0
        return total_gaps + self.seq_length

    def from_seq_to_align_index(self, seq_index: int) -> int:
        """convert a sequence index into an alignment index"""
        if seq_index < 0:
            raise NotImplementedError(f"{seq_index} negative align_index not supported")
        # TODO convert this to numba function

        total_gaps = 0
        for gap_index, gap_length in self.gaps:
            if seq_index < gap_index:
                break

            total_gaps += gap_length

        return seq_index + total_gaps

    def from_align_to_seq_index(self, align_index: int) -> int:
        """converts alignment index to sequence index"""
        if align_index < 0:
            raise NotImplementedError(
                f"{align_index} negative align_index not supported"
            )

        # TODO convert this to numba function
        total_gaps = 0
        for seq_pos, gap_length in self.gaps:
            aln_start = seq_pos + total_gaps
            aln_end = aln_start + gap_length
            if align_index < aln_start:
                seq_index = align_index - total_gaps
                break
            if aln_start <= align_index <= aln_end:
                # align_index between gaps
                seq_index = seq_pos
                break
            total_gaps += gap_length
        else:
            # align_index is after the last gap
            seq_index = align_index - total_gaps
        return seq_index
