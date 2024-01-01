from __future__ import annotations

import typing

from collections import defaultdict
from dataclasses import dataclass

import numpy

from cogent3 import make_seq
from cogent3.core.alignment import Alignment

from ensembl_lite._db_base import SqliteDbMixin, _compressed_array_proxy


class AlignRecordType(typing.TypedDict):
    source: str
    block_id: int
    species: str
    seqid: str
    start: int
    end: int
    strand: str
    gap_spans: numpy.ndarray


ReturnType = typing.Tuple[str, tuple]  # the sql statement and corresponding values


class AlignDb(SqliteDbMixin):
    # table schema for user provided annotations
    table_name = "align"
    _align_schema = {
        "source": "TEXT",  # the file path
        "block_id": "INT",  # an alignment number
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

    def add_records(self, records: typing.Sequence[AlignRecordType]):
        # bulk insert
        col_order = [
            row[1]
            for row in self.db.execute(
                f"PRAGMA table_info({self.table_name})"
            ).fetchall()
        ]
        for i in range(len(records)):
            records[i]["gap_spans"] = _compressed_array_proxy(records[i]["gap_spans"])
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
    ) -> list[int]:
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
    ):
        # We need the block IDs for all records for a species whose coordinates
        # lie in the range (start, end).
        # We then search for all records with each block id.
        # We return full records.
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
            results[record["block_id"]].append(AlignRecordType(**record))
        return results.values()


def get_alignment(
    align_db: AlignDb,
    genomes: dict,
    species: str,
    seqid: str,
    start: int | None = None,
    end: int | None = None,
) -> typing.Generator[Alignment]:
    """yields cogent3 Alignments"""
    from ensembl_lite.convert import gap_coords_to_seq

    if species not in genomes:
        raise ValueError(f"unknown species {species!r}")

    if start is not None:  # deal with case where numpy scalars are input
        start = int(start)
    if end is not None:
        end = int(end)

    # todo connect to annotation db that has been filtered for all records for
    #  each species that fall within the coordinates of the records
    align_records = align_db.get_records_matching(
        species=species, seqid=seqid, start=start, end=end
    )

    # sample the sequences
    for block in align_records:
        seqs = []
        # we get the gaps corresponding to the reference sequence
        # and convert them to a GapPosition instance. We then convert
        # the start, end into align_start, align_end. Those values are
        # used for all other species -- they are converted into sequence
        # coordinates for each species -- selecting their sequence and,
        # building the aligned instance, selecting the annotation subset.
        for align_record in block:
            if align_record["species"] == species and align_record["seqid"] == seqid:
                # start, end are genomic positions and the align_record
                # start / end are also genomic positions
                genome_start = align_record["start"]
                genome_end = align_record["end"]
                gaps = GapPositions(
                    align_record["gap_spans"], seq_length=genome_end - genome_start
                )

                # We use the GapPosition object to identify the alignment
                # positions the start / end correspond to. The alignment
                # positions are used below for slicing each sequence in the
                # alignment.

                # make sure the sequence start and end are within this
                # aligned block
                seq_start = max(start or genome_start, genome_start)
                seq_end = min(end or genome_end, genome_end)
                # make these coordinates relative to the aligned segment
                if align_record["strand"] == "-":
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
            raise ValueError(f"no matching alignment record for {species!r}")

        for align_record in block:
            species = align_record["species"]
            genome = genomes[species]
            seqid = align_record["seqid"]
            # We need to convert the alignment coordinates into sequence
            # coordinates for this species.
            genome_start = align_record["start"]
            genome_end = align_record["end"]
            gaps = GapPositions(
                align_record["gap_spans"], seq_length=genome_end - genome_start
            )

            # We use the alignment indices derived for the reference sequence
            # above
            seq_start = gaps.from_align_to_seq_index(align_start)
            seq_end = gaps.from_align_to_seq_index(align_end)
            seq_length = seq_end - seq_start
            if align_record["strand"] == "-":
                # if it's neg strand, the alignment start is the genome end
                seq_start = gaps.seq_length - seq_end
            s = genome.get_seq(
                seqid=seqid,
                start=genome_start + seq_start,
                end=genome_start + seq_start + seq_length,
            )
            # we now trim the gaps for this sequence to the sub-alignment
            gaps = gaps[align_start:align_end]

            s = make_seq(s, name=seqid, moltype="dna")
            if align_record["strand"] == "-":
                s = s.rc()

            aligned = gap_coords_to_seq(gaps.gaps, s)
            seqs.append(aligned)

        # todo need to add annotation_db
        yield Alignment(seqs)


def _adjust_gap_starts(gaps: numpy.ndarray, new_start: int) -> numpy.ndarray:
    return numpy.array([gaps.T[0] - new_start, gaps.T[1]], dtype=gaps.dtype).T


def _ends_within_gap(gaps: numpy.ndarray, align_index: int) -> numpy.ndarray:
    """
    return the gaps array ending at align_index

    Parameters
    ----------
    gaps
        gaps are [[seq index, gap length], ...]
    align_index
        position to seek
    """
    assert gaps[0, 0] <= align_index  # must fall within the gaps
    total_gaps = 0
    for i, (gap_index, gap_length) in enumerate(gaps):
        gap_start = gap_index + total_gaps
        gap_end = gap_start + gap_length
        # align index can fall between gaps or within a gap
        if align_index <= gap_start:
            new_gaps = gaps[:i]
            break
        elif gap_start < align_index <= gap_end:
            # we end within a gap, so adjust the length of that gap
            new_gaps = gaps[: i + 1]
            new_gaps[-1, 1] = align_index - gap_start
            break

        total_gaps += gap_length
    else:
        # align_index is after the last gap, so result has all gaps
        raise RuntimeError(f"{gaps=}  {align_index=}")  # not copying! bad idea?

    return new_gaps


def _starts_within_gap(gaps: numpy.ndarray, align_index: int) -> numpy.ndarray:
    """
    return the gaps array starting with align_index

    Parameters
    ----------
    gaps
        gaps are [[seq index, gap length], ...]
    align_index
        position to seek
    """
    assert gaps[0, 0] <= align_index  # must fall within the gaps
    total_gaps = 0
    for i, (gap_index, gap_length) in enumerate(gaps):
        gap_start = gap_index + total_gaps
        gap_end = gap_start + gap_length
        # align index can fall between gaps or within a gap
        if gap_start <= align_index < gap_end:
            new_gaps = gaps[i:]
            gap_start_diff = align_index - gap_start
            new_gaps[0, 1] = gap_length - gap_start_diff
            break
        elif align_index < gap_start:
            new_gaps = gaps[i:]
            break
        elif align_index == gap_end:
            new_gaps = gaps[i + 1 :]
            break
        total_gaps += gap_length
    else:
        # align_index is after the last gap, so result has all gaps
        raise RuntimeError(f"{gaps=}  {align_index=}")

    return new_gaps


def _within_a_gap(gaps: numpy.ndarray, start: int, stop: int) -> bool:
    """return True if start/stop fall within a gap

    Parameters
    ----------
    gaps
        numpy 2D array
    start, stop
        start and stop are align indices
    """
    # todo convert to numba
    cumsum_gap_length = 0
    for gap_start, gap_length in gaps:
        gap_start += cumsum_gap_length
        if gap_start <= start < stop <= gap_start + gap_length:
            return True
        cumsum_gap_length += gap_length
    return False


@dataclass
class GapPositions:
    # 2D numpy int array,
    # each row is a gap
    # column 0 is sequence index of gap **relative to the alignment**
    # column 1 is gap length
    gaps: numpy.ndarray
    seq_length: int

    def __post_init__(self):
        if not len(self.gaps):
            # can get have a zero length array with shape != (0, 0)
            # e.g. by slicing gaps[:0], but since there's no data
            # we force it to have zero elements on both dimensions
            self.gaps = self.gaps.reshape((0, 0))

        # make gap array immutable
        self.gaps.flags.writeable = False

    def __getitem__(self, item: slice) -> typing.Self:
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
        total_gap_length = gaps[:, 1].sum() if len(gaps) else 0
        if (
            not total_gap_length
            or stop <= gaps[0, 0]
            or start > gaps[-1][0] + total_gap_length
        ):
            # no gaps or slice is before first gap or after last gap
            gaps = numpy.empty(shape=(0, 0), dtype=gaps.dtype)
            return type(self)(gaps=gaps, seq_length=stop - start)

        if start < gaps[0, 0] and stop > gaps[-1, 0] + total_gap_length:
            # slice result contains all gaps and shift gap left
            gaps = _adjust_gap_starts(gaps, start)
            seq_length = self.from_align_to_seq_index(stop) - start
            return type(self)(gaps=gaps, seq_length=seq_length)
        elif start < gaps[0, 0]:
            # start is in seq coords, ends within gaps
            seq_length = self.from_align_to_seq_index(stop) - start
            gaps = _ends_within_gap(gaps, stop)
            gaps = _adjust_gap_starts(gaps, start)
        elif stop > total_gap_length + gaps[-1, 0]:
            # slice starts within gaps
            gaps = _starts_within_gap(gaps, start)
            seq_start = self.from_align_to_seq_index(start)
            gaps = _adjust_gap_starts(gaps, seq_start)
            seq_length = self.from_align_to_seq_index(stop) - seq_start
        elif _within_a_gap(gaps, start, stop):
            # falls within a gap
            gaps = numpy.array([[0, stop - start]], dtype=gaps.dtype)
            return type(self)(gaps=gaps, seq_length=0)
        else:
            # slice is within the gaps
            gaps = _ends_within_gap(gaps, stop)
            gaps = _starts_within_gap(gaps, start)
            seq_start = self.from_align_to_seq_index(start)
            gaps = _adjust_gap_starts(gaps, seq_start)
            seq_length = self.from_align_to_seq_index(stop) - seq_start

        return type(self)(gaps=gaps, seq_length=seq_length)

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
        gaps = self.gaps.copy()
        total_gaps = 0
        for gap_index, gap_length in gaps:
            gap_start = gap_index + total_gaps
            gap_end = gap_start + gap_length
            if align_index < gap_start:
                seq_index = align_index - total_gaps
                break
            if gap_start <= align_index <= gap_end:
                # align_index between gaps
                seq_index = gap_index
                break
            total_gaps += gap_length
        else:
            # align_index is after the last gap
            seq_index = align_index - total_gaps
        return seq_index
