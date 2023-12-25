from __future__ import annotations

import typing

from collections import defaultdict
from dataclasses import dataclass

import numpy

from cogent3 import make_aligned_seqs, make_seq
from cogent3.core.alignment import Aligned, Alignment

from ensembl_lite._db_base import SqliteDbMixin, _compressed_array_proxy


class AlignRecordType(typing.TypedDict):
    source: str
    block_id: int
    species: str
    coord_name: str
    start: int
    end: int
    strand: str
    gap_spans: numpy.ndarray


ReturnType = typing.Tuple[str, tuple]  # the sql statement and corresponding values
OptInt = typing.Optional[int]
OptStr = typing.Optional[int]


class AlignDb(SqliteDbMixin):
    # table schema for user provided annotations
    table_name = "align"
    _align_schema = {
        "source": "TEXT",  # the file path
        "block_id": "INT",  # an alignment number
        "species": "TEXT",
        "coord_name": "TEXT",
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
        coord_name: str,
        start: int | None,
        end: int | None,
    ) -> list[int]:
        sql = f"SELECT block_id from {self.table_name} WHERE species = ? AND coord_name = ?"
        values = species, coord_name
        if start and end:
            sql = f"{sql} AND start > ? AND end < ?"
            values += (start, end)
        elif start:
            sql = f"{sql} AND start > ?"
            values += (start,)
        elif end:
            sql = f"{sql} AND end < ?"
            values += (end,)

        return self.db.execute(sql, values).fetchall()

    def get_records_matching(
        self,
        *,
        species,
        coord_name: str,
        start: OptInt = None,
        end: OptInt = None,
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
                species=species, coord_name=coord_name, start=start, end=end
            )
        ]

        values = " ".join("?" * len(block_ids))
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
    coord_name: str,
    start=None,
    end=None,
) -> Alignment:
    """return a cogent3 Alignment"""
    from ensembl_lite.convert import gap_coords_to_seq

    # todo connect to annotation db has been filtered for all records for
    #  each species that fall within
    # the coordinates of the records
    align_records = align_db.get_records_matching(
        species=species, coord_name=coord_name, start=start, end=end
    )
    # sample the sequences
    seqs = []
    for block in align_records:
        for record in block:
            species = record["species"]
            genome = genomes[species]
            s = genome.get_seq(
                coord_name=record["coord_name"],
                start=record["start"],
                end=record["end"],
            )
            s = make_seq(s, name=record["coord_name"], moltype="dna")
            aligned = gap_coords_to_seq(record["gap_spans"], s)
            seqs.append(aligned)

    return Alignment(seqs)


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
        if gap_start == align_index:
            new_gaps = gaps[:i]
            break
        elif align_index == gap_end:
            new_gaps = gaps[: i + 1]
            break
        elif align_index < gap_start:
            # align_index is before this gap
            new_gaps = gaps[:i]
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
            gap_start_diff = align_index - gap_index
            new_gaps[0, 1] = gap_length - gap_start_diff
            break
        elif align_index == gap_end or align_index < gap_start:
            new_gaps = gaps[i + 1 :]
            break
        total_gaps += gap_length
    else:
        # align_index is after the last gap, so result has all gaps
        raise RuntimeError(f"{gaps=}  {align_index=}")

    return new_gaps


@dataclass
class GapPositions:
    # 2D numpy int array,
    # each row is a gap
    # column 0 is sequence index of gap
    # column 1 is gap length
    gaps: numpy.ndarray
    seq_length: int

    def __post_init__(self):
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
        # slice is before first gap or after last gap
        if stop < gaps[0, 0] or start > gaps[-1].sum():
            gaps = numpy.empty(shape=(0, 0), dtype=gaps.dtype)
            return type(self)(gaps=gaps, seq_length=stop - start)

        total_gaps = gaps[:, 1].sum()
        if start < gaps[0, 0] and stop > gaps[-1, 0] + total_gaps:
            # slice result contains all gaps, so we shift gaps left
            gaps = _adjust_gap_starts(gaps, start)
            seq_length = self.from_align_to_seq_index(stop) - start
            return type(self)(gaps=gaps, seq_length=seq_length)

        if start < gaps[0, 0]:
            # start is in seq coords, ends within gaps
            seq_length = self.from_align_to_seq_index(stop) - start
            gaps = _ends_within_gap(gaps, stop)
            gaps = _adjust_gap_starts(gaps, start)
        elif stop > total_gaps + gaps[-1, 0]:
            # slice starts within gaps
            gaps = _starts_within_gap(gaps, start)
            seq_start = self.from_align_to_seq_index(start)
            gaps = _adjust_gap_starts(gaps, seq_start)
            seq_length = self.from_align_to_seq_index(stop) - seq_start
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
