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


@dataclass
class GapPositions:
    # 2D numpy int array,
    # each row is a gap
    # column 0 is gap index
    # column 1 is gap length
    gaps: numpy.ndarray
    seq_length: int

    def from_seq_to_align_index(self, index: int) -> int:
        """convert a sequence index into an alignment index"""
        # TODO edge cases that should raise an IndexError
        #  when the index is invalid given gaps and align length
        #  when the result lies outside the alignment
        #  negative indices

        first = numpy.argmax(self.gaps[:, 0] > index)
        return self.gaps[:first, 1].sum(axis=0) + index

    def from_align_to_seq_index(self, index: int) -> int:
        """converts alignment index to sequence index"""
        # TODO edge cases that should raise an IndexError
        #  if the result falls within a gap run at the sequence start
        #  OR end
        #  negative indices

        new_index = numpy.argmax(self.gaps.sum(axis=1).cumsum() > index)
        sub = self.gaps[:new_index, 1].sum()
        return index - sub
