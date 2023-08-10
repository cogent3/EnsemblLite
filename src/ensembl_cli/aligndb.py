import typing

import numpy

from ensembl_cli._db_base import SqliteDbMixin


class AlignRecordType(typing.TypedDict):
    source: str
    block_id: str
    species: str
    coord_name: str
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
        "coord_name": "TEXT",
        "start": "INTEGER",
        "end": "INTEGER",
        "strand": "TEXT",
        "gap_spans": "array",
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

    def add_records(self, records: typing.Iterable[AlignRecordType]):
        # bulk insert
        col_order = [
            row[1]
            for row in self.db.execute(
                f"PRAGMA table_info({self.table_name})"
            ).fetchall()
        ]
        for i in range(len(records)):
            records[i] = [records[i][c] for c in col_order]

        val_placeholder = ", ".join("?" * len(col_order))
        sql = f"INSERT INTO {self.table_name} ({', '.join(col_order)}) VALUES ({val_placeholder})"
        self.db.executemany(sql, records)
