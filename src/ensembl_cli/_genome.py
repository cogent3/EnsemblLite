import typing

from ensembl_cli._db_base import SqliteDbMixin


OptInt = typing.Optional[int]


# make a plugin hook that queries for sequence strings, this is our distributed version
class Genome(SqliteDbMixin):
    table_name = "genome"
    _genome_schema = {"coord_name": "TEXT PRIMARY KEY", "seq": "TEXT", "length": "INT"}
    _metadata_schema = {"species": "TEXT"}

    def __init__(self, *, source: str = ":memory:", species: str = None):
        self.source = source
        self._init_tables()
        # the metadata table stores species info
        self._execute_sql("INSERT INTO metadata(species) VALUES (?)", (species,))
        self.db.commit()

    def add_record(self, *, coord_name: str, seq: str):
        sql = f"INSERT INTO {self.table_name}(coord_name, seq, length) VALUES (?, ?, ?)"
        self._execute_sql(sql, (coord_name, seq, len(seq)))
        self.db.commit()

    def add_records(self, *, records: typing.Iterable[list[str, str]]):
        sql = f"INSERT INTO {self.table_name}(coord_name, seq, length) VALUES (?, ?, ?)"
        self.db.executemany(sql, [(n, s, len(s)) for n, s in records])

    def get_seq(
        self, *, coord_name: str, start: OptInt = None, end: OptInt = None
    ) -> str:
        """

        Parameters
        ----------
        coord_name
            name of chromosome etc..
        start
            starting position of slice in python coordinates, defaults
            to 0
        end
            ending position of slice in python coordinates, defaults
            to length of coordinate
        """
        if start is not None:
            start += 1  # SQLite counts from 1
        else:
            start = 1

        if end is None:
            sql = f"SELECT SUBSTR(seq, ?, length) FROM {self.table_name} where coord_name = ?"
            values = start, coord_name
        else:
            end -= start - 1
            sql = (
                f"SELECT SUBSTR(seq, ?, ?) FROM {self.table_name} where coord_name = ?"
            )
            values = start, end, coord_name

        return self._execute_sql(sql, values).fetchone()[0]
