import gzip
import inspect
import sqlite3
import typing

import numpy


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


def array_to_sqlite(data):
    return gzip.compress(data.tobytes())


def sqlite_to_array(data):
    result = numpy.frombuffer(gzip.decompress(data), dtype=numpy.int32)
    dim = result.shape[0] // 2
    return result.reshape((dim, 2))


# registering the conversion functions with sqlite
sqlite3.register_adapter(numpy.ndarray, array_to_sqlite)
sqlite3.register_converter("array", sqlite_to_array)


def _make_table_sql(
    table_name: str,
    columns: dict,
) -> str:
    """makes the SQL for creating a table

    Parameters
    ----------
    table_name : str
        name of the table
    columns : dict
        {<column name>: <column SQL type>, ...}

    Returns
    -------
    str
    """
    columns_types = ", ".join([f"{name} {ctype}" for name, ctype in columns.items()])
    return f"CREATE TABLE IF NOT EXISTS {table_name} ({columns_types});"


def _add_record_sql(
    table_name: str,
    data: dict,
) -> ReturnType:
    """creates SQL defining the table

    Parameters
    ----------
    table_name : str
        name of the table
    data : dict
        {<column name>: column type}

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """
    cols = ", ".join(k.lower() for k in data)
    pos = ", ".join("?" * len(data))
    sql = f"INSERT INTO {table_name} ({cols}) VALUES ({pos});"
    return sql, tuple(data.values())


class SqliteDbMixin:
    table_name = None
    _db = None
    source = None

    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        init_sig = inspect.signature(cls.__init__)
        bargs = init_sig.bind_partial(cls, *args, **kwargs)
        bargs.apply_defaults()
        init_vals = bargs.arguments
        init_vals.pop("self", None)
        obj._serialisable = init_vals
        return obj

    def __repr__(self):
        name = self.__class__.__name__
        total_records = len(self)
        args = ", ".join(
            f"{k}={repr(v) if isinstance(v, str) else v}"
            for k, v in self._serialisable.items()
            if k != "data"
        )
        return f"{name}({args}, total_records={total_records})"

    def __len__(self):
        return self.num_records()

    def __eq__(self, other):
        return isinstance(other, self.__class__) and other.db is self.db

    def _init_tables(self) -> None:
        # is source an existing db
        self._db = sqlite3.connect(
            self.source,
            detect_types=sqlite3.PARSE_DECLTYPES,
            check_same_thread=False,
        )
        self._db.row_factory = sqlite3.Row

        # A bit of magic.
        # Assumes schema attributes named as `_<table name>_schema`
        for attr in dir(self):
            if attr.endswith("_schema"):
                table_name = attr.split("_")[1]
                attr = getattr(self, attr)
                sql = _make_table_sql(table_name, attr)
                self._execute_sql(sql)

    @property
    def db(self) -> sqlite3.Connection:
        if self._db is None:
            self._db = sqlite3.connect(
                self.source,
                detect_types=sqlite3.PARSE_DECLTYPES,
                check_same_thread=False,
            )
            self._db.row_factory = sqlite3.Row

        return self._db

    def _execute_sql(self, cmnd: str, values=None) -> sqlite3.Cursor:
        with self.db:
            # context manager ensures safe transactions
            cursor = self.db.cursor()
            cursor.execute(cmnd, values or [])
            return cursor

    def num_records(self):
        sql = f"SELECT COUNT(*) as count FROM {self.table_name}"
        return list(self._execute_sql(sql).fetchone())[0]

    def close(self):
        self.db.commit()
        self.db.close()
