import dataclasses
import inspect
import sqlite3
import typing

import numpy

from ensembl_lite.util import blosc_compress_it, blosc_decompress_it


@dataclasses.dataclass
class _compressed_array_proxy:
    """this exists only to automate conversion of a customised sqlite type"""

    array: numpy.ndarray


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

_compressor = blosc_compress_it()
_decompressor = blosc_decompress_it()


def compressed_array_to_sqlite(data):
    return _compressor(data.array.astype(numpy.int32).tobytes())


def decompressed_sqlite_to_array(data):
    result = numpy.frombuffer(_decompressor(data), dtype=numpy.int32)
    if len(result):
        dim = result.shape[0] // 2
        result = result.reshape((dim, 2))
    return result


# registering the conversion functions with sqlite
# since these conversion functions are tied to a type, need to ensure the
# type will be unique to this tool, best way is to use <libname_type> and
# wrap a fundamental type with a proxy
sqlite3.register_adapter(_compressed_array_proxy, compressed_array_to_sqlite)
sqlite3.register_converter("compressed_array", decompressed_sqlite_to_array)


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
        obj._init_vals = init_vals
        return obj

    def __getstate__(self):
        return {**self._init_vals}

    def __setstate__(self, state):
        # this will reset connections to read only db's
        obj = self.__class__(**state)
        self.__dict__.update(obj.__dict__)

    def __repr__(self):
        name = self.__class__.__name__
        total_records = len(self)
        args = ", ".join(
            f"{k}={repr(v) if isinstance(v, str) else v}"
            for k, v in self._init_vals.items()
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
                table_name = "_".join(attr.split("_")[1:-1])
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

    def get_distinct(self, column: str) -> set[str]:
        sql = f"SELECT DISTINCT {column} from {self.table_name}"
        return {r[column] for r in self._execute_sql(sql).fetchall()}
