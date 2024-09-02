import contextlib
import dataclasses
import functools
import io
import os
import sqlite3

import numpy

from ensembl_tui import _util as elt_util

ReturnType = tuple[str, tuple]  # the sql statement and corresponding values


@functools.singledispatch
def array_to_blob(data: numpy.ndarray) -> bytes:
    with io.BytesIO() as out:
        numpy.save(out, data)
        out.seek(0)
        output = out.read()
    return output


@array_to_blob.register
def _(data: bytes) -> bytes:
    # already a blob
    return data


@functools.singledispatch
def blob_to_array(data: bytes) -> numpy.ndarray:
    with io.BytesIO(data) as out:
        out.seek(0)
        result = numpy.load(out)
    return result


@blob_to_array.register
def _(data: numpy.ndarray) -> numpy.ndarray:
    return data


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
    primary_key = columns.pop("PRIMARY KEY", None)
    columns_types = ", ".join([f"{name} {ctype}" for name, ctype in columns.items()])
    if primary_key:
        columns_types = f"{columns_types}, PRIMARY KEY ({','.join(primary_key)})"
    sql = f"CREATE TABLE IF NOT EXISTS {table_name} ({columns_types})"
    return sql


class SqliteDbMixin(elt_util.SerialisableMixin):
    table_name = None
    _db = None
    source: elt_util.PathType = ":memory:"

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
        self._db = self._db or sqlite3.connect(
            self.source,
            detect_types=sqlite3.PARSE_DECLTYPES,
            check_same_thread=False,
        )
        self._db.row_factory = sqlite3.Row

        # try and reduce memory usage
        cursor = self._db.cursor()
        cursor.execute("PRAGMA cache_size = -2048;")

        # A bit of magic.
        # Assumes schema attributes named as `_<table name>_schema`
        for attr in dir(self):
            if attr.endswith("_schema"):
                table_name = "_".join(attr.split("_")[1:-1])
                attr = getattr(self, attr)
                sql = _make_table_sql(table_name, attr)
                cursor.execute(sql)

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

    def make_indexes(self):
        """adds db indexes for core attributes"""
        sql = "CREATE INDEX IF NOT EXISTS %(index)s on %(table)s(%(col)s)"
        for table_name, columns in self._index_columns.items():
            for col in columns:
                index = f"{col}_index"
                self._execute_sql(sql % dict(table=table_name, index=index, col=col))


# HDF5 base class
@dataclasses.dataclass
class Hdf5Mixin(elt_util.SerialisableMixin):
    """HDF5 sequence data storage"""

    _file = None
    _is_open = False

    def __getstate__(self):
        if set(self.mode) & {"w", "a"}:
            raise NotImplementedError(f"pickling not supported for mode={self.mode!r}")
        return self._init_vals.copy()

    def __setstate__(self, state):
        obj = self.__class__(**state)
        self.__dict__.update(obj.__dict__)
        # because we have a __del__ method, and self attributes point to
        # attributes on obj, we need to modify obj state so that garbage
        # collection does not screw up self
        obj._is_open = False
        obj._file = None

    def __del__(self):
        self.close()

    def close(self):
        """closes the hdf5 file"""
        # hdf5 dumps content to stdout if resource already closed, so
        # we trap that here, and capture expected exceptions raised in the process
        with open(os.devnull, "w") as devnull:
            with (
                contextlib.redirect_stderr(devnull),
                contextlib.redirect_stdout(devnull),
            ):
                with contextlib.suppress(ValueError, AttributeError):
                    if self._is_open:
                        self._file.flush()

                with contextlib.suppress(AttributeError):
                    self._file.close()

        self._is_open = False
