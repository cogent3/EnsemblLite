import contextlib
import dataclasses
import functools
import io
import os
import pathlib

import duckdb
import numpy
import typing_extensions

from ensembl_tui import _util as eti_util

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
    return f"CREATE TABLE IF NOT EXISTS {table_name} ({columns_types})"


# HDF5 base class
@dataclasses.dataclass
class Hdf5Mixin(eti_util.SerialisableMixin):
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


class ViewMixin:
    _source: pathlib.Path  # override in subclass

    @property
    def species(self) -> str:
        return self._source.name


@dataclasses.dataclass(slots=True)
class DuckdbParquetBase:
    source: dataclasses.InitVar[pathlib.Path]
    # db is for testing purposes
    db: dataclasses.InitVar[duckdb.DuckDBPyConnection | None] = None
    _source: pathlib.Path = dataclasses.field(init=False)
    _conn: duckdb.DuckDBPyConnection = dataclasses.field(init=False, default=None)  # type: ignore
    _tables: tuple[str, ...] | tuple = ()

    def __post_init__(
        self,
        source: pathlib.Path,
        db: duckdb.DuckDBPyConnection | None,
    ) -> None:
        source = pathlib.Path(source)
        self._source = source
        if db:
            self._conn = db
            return

        self._conn = None

        if not source.is_dir():
            msg = f"{self._source} is not a directory"
            raise OSError(msg)

        if hasattr(self, "_post_init"):
            self._post_init()

    @property
    def conn(self) -> duckdb.DuckDBPyConnection:
        if self._conn is None:
            self._conn = duckdb.connect(":memory:")
            for table in self._tables:
                parquet_file = self._source / f"{table}.parquet"
                if not parquet_file.exists():
                    msg = f"{parquet_file} does not exist"
                    raise FileNotFoundError(msg)

                sql = f"CREATE TABLE {table} AS SELECT * FROM read_parquet('{parquet_file}')"
                self._conn.sql(sql)

        return self._conn

    def __len__(self) -> int:
        return self.num_records()

    def __eq__(self, other: typing_extensions.Self) -> bool:
        return other.conn is self.conn

    @property
    def source(self) -> pathlib.Path:
        return self._source

    def close(self) -> None:
        self.conn.close()

    def num_records(self) -> int:  # pragma: no cover
        # override in subclass
        raise NotImplementedError
