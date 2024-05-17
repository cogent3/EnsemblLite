import dataclasses
import sqlite3

from ensembl_lite._util import SerialisableMixin


ReturnType = tuple[str, tuple]  # the sql statement and corresponding values


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


class SqliteDbMixin(SerialisableMixin):
    table_name = None
    _db = None
    source = None

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
        sql = f"CREATE INDEX IF NOT EXISTS %s on {self.table_name}(%s)"
        for col in self._index_columns:
            self._execute_sql(sql % (col, col))


# HDF5 base class
@dataclasses.dataclass
class Hdf5Mixin(SerialisableMixin):
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
        if self._is_open and self._file is not None:
            self._file.flush()
        if self._file is not None:
            self._file.close()
        self._is_open = False

    def close(self):
        if self._is_open:
            self._file.flush()
        self._file.close()
        self._is_open = False
