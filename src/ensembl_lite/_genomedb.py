import typing

from cogent3 import get_app, make_seq
from cogent3.app.composable import define_app
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.core.sequence import Sequence

from ensembl_lite._db_base import SqliteDbMixin


OptInt = typing.Optional[int]
OptionalStr = typing.Optional[str]

_SEQDB_NAME = "genome_sequence.seqdb"
_ANNOTDB_NAME = "features.gff3db"
# todo: make a variant that wraps a directory of compressed sequence files
# todo: or compresses sequence records on write into db, and just inflates then returns substring


class GenomeSeqsDb(SqliteDbMixin):
    table_name = "genome"
    _genome_schema = {"coord_name": "TEXT PRIMARY KEY", "seq": "TEXT", "length": "INT"}
    _metadata_schema = {"species": "TEXT"}

    def __init__(self, *, source: str = ":memory:", species: str = None):
        self.source = source
        self._init_tables()
        # the metadata table stores species info
        self._execute_sql("INSERT INTO metadata(species) VALUES (?)", (species,))
        self.db.commit()

    def __hash__(self):
        return id(self)

    def add_record(self, *, coord_name: str, seq: str):
        sql = f"INSERT INTO {self.table_name}(coord_name, seq, length) VALUES (?, ?, ?)"
        self._execute_sql(sql, (coord_name, seq, len(seq)))
        self.db.commit()

    def add_records(self, *, records: typing.Iterable[list[str, str]]):
        sql = f"INSERT INTO {self.table_name}(coord_name, seq, length) VALUES (?, ?, ?)"
        self.db.executemany(sql, [(n, s, len(s)) for n, s in records])
        self.db.commit()

    def get_seq(
        self, *, coord_name: str, start: OptInt = None, stop: OptInt = None
    ) -> str:
        """

        Parameters
        ----------
        coord_name
            name of chromosome etc..
        start
            starting position of slice in python coordinates, defaults
            to 0
        stop
            ending position of slice in python coordinates, defaults
            to length of coordinate
        """
        if start is not None:
            start += 1  # SQLite counts from 1
        else:
            start = 1

        if stop is None:
            sql = f"SELECT SUBSTR(seq, ?, length) FROM {self.table_name} where coord_name = ?"
            values = start, coord_name
        else:
            stop -= start - 1
            sql = (
                f"SELECT SUBSTR(seq, ?, ?) FROM {self.table_name} where coord_name = ?"
            )
            values = start, stop, coord_name

        return self._execute_sql(sql, values).fetchone()[0]


@define_app
def _str_to_bytes(data: str) -> bytes:
    """converts string to bytes"""
    return data.encode("utf8")


@define_app
def _bytes_to_str(data: bytes) -> str:
    """converts bytes into string"""
    return data.decode("utf8")


compress_it = _str_to_bytes() + get_app("compress")
decompress_it = get_app("decompress") + _bytes_to_str()


class CompressedGenomeSeqsDb(GenomeSeqsDb):
    _genome_schema = {"coord_name": "TEXT PRIMARY KEY", "seq": "BLOB", "length": "INT"}

    def __hash__(self):
        return id(self)

    def add_record(self, *, coord_name: str, seq: str):
        sql = f"INSERT INTO {self.table_name}(coord_name, seq, length) VALUES (?, ?, ?)"
        self._execute_sql(sql, (coord_name, compress_it(seq), len(seq)))
        self.db.commit()

    def add_records(self, *, records: typing.Iterable[list[str, str]]):
        sql = f"INSERT INTO {self.table_name}(coord_name, seq, length) VALUES (?, ?, ?)"
        self.add_compressed_records(records=[(n, compress_it(s)) for n, s in records])

    def add_compressed_records(self, *, records: typing.Iterable[list[str, bytes]]):
        """sequences already compressed"""
        sql = f"INSERT INTO {self.table_name}(coord_name, seq, length) VALUES (?, ?, ?)"
        self.db.executemany(sql, [(n, s, len(s)) for n, s in records])
        self.db.commit()

    def get_seq(
        self, *, coord_name: str, start: OptInt = None, stop: OptInt = None
    ) -> str:
        """

        Parameters
        ----------
        coord_name
            name of chromosome etc..
        start
            starting position of slice in python coordinates, defaults
            to 0
        stop
            ending position of slice in python coordinates, defaults
            to length of coordinate
        """
        sql = f"SELECT seq FROM {self.table_name} where coord_name = ?"

        seq = decompress_it(self._execute_sql(sql, (coord_name,)).fetchone()[0])
        return seq[start:stop] if start or stop else seq


# todo: this wrapping class is required for memory efficiency because
# the cogent3 SequeceCollection class is not designed for large sequence
# collections, either large sequences or large numbers of sequences. The
# correct solution is to improve that.
class Genome:
    """connecting sequences and their annotations"""

    def __init__(
        self,
        *,
        species: str,
        seqs: GenomeSeqsDb | CompressedGenomeSeqsDb,
        annots: GffAnnotationDb,
    ) -> None:
        self.species = species
        self._seqs = seqs
        self._annotdb = annots

    def get_seq(self, *, seqid: str, start: OptInt = None, stop: OptInt = None) -> str:
        """

        Parameters
        ----------
        seqid
            name of chromosome etc..
        start
            starting position of slice in python coordinates, defaults
            to 0
        stop
            ending position of slice in python coordinates, defaults
            to length of coordinate
        """
        seq = self._seqs.get_seq(coord_name=seqid, start=start, stop=stop)
        seq = make_seq(seq, name=seqid, moltype="dna")
        seq.annotation_offset = start or 0
        seq.annotation_db = self._annotdb
        return seq

    def get_features(
        self,
        *,
        biotype: str = None,
        seqid: str = None,
        name: str = None,
        start: int = None,
        stop: int = None,
    ) -> typing.Iterable[Feature]:
        """yields features in blocks of seqid"""
        kwargs = {k: v for k, v in locals().items() if k not in ("self", "seqid") and v}
        if seqid:
            seqids = [seqid]
        else:
            seqids = {
                ft["seqid"] for ft in self._annotdb.get_features_matching(**kwargs)
            }

        for seqid in seqids:
            seq = self.get_seq(seqid=seqid)
            yield from seq.get_features(**kwargs)


def get_feature_table(genome: Genome) -> Table:
    ...
