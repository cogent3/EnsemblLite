import typing

from cogent3 import make_seq
from cogent3.app.composable import define_app
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.core.sequence import Sequence

from ensembl_lite._config import InstalledConfig
from ensembl_lite._db_base import SqliteDbMixin
from ensembl_lite._homologydb import species_genes
from ensembl_lite.util import elt_compress_it, elt_decompress_it


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


class CompressedGenomeSeqsDb(GenomeSeqsDb):
    _genome_schema = {"coord_name": "TEXT PRIMARY KEY", "seq": "BLOB", "length": "INT"}

    def __hash__(self):
        return id(self)

    def add_record(self, *, coord_name: str, seq: str):
        sql = f"INSERT INTO {self.table_name}(coord_name, seq, length) VALUES (?, ?, ?)"
        self._execute_sql(sql, (coord_name, elt_compress_it(seq), len(seq)))
        self.db.commit()

    def add_records(self, *, records: typing.Iterable[list[str, str]]):
        self.add_compressed_records(
            records=[(n, elt_compress_it(s)) for n, s in records]
        )

    def add_compressed_records(self, *, records: typing.Iterable[list[str, bytes]]):
        """sequences already compressed"""

        sql = f"INSERT INTO {self.table_name}(coord_name, seq, length) VALUES (?, ?, ?)"

        self.db.executemany(sql, [(n, s, len(s)) for n, s in records])
        self.db.commit()

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
        sql = f"SELECT seq FROM {self.table_name} where coord_name = ?"

        seq = elt_decompress_it(self._execute_sql(sql, (coord_name,)).fetchone()[0])
        return seq[start:end] if start or end else seq


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

    def get_seq(self, *, seqid: str, start: OptInt = None, end: OptInt = None) -> str:
        """

        Parameters
        ----------
        seqid
            name of chromosome etc..
        start
            starting position of slice in python coordinates, defaults
            to 0
        end
            ending position of slice in python coordinates, defaults
            to length of coordinate
        """
        seq = self._seqs.get_seq(coord_name=seqid, start=start, end=end)
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
            try:
                seq = self.get_seq(seqid=seqid)
            except TypeError:
                msg = f"ERROR (report me): {self.species!r}, {seqid!r}"
                raise TypeError(msg)
            yield from seq.get_features(**kwargs)

    def close(self):
        self._seqs.close()
        self._annotdb.db.close()


def load_genome(*, cfg: InstalledConfig, species: str):
    """returns the Genome with bound seqs and features"""
    genome_path = cfg.installed_genome(species) / _SEQDB_NAME
    seqs = CompressedGenomeSeqsDb(source=genome_path, species=species)
    ann_path = cfg.installed_genome(species) / _ANNOTDB_NAME
    ann = GffAnnotationDb(source=ann_path)
    return Genome(species=species, seqs=seqs, annots=ann)


def get_seqs_for_ids(
    *,
    cfg: InstalledConfig,
    species: str,
    names: list[str],
    make_seq_name: typing.Callable = None,
) -> typing.Iterable[Sequence]:
    genome = load_genome(cfg=cfg, species=species)
    # is it possible to do batch query for all names?
    for name in names:
        feature = list(genome.get_features(name=f"%{name}"))[0]
        transcripts = list(feature.get_children(biotype="mRNA"))
        if not transcripts:
            continue

        longest = max(transcripts, key=lambda x: len(x))
        cds = list(longest.get_children(biotype="CDS"))
        if not cds:
            continue

        feature = cds[0]
        seq = feature.get_slice()
        if callable(make_seq_name):
            seq.name = make_seq_name(feature)
        else:
            seq.name = f"{species}-{name}"
        seq.info["species"] = species
        seq.info["name"] = name
        # disconnect from annotation so the closure of the genome
        # does not cause issues when run in parallel
        seq.annotation_db = None
        yield seq

    genome.close()
    del genome


@define_app
def get_selected_seqs(species_gene_ids: species_genes, config: InstalledConfig) -> list:
    """return gene sequences when given a species_gene_id instance

    Notes
    -----
    This function becomes a class, created using config. Calling the class
    instance with a species_genes instance is used to extract the list of gene
    ID's from the species.
    """
    species = species_gene_ids.species
    gene_ids = species_gene_ids.gene_ids
    return list(get_seqs_for_ids(cfg=config, species=species, names=gene_ids))
