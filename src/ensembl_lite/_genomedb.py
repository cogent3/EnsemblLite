import typing

import click

from cogent3 import make_seq, make_table
from cogent3.app.composable import define_app
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.core.sequence import Sequence
from cogent3.util.table import Table

from ensembl_lite._config import InstalledConfig
from ensembl_lite._db_base import SqliteDbMixin
from ensembl_lite._homologydb import species_genes
from ensembl_lite.util import elt_compress_it, elt_decompress_it


_SEQDB_NAME = "genome_sequence.seqdb"
_ANNOTDB_NAME = "features.gff3db"


class GenomeSeqsDb(SqliteDbMixin):
    """class to be replaced by cogent3 sequence collection when that
    has been modernised"""

    table_name = "genome"
    _genome_schema = {"seqid": "TEXT PRIMARY KEY", "seq": "TEXT", "length": "INT"}
    _metadata_schema = {"species": "TEXT"}

    def __init__(self, *, source: str = ":memory:", species: str = None):
        self.source = source
        self._init_tables()
        # the metadata table stores species info
        self._execute_sql("INSERT INTO metadata(species) VALUES (?)", (species,))
        self.db.commit()

    def __hash__(self):
        return id(self)

    def add_record(self, *, seqid: str, seq: str):
        sql = f"INSERT INTO {self.table_name}(seqid, seq, length) VALUES (?, ?, ?)"
        self._execute_sql(sql, (seqid, seq, len(seq)))
        self.db.commit()

    def add_records(self, *, records: typing.Iterable[list[str, str]]):
        sql = f"INSERT INTO {self.table_name}(seqid, seq, length) VALUES (?, ?, ?)"
        self.db.executemany(sql, [(n, s, len(s)) for n, s in records])
        self.db.commit()

    def get_seq(
        self, *, seqid: str, start: int | None = None, end: int | None = None
    ) -> str:
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
        if start is not None:
            start += 1  # SQLite counts from 1
        else:
            start = 1

        if end is None:
            sql = (
                f"SELECT SUBSTR(seq, ?, length) FROM {self.table_name} where seqid = ?"
            )
            values = start, seqid
        else:
            end -= start - 1
            sql = f"SELECT SUBSTR(seq, ?, ?) FROM {self.table_name} where seqid = ?"
            values = start, end, seqid

        return self._execute_sql(sql, values).fetchone()[0]


class CompressedGenomeSeqsDb(GenomeSeqsDb):
    """class to be replaced by cogent3 sequence collection when that
    has been modernised"""

    _genome_schema = {"seqid": "TEXT PRIMARY KEY", "seq": "BLOB", "length": "INT"}

    def __hash__(self):
        return id(self)

    def add_record(self, *, seqid: str, seq: str):
        sql = f"INSERT INTO {self.table_name}(seqid, seq, length) VALUES (?, ?, ?)"
        self._execute_sql(sql, (seqid, elt_compress_it(seq), len(seq)))
        self.db.commit()

    def add_records(self, *, records: typing.Iterable[list[str, str]]):
        self.add_compressed_records(
            records=[(n, elt_compress_it(s)) for n, s in records]
        )

    def add_compressed_records(self, *, records: typing.Iterable[list[str, bytes]]):
        """sequences already compressed"""

        sql = f"INSERT INTO {self.table_name}(seqid, seq, length) VALUES (?, ?, ?)"

        self.db.executemany(sql, [(n, s, len(s)) for n, s in records])
        self.db.commit()

    def get_seq(
        self, *, seqid: str, start: int | None = None, end: int | None = None
    ) -> str:
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
        sql = f"SELECT seq FROM {self.table_name} where seqid = ?"

        seq = elt_decompress_it(self._execute_sql(sql, (seqid,)).fetchone()[0])
        return seq[start:end] if start or end else seq


# todo: this wrapping class is required for memory efficiency because
# the cogent3 SequeceCollection class is not designed for large sequence
# collections, either large sequences or large numbers of sequences. The
# correct solution is to improve that.
class Genome:
    """class to be replaced by cogent3 sequence collection when that
    has been modernised"""

    def __init__(
        self,
        *,
        species: str,
        seqs: GenomeSeqsDb | CompressedGenomeSeqsDb,
        annots: GffAnnotationDb,
    ) -> None:
        self.species = species
        self._seqs = seqs
        self.annotation_db = annots

    def get_seq(
        self,
        *,
        seqid: str,
        start: int | None = None,
        end: int | None = None,
        namer: typing.Callable | None = None,
    ) -> str:
        """returns annotated sequence

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
        namer
            callback for naming the sequence. Callback must take four
            arguments: species, seqid,start, end. Default is
            species:seqid:start-end.
        Notes
        -----
        Annotations partially within region are included.
        """
        seq = self._seqs.get_seq(seqid=seqid, start=start, end=end)
        if namer:
            name = namer(self.species, seqid, start, end)
        else:
            name = f"{self.species}:{seqid}:{start}-{end}"
        seq = make_seq(seq, name=name, moltype="dna")
        if self.annotation_db:
            seq.annotation_offset = start or 0
            seq.annotation_db = self.annotation_db.subset(
                seqid=seqid, start=start, end=end, allow_partial=True
            )
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
                ft["seqid"] for ft in self.annotation_db.get_features_matching(**kwargs)
            }
        for seqid in seqids:
            try:
                seq = self.get_seq(seqid=seqid)
            except TypeError:
                msg = f"ERROR (report me): {self.species!r}, {seqid!r}"
                raise TypeError(msg)
            # because self.get_seq() automatically names seqs differently
            seq.name = seqid
            yield from seq.get_features(**kwargs)

    def close(self):
        self._seqs.close()
        self.annotation_db.db.close()


def load_genome(*, config: InstalledConfig, species: str):
    """returns the Genome with bound seqs and features"""
    genome_path = config.installed_genome(species) / _SEQDB_NAME
    seqs = CompressedGenomeSeqsDb(source=genome_path, species=species)
    ann_path = config.installed_genome(species) / _ANNOTDB_NAME
    ann = GffAnnotationDb(source=ann_path)
    return Genome(species=species, seqs=seqs, annots=ann)


def get_seqs_for_ids(
    *,
    config: InstalledConfig,
    species: str,
    names: list[str],
    make_seq_name: typing.Callable = None,
) -> typing.Iterable[Sequence]:
    genome = load_genome(config=config, species=species)
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
    return list(get_seqs_for_ids(config=config, species=species, names=gene_ids))


def get_annotations_for_species(
    *, config: InstalledConfig, species: str
) -> GffAnnotationDb:
    """returns the annotation Db for species"""
    path = config.installed_genome(species=species)
    if not path.exists():
        click.secho(f"{species!r} not in {str(config.install_path.parent)!r}", fg="red")
        exit(1)
    # TODO: this filename should be defined in one place
    path = path / "features.gff3db"
    if not path.exists():
        click.secho(f"{path.name!r} is missing", fg="red")
        exit(1)
    return GffAnnotationDb(source=path)


def get_gene_table_for_species(
    *, annot_db: GffAnnotationDb, limit: int | None, species: str | None = None
) -> Table:
    """
    returns gene data from a GffDb

    Parameters
    ----------
    annot_db
        feature db
    limit
        limit number of records to
    species
        species name, overrides inference from annot_db.source
    """
    species = species or annot_db.source.parent.name

    columns = (
        "species",
        "name",
        "seqid",
        "source",
        "biotype",
        "start",
        "end",
        "score",
        "strand",
        "phase",
    )
    rows = []
    for i, record in enumerate(annot_db.get_records_matching(biotype="gene")):
        rows.append([species] + [record.get(c, None) for c in columns[1:]])
        if i == limit:
            break

    return make_table(header=columns, data=rows)


def get_species_summary(
    *, annot_db: GffAnnotationDb, species: str | None = None
) -> Table:
    """
    returns the Table summarising data for species_name

    Parameters
    ----------
    annot_db
        feature db
    species
        species name, overrides inference from annot_db.source
    """
    from .species import Species

    # for now, just biotype
    species = species or annot_db.source.parent.name
    counts = annot_db.biotype_counts()
    try:
        common_name = Species.get_common_name(species)
    except ValueError:
        common_name = species

    return Table(
        header=("biotype", "count"),
        data=list(counts.items()),
        title=f"{common_name} features",
    )
