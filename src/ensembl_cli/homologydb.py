import typing

from ensembl_cli._db_base import SqliteDbMixin


class HomologyRecordType(typing.TypedDict):
    source: str
    species_1: str
    gene_id_1: str
    prot_id_1: str
    species_2: str
    gene_id_2: str
    prot_id_2: str
    relationship: str


# the homology db stores pairwise relationship information
class HomologyDb(SqliteDbMixin):
    table_name = "homology"

    _homology_schema = {
        "source": "TEXT",  # the file path
        "species_1": "TEXT",
        "gene_id_1": "TEXT",
        "prot_id_1": "TEXT",
        "species_2": "TEXT",
        "gene_id_2": "TEXT",
        "prot_id_2": "TEXT",
        "relationship": "TEXT",  # defined by Ensembl
    }

    def __init__(self, source=":memory:"):
        self.source = source
        self._init_tables()

    def add_records(
        self, records: typing.Sequence[HomologyRecordType], col_order: list[str]
    ):
        # bulk insert
        val_placeholder = ", ".join("?" * len(col_order))
        sql = f"INSERT INTO {self.table_name} ({', '.join(col_order)}) VALUES ({val_placeholder})"
        self.db.executemany(sql, records)
