from __future__ import annotations

import dataclasses
import typing

from rich.progress import track

from ensembl_lite._config import InstalledConfig
from ensembl_lite._db_base import SqliteDbMixin


_HOMOLOGYDB_NAME = "homologies.sqlitedb"


@dataclasses.dataclass(slots=True, eq=True)
class HomologyRecord:
    source: str | None = None
    species_1: str | None = None
    gene_id_1: str | None = None
    prot_id_1: str | None = None
    species_2: str | None = None
    gene_id_2: str | None = None
    prot_id_2: str | None = None
    relationship: str | None = None

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, item, value):
        setattr(self, item, value)


def grouped_related(
    data: list[HomologyRecord],
) -> set[frozenset[tuple[str, str]]]:
    """determines related groups of genes

    Parameters
    ----------
    data
        list of full records from the HomologyDb

    Returns
    -------
    a data structure that can be json serialised

    Notes
    -----
    I assume that for a specific relationship type, a gene can only belong
    to one group.
    """
    # grouped is {gene id: set(group)}. So gene's that belong to the same
    # group have the same value
    grouped = {}
    for record in track(data, description="Grouping related...", transient=True):
        pair = [
            (record.species_1, record.gene_id_1),
            (record.species_2, record.gene_id_2),
        ]
        if record.gene_id_1 in grouped:
            val = grouped[record.gene_id_1]
        elif record.gene_id_2 in grouped:
            val = grouped[record.gene_id_2]
        else:
            val = set()
        val.update(pair)
        grouped[record.gene_id_1] = grouped[record.gene_id_2] = val
    return {frozenset(v) for v in grouped.values()}


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
        self,
        *,
        records: typing.Sequence,
        col_order: typing.Sized[str],
    ) -> None:
        # bulk insert
        val_placeholder = ", ".join("?" * len(col_order))
        sql = f"INSERT INTO {self.table_name} ({', '.join(col_order)}) VALUES ({val_placeholder})"
        self.db.executemany(sql, records)
        self.db.commit()

    def get_related_to(
        self, *, gene_id: str, relationship_type: str
    ) -> set[tuple[str, str]]:
        """return genes with relationship type to gene_id"""
        sql = (
            f"SELECT species_1, gene_id_1, species_2, gene_id_2 from {self.table_name}"
            f" WHERE relationship = ? AND (gene_id_1=? OR gene_id_2=?)"
        )
        result = set()
        for r in self._execute_sql(
            sql, (relationship_type, gene_id, gene_id)
        ).fetchall():
            for num in (1, 2):
                species_key = f"species_{num}"
                gene_id_key = f"gene_id_{num}"
                result.add((r[species_key], r[gene_id_key]))
        return result

    def get_related_groups(
        self, relationship_type: str
    ) -> set[frozenset[tuple[str, str]]]:
        """returns all groups of relationship type"""
        # get all gene ID's first
        sql = f"SELECT * from {self.table_name} WHERE relationship=?"
        results = [
            HomologyRecord(**dict(zip(r.keys(), r)))
            for r in self._execute_sql(sql, (relationship_type,)).fetchall()
        ]
        return grouped_related(results)


def load_homology_db(
    *,
    config: InstalledConfig,
) -> HomologyDb:
    return HomologyDb(source=config.homologies_path / _HOMOLOGYDB_NAME)


@dataclasses.dataclass(slots=True)
class species_genes:
    """contains gene IDs for species"""

    species: str
    gene_ids: list[str] = None

    def __hash__(self):
        return hash(self.species)

    def __post_init__(self):
        self.gene_ids = []


def id_by_species_group(related) -> tuple[list[species_genes], dict[str, int]]:
    """returns species gene sets and relationship index"""
    sp_groups = {}
    id_group_map = {}
    for group_num, group in enumerate(related):
        for sp, gene_id in group:
            val = sp_groups[sp] if sp in sp_groups else species_genes(species=sp)
            val.gene_ids.append(gene_id)
            sp_groups[sp] = val
            id_group_map[gene_id] = group_num
    return list(sp_groups.values()), id_group_map
