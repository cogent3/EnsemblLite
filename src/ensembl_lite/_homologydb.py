from __future__ import annotations

import typing

from typing import Iterable, Sized

from rich.progress import track

from ensembl_lite._db_base import SqliteDbMixin


class HomologyRecordType(typing.TypedDict):
    source: str
    species_1: str
    gene_id_1: str
    prot_id_1: str
    species_2: str
    gene_id_2: str
    prot_id_2: str
    relationship: str


def grouped_related(
    data: list[HomologyRecordType],
) -> typing.Iterable[tuple[tuple[str, str]]]:
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
    result = {}
    for record in track(data, description="Grouping related...", transient=True):
        pair = [
            (record["species_1"], record["gene_id_1"]),
            (record["species_2"], record["gene_id_2"]),
        ]

        for member in pair:
            if member in result:
                # one member, rel has already been encountered
                # so we get this value and update it with the new pair
                val = result[member]
                break
        else:
            val = set()

        val.update(pair)
        for member in pair:
            result[member] = val
    return [tuple(v) for v in result.values()]


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
        self, *, records: typing.Sequence[HomologyRecordType], col_order: Sized[str]
    ):
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
    ) -> Iterable[tuple[tuple[str, str]]]:
        """returns all groups of relationship type"""
        # get all gene ID's first
        sql = f"SELECT * from {self.table_name} WHERE relationship=?"
        results = [
            HomologyRecordType(iter(zip(r.keys(), r)))
            for r in self._execute_sql(sql, (relationship_type,)).fetchall()
        ]
        return grouped_related(results)

    def get_distinct(self, column: str) -> set[str]:
        sql = f"SELECT DISTINCT {column} from {self.table_name}"
        return {r[column] for r in self._execute_sql(sql).fetchall()}
