from __future__ import annotations

import dataclasses
import typing

from cogent3.parse.table import FilteringParser
from cogent3.util.io import iter_splitlines
from rich.progress import track

from ensembl_lite._config import InstalledConfig
from ensembl_lite._db_base import SqliteDbMixin
from ensembl_lite._util import PathType


_HOMOLOGYDB_NAME = "homologies.sqlitedb"


@dataclasses.dataclass(slots=True, eq=True)
class HomologyRecord:
    source: str | None = None
    gene_id_1: str | None = None
    gene_id_2: str | None = None
    relationship: str | None = None

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, item, value):
        setattr(self, item, value)


@dataclasses.dataclass(slots=True)
class species_genes:
    """contains gene IDs for species"""

    species: str
    gene_ids: list[str] = None

    def __hash__(self):
        return hash(self.species)

    def __eq__(self, other):
        return self.species == other.species and self.gene_ids == other.gene_ids

    def __post_init__(self):
        if self.gene_ids is None:
            self.gene_ids = []


@dataclasses.dataclass
class homolog_group:
    """has species_genes instances belonging to the same ortholog group"""

    relationship: str
    gene_ids: typing.Optional[set[str, ...]] = None

    def __post_init__(self):
        self.gene_ids = self.gene_ids if self.gene_ids else set()

    def __hash__(self):
        # allow hashing, but bearing in mind we are updating
        # gene values
        return hash((hash(self.relationship), id(self.gene_ids)))

    def __eq__(self, other):
        return (
            self.relationship == other.relationship and self.gene_ids == other.gene_ids
        )

    def __getstate__(self) -> tuple[str, set[str]]:
        return self.relationship, self.gene_ids

    def __setstate__(self, state):
        relationship, gene_ids = state
        self.relationship = relationship
        self.gene_ids = gene_ids

    def __len__(self):
        return len(self.gene_ids)

    def __or__(self, other):
        if other.relationship != self.relationship:
            raise ValueError(
                f"relationship type {self.relationship!r} != {other.relationship!r}"
            )
        return self.__class__(
            relationship=self.relationship, gene_ids=self.gene_ids | other.gene_ids
        )


T = dict[str, tuple[homolog_group, ...]]


def grouped_related(
    data: list[HomologyRecord],
) -> T:
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
    # grouped is {<relationship type>: {gene id: homolog_group}. So gene's
    # that belong to the same group have the same value
    grouped = {}
    for record in data:
        rel_type = record.relationship
        relationship = grouped.get(rel_type, {})
        pair = {record.gene_id_1, record.gene_id_2}

        if record.gene_id_1 in relationship:
            val = relationship[record.gene_id_1]
        elif record.gene_id_2 in relationship:
            val = relationship[record.gene_id_2]
        else:
            val = homolog_group(relationship=record.relationship)
        val.gene_ids |= pair

        relationship[record.gene_id_1] = relationship[record.gene_id_2] = val
        grouped[rel_type] = relationship

    reduced = {}
    for rel_type, groups in grouped.items():
        reduced[rel_type] = tuple(set(groups.values()))

    return reduced


def _gene_id_to_group(series: tuple[homolog_group, ...]) -> dict[str:homolog_group]:
    """converts series of homolog_group instances to {geneid: groupl, ..}"""
    result = {}
    for group in series:
        result.update({gene_id: group for gene_id in group.gene_ids})
    return result


def _add_unique(
    a: dict[str, homolog_group],
    b: dict[str, homolog_group],
    combined: dict[str, homolog_group],
) -> dict[str, homolog_group]:
    unique = a.keys() - b.keys()
    combined.update(**{gene_id: a[gene_id] for gene_id in unique})
    return combined


def merge_grouped(group1: T, group2: T) -> T:
    """merges homolog_group with overlapping members"""
    joint = {}
    groups = group1, group2
    rel_types = group1.keys() | group2.keys()
    for rel_type in rel_types:
        if any(rel_type not in grp for grp in groups):
            joint[rel_type] = group1.get(rel_type, group2.get(rel_type))
            continue

        # expand values to dicts
        grp1 = _gene_id_to_group(group1[rel_type])
        grp2 = _gene_id_to_group(group2[rel_type])

        # if a group is unique for a relationship type, not one member
        # will be present in the other group
        # add groups that are truly unique to each
        rel_type_group = {}
        # unique to grp 1
        rel_type_group = _add_unique(grp1, grp2, rel_type_group)
        # unique to grp 2
        rel_type_group = _add_unique(grp2, grp1, rel_type_group)

        shared_ids = grp1.keys() & grp2.keys()
        skip = set()  # id's for groups already processed
        for gene_id in shared_ids:
            if gene_id in skip:
                continue
            merged = grp1[gene_id] | grp2[gene_id]
            rel_type_group.update({gene_id: merged for gene_id in merged.gene_ids})
            skip.update(merged.gene_ids)

        joint[rel_type] = tuple(set(rel_type_group.values()))

    return joint


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
    _index_columns = "relationship", "gene_id_1", "gene_id_2"

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
    ) -> typing.Sequence[homolog_group]:
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


class LoadHomologies:
    def __init__(self, allowed_species: set):
        self._allowed_species = allowed_species
        # map the Ensembl columns to HomologyDb columns

        self.src_cols = (
            "homology_type",
            "species",
            "gene_stable_id",
            "homology_species",
            "homology_gene_stable_id",
        )
        self.dest_col = (
            "relationship",
            "species_1",
            "gene_id_1",
            "species_2",
            "gene_id_2",
        )
        self._reader = FilteringParser(
            row_condition=self._matching_species, columns=self.src_cols, sep="\t"
        )

    def _matching_species(self, row):
        return {row[1], row[3]} <= self._allowed_species

    def __call__(self, paths: typing.Iterable[PathType]) -> list[HomologyRecord]:
        final = []
        for path in paths:
            rows = list(self._reader(iter_splitlines(path)))
            header = rows.pop(0)
            assert list(header) == list(self.src_cols), (header, self.src_cols)
            # we select only the relationship type and gene IDs
            final.extend(
                [
                    HomologyRecord(
                        relationship=row[0],
                        gene_id_1=row[2],
                        gene_id_2=row[4],
                        source=path,
                    )
                    for row in rows
                ]
            )

        return final
