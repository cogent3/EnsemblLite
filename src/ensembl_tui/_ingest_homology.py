# consumes homology data from Ensembl tsv files and writes a parquet file
# for representation by _homology.HomologyDb
import pathlib
import shutil
import typing

import duckdb
from cogent3.app import typing as c3_types
from cogent3.app.composable import LOADER, define_app

from ensembl_tui import _config as eti_config
from ensembl_tui import _ingest_annotation as eti_annotation
from ensembl_tui import _util as eti_util
from ensembl_tui._homology import homolog_group

T = dict[str, tuple[homolog_group, ...]]


def grouped_related(
    data: typing.Iterable[tuple[str, dict[str, str]]],
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
    for rel_type, gene_species in data:
        relationship = grouped.get(rel_type, {})
        pair = gene_species.keys()
        gene_id_1, gene_id_2 = pair
        if gene_id_1 in relationship:
            val = relationship[gene_id_1]
        elif gene_id_2 in relationship:
            val = relationship[gene_id_2]
        else:
            val = homolog_group(relationship=rel_type)
        val.gene_ids |= gene_species

        relationship[gene_id_1] = relationship[gene_id_2] = val
        grouped[rel_type] = relationship

    return {
        rel_type: tuple(set(groups.values())) for rel_type, groups in grouped.items()
    }


@define_app(app_type=LOADER)
class load_homologies:  # noqa: N801
    """app to load homology groups from a single Ensembl tsv file"""

    def __init__(self, allowed_species: set[str]) -> None:
        # map the Ensembl columns to HomologyDb columns
        self._src_cols = (
            "homology_type",
            "species",
            "gene_stable_id",
            "homology_species",
            "homology_gene_stable_id",
        )
        self._create_sql = (
            f"CREATE TABLE my_table AS SELECT {','.join(self._src_cols)} FROM"
            " read_csv_auto('{}', delim='\t', header=True)"
        )
        allowed = ", ".join(f"{sp!r}" for sp in allowed_species)
        self._select_sql = (
            f"SELECT {','.join(self._src_cols)} FROM my_table"
            f" WHERE species IN ({allowed}) AND homology_species IN ({allowed})"
        )

    def main(self, path: c3_types.IdentifierType) -> c3_types.SerialisableType:
        conn = duckdb.connect(":memory:")
        sql = self._create_sql.format(path)
        conn.sql(sql)
        return grouped_related(
            (rel_type, {gene_1: sp_1, gene_2: sp_2})
            for rel_type, sp_1, gene_1, sp_2, gene_2 in conn.sql(
                self._select_sql,
            ).fetchall()
        )  # type: ignore


class HomologyAggregator:
    def __init__(self, conn: duckdb.DuckDBPyConnection) -> None:
        self._relationship_indexer = eti_util.unique_value_indexer()
        self._species_indexer = eti_util.unique_value_indexer()
        self._stableids_indexer = eti_util.unique_value_indexer()
        self._homology_indexer = eti_util.category_indexer()
        self.conn = conn

    def _make_stableid_id(self, *, stableid: str, species: str) -> int:
        """returns the stableid.id value for (species,stableid)"""
        species_id = self._species_indexer(species)
        return self._stableids_indexer((species_id, stableid))

    def add_records(
        self,
        *,
        records: typing.Sequence[homolog_group],
        relationship_type: str,
    ) -> None:
        """inserts homology data from records

        Parameters
        ----------
        records
            a sequence of homolog group instances, all with the same
            relationship type
        relationship_type
            the relationship type
        """
        if not relationship_type:
            msg = f"invalid {relationship_type=!r}"
            raise ValueError(msg)

        relationship_id = self._relationship_indexer(relationship_type)
        homology_values = []
        values = []
        for group in records:
            if group.relationship != relationship_type:
                msg = f"{group.relationship=} != {relationship_type=}"
                raise ValueError(msg)

            # get geneids and species for this group, storing the
            # geneid id for each record
            gene_ids = {
                self._make_stableid_id(stableid=gene_id, species=species)
                for gene_id, species in group.gene_ids.items()
            }

            # now get the homology id for this group
            homology_id = self._homology_indexer(relationship_id, gene_ids)
            homology_values.append((homology_id, relationship_id))
            values.extend([(int(gene_id), int(homology_id)) for gene_id in gene_ids])

        # create the homology table entries
        sql = (
            "INSERT OR IGNORE INTO homology(homology_id, relationship_id) VALUES (?, ?)"
        )
        self.conn.executemany(sql, parameters=homology_values)
        sql = "INSERT OR IGNORE INTO member(stableid_id, homology_id) VALUES (?, ?)"
        self.conn.executemany(sql, parameters=values)

    def commit(self):
        # add the indexes for species and relationships
        relationships = list(self._relationship_indexer)
        sql = "INSERT OR IGNORE INTO relationship(relationship_id, homology_type) VALUES (?,?)"
        self.conn.executemany(sql, parameters=relationships)
        species = list(self._species_indexer)
        sql = "INSERT OR IGNORE INTO species(species_id, species_db) VALUES (?,?)"
        self.conn.executemany(sql, parameters=species)
        # now the same for stableids
        stableids = [(index, *attr) for index, attr in self._stableids_indexer]
        sql = "INSERT OR IGNORE INTO stableid(stableid_id, species_id, stableid) VALUES (?,?,?)"
        self.conn.executemany(sql, parameters=stableids)
        self.conn.commit()


def make_homology_aggregator_db() -> HomologyAggregator:
    defns = {
        "relationship": {
            "homology_type": "TEXT",
            "relationship_id": "INTEGER PRIMARY KEY",
        },
        "homology": {
            "homology_id": "INTEGER PRIMARY KEY DEFAULT nextval('homology_id_seq')",
            "relationship_id": "INTEGER",
        },
        "species": {
            "species_id": "INTEGER PRIMARY KEY",
            "species_db": "TEXT",
        },
        "stableid": {
            "stableid_id": "INTEGER PRIMARY KEY",
            "stableid": "TEXT",
            "species_id": "INTEGER",
        },
        "member": {  # gene membership of a specific homolog group
            "member_id": "INTEGER PRIMARY KEY DEFAULT nextval('member_id_seq')",
            "stableid_id": "INTEGER",
            "homology_id": "INTEGER",
        },
    }

    conn = duckdb.connect(":memory:")
    conn.sql("CREATE SEQUENCE homology_id_seq")
    conn.sql("CREATE SEQUENCE member_id_seq")

    for table_name, schema in defns.items():
        columns = ", ".join(f"{k} {v}" for k, v in schema.items())
        conn.sql(f"CREATE TABLE {table_name} ({columns})")

    sql = """CREATE VIEW IF NOT EXISTS homology_member AS
        SELECT h.homology_id as homology_id,
               h.relationship_id as relationship_id,
               r.homology_type as homology_type,
               m.stableid_id as stableid_id
        FROM homology h
        JOIN relationship r ON h.relationship_id = r.relationship_id
        JOIN member m ON m.homology_id = h.homology_id
        """
    conn.sql(sql)
    # this view relates stabelids to species db names
    # since that cannot be reliably inferred from the stableid
    # naming scheme. Possibly just need the species_db and stableid columns?
    sql = """CREATE VIEW IF NOT EXISTS gene_species_attr AS
        SELECT sp.species_db as species_db,
               sp.species_id as species_id,
               st.stableid as stableid,
               st.stableid_id as stableid_id
        FROM species sp
        JOIN stableid st ON st.species_id = sp.species_id
        """
    conn.sql(sql)
    # this view condenses the homology member data so that we can
    # query for a relationship type and get the stableid and species
    # bear in mind that for a given relationship type, I'm assuming
    # that a stableid can only belong to one homology group (homology_id)
    sql = """
        CREATE VIEW IF NOT EXISTS homology_groups_attr AS
        SELECT gs.stableid as stableid,
               gs.species_db as species_db,
               hm.homology_id as homology_id,
               hm.homology_type as homology_type
        FROM gene_species_attr gs
        JOIN homology_member hm ON hm.stableid_id = gs.stableid_id
        """
    conn.sql(sql)
    return HomologyAggregator(conn)


def write_homology_views(agg: HomologyAggregator, outdir: pathlib.Path) -> None:
    # we write out the gene_species_attr and homology_groups_attr views
    # to the output directory
    outdir.mkdir(parents=True, exist_ok=True)
    _ = eti_annotation.export_parquet(
        con=agg.conn,
        table_name="homology_groups_attr",
        dest_dir=outdir,
    )


def local_install_homology(
    config: eti_config.Config,
    force_overwrite: bool = False,
    max_workers: int = 1,
    verbose: bool = False,
):
    from rich.progress import track

    if force_overwrite:
        shutil.rmtree(config.install_homologies, ignore_errors=True)

    # config.install_homologies.mkdir(parents=True, exist_ok=True)

    tsv_paths = []
    for sp in config.db_names:
        path = config.staging_homologies / sp
        tsv_paths.extend(list(path.glob("*.tsv*")))

    if verbose:
        eti_util.print_colour(f"homologies {max_workers=}", "yellow")

    loader = load_homologies(allowed_species=set(config.db_names))

    tasks = eti_util.get_iterable_tasks(
        func=loader,
        series=tsv_paths,
        max_workers=max_workers,
    )
    agg = make_homology_aggregator_db()
    for grouped in track(tasks, total=len(tsv_paths), description="aggregating"):
        for rel_type, records in grouped.items():
            agg.add_records(records=records, relationship_type=rel_type)

    agg.commit()
