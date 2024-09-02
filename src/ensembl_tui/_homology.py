import dataclasses
import typing

import blosc2
from cogent3 import make_unaligned_seqs
from cogent3.app.composable import LOADER, NotCompleted, define_app
from cogent3.app.io import compress, decompress, pickle_it, unpickle_it
from cogent3.app.typing import (
    IdentifierType,
    SeqsCollectionType,
    SerialisableType,
)
from cogent3.parse.table import FilteringParser
from cogent3.util.io import PathType, iter_splitlines

from ensembl_tui import _config as elt_config
from ensembl_tui import _genome as elt_genome
from ensembl_tui import _storage_mixin as elt_mixin

HOMOLOGY_STORE_NAME = "homologies.homology-sqlitedb"

compressor = compress(compressor=blosc2.compress2)
decompressor = decompress(decompressor=blosc2.decompress2)
pickler = pickle_it()
unpickler = unpickle_it()
inflate = decompressor + unpickler


@dataclasses.dataclass(slots=True)
class species_genes:
    """contains gene IDs for species"""

    species: str
    gene_ids: list[str] | None = None

    def __hash__(self):
        return hash(self.species)

    def __eq__(self, other):
        return self.species == other.species and self.gene_ids == other.gene_ids

    def __post_init__(self):
        self.gene_ids = list(self.gene_ids) if self.gene_ids else []

    def __getstate__(self) -> tuple[str, list[str]]:
        return self.species, self.gene_ids

    def __setstate__(self, args):
        species, gene_ids = args
        self.species = species
        self.gene_ids = gene_ids


@dataclasses.dataclass
class homolog_group:
    """has species_genes instances belonging to the same ortholog group"""

    relationship: str
    # gene id -> species
    gene_ids: dict[str, str] | None = None
    source: str | None = None

    def __post_init__(self):
        self.gene_ids = self.gene_ids or {}
        if self.source is None:
            self.source = next(iter(self.gene_ids), None)

    def __hash__(self):
        # allow hashing, but bearing in mind we are updating
        # gene values
        return hash((hash(self.relationship), id(self.gene_ids)))

    def __eq__(self, other):
        return (
            self.relationship == other.relationship and self.gene_ids == other.gene_ids
        )

    def __getstate__(self) -> tuple[str, dict[str, str] | None, str | None]:
        return self.relationship, self.gene_ids, self.source

    def __setstate__(self, state: tuple[str, dict[str, str] | None, str | None]):
        relationship, gene_ids, source = state
        self.relationship = relationship
        self.gene_ids = gene_ids
        self.source = source

    def __len__(self):
        return len(self.gene_ids or ())

    def __or__(self, other):
        if other.relationship != self.relationship:
            raise ValueError(
                f"relationship type {self.relationship!r} != {other.relationship!r}",
            )
        gene_ids = {**(self.gene_ids or {}), **(other.gene_ids or {})}
        return self.__class__(relationship=self.relationship, gene_ids=gene_ids)

    def species_ids(self) -> dict[str, tuple[str, ...]]:
        """returns {species: gene_ids, ...}"""
        result = {}
        gene_ids = self.gene_ids or {}
        for gene_id, sp in gene_ids.items():
            ids = result.get(sp, [])
            ids.append(gene_id)
            result[sp] = ids
        return result


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

    reduced = {}
    for rel_type, groups in grouped.items():
        reduced[rel_type] = tuple(set(groups.values()))

    return reduced


def _gene_id_to_group(series: tuple[homolog_group, ...]) -> dict[str, homolog_group]:
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
class HomologyDb(elt_mixin.SqliteDbMixin):
    table_names = "homology", "relationship", "member", "species", "stableid"

    _relationship_schema = {  # e.g. ortholog_one2one
        "homology_type": "TEXT",
    }

    _homology_schema = {  # e.g. an individual homolog group of homology_type
        "relationship_id": "INTEGER",
    }

    _species_schema = {"species_db": "TEXT"}

    _stableid_schema = {
        "stableid": "TEXT PRIMARY KEY",
        "species_id": "INTEGER",
    }

    _member_schema = {  # gene membership of a specific homolog group
        "stableid_id": "INTEGER",  # from stableid table
        "homology_id": "INTEGER",  # from homology table
        "PRIMARY KEY": ("stableid_id", "homology_id"),
    }

    _index_columns = {
        "homology": ("relationship_id",),
        "relationship": ("homology_type",),
        "member": ("stableid_id", "homology_id"),
        "species": ("species_db",),
        "stableid": ("stableid", "species_id"),
    }

    def __init__(self, source: PathType = ":memory:"):
        self.source = source
        self._relationship_types = {}
        self._species_ids = {}
        self._init_tables()
        self._create_views()

    def _create_views(self):
        """define views to simplify queries"""
        sql = """CREATE VIEW IF NOT EXISTS homology_member AS
        SELECT h.rowid as homology_id,
               h.relationship_id as relationship_id,
               r.homology_type as homology_type,
               m.stableid_id as stableid_id
        FROM homology h
        JOIN relationship r ON h.relationship_id = r.rowid
        JOIN member m ON m.homology_id = h.rowid
        """
        self._execute_sql(sql)
        sql = """
        CREATE VIEW IF NOT EXISTS gene_species AS
        SELECT sp.species_db as species_db,
               sp.rowid as species_id,
               st.stableid as stableid,
               st.rowid as stableid_id
        FROM species sp
        JOIN stableid st ON st.species_id = sp.rowid
        """
        self._execute_sql(sql)
        sql = """
        CREATE VIEW IF NOT EXISTS related_groups AS
        SELECT gs.stableid as stableid,
               gs.species_db as species_db,
               hm.homology_id as homology_id,
               hm.homology_type as homology_type
        FROM gene_species gs
        JOIN homology_member hm ON hm.stableid_id = gs.stableid_id
        """
        self._execute_sql(sql)

    def _make_species_id(self, species: str) -> int:
        """returns the species.id value for species"""
        if species not in self._species_ids:
            sql = "INSERT INTO species(species_db) VALUES (?) RETURNING rowid"
            result = self.db.execute(sql, (species,)).fetchone()[0]
            self._species_ids[species] = result
        return self._species_ids[species]

    def _make_stableid_id(self, *, stableid: str, species: str) -> int:
        """returns the stableid.id value for gene_id"""
        species_id = self._make_species_id(species)
        sql = """
        INSERT OR IGNORE INTO stableid(stableid,species_id) VALUES (?,?) RETURNING rowid
        """
        r = self.db.execute(sql, (stableid, species_id)).fetchone()
        if r is None:
            sql = "SELECT rowid FROM stableid WHERE stableid = ? AND species_id = ?"
            r = self.db.execute(sql, (stableid, species_id)).fetchone()
        return r[0]

    def _make_relationship_type_id(self, rel_type: str) -> int:
        """returns the relationship.id value for relationship_type"""
        if rel_type not in self._relationship_types:
            sql = "INSERT INTO relationship(homology_type) VALUES (?) RETURNING rowid"
            result = self.db.execute(sql, (rel_type,)).fetchone()[0]
            self._relationship_types[rel_type] = result
        return self._relationship_types[rel_type]

    def _get_homology_group_id(
        self,
        *,
        relationship_id: int,
        gene_ids: tuple[str, ...],
    ) -> int:
        """creates a new homolog table entry for this relationship id / group"""
        # need a join of homology by relationship
        # I want to get the homology_id from member table where
        # member.stableid_id corresponds to gene_ids
        id_placeholders = ",".join("?" * len(gene_ids))
        sql = f"""
        SELECT hm.homology_id
        FROM homology_member hm
        WHERE hm.relationship_id = ? AND hm.stableid_id IN ({id_placeholders})
        """
        result = self.db.execute(sql, (relationship_id,) + gene_ids).fetchone()
        if result:
            return result[0]
        # this group not seen before, so we just create a homology entry
        sql = "INSERT INTO homology(relationship_id) VALUES (?) RETURNING rowid"
        result = self.db.execute(sql, (relationship_id,)).fetchone()

        return result[0]

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
        assert relationship_type is not None
        rel_type_id = self._make_relationship_type_id(relationship_type)
        # we now iterate over the homology groups
        # we get a new homology id, then add all genes for that group
        # using the IGNORE to skip duplicates
        sql = "INSERT OR IGNORE INTO member(stableid_id,homology_id) VALUES (?, ?)"
        for group in records:
            if group.relationship != relationship_type:
                raise ValueError(f"{group.relationship=} != {relationship_type=}")

            # get geneids and species for this group, storing the
            # geneid id for each record
            gene_ids = tuple(
                self._make_stableid_id(stableid=gene_id, species=species)
                for gene_id, species in group.gene_ids.items()
            )
            # now get the homology id for this group
            homology_id = self._get_homology_group_id(
                relationship_id=rel_type_id,
                gene_ids=gene_ids,
            )
            values = [(int(gene_id), int(homology_id)) for gene_id in gene_ids]
            self.db.executemany(sql, values)

        self.db.commit()

    def get_related_to(self, *, gene_id: str, relationship_type: str) -> homolog_group:
        """return genes with relationship type to gene_id"""
        result = homolog_group(relationship=relationship_type, source=gene_id)
        sql = """
        SELECT st.rowid as stableid_id
        FROM stableid st
        WHERE st.stableid = ?
        """
        stableid_id = self.db.execute(sql, (gene_id,)).fetchone()
        if not stableid_id:
            return result

        sql = """
        SELECT hm.homology_id as homology_id
        FROM homology_member hm
        WHERE hm.homology_type = ? AND hm.stableid_id = ?
        """
        homology_id = self._execute_sql(
            sql,
            (relationship_type, stableid_id[0]),
        ).fetchone()

        if not homology_id:
            return result

        homology_id = homology_id["homology_id"]
        sql = """
        SELECT gs.stableid as stableid,
               gs.species_db as species_db
        FROM gene_species gs
        JOIN homology_member hm ON hm.stableid_id = gs.stableid_id
        WHERE hm.homology_id = ?
        """
        for record in self._execute_sql(sql, (homology_id,)).fetchall():
            result.gene_ids[record["stableid"]] = record["species_db"]

        return result

    def get_related_groups(self, relationship_type: str) -> list[homolog_group]:
        """returns all groups of relationship type"""
        sql = """
        SELECT stableid, species_db, homology_id
        FROM related_groups
        WHERE homology_type = ?
        """
        results = {}
        for result in self._execute_sql(sql, (relationship_type,)).fetchall():
            record = results.get(
                result["homology_id"],
                homolog_group(relationship=relationship_type),
            )
            record.gene_ids |= {result["stableid"]: result["species_db"]}
            results[result["homology_id"]] = record
        return list(results.values())

    def num_records(self):
        return list(
            self._execute_sql("SELECT COUNT(*) as count FROM member").fetchone(),
        )[0]


def load_homology_db(
    *,
    path: PathType,
) -> HomologyDb:
    return HomologyDb(source=path)


@define_app(app_type=LOADER)
class load_homologies:
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
            row_condition=self._matching_species,
            columns=self.src_cols,
            sep="\t",
        )

    def _matching_species(self, row):
        return {row[1], row[3]} <= self._allowed_species

    def main(self, path: IdentifierType) -> SerialisableType:
        parser = self._reader(iter_splitlines(path, chunk_size=500_000))
        header = next(parser)
        assert list(header) == list(self.src_cols), (header, self.src_cols)
        return grouped_related(
            (row[0], {row[2]: row[1], row[4]: row[3]}) for row in parser
        )


@define_app
class collect_seqs:
    """given a config and homolog group, loads genome instances on demand
    and extracts sequences"""

    def __init__(
        self,
        config: elt_config.InstalledConfig,
        make_seq_name: typing.Callable | None = None,
        verbose: bool = False,
    ):
        self._config = config
        self._genomes = {}
        self._namer = make_seq_name
        self._verbose = verbose

    def main(self, homologs: homolog_group) -> SeqsCollectionType:
        namer = self._namer
        seqs = []
        for species, sp_genes in homologs.species_ids().items():
            if species not in self._genomes:
                self._genomes[species] = elt_genome.load_genome(
                    config=self._config,
                    species=species,
                )
            genome = self._genomes[species]
            for name in sp_genes:
                cds = list(genome.get_gene_cds(name=name, is_canonical=True))
                if not cds:
                    if self._verbose:
                        print(f"no cds for {name=} {type(name)=}")
                    continue

                feature = cds[0]
                seq = feature.get_slice()
                seq.name = f"{species}-{name}" if namer is None else namer(feature)
                seq.info["species"] = species
                seq.info["name"] = name
                # disconnect from annotation so the closure of the genome
                # does not cause issues when run in parallel
                seq.annotation_db = None
                seqs.append(seq)

        if not seqs:
            return NotCompleted(
                type="FAIL",
                origin=self,
                message=f"no CDS for {homologs.source=}",
            )

        return make_unaligned_seqs(data=seqs, moltype="dna")
