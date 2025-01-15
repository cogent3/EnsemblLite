import dataclasses
import typing

import blosc2
import typing_extensions
from cogent3 import make_table, make_unaligned_seqs
from cogent3.app.composable import NotCompleted, define_app
from cogent3.app.io import compress, decompress, pickle_it, unpickle_it
from cogent3.app.typing import (
    SeqsCollectionType,
)
from cogent3.util.io import PathType

from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome
from ensembl_tui import _storage_mixin as eti_storage

compressor = compress(compressor=blosc2.compress2)
decompressor = decompress(decompressor=blosc2.decompress2)
pickler = pickle_it()
unpickler = unpickle_it()
inflate = decompressor + unpickler


@dataclasses.dataclass(slots=True)
class homolog_group:  # noqa: N801
    """encloses {stableid: species, ...} of homologous sequences"""

    relationship: str
    # gene id -> species
    gene_ids: dict[str, str] = None  # type: ignore
    source: str | None = None

    def __post_init__(self) -> None:
        self.gene_ids = self.gene_ids or {}
        if self.source is None:
            self.source = next(iter(self.gene_ids), None)

    def __hash__(self) -> int:
        # allow hashing, but bearing in mind we are updating
        # gene values
        return hash((hash(self.relationship), id(self.gene_ids)))

    def __eq__(self, other: typing_extensions.Self) -> bool:
        return (
            self.relationship == other.relationship and self.gene_ids == other.gene_ids
        )

    def __getstate__(self) -> tuple[str, dict[str, str] | None, str | None]:
        return self.relationship, self.gene_ids, self.source

    def __setstate__(
        self,
        state: tuple[str, dict[str, str] | None, str | None],
    ) -> None:
        relationship, gene_ids, source = state
        self.relationship = relationship
        self.gene_ids = gene_ids
        self.source = source

    def __len__(self) -> int:
        return len(self.gene_ids or ())

    def __or__(self, other: typing_extensions.Self) -> typing_extensions.Self:
        if other.relationship != self.relationship:
            msg = f"relationship type {self.relationship!r} != {other.relationship!r}"
            raise ValueError(msg)
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


# the homology db stores pairwise relationship information
@dataclasses.dataclass(slots=True)
class HomologyDb(eti_storage.DuckdbParquetBase):
    _tables: tuple[str, str] = ("homology_groups_attr",)

    def get_related_to(self, *, gene_id: str, relationship_type: str) -> homolog_group:
        """return genes with relationship type to gene_id"""
        result = homolog_group(relationship=relationship_type, source=gene_id)
        sql = """
        SELECT homology_id
        FROM homology_groups_attr
        WHERE homology_type = ? AND stableid = ?
        """
        homology_id = self.conn.sql(
            sql,
            params=(relationship_type, gene_id),
        ).fetchone()

        if not homology_id:
            return result

        homology_id = homology_id[0]
        sql = """SELECT
        homology_id,
        STRING_AGG(stableid, ' ') AS agg_stableid,
        STRING_AGG(species_db, ' ') AS agg_species_db
        FROM homology_groups_attr
        WHERE homology_id = ?
        GROUP BY homology_id
        """
        for _, stableids, species_dbs in self.conn.sql(
            sql,
            params=(homology_id,),
        ).fetchall():
            result.gene_ids |= dict(
                zip(stableids.split(), species_dbs.split(), strict=False),
            )

        return result

    def get_related_groups(self, relationship_type: str) -> list[homolog_group]:
        """returns all groups of relationship type"""
        sql = """
        SELECT 
        STRING_AGG(stableid, ' ') AS agg_stableid,
        STRING_AGG(species_db, ' ') AS agg_species_db
        FROM homology_groups_attr
        WHERE homology_type = ?
        GROUP BY homology_id
        """
        results = []
        for stableids, species_dbs in self.conn.sql(
            sql,
            params=(relationship_type,),
        ).fetchall():
            result = homolog_group(
                relationship=relationship_type,
                gene_ids=dict(
                    zip(stableids.split(), species_dbs.split(), strict=False),
                ),
            )
            results.append(result)
        return results

    def num_records(self) -> int:
        """number of distinct homology_id's"""
        sql = "SELECT COUNT(DISTINCT homology_id) FROM homology_groups_attr"
        return self.conn.sql(sql).fetchone()[0]

    def count_distinct(self, species: bool = False, homology_type: bool = False) -> int:
        columns = []
        if species:
            columns.append("species_db")
        if homology_type:
            columns.append("homology_type")

        if not columns:
            msg = "must specify at least one of species or homology_type"
            raise ValueError(msg)

        cols = ",".join(columns)
        sql = f"SELECT {cols}, COUNT(*) FROM homology_groups_attr GROUP BY {cols}"
        return make_table(
            header=[*columns, "count"],
            data=self.conn.sql(sql).fetchall(),
        )


def load_homology_db(
    *,
    path: PathType,
) -> HomologyDb:
    return HomologyDb(source=path)


@define_app
class collect_seqs:
    """given a config and homolog group, loads genome instances on demand
    and extracts sequences"""

    def __init__(
        self,
        config: eti_config.InstalledConfig,
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
                self._genomes[species] = eti_genome.load_genome(
                    config=self._config,
                    species=species,
                )
            genome = self._genomes[species]
            for name in sp_genes:
                cds = list(genome.get_cds(stable_id=name))
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
                source=homologs.source,
            )

        return make_unaligned_seqs(
            data=seqs,
            moltype="dna",
            info=dict(source=homologs.source),
        )
