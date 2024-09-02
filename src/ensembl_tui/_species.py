import os
import pathlib
import typing

from cogent3 import load_table
from cogent3.core.tree import TreeNode
from cogent3.util.table import Table

from ensembl_tui import _util as elt_util

SPECIES_NAME = "species.tsv"
StrOrNone = typing.Union[str, type(None)]


def load_species(species_path: elt_util.PathType) -> list[list[str]]:
    """returns [[latin_name, common_name, stableid prefix],..] from species_path

    if species_path does not exist, defaults to default one"""
    if not os.path.exists(species_path):
        species_path = elt_util.get_resource_path("species.tsv")

    table = load_table(species_path)
    return table.to_list()


_species_common_map = load_species(os.path.join(elt_util.ENSEMBLDBRC, SPECIES_NAME))


class SpeciesNameMap:
    """mapping between common names and latin names"""

    def __init__(self, species_common=_species_common_map):
        """provides latin name:common name mappings"""
        self._species_common = {}
        self._common_species = {}
        self._species_ensembl = {}
        self._ensembl_species = {}
        self._stableid_species = {}  # stable id prefix to species map
        for names in species_common:
            names = list(map(elt_util.CaseInsensitiveString, names))
            self.amend_species(*names)

    def __str__(self) -> str:
        return str(self.to_table())

    def __repr__(self) -> str:
        return repr(self.to_table())

    def __contains__(self, item) -> bool:
        item = elt_util.CaseInsensitiveString(item)
        return any(
            item in attr
            for attr in (
                self._species_common,
                self._common_species,
                self._ensembl_species,
            )
        )

    def _repr_html_(self) -> str:
        table = self.to_table()
        return table._repr_html_()

    def get_common_name(self, name: str, level="raise") -> StrOrNone:
        """returns the common name for the given name (which can be either a
        species name or the ensembl version)"""
        name = elt_util.CaseInsensitiveString(name)
        if name in self._ensembl_species:
            name = self._ensembl_species[name]

        if name in self._species_common:
            common_name = self._species_common[name]
        elif name in self._common_species:
            common_name = name
        else:
            common_name = None

        if common_name is None:
            msg = f"Unknown species name: {name}"
            if level == "raise":
                raise ValueError(msg)
            if level == "warn":
                print(f"WARN: {msg}")

        return common_name

    def get_species_name(self, name: str, level="ignore") -> StrOrNone:
        """returns the species name for the given common name"""
        name = elt_util.CaseInsensitiveString(name)
        if name in self._species_common:
            return name

        species_name = None
        level = level.lower().strip()
        for data in [self._common_species, self._ensembl_species]:
            if name in data:
                species_name = data[name]
        if species_name is None:
            msg = f"Unknown common name: {name}"
            if level == "raise":
                raise ValueError(msg)
            if level == "warn":
                print(f"WARN: {msg}")

        return species_name

    def get_species_names(self) -> typing.Sequence[StrOrNone]:
        """returns the list of species names"""
        return sorted(self._species_common.keys())

    def get_ensembl_db_prefix(self, name: str) -> str:
        """returns a string of the species name in the format used by
        ensembl"""
        name = elt_util.CaseInsensitiveString(name)
        if name in self._common_species:
            name = self._common_species[name]
        try:
            species_name = self.get_species_name(name, level="raise")
        except ValueError as e:
            if name not in self._species_common:
                raise ValueError(f"Unknown name {name}") from e
            species_name = name

        return str(species_name.lower().replace(" ", "_"))

    def get_db_prefix_from_stableid(self, stableid: str) -> str:
        """returns the db name from a stableid"""
        prefix = elt_util.get_stableid_prefix(stableid)
        species = self._stableid_species[prefix]
        return species.replace(" ", "_").lower()

    def _purge_species(self, species_name):
        """removes a species record"""
        species_name = elt_util.CaseInsensitiveString(species_name)
        if species_name not in self._species_common:
            return
        common_name = self._species_common.pop(species_name)
        ensembl_name = self._species_ensembl.pop(species_name)
        self._ensembl_species.pop(ensembl_name)
        self._common_species.pop(common_name)

    def amend_species(self, species_name, common_name, stableid_prefix=None):
        """add a new species, and common name"""
        species_name = elt_util.CaseInsensitiveString(species_name)
        common_name = elt_util.CaseInsensitiveString(common_name)
        assert "_" not in species_name, "'_' in species_name, not a Latin name?"
        self._purge_species(species_name)  # remove if existing
        self._species_common[species_name] = common_name
        self._common_species[common_name] = species_name
        ensembl_name = species_name.lower().replace(" ", "_")
        self._species_ensembl[species_name] = ensembl_name
        self._ensembl_species[ensembl_name] = species_name
        if stableid_prefix:
            # make sure stableid just a string
            for prefix in stableid_prefix.split(","):
                self._stableid_species[prefix] = ensembl_name

    def add_stableid_prefix(
        self,
        species_name: str,
        stableid_prefix: str | elt_util.CaseInsensitiveString,
    ):
        self._stableid_species[str(stableid_prefix)] = self.get_species_name(
            species_name,
        )

    def to_table(self):
        """returns cogent3 Table"""
        rows = []
        for common in self._common_species:
            species = self._common_species[common]
            ensembl = self._species_ensembl[species]
            # all prefixes for this species
            stableids = ",".join(
                [k for k, v in self._stableid_species.items() if v == ensembl],
            )
            rows += [[species, common, stableids]]
        return Table(
            [
                "Species name",
                "Common name",
                "Ensembl stableid Prefix",
            ],
            data=rows,
            space=2,
        ).sorted()

    def update_from_file(self, species_path: pathlib.Path) -> None:
        """updates instance from tab delimited table at species_path"""
        table = load_table(species_path)
        columns = "Species name", "Ensembl stableid Prefix"
        for sp, prefixes in table.to_list(columns=columns):
            db_name = self.get_ensembl_db_prefix(sp)
            for prefix in prefixes.split(","):
                self._stableid_species[prefix] = db_name


Species = SpeciesNameMap()


def species_from_ensembl_tree(tree: TreeNode) -> dict[str, str]:
    """get species identifiers from an Ensembl tree"""
    tip_names = tree.get_tip_names()
    selected_species = {}
    for tip_name in tip_names:
        name_fields = tip_name.lower().split("_")
        # produce parts of name starting with highly specific to
        # more general and look for matches
        for j in range(len(name_fields) + 1, 1, -1):
            n = "_".join(name_fields[:j])
            if n in Species:
                selected_species[Species.get_common_name(n)] = n
                break
        else:
            raise ValueError(f"cannot establish species for {'_'.join(name_fields)}")

    return selected_species
