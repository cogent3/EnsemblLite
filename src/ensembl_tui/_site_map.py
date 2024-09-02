import typing
from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import cache

from cogent3.util.misc import extend_docstring_from

_ensembl_site_map = {}


class register_ensembl_site_map:
    """
    registration decorator for Ensembl site-map classes

    The registration key must be a string that of the domain name.

    Parameters
    ----------
    domain: str of domain name, must be unique
    """

    def __init__(self, domain: str):
        if not isinstance(domain, str):
            raise TypeError(f"{domain!r} is not a string")

        domain = domain.strip()
        if not domain:
            raise ValueError("cannot have empty string domain")

        assert (
            domain not in _ensembl_site_map
        ), f"{domain!r} already in {list(_ensembl_site_map)}"

        self._domain = domain

    def __call__(self, func):
        # pass through
        _ensembl_site_map[self._domain] = func
        return func


StrOrNone = typing.Union[str, type(None)]


class SiteMapABC(ABC):
    @abstractmethod
    def get_seqs_path(self, ensembl_name: str) -> str:
        """returns the path to genome sequences for species_db_name"""
        ...

    @abstractmethod
    def get_annotations_path(self, ensembl_name: str) -> str: ...

    @property
    def alignments_path(self) -> StrOrNone:
        return self._alignments_path

    @property
    def homologies_path(self) -> StrOrNone:
        return self._homologies_path

    @property
    def trees_path(self) -> StrOrNone:
        return self._trees_path


@dataclass(slots=True)
class SiteMap:
    """records the locations of specific attributes relative to an Ensembl release"""

    site: str
    _seqs_path: str = "fasta"
    _annotations_path: str = "gff3"
    _alignments_path: str | None = None
    _homologies_path: str | None = None
    _trees_path: str | None = None


class EnsemblPrimary(SiteMapABC, SiteMap):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_seqs_path(self, ensembl_name: str) -> str:
        """path to unmasked genome sequences"""
        return f"{self._seqs_path}/{ensembl_name}/dna"

    def get_annotations_path(self, ensembl_name: str) -> str:
        return f"{self._annotations_path}/{ensembl_name}"


@extend_docstring_from(SiteMap)
@register_ensembl_site_map("ftp.ensembl.org")
def ensembl_main_sitemap():
    """the main Ensembl site map"""
    return EnsemblPrimary(
        site="ftp.ensembl.org",
        _alignments_path="maf/ensembl-compara/multiple_alignments",
        _homologies_path="tsv/ensembl-compara/homologies",
        _trees_path="compara/species_trees",
    )


# for bacteria we have, but complexities related to the bacterial collection
# a species belongs to. For example
# https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-57/fasta/bacteria_15_collection/_butyribacterium_methylotrophicum_gca_001753695/dna/
# so to address this, the sitemap class needs to download the species table from ensembl bacteria
# and cache the collection/species mapping
# site = "ftp.ensemblgenomes.ebi.ac.uk",
# _genomes_path = "fasta",
# _annotations_path = "gff3",
# _homologies_path = "pan_ensembl/tsv/ensembl-compara/homologies",


@cache
def get_site_map(domain: str) -> SiteMapABC:
    """returns a site map instance"""
    return _ensembl_site_map[domain]()
