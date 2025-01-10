from __future__ import annotations

import re
from dataclasses import dataclass

import typing_extensions

from ._species import Species

_release = re.compile(r"\d+")


def get_version_from_name(name):
    """returns the release and build identifiers from an ensembl db_name"""
    r = _release.search(name)
    if r is None:
        return None, None

    # first number run is release, followed by build
    # note, for the ensemblgenomes naming system, the second digit run is the
    # standard Ensembl release and the first is for the specified genome
    release = name[r.start() : r.end()]
    b = [s for s in _name_delim.split(name[r.end() :]) if s]

    return release, b


_name_delim = re.compile("_")


def get_dbtype_from_name(name: str) -> str:
    """returns the data base type from the name"""
    parts = _release.split(name)
    parts = [s for s in _name_delim.split(parts[0]) if s]
    return parts[1] if parts[0] == "ensembl" else parts[-1]


def get_db_prefix(name: str) -> str:
    """returns the db prefix, typically an organism or `ensembl'"""
    parts = _release.split(name)
    parts = [s for s in _name_delim.split(parts[0]) if s]
    if parts[0] == "ensembl":
        prefix = "ensembl"
    elif len(parts) > 2:
        prefix = "_".join(parts[:-1])
    else:
        msg = f"Unknown name structure: {'_'.join(parts)}"
        raise ValueError(msg)
    return prefix


class EnsemblDbName:
    """container for a db name, inferring different attributes from the name,
    such as species, version, build"""

    def __init__(self, db_name: str) -> None:
        """db_name: and Emsembl database name"""
        self.name = db_name
        self.type = get_dbtype_from_name(db_name)
        self.prefix = get_db_prefix(db_name)

        release, build = get_version_from_name(db_name)
        self.release = release
        self.general_release = self.release

        self.build = None
        if build and len(build) == 1:
            if self.type != "compara":
                self.build = build[0]
            else:
                self.general_release = build[0]
        elif build:
            self.build = build[1]
            self.general_release = build[0]

        self.species = Species.get_species_name(self.prefix)

    def __repr__(self) -> str:
        build = f"; build='{self.build}'" if self.build is not None else ""
        return f"db(prefix='{self.prefix}'; type='{self.type}'; release='{self.release}'{build})"

    def __str__(self) -> str:
        return self.name

    def __lt__(self, other: typing_extensions.Self | str) -> bool:
        if isinstance(other, type(self)):
            return self.name < other.name
        return self.name < other

    def __eq__(self, other: typing_extensions.Self) -> bool:
        if isinstance(other, type(self)):
            return self.name == other.name
        return self.name == other

    def __ne__(self, other: typing_extensions.Self | str) -> bool:
        if isinstance(other, type(self)):
            return self.name != other.name
        return self.name != other

    def __hash__(self) -> int:
        return hash(self.name)


@dataclass(slots=True)
class EmfName:
    """stores information from EMF SEQ records"""

    species: str
    seqid: str
    start: int
    stop: int
    strand: str
    coord_length: str

    def __post_init__(self) -> None:
        # adjust the lengths to be ints and put into python coord
        self.start = int(self.start) - 1
        self.stop = int(self.stop)

    def __str__(self) -> str:
        attrs = "species", "seqid", "start", "stop", "strand"
        n = [str(getattr(self, attr)) for attr in attrs]
        return ":".join(n)

    def __hash__(self) -> int:
        return hash(str(self))

    def to_dict(self) -> dict:
        attrs = "species", "seqid", "start", "stop", "strand"
        return {attr: getattr(self, attr) for attr in attrs}


@dataclass(slots=True)
class MafName:
    """stores source information from Maf records"""

    species: str
    seqid: str
    start: int
    stop: int
    strand: str
    coord_length: str | int | None

    def __post_init__(self) -> None:
        # adjust the lengths to be ints
        self.start = int(self.start)
        self.stop = int(self.stop)
        self.coord_length = int(self.coord_length) if self.coord_length else None

    def __str__(self) -> str:
        attrs = "species", "seqid", "start", "stop", "strand"
        n = [str(getattr(self, attr)) for attr in attrs]
        return ":".join(n)

    def __hash__(self) -> int:
        return hash(str(self))

    def to_dict(self) -> dict:
        attrs = "species", "seqid", "start", "stop", "strand"
        return {attr: getattr(self, attr) for attr in attrs}
