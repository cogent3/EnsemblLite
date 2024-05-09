import configparser
import fnmatch
import pathlib

from dataclasses import dataclass
from typing import Iterable

from ensembl_lite._species import Species, species_from_ensembl_tree
from ensembl_lite._util import PathType


INSTALLED_CONFIG_NAME = "installed.cfg"
DOWNLOADED_CONFIG_NAME = "downloaded.cfg"

_COMPARA_NAME: str = "compara"
_ALIGNS_NAME: str = "aligns"
_HOMOLOGIES_NAME: str = "homologies"
_GENOMES_NAME: str = "genomes"


@dataclass
class Config:
    host: str
    remote_path: str
    release: str
    staging_path: PathType
    install_path: PathType
    species_dbs: dict[str, list[str]]
    align_names: Iterable[str]
    tree_names: Iterable[str]
    homologies: bool

    def update_species(self, species: dict[str, list[str]]):
        if not species:
            return
        for k in species:
            if k not in Species:
                raise ValueError(f"Unknown species {k}")
        self.species_dbs |= species

    @property
    def db_names(self) -> Iterable[str]:
        for species in self.species_dbs:
            yield Species.get_ensembl_db_prefix(species)

    @property
    def staging_genomes(self):
        return self.staging_path / _GENOMES_NAME

    @property
    def install_genomes(self):
        return self.install_path / _GENOMES_NAME

    @property
    def staging_homologies(self):
        return self.staging_path / _COMPARA_NAME / _HOMOLOGIES_NAME

    @property
    def install_homologies(self):
        return self.install_path / _COMPARA_NAME / _HOMOLOGIES_NAME

    @property
    def staging_aligns(self):
        return self.staging_path / _COMPARA_NAME / _ALIGNS_NAME

    @property
    def install_aligns(self):
        return self.install_path / _COMPARA_NAME / _ALIGNS_NAME

    def to_dict(self):
        """returns cfg as a dict"""
        if not self.db_names:
            raise ValueError("no db names")

        data = {
            "remote path": {"path": str(self.remote_path), "host": str(self.host)},
            "local path": {
                "staging_path": str(self.staging_path),
                "install_path": str(self.install_path),
            },
            "release": {"release": self.release},
        }

        if self.align_names or self.tree_names:
            data["compara"] = {}

        if self.align_names:
            data["compara"]["align_names"] = "".join(self.align_names)
        if self.tree_names:
            data["compara"]["tree_names"] = "".join(self.tree_names)

        for db_name in self.db_names:
            data[db_name] = {"db": "core"}

        return data

    def write(self):
        """writes a ini to staging_path/DOWNLOADED_CONFIG_NAME"""
        parser = configparser.ConfigParser()
        cfg = self.to_dict()
        for section, settings in cfg.items():
            parser.add_section(section)
            for option, val in settings.items():
                parser.set(section, option=option, value=val)
        self.staging_path.mkdir(parents=True, exist_ok=True)
        with (self.staging_path / DOWNLOADED_CONFIG_NAME).open(mode="w") as out:
            parser.write(out, space_around_delimiters=True)


@dataclass
class InstalledConfig:
    release: str
    install_path: PathType

    def __hash__(self):
        return id(self)

    def __post_init__(self):
        self.install_path = pathlib.Path(self.install_path)

    @property
    def compara_path(self):
        return self.install_path / _COMPARA_NAME

    @property
    def homologies_path(self):
        return self.compara_path / _HOMOLOGIES_NAME

    @property
    def aligns_path(self):
        return self.compara_path / _ALIGNS_NAME

    @property
    def genomes_path(self):
        return self.install_path / _GENOMES_NAME

    def installed_genome(self, species: str) -> PathType:
        db_name = Species.get_ensembl_db_prefix(species)
        return self.genomes_path / db_name

    def list_genomes(self):
        """returns list of installed genomes"""
        return [p.name for p in self.genomes_path.glob("*") if p.name in Species]

    def path_to_alignment(self, pattern: str) -> PathType | None:
        """returns the full path to alignment matching the name

        Parameters
        ----------
        pattern
            glob pattern for the Ensembl alignment name
        """
        align_dirs = [
            d for d in self.aligns_path.glob("*") if fnmatch.fnmatch(d.name, pattern)
        ]
        if not align_dirs:
            return None

        if len(align_dirs) > 1:
            raise ValueError(
                f"{pattern!r} matches too many directories in {self.aligns_path}"
            )

        return align_dirs[0]


def write_installed_cfg(config: Config) -> PathType:
    """writes an ini file under config.installed_path"""
    parser = configparser.ConfigParser()
    parser.add_section("release")
    parser.set("release", "release", config.release)
    # create all the genome
    outpath = config.install_path / INSTALLED_CONFIG_NAME
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with outpath.open(mode="w") as out:
        parser.write(out)
    return outpath


def read_installed_cfg(path: PathType) -> InstalledConfig:
    """reads an ini file under config.installed_path"""
    parser = configparser.ConfigParser()
    path = (
        path if path.name == INSTALLED_CONFIG_NAME else (path / INSTALLED_CONFIG_NAME)
    )
    if not path.exists():
        print(f"{str(path)} does not exist, exiting")
        exit(1)

    parser.read(path)
    release = parser.get("release", "release")
    return InstalledConfig(release=release, install_path=path.parent)


def read_config(config_path) -> Config:
    """returns ensembl release, local path, and db specifics from the provided
    config path"""
    from ensembl_lite._download import download_ensembl_tree

    parser = configparser.ConfigParser()

    with config_path.expanduser().open() as f:
        parser.read_file(f)

    release = parser.get("release", "release")
    host = parser.get("remote path", "host")
    remote_path = parser.get("remote path", "path")
    remote_path = remote_path[:-1] if remote_path.endswith("/") else remote_path
    staging_path = (
        pathlib.Path(parser.get("local path", "staging_path")).expanduser().absolute()
    )
    install_path = (
        pathlib.Path(parser.get("local path", "install_path")).expanduser().absolute()
    )
    homologies = parser.has_option("compara", "homologies")
    species_dbs = {}
    get_option = parser.get
    align_names = []
    tree_names = []
    for section in parser.sections():
        if section in ("release", "remote path", "local path"):
            continue

        if section == "compara":
            value = get_option(section, "align_names", fallback=None)
            align_names = [] if value is None else [n.strip() for n in value.split(",")]
            value = get_option(section, "tree_names", fallback=None)
            tree_names = [] if value is None else [n.strip() for n in value.split(",")]
            continue

        dbs = [db.strip() for db in get_option(section, "db").split(",")]

        # handle synonyms
        species = Species.get_species_name(section, level="raise")
        species_dbs[species] = dbs

    # we also want homologies if we want alignments
    homologies = homologies or bool(align_names)

    if tree_names:
        # add all species in the tree to species_dbs
        for tree_name in tree_names:
            tree = download_ensembl_tree(host, remote_path, release, tree_name)
            sp = species_from_ensembl_tree(tree)
            species_dbs.update(sp)

    return Config(
        host=host,
        remote_path=remote_path,
        release=release,
        staging_path=staging_path,
        install_path=install_path,
        species_dbs=species_dbs,
        align_names=align_names,
        tree_names=tree_names,
        homologies=homologies,
    )
