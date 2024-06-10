import dataclasses
import pathlib
import typing

from abc import ABC, abstractmethod
from typing import Any, Optional

import click
import h5py
import numpy

from cogent3 import get_moltype, load_annotations, make_seq, make_table, open_
from cogent3.app.composable import define_app
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.core.sequence import Sequence
from cogent3.parse.fasta import MinimalFastaParser
from cogent3.parse.gff import GffRecord, gff_parser, is_gff3
from cogent3.util.io import iter_splitlines
from cogent3.util.table import Table
from numpy.typing import NDArray

from ensembl_lite._config import Config, InstalledConfig
from ensembl_lite._db_base import Hdf5Mixin
from ensembl_lite._species import Species
from ensembl_lite._util import _HDF5_BLOSC2_KWARGS, PathType


_SEQDB_NAME = "genome_sequence.hdf5_blosc2"
_ANNOTDB_NAME = "features.gff3db"
_typed_id = re.compile(
    r"\b[a-z]+:", flags=re.IGNORECASE
)  # ensembl stableid's prefixed by the type
_feature_id = re.compile(r"(?<=\bID=)[^;]+")
_exon_id = re.compile(r"(?<=\bexon_id=)[^;]+")
_parent_id = re.compile(r"(?<=\bParent=)[^;]+")


def _lower_case_match(match) -> str:
    return match.group(0).lower()


def tidy_gff3_stableids(attrs: str) -> str:
    """makes the feature type prefix lowercase in gff3 attribute fields"""
    return _typed_id.sub(_lower_case_match, attrs)


class EnsemblGffRecord(GffRecord):
    __slots__ = GffRecord.__slots__ + ("feature_id",)

    def __init__(self, feature_id: Optional[int] = None, **kwargs):
        is_canonical = kwargs.pop("is_canonical", None)
        super().__init__(**kwargs)
        self.feature_id = feature_id
        if is_canonical:
            self.attrs = "Ensembl_canonical;" + self.attrs or ""

    def __hash__(self) -> int:
        return hash(self.name)

    def __eq__(self, other):
        return self.name == getattr(other, "name", other)

    @property
    def stableid(self):
        return _typed_id.sub("", self.name or "")

    @property
    def is_canonical(self):
        attrs = self.attrs or ""
        return "Ensembl_canonical" in attrs

    def update_from_attrs(self) -> None:
        """updates attributes from the attrs string

        Notes
        -----
        also updates biotype from the prefix in the name
        """
        attrs = self.attrs
        id_regex = _feature_id if "ID=" in attrs else _exon_id
        attr = tidy_gff3_stableids(attrs)
        if feature_id := id_regex.search(attr):
            self.name = feature_id.group()

        if pid := _parent_id.search(attr):
            parents = pid.group().split(",")
            # now sure how to handle multiple-parent features
            # so taking first ID as the parent for now
            self.parent_id = parents[0]

        if ":" in (self.name or ""):
            biotype = self.name.split(":")[0]
            self.biotype = "mrna" if biotype == "transcript" else biotype

    @property
    def size(self) -> int:
        """the sum of span segments"""
        return 0 if self.spans is None else sum(abs(s - e) for s, e in self.spans)


def custom_gff_parser(
    path: PathType, num_fake_ids: int
) -> tuple[dict[str, EnsemblGffRecord], int]:
    """replacement for cogent3 merged_gff_records"""
    reduced = {}
    gff3 = is_gff3(path)
    for record in gff_parser(
        iter_splitlines(path),
        gff3=gff3,
        make_record=EnsemblGffRecord,
    ):
        record.update_from_attrs()
        if not record.name:
            record.name = f"unknown-{num_fake_ids}"
            num_fake_ids += 1

        if record.name not in reduced:
            record.spans = record.spans or []
            reduced[record] = record

        reduced[record].spans.append([record.start, record.stop])
        reduced[record].start = min(reduced[record].start, record.start)
        reduced[record].stop = max(reduced[record].stop, record.stop)

    return reduced, num_fake_ids


def make_annotation_db(src_dest: tuple[PathType, PathType]) -> bool:
    """convert gff3 file into aan AnnotationDb

    Parameters
    ----------
    src_dest
        path to gff3 file, path to write AnnotationDb
    """
    src, dest = src_dest
    if dest.exists():
        return True

    db = load_annotations(path=src, write_path=dest)
    db.close()
    del db
    return True


def _rename(label: str) -> str:
    return label.split()[0]


@define_app
class fasta_to_hdf5:
    def __init__(self, config: Config, label_to_name=_rename):
        self.config = config
        self.label_to_name = label_to_name

    def main(self, db_name: str) -> bool:
        src_dir = self.config.staging_genomes / db_name
        dest_dir = self.config.install_genomes / db_name

        seq_store = SeqsDataHdf5(
            source=dest_dir / _SEQDB_NAME,
            species=Species.get_species_name(db_name),
            mode="w",
        )

        src_dir = src_dir / "fasta"
        for path in src_dir.glob("*.fa.gz"):
            for label, seq in MinimalFastaParser(iter_splitlines(path)):
                seqid = self.label_to_name(label)
                seq_store.add_record(seqid=seqid, seq=seq)
                del seq

        seq_store.close()

        return True


def _get_seqs(src: PathType) -> list[tuple[str, str]]:
    with open_(src) as infile:
        data = infile.read().splitlines()
    name_seqs = list(MinimalFastaParser(data))
    return [(_rename(name), seq) for name, seq in name_seqs]


T = tuple[PathType, list[tuple[str, str]]]


class SeqsDataABC(ABC):
    """interface for genome sequence storage"""

    # the storage reference, e.g. path to file
    source: PathType
    species: str
    mode: str  # as per standard file opening modes, r, w, a
    _is_open = False
    _file: Optional[Any] = None

    @abstractmethod
    def __hash__(self): ...

    @abstractmethod
    def add_record(self, *, seqid: str, seq: str): ...

    @abstractmethod
    def add_records(self, *, records: typing.Iterable[list[str, str]]): ...

    @abstractmethod
    def get_seq_str(
        self, *, seqid: str, start: Optional[int] = None, stop: Optional[int] = None
    ) -> str: ...

    @abstractmethod
    def get_seq_arr(
        self, *, seqid: str, start: Optional[int] = None, stop: Optional[int] = None
    ) -> NDArray[numpy.uint8]: ...

    @abstractmethod
    def get_coord_names(self) -> tuple[str]: ...

    @abstractmethod
    def close(self): ...


@define_app
class str2arr:
    """convert string to array of uint8"""

    def __init__(self, moltype: str = "dna", max_length=None):
        moltype = get_moltype(moltype)
        self.canonical = "".join(moltype)
        self.max_length = max_length
        extended = "".join(list(moltype.alphabets.degen))
        self.translation = b"".maketrans(
            extended.encode("utf8"),
            "".join(chr(i) for i in range(len(extended))).encode("utf8"),
        )

    def main(self, data: str) -> numpy.ndarray:
        if self.max_length:
            data = data[: self.max_length]

        b = data.encode("utf8").translate(self.translation)
        return numpy.array(memoryview(b), dtype=numpy.uint8)


@define_app
class arr2str:
    """convert array of uint8 to str"""

    def __init__(self, moltype: str = "dna", max_length=None):
        moltype = get_moltype(moltype)
        self.canonical = "".join(moltype)
        self.max_length = max_length
        extended = "".join(list(moltype.alphabets.degen))
        self.translation = b"".maketrans(
            "".join(chr(i) for i in range(len(extended))).encode("utf8"),
            extended.encode("utf8"),
        )

    def main(self, data: numpy.ndarray) -> str:
        if self.max_length:
            data = data[: self.max_length]

        b = data.tobytes().translate(self.translation)
        return bytearray(b).decode("utf8")


@dataclasses.dataclass
class SeqsDataHdf5(Hdf5Mixin, SeqsDataABC):
    """HDF5 sequence data storage"""

    def __init__(
        self,
        source: PathType,
        species: Optional[str] = None,
        mode: str = "r",
        in_memory: bool = False,
    ):
        # note that species are converted into the Ensembl db prefix

        source = pathlib.Path(source)
        self.source = source

        if mode == "r" and not source.exists():
            raise OSError(f"{self.source!s} not found")

        species = Species.get_ensembl_db_prefix(species) if species else None
        self.mode = "w-" if mode == "w" else mode
        if in_memory:
            h5_kwargs = dict(
                driver="core",
                backing_store=False,
            )
        else:
            h5_kwargs = {}

        try:
            self._file = h5py.File(source, mode=self.mode, **h5_kwargs)
        except OSError:
            print(source)
            raise
        self._str2arr = str2arr(moltype="dna")
        self._arr2str = arr2str(moltype="dna")
        self._is_open = True
        if "r" not in self.mode and "species" not in self._file.attrs:
            assert species
            self._file.attrs["species"] = species

        if (
            species
            and (file_species := self._file.attrs.get("species", None)) != species
        ):
            raise ValueError(f"{self.source.name!r} {file_species!r} != {species}")
        self.species = self._file.attrs["species"]

    def __hash__(self):
        return id(self)

    def add_record(self, *, seqid: str, seq: str):
        seq = self._str2arr(seq)
        if seqid in self._file:
            stored = self._file[seqid]
            if (seq == stored).all():
                # already seen this seq
                return
            # but it's different, which is a problem
            raise ValueError(f"{seqid!r} already present but with different seq")
        self._file.create_dataset(
            name=seqid, data=seq, chunks=True, **_HDF5_BLOSC2_KWARGS
        )

    def add_records(self, *, records: typing.Iterable[list[str, str]]):
        for seqid, seq in records:
            self.add_record(seqid=seqid, seq=seq)

    def get_seq_str(
        self, *, seqid: str, start: Optional[int] = None, stop: Optional[int] = None
    ) -> str:
        return self._arr2str(self.get_seq_arr(seqid=seqid, start=start, stop=stop))

    def get_seq_arr(
        self, *, seqid: str, start: Optional[int] = None, stop: Optional[int] = None
    ) -> NDArray[numpy.uint8]:
        if not self._is_open:
            raise OSError(f"{self.source.name!r} is closed")

        return self._file[seqid][start:stop]

    def get_coord_names(self) -> tuple[str]:
        """names of chromosomes / contig"""
        return tuple(self._file)


# todo: this wrapping class is required for memory efficiency because
#  the cogent3 SequenceCollection class is not designed for large sequence
#  collections, either large sequences or large numbers of sequences. The
#  longer term solution is improving SequenceCollections,
#  which is underway 🎉
class Genome:
    """class to be replaced by cogent3 sequence collection when that
    has been modernised"""

    def __init__(
        self,
        *,
        species: str,
        seqs: SeqsDataABC,
        annots: GffAnnotationDb,
    ) -> None:
        self.species = species
        self._seqs = seqs
        self.annotation_db = annots

    def get_seq(
        self,
        *,
        seqid: str,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        namer: typing.Callable | None = None,
    ) -> str:
        """returns annotated sequence

        Parameters
        ----------
        seqid
            name of chromosome etc..
        start
            starting position of slice in python coordinates, defaults
            to 0
        stop
            ending position of slice in python coordinates, defaults
            to length of coordinate
        namer
            callback for naming the sequence. Callback must take four
            arguments: species, seqid,start, stop. Default is
            species:seqid:start-stop.
        Notes
        -----
        Annotations partially within region are included.
        """
        seq = self._seqs.get_seq_str(seqid=seqid, start=start, stop=stop)
        if namer:
            name = namer(self.species, seqid, start, stop)
        else:
            name = f"{self.species}:{seqid}:{start}-{stop}"
        # we use seqid to make the sequence here because that identifies the
        # parent seq identity, required for querying annotations
        seq = make_seq(seq, name=seqid, moltype="dna")
        seq.name = name
        if self.annotation_db:
            seq.annotation_offset = start or 0
            seq.annotation_db = self.annotation_db.subset(
                seqid=seqid, start=start, stop=stop, allow_partial=True
            )
        return seq

    def get_features(
        self,
        *,
        biotype: str = None,
        seqid: str = None,
        name: str = None,
        start: int = None,
        stop: int = None,
    ) -> typing.Iterable[Feature]:
        """yields features in blocks of seqid"""
        kwargs = {k: v for k, v in locals().items() if k not in ("self", "seqid") and v}
        if seqid:
            seqids = [seqid]
        else:
            seqids = {
                ft["seqid"] for ft in self.annotation_db.get_features_matching(**kwargs)
            }
        for seqid in seqids:
            try:
                seq = self.get_seq(seqid=seqid)
            except TypeError:
                msg = f"ERROR (report me): {self.species!r}, {seqid!r}"
                raise TypeError(msg)
            # because self.get_seq() automatically names seqs differently
            seq.name = seqid
            yield from seq.get_features(**kwargs)

    def get_ids_for_biotype(self, biotype, limit=None):
        annot_db = self.annotation_db
        sql = "SELECT name from gff WHERE biotype=?"
        if limit:
            sql += " LIMIT ?"
        for result in annot_db._execute_sql(sql, (biotype, limit)):
            yield result["name"].split(":")[-1]

    def close(self):
        self._seqs.close()
        self.annotation_db.db.close()


def load_genome(*, config: InstalledConfig, species: str):
    """returns the Genome with bound seqs and features"""
    genome_path = config.installed_genome(species) / _SEQDB_NAME
    seqs = SeqsDataHdf5(source=genome_path, species=species, mode="r")
    ann_path = config.installed_genome(species) / _ANNOTDB_NAME
    ann = GffAnnotationDb(source=ann_path)
    return Genome(species=species, seqs=seqs, annots=ann)


def get_seqs_for_ids(
    *,
    config: InstalledConfig,
    species: str,
    names: list[str],
    make_seq_name: typing.Callable = None,
) -> typing.Iterable[Sequence]:
    genome = load_genome(config=config, species=species)
    # is it possible to do batch query for all names?
    for name in names:
        feature = list(genome.get_features(name=f"%{name}"))[0]
        transcripts = list(feature.get_children(biotype="mRNA"))
        if not transcripts:
            continue

        longest = max(transcripts, key=lambda x: len(x))
        cds = list(longest.get_children(biotype="CDS"))
        if not cds:
            continue

        feature = cds[0]
        seq = feature.get_slice()
        if callable(make_seq_name):
            seq.name = make_seq_name(feature)
        else:
            seq.name = f"{species}-{name}"
        seq.info["species"] = species
        seq.info["name"] = name
        # disconnect from annotation so the closure of the genome
        # does not cause issues when run in parallel
        seq.annotation_db = None
        yield seq

    genome.close()
    del genome


def get_annotations_for_species(
    *, config: InstalledConfig, species: str
) -> GffAnnotationDb:
    """returns the annotation Db for species"""
    path = config.installed_genome(species=species)
    if not path.exists():
        click.secho(f"{species!r} not in {str(config.install_path.parent)!r}", fg="red")
        exit(1)
    # TODO: this filename should be defined in one place
    path = path / "features.gff3db"
    if not path.exists():
        click.secho(f"{path.name!r} is missing", fg="red")
        exit(1)
    return GffAnnotationDb(source=path)


def get_gene_table_for_species(
    *, annot_db: GffAnnotationDb, limit: Optional[int], species: Optional[str] = None
) -> Table:
    """
    returns gene data from a GffDb

    Parameters
    ----------
    annot_db
        feature db
    limit
        limit number of records to
    species
        species name, overrides inference from annot_db.source
    """
    species = species or annot_db.source.parent.name

    columns = (
        "species",
        "name",
        "seqid",
        "source",
        "biotype",
        "start",
        "stop",
        "score",
        "strand",
        "phase",
    )
    rows = []
    for i, record in enumerate(annot_db.get_records_matching(biotype="gene")):
        rows.append([species] + [record.get(c, None) for c in columns[1:]])
        if i == limit:
            break

    return make_table(header=columns, data=rows)


def get_species_summary(
    *, annot_db: GffAnnotationDb, species: Optional[str] = None
) -> Table:
    """
    returns the Table summarising data for species_name

    Parameters
    ----------
    annot_db
        feature db
    species
        species name, overrides inference from annot_db.source
    """
    from ._species import Species

    # for now, just biotype
    species = species or annot_db.source.parent.name
    counts = annot_db.biotype_counts()
    try:
        common_name = Species.get_common_name(species)
    except ValueError:
        common_name = species

    return Table(
        header=("biotype", "count"),
        data=list(counts.items()),
        title=f"{common_name} features",
    )
