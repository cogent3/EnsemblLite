import dataclasses
import pathlib
import typing

from abc import ABC, abstractmethod
from typing import Any, Optional

import click
import h5py
import numpy

from cogent3 import get_moltype, make_seq, make_table, make_unaligned_seqs
from cogent3.app.composable import NotCompleted, define_app
from cogent3.app.typing import UnalignedSeqsType
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.core.sequence import Sequence
from cogent3.util.table import Table
from numpy.typing import NDArray

from ensembl_lite._config import InstalledConfig
from ensembl_lite._db_base import Hdf5Mixin
from ensembl_lite._homologydb import homolog_group
from ensembl_lite._species import Species
from ensembl_lite._util import _HDF5_BLOSC2_KWARGS, PathType


_SEQDB_NAME = "genome_sequence.hdf5_blosc2"
_ANNOTDB_NAME = "features.gff3db"


class SeqsDataABC(ABC):
    """interface for genome sequence storage"""

    # the storage reference, e.g. path to file
    source: PathType
    species: str
    mode: str  # as per standard file opening modes, r, w, a
    _is_open = False
    _file: Optional[Any] = None

    def __hash__(self):
        return id(self)

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
    def get_coord_names(self) -> tuple[str]:
        """names of chromosomes / contig"""
        ...

    @abstractmethod
    def close(self):
        """closes the resource"""
        ...


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
        return numpy.array(memoryview(bytearray(b)), dtype=numpy.uint8)


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
        self._file.flush()

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


@define_app
class collect_seqs:
    """given a config and homolog group, loads genome instances on demand
    and extracts sequences"""

    def __init__(self, config: InstalledConfig, make_seq_name: typing.Callable = None):
        self._config = config
        self._genomes = {}
        self._namer = make_seq_name

    def main(self, homologs: homolog_group) -> UnalignedSeqsType:
        namer = self._namer
        seqs = []
        for species, sp_genes in homologs.items():
            if species not in self._genomes:
                self._genomes[species] = load_genome(
                    config=self._config, species=species
                )
            genome = self._genomes[species]
            for name in sp_genes.gene_ids:
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
                seq.name = f"{species}-{name}" if namer is None else namer(feature)
                seq.info["species"] = species
                seq.info["name"] = name
                # disconnect from annotation so the closure of the genome
                # does not cause issues when run in parallel
                seq.annotation_db = None
                seqs.append(seq)

        if not seqs:
            return NotCompleted(
                type="FAIL", origin=self, message=f"no CDS for {homologs}"
            )

        return make_unaligned_seqs(data=seqs, moltype="dna")
