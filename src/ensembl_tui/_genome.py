import dataclasses
import functools
import pathlib
import re
import sys
import typing
import uuid
from abc import ABC, abstractmethod
from typing import Any

import h5py
import numpy
from cogent3 import get_moltype, make_seq
from cogent3.app.composable import define_app
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import (
    OptionalInt,
    OptionalStr,
)
from cogent3.core.sequence import Sequence
from cogent3.util.table import Table
from numpy.typing import NDArray

import ensembl_tui._annotation as eti_annots
from ensembl_tui import _config as eti_config
from ensembl_tui import _species as eti_species
from ensembl_tui import _storage_mixin as eti_storage
from ensembl_tui import _util as eti_util
from ensembl_tui._faster_fasta import quicka_parser

SEQ_STORE_NAME = "genome.seqs-hdf5_blosc2"

_typed_id = re.compile(
    r"\b[a-z]+:",
    flags=re.IGNORECASE,
)  # ensembl stableid's prefixed by the type
_feature_id = re.compile(r"(?<=\bID=)[^;]+")
_exon_id = re.compile(r"(?<=\bexon_id=)[^;]+")
_parent_id = re.compile(r"(?<=\bParent=)[^;]+")
_symbol = re.compile(r"(?<=\bName=)[^;]+")
_description = re.compile(r"(?<=\bdescription=)[^;]+")


def _rename(label: str) -> str:
    return label.split()[0]


@define_app
class fasta_to_hdf5:  # noqa: N801
    def __init__(self, config: eti_config.Config, label_to_name=_rename) -> None:
        self.config = config
        self.label_to_name = label_to_name

    def main(self, db_name: str) -> bool:
        src_dir = self.config.staging_genomes / db_name
        dest_dir = self.config.install_genomes / db_name

        seq_store = SeqsDataHdf5(
            source=dest_dir / SEQ_STORE_NAME,
            species=eti_species.Species.get_species_name(db_name),
            mode="w",
        )

        src_dir = src_dir / "fasta"
        for path in src_dir.glob("*.fa.gz"):
            for label, seq in quicka_parser(path):
                seqid = self.label_to_name(label)
                seq_store.add_record(seq, seqid)
                del seq

        seq_store.close()

        return True


T = tuple[eti_util.PathType, list[tuple[str, str]]]


class SeqsDataABC(ABC):
    """interface for genome sequence storage"""

    # the storage reference, e.g. path to file
    source: eti_util.PathType
    species: str
    mode: str  # as per standard file opening modes, r, w, a
    _is_open = False
    _file: Any | None = None

    @abstractmethod
    def __hash__(self) -> int: ...

    @abstractmethod
    def add_record(self, seq: str, seqid: str) -> None: ...

    @abstractmethod
    def add_records(self, *, records: typing.Iterable[list[str, str]]) -> None: ...

    @abstractmethod
    def get_seq_str(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
    ) -> str: ...

    @abstractmethod
    def get_seq_arr(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
    ) -> NDArray[numpy.uint8]: ...

    @abstractmethod
    def get_coord_names(self) -> tuple[str]: ...

    @abstractmethod
    def close(self) -> None: ...


@define_app
class str2arr:  # noqa: N801
    """convert string to array of uint8"""

    def __init__(self, moltype: str = "dna", max_length: int | None = None) -> None:
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
class arr2str:  # noqa: N801
    """convert array of uint8 to str"""

    def __init__(self, moltype: str = "dna", max_length: int | None = None) -> None:
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
class SeqsDataHdf5(eti_storage.Hdf5Mixin, SeqsDataABC):
    """HDF5 sequence data storage"""

    def __init__(
        self,
        source: eti_util.PathType,
        species: str | None = None,
        mode: str = "r",
        in_memory: bool = False,
    ) -> None:
        # note that species are converted into the Ensembl db prefix
        in_memory = in_memory or "memory" in str(source)
        source = uuid.uuid4().hex if in_memory else source
        self.source = pathlib.Path(source)

        if not in_memory and mode == "r" and not self.source.exists():
            msg = f"{self.source!s} not found"
            raise OSError(msg)

        species = (
            eti_species.Species.get_ensembl_db_prefix(species) if species else None
        )
        self.mode = "w-" if mode == "w" else mode
        h5_kwargs = (
            {
                "driver": "core",
                "backing_store": False,
            }
            if in_memory
            else {}
        )
        try:
            self._file: h5py.File = h5py.File(source, mode=self.mode, **h5_kwargs)
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
            msg = f"{self.source.name!r} {file_species!r} != {species}"
            raise ValueError(msg)
        self.species = self._file.attrs["species"]

    def __hash__(self) -> int:
        return id(self)

    @functools.singledispatchmethod
    def add_record(self, seq: str, seqid: str) -> None:
        seq = self._str2arr(seq)
        self.add_record(seq, seqid)

    @add_record.register
    def _(self, seq: numpy.ndarray, seqid: str) -> None:
        if seqid in self._file:
            stored = self._file[seqid]
            if (seq == stored).all():
                # already seen this seq
                return
            # but it's different, which is a problem
            msg = f"{seqid!r} already present but with different seq"
            raise ValueError(msg)

        self._file.create_dataset(
            name=seqid,
            data=seq,
            chunks=True,
            **eti_util.HDF5_BLOSC2_KWARGS,
        )

    def add_records(self, *, records: typing.Iterable[list[str, str]]) -> None:
        for seqid, seq in records:
            self.add_record(seq, seqid)

    def get_seq_str(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
    ) -> str:
        return self._arr2str(self.get_seq_arr(seqid=seqid, start=start, stop=stop))

    def get_seq_arr(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
    ) -> NDArray[numpy.uint8]:
        if not self._is_open:
            msg = f"{self.source.name!r} is closed"
            raise OSError(msg)

        return self._file[seqid][start:stop]

    def get_coord_names(self) -> tuple[str]:
        """names of chromosomes / contig"""
        return tuple(self._file)


@dataclasses.dataclass(slots=True)
class genome_segment:  # noqa: N801
    species: str
    seqid: str
    start: int
    stop: int
    strand: str
    unique_id: str | None = None

    def __post_init__(self) -> None:
        self.unique_id = (
            eti_util.sanitise_stableid(self.unique_id)
            if self.unique_id
            else f"{self.species}-{self.seqid}-{self.start}-{self.stop}"
        )

    @property
    def source(self) -> str | None:
        return self.unique_id


# TODO: this wrapping class is required for memory efficiency because
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
        annots: eti_annots.Annotations,
    ) -> None:
        self.species = species
        self._seqs = seqs
        self.annotation_db = annots

    @property
    def seqs(self) -> SeqsDataABC:
        return self._seqs

    def get_seq(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        namer: typing.Callable | None = None,
        with_annotations: bool = True,
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
        with_annotations
            assign annotation_db to seq

        Notes
        -----
        Full annotations are bound to the instance.
        """
        seq = self._seqs.get_seq_str(seqid=seqid, start=start, stop=stop)
        if namer:
            name = namer(self.species, seqid, start, stop)
        else:
            name = f"{self.species}:{seqid}:{start}-{stop}"
        # we use seqid to make the sequence here because that identifies the
        # parent seq identity, required for querying annotations
        try:
            seq = make_seq(seq, name=seqid, moltype="dna", annotation_offset=start or 0)
        except TypeError:
            # older version of cogent3
            seq = make_seq(seq, name=seqid, moltype="dna")
            seq.annotation_offset = start or 0

        seq.name = name
        seq.annotation_db = self.annotation_db if with_annotations else None
        return seq

    def get_features(
        self,
        *,
        biotype: OptionalStr = None,
        seqid: OptionalStr = None,
        name: OptionalStr = None,
        start: OptionalInt = None,
        stop: OptionalInt = None,
        limit: OptionalInt = None,
    ) -> typing.Iterable[Feature]:
        for ft in self.annotation_db.get_features_matching(
            biotype=biotype,
            seqid=seqid,
            name=name,
            start=start,
            stop=stop,
            limit=limit,
        ):
            seqid = ft.seqid
            ft.spans = numpy.array(ft.spans)
            start = int(ft.spans.min())
            stop = int(ft.spans.max())
            ft.spans = ft.spans - start
            seq = self.get_seq(
                seqid=seqid,
                start=start,
                stop=stop,
                with_annotations=True,
            )
            # because self.get_seq() automatically names seqs differently
            seq.name = seqid
            yield seq.make_feature(ft)

    def get_cds(
        self,
        *,
        stable_id: str,
        biotype: str = "protein_coding",
    ) -> typing.Iterable[Feature]:
        gene = next(
            iter(
                self.annotation_db.get_features_matching(
                    stable_id=stable_id,
                    biotype=biotype,
                ),
            ),
        )
        cds = self.annotation_db.get_cds(gene=gene)
        seq = self.get_seq(seqid=cds.seqid, start=cds.start, stop=cds.stop)
        cds["spans"] = cds["spans"] - cds.start
        data = dict(cds)
        data.pop("xattr", None)
        # cogent3 currently only handles character data for strand
        data["strand"] = "-" if data.get("strand") == -1 else "+"

        try:
            yield seq.make_feature(feature=data)
        except ValueError:
            msg = f"invalid location data for {cds!r}"
            raise ValueError(msg)

    def get_ids_for_biotype(
        self,
        biotype: str,
        limit: OptionalInt = None,
    ) -> typing.Iterable[str]:
        genes = self.annotation_db.genes
        return genes.get_ids_for_biotype(biotype=biotype, limit=limit)

    def close(self) -> None:
        self._seqs.close()
        self.annotation_db.close()


def load_genome(*, config: eti_config.InstalledConfig, species: str) -> Genome:
    """returns the Genome with bound seqs and features"""
    genome_path = config.installed_genome(species) / SEQ_STORE_NAME
    seqs = SeqsDataHdf5(source=genome_path, species=species, mode="r")
    ann = eti_annots.Annotations(source=config.installed_genome(species))
    return Genome(species=species, seqs=seqs, annots=ann)


def get_seqs_for_ids(
    *,
    config: eti_config.InstalledConfig,
    species: str,
    names: list[str],
    make_seq_name: typing.Callable | None = None,
) -> typing.Iterable[Sequence]:
    genome = load_genome(config=config, species=species)
    # is it possible to do batch query for all names?
    for name in names:
        cds = list(
            genome.get_cds(
                stable_id=name,
                biotype="protein_coding",
            ),
        )
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


def load_annotations_for_species(*, path: pathlib.Path) -> eti_annots.Annotations:
    """returns the annotation Db for species"""
    if not path.exists():
        eti_util.print_colour(
            text=f"{path.name!r} is missing",
            colour="red",
        )
        sys.exit(1)
    try:
        return eti_annots.Annotations(source=path)
    except FileNotFoundError:
        eti_util.print_colour(
            text=f"expected files not in {str(path)!r}",
            colour="red",
        )
        sys.exit(1)


def _get_all_gene_segments(
    *,
    annot_db: eti_annots.Annotations,
    limit: int | None,
    biotype: str | None,
) -> list[eti_annots.GeneData]:
    return list(annot_db.get_features_matching(biotype=biotype, limit=limit))


def _get_selected_gene_segments(
    *,
    annot_db: eti_annots.Annotations,
    limit: int | None,
    stableids: list[str],
    biotype: str | None,
) -> list[eti_annots.GeneData]:
    result = []
    for stable_id in stableids:
        record = list(
            annot_db.get_features_matching(
                biotype=biotype,
                stable_id=stable_id,
                limit=limit,
            ),
        )
        result.extend(record)
    return result


def get_gene_segments(
    *,
    annot_db: eti_annots.Annotations,
    limit: int | None = None,
    species: str | None = None,
    stableids: list[str] | None = None,
    biotype: str = "protein_coding",
) -> list[genome_segment]:
    """return genome segment information for genes

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
    records = (
        _get_selected_gene_segments(
            annot_db=annot_db,
            limit=limit,
            stableids=stableids,
            biotype=biotype,
        )
        if stableids
        else _get_all_gene_segments(annot_db=annot_db, limit=limit, biotype=biotype)
    )
    for i, record in enumerate(records):
        segment = genome_segment(
            species=species,
            start=record["start"],
            stop=record["stop"],
            strand=record["strand"],
            seqid=record["seqid"],
            unique_id=record["name"],
        )
        records[i] = segment
    return records


def get_gene_table_for_species(
    *,
    annot_db: eti_annots.Annotations,
    limit: int | None = None,
) -> Table:
    """
    returns gene data from Annotations

    Parameters
    ----------
    annot_db
        feature db
    limit
        limit number of records to
    """
    table = annot_db.genes.gene_table
    if limit:
        table = table[:limit]
    return table


def get_species_summary(
    *,
    annot_db: eti_annots.Annotations,
    species: str | None = None,
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
    # for now, just biotype
    species = species or annot_db.source.parent.name
    counts = annot_db.biotypes.count_distinct
    try:
        common_name = eti_species.Species.get_common_name(species)
    except ValueError:
        common_name = species

    counts.title = f"{common_name} features"
    counts.format_column("count", lambda x: f"{x:,}")
    return counts
