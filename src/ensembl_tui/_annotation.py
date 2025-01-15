import dataclasses
import functools
import pathlib
import typing

import cogent3
import duckdb
import numpy
from cogent3.core.annotation_db import (
    AnnotationDbABC,
    FeatureDataType,
)

import ensembl_tui._mysql_core_attr as core_tables
from ensembl_tui import _storage_mixin as eti_storage

if typing.TYPE_CHECKING:
    from cogent3.util.table import Table

OptInt = int | None
OptStr = str | None
OptBool = bool | None
StrOrBool = str | bool
FeatureDictVals = str | int | numpy.ndarray

GENE_ATTR_COLUMNS = (
    "stable_id",
    "biotype",
    "seqid",
    "start",
    "stop",
    "strand",
    "canonical_transcript_id",
    "symbol",
    "gene_id",
    "description",
)


class FeatureDataMixin:  # supports getitem as a dict on properties
    def __getitem__(self, key: str) -> FeatureDictVals:
        return getattr(self, key)

    def __setitem__(self, key: str, value: FeatureDictVals) -> None:
        setattr(self, key, value)

    def pop(self, key: str) -> FeatureDictVals:
        value = getattr(self, key)
        delattr(self, key)
        return value

    def get(self, key: str, default: typing.Any = None) -> FeatureDictVals:
        return getattr(self, key, default)


@dataclasses.dataclass(slots=True)
class FeatureDataBase(FeatureDataMixin):
    seqid: str = dataclasses.field(kw_only=True)
    start: int = dataclasses.field(kw_only=True)
    stop: int = dataclasses.field(kw_only=True)
    spans: numpy.ndarray[numpy.int32] = dataclasses.field(kw_only=True)
    strand: int = dataclasses.field(kw_only=True)
    name: str = dataclasses.field(kw_only=True)
    biotype: str | None = dataclasses.field(kw_only=True, default=None)

    def __dict__(self) -> dict:
        return dataclasses.asdict(self)

    def __iter__(self):  # noqa: ANN204
        feature_fields = {"seqid", "biotype", "name", "spans", "strand"}
        xattr = {}
        seen_name = False
        for field in dataclasses.fields(self):
            if field.name not in feature_fields:
                xattr[field.name] = getattr(self, field.name)
                continue
            yield field.name, getattr(self, field.name)
            if field.name == "name":
                seen_name = True

        if not seen_name:
            yield "name", self.name

        # for now, we don't return xattr until cogent3 Feature
        # class supports it
        # yield "xattr", xattr


@dataclasses.dataclass(slots=True)
class GeneData(FeatureDataBase):
    canonical_transcript_id: int
    stable_id: str
    gene_id: str
    symbol: str
    name: str | None = dataclasses.field(kw_only=True, default=None)
    description: str | None = dataclasses.field(kw_only=True, default=None)

    def __post_init__(self) -> None:
        self.name = self.stable_id


@dataclasses.dataclass(slots=True)
class TranscriptData(FeatureDataBase):
    transcript_id: int
    stable_id: str
    gene_stable_id: str
    name: str | None = dataclasses.field(kw_only=True, default=None)
    gene_id: int | None = dataclasses.field(kw_only=True, default=None)
    symbol: str | None = dataclasses.field(kw_only=True, default=None)

    def __post_init__(self) -> None:
        self.name = self.stable_id


@dataclasses.dataclass(slots=True)
class CdsData(FeatureDataBase):
    stable_id: str
    gene_stable_id: str
    name: str | None = dataclasses.field(kw_only=True, default=None)
    transcript_id: int | None = dataclasses.field(kw_only=True, default=None)

    def __post_init__(self) -> None:
        self.name = self.stable_id


@dataclasses.dataclass(slots=True)
class RepeatData(FeatureDataBase):
    repeat_type: str
    repeat_class: str
    repeat_name: str


def _matching_conditions(
    equals_conds: dict[str, str | int] | None = None,
    like_conds: dict[str, str] | None = None,
    allow_partial: bool = True,
) -> str:
    """creates WHERE clause

    Parameters
    ----------
    equals_conds
        column name and values to be matched by equals
    like_conds
        column name and values to be matched by ILIKE (case-insensitive)
    allow_partial
        if False, only records within start, stop are included. If True,
        all records that overlap the segment defined by start, stop are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """
    equals_conds = equals_conds or {}
    start = equals_conds.pop("start", None)
    stop = equals_conds.pop("stop", None)

    sql = []
    if equals_conds:
        conds = []
        for col, val in equals_conds.items():
            # conditions are filtered for None before here, so we should add
            # an else where the op is assigned !=
            if isinstance(val, tuple | set | list):
                vals = ",".join(f"{v!r}" for v in val)
                conds.append(f"{col} IN ({vals})")
            elif val is not None:
                conds.append(f"{col} = {val!r}")

        sql.append(" AND ".join(conds))
    if like_conds:
        sql.extend(f"{col} ILIKE '%{val}%'" for col, val in like_conds.items())
    if start is not None and stop is not None:
        if allow_partial:
            # allow matches that overlap the segment
            cond = [
                f"(start >= {start} AND stop <= {stop})",  # lies within the segment
                f"(start <= {start} AND stop > {start})",  # straddles beginning of segment
                f"(start < {stop} AND stop >= {stop})",  # straddles stop of segment
                f"(start <= {start} AND stop >= {stop})",  # includes segment
            ]
            cond = " OR ".join(cond)
        else:
            # only matches within bounds
            cond = f"start >= {start} AND stop <= {stop}"
        sql.append(f"({cond})")
    elif start is not None:
        # if query has no stop, then any feature containing start
        cond = f"(start <= {start} AND {start} < stop)"
        sql.append(f"({cond})")
    elif stop is not None:
        # if query has no start, then any feature containing stop
        cond = f"(start <= {stop} AND {stop} < stop)"
        sql.append(f"({cond})")

    return " AND ".join(sql)


def _select_records_sql(
    *,
    table_name: str,
    equals_conds: dict[str, str | int] | None = None,
    like_conds: dict[str, str] | None = None,
    columns: typing.Sequence[str] | None = None,
    allow_partial: bool = True,
) -> str:
    """create SQL select statement and values

    Parameters
    ----------
    table_name
        containing the data to be selected from
    columns
        values to select
    equals_conds
        the WHERE condition = value
    like_conds
        the WHERE condition ILIKE '%value%'
    start, stop
        select records whose (start, stop) values lie between start and stop,
        or overlap them if (allow_partial is True)
    allow_partial
        if False, only records within start, stop are included. If True,
        all records that overlap the segment defined by start, stop are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """

    conditions = _matching_conditions(
        equals_conds=equals_conds,
        like_conds=like_conds,
        allow_partial=allow_partial,
    )
    cols = "*" if columns is None else f"{', '.join(columns)}"
    sql = f"SELECT {cols} FROM {table_name}"
    return f"{sql} WHERE {conditions}" if conditions else sql


@dataclasses.dataclass
class BiotypeView(eti_storage.DuckdbParquetBase, eti_storage.ViewMixin):
    _tables: tuple[str] = ("gene_attr",)

    def num_records(self) -> int:
        """returns the number of distinct biotypes"""
        return len(self.distinct)

    @functools.cached_property
    def distinct(self) -> tuple[str, ...]:
        sql = f"SELECT DISTINCT biotype FROM {self._tables[0]}"
        return tuple(r[0] for r in self.conn.sql(sql).fetchall())

    @functools.cached_property
    def count_distinct(self) -> "Table":
        sql = (
            f"SELECT biotype, COUNT(*) AS freq FROM {self._tables[0]} GROUP BY biotype"
        )
        got = self.conn.sql(sql).fetchall()
        return cogent3.make_table(
            header=["biotype", "count"],
            data=got,
            index_name="biotype",
        )


@dataclasses.dataclass
class GeneView(eti_storage.DuckdbParquetBase, eti_storage.ViewMixin):
    _tables: tuple[str, str] = ("gene_attr", "transcript_attr")

    def num_records(self) -> int:
        """returns the number of distinct genes as identified by stable_id"""
        sql = "SELECT DISTINCT stable_id FROM gene_attr"
        return len(self.conn.sql(sql).fetchall())

    def get_features_matching(
        self,
        *,
        seqid: OptStr = None,
        biotype: OptStr = None,
        stable_id: OptStr = None,
        start: OptInt = None,
        stop: OptInt = None,
        strand: OptStr = None,
        symbol: OptStr = None,
        description: OptStr = None,
        **kwargs,  # noqa: ANN003
    ) -> typing.Iterator[FeatureDataBase]:
        # add supoport for querying by symbol and description
        limit = kwargs.pop("limit", None)
        local_vars = locals()
        if kwargs := {
            k: v
            for k, v in local_vars.items()
            if k not in ("self", "kwargs", "columns", "limit") and v is not None
        }:
            like_conds = (
                {"description": kwargs.pop("description")} if description else None
            )
            sql = _select_records_sql(
                equals_conds=kwargs,
                like_conds=like_conds,
                table_name="gene_attr",
                columns=GENE_ATTR_COLUMNS,
            )
        else:
            sql = f"SELECT {','.join(GENE_ATTR_COLUMNS)} FROM gene_attr"

        sql += f" LIMIT {limit}" if limit else ""

        for record in self.conn.sql(sql).fetchall():
            data = dict(zip(GENE_ATTR_COLUMNS, record, strict=True))
            start, stop = (
                data.get("start"),
                data.get("stop"),
            )
            data["spans"] = numpy.array([sorted([start, stop])], dtype=numpy.int32)  # type: ignore
            yield GeneData(**data)

    def get_by_stable_id(self, stable_id: str) -> typing.Iterator[GeneData]:
        yield from self.get_features_matching(stable_id=stable_id)

    def get_by_symbol(self, symbol: str) -> typing.Iterator[GeneData]:
        yield from self.get_features_matching(symbol=symbol)

    def get_by_description(self, description: str) -> typing.Iterator[GeneData]:
        yield from self.get_features_matching(description=description)

    def get_cds(self, *, gene: GeneData) -> CdsData:
        # for now, we only support getting the canonical transcript
        transcript_id = gene["canonical_transcript_id"]
        columns = (
            "transcript_id",
            "seqid",
            "start",
            "stop",
            "strand",
            "cds_spans",
            "cds_stable_id",
        )
        sql = f"SELECT {','.join(columns)} FROM transcript_attr WHERE transcript_id = {transcript_id}"
        if not (record := self.conn.sql(sql).fetchone()):
            msg = f"No CDS spans found for {gene=}"
            raise ValueError(msg)

        transcript = dict(zip(columns, record, strict=True))
        if not (spans := transcript.pop("cds_spans", None)):
            msg = f"No CDS spans found for {gene=}"
            raise ValueError(msg)

        spans = core_tables.blob_to_array(spans)
        stable_id = transcript.pop("cds_stable_id")
        return CdsData(
            **{
                **transcript,
                "spans": spans,
                "stable_id": stable_id,
                "gene_stable_id": gene.stable_id,
            },
        )

    @functools.singledispatchmethod
    def get_feature_children(
        self,
        feature: FeatureDataBase,
    ) -> typing.Iterator[FeatureDataBase]:
        msg = f"{type(feature)=} not supported"
        raise NotImplementedError(msg)

    @get_feature_children.register(GeneData)
    def _(self, gene: GeneData) -> typing.Iterator[TranscriptData]:
        columns = (
            "transcript_id",
            "seqid",
            "start",
            "stop",
            "strand",
            "transcript_spans",
            "transcript_stable_id",
            "gene_id",
        )
        sql = f"SELECT {','.join(columns)} FROM transcript_attr WHERE transcript_id = {gene['canonical_transcript_id']}"
        for record in self.conn.sql(sql).fetchall():
            transcript = dict(zip(columns, record, strict=True))
            spans = transcript.pop("transcript_spans")
            spans = core_tables.blob_to_array(spans)
            stable_id = transcript.pop("transcript_stable_id")
            yield TranscriptData(
                **{
                    **transcript,
                    "spans": spans,
                    "biotype": gene.biotype,
                    "stable_id": stable_id,
                    "gene_stable_id": gene.stable_id,
                    "symbol": gene.symbol,
                },
            )

    @get_feature_children.register
    def _(self, transcript: TranscriptData) -> typing.Iterator[CdsData]:
        columns = (
            "transcript_id",
            "seqid",
            "start",
            "stop",
            "strand",
            "cds_spans",
            "cds_stable_id",
        )
        sql = f"SELECT {','.join(columns)} FROM transcript_attr WHERE transcript_id = {transcript['transcript_id']}"
        if not (record := self.conn.sql(sql).fetchone()):
            msg = f"No CDS spans found for {transcript=}"
            raise ValueError(msg)

        data = dict(zip(columns, record, strict=True))
        if not (spans := data.pop("cds_spans", None)):
            msg = f"No CDS spans found for {transcript=}"
            raise ValueError(msg)

        spans = core_tables.blob_to_array(spans)
        stable_id = data.pop("cds_stable_id")
        yield CdsData(
            **{
                **data,
                "spans": spans,
                "stable_id": stable_id,
                "gene_stable_id": transcript.gene_stable_id,
            },
        )

    @functools.singledispatchmethod
    def get_feature_parent(self, feature: FeatureDataType) -> FeatureDataBase:
        msg = f"type {feature=} has no parents"
        raise ValueError(msg)

    @get_feature_parent.register
    def _(self, transcript: TranscriptData) -> GeneData:
        sql = f"SELECT {','.join(GENE_ATTR_COLUMNS)} FROM gene_attr WHERE gene_id = {transcript['gene_id']}"
        if not (record := self.conn.sql(sql).fetchone()):
            msg = f"No gene spans found for {transcript=}"
            raise ValueError(msg)

        gene = dict(zip(GENE_ATTR_COLUMNS, record, strict=True))
        spans = numpy.array([sorted([gene["start"], gene["stop"]])], dtype=numpy.int32)
        stable_id = gene.pop("stable_id")
        return GeneData(**{**gene, "spans": spans, "stable_id": stable_id})

    @get_feature_parent.register
    def _(self, cds: CdsData) -> TranscriptData:
        columns = (
            "transcript_id",
            "seqid",
            "start",
            "stop",
            "strand",
            "transcript_spans",
            "transcript_stable_id",
        )
        sql = f"SELECT {','.join(columns)} FROM transcript_attr WHERE transcript_id = {cds['transcript_id']}"
        if not (record := self.conn.sql(sql).fetchone()):
            msg = f"No transcript spans found for {cds=}"
            raise ValueError(msg)

        transcript = dict(zip(columns, record, strict=True))
        spans = transcript.pop("transcript_spans")
        spans = core_tables.blob_to_array(spans)
        stable_id = transcript.pop("transcript_stable_id")
        return TranscriptData(
            **{
                **transcript,
                "spans": spans,
                "biotype": cds["biotype"],
                "stable_id": stable_id,
                "gene_stable_id": cds.gene_stable_id,
            },
        )

    def count_distinct(
        self,
        *,
        seqid: StrOrBool = False,
        biotype: OptBool = False,
    ) -> "Table | None":
        if not any((seqid, biotype)):
            return None

        local_vars = locals()
        if constraints := {k: v for k, v in local_vars.items() if isinstance(v, str)}:
            where_clause = f"WHERE {_matching_conditions(equals_conds=constraints)}"
        else:
            where_clause = ""

        header = [c for c in ("biotype", "seqid") if local_vars[c]]
        sql = (
            f"SELECT {', '.join(header)}, COUNT(*) as count FROM gene_attr"
            f" {where_clause} GROUP BY {', '.join(header)};"
        )
        return cogent3.make_table(
            header=[*header, "count"],
            data=self.conn.sql(sql).fetchall(),
            column_templates={"count": lambda x: f"{x:,}"},
        )

    def get_ids_for_biotype(
        self,
        biotype: str,
        seqid: str | None = None,
        limit: int | None = None,
    ) -> list[str]:
        sql = "SELECT stable_id from gene_attr WHERE biotype=?"
        val = (biotype,)
        if seqid:
            sql += " AND seqid=?"
            val = (*val, seqid)
        if limit:
            sql += " LIMIT ?"
            val = (*val, limit)
        return [r[0] for r in self.conn.sql(sql, params=val).fetchall()]

    @functools.cached_property
    def gene_table(self) -> "Table":
        """return a Table with all gene data"""
        columns = (
            "species",
            "name",
            "seqid",
            "source",
            "biotype",
            "start",
            "stop",
            "strand",
            "symbol",
            "description",
        )
        rows = []
        rows.extend(
            [self.species] + [record.get(c, None) for c in columns[1:]]
            for record in self.get_features_matching()
        )
        header = ["stableid" if c == "name" else c for c in columns]
        return cogent3.make_table(header=header, data=rows)


@dataclasses.dataclass
class RepeatView(eti_storage.DuckdbParquetBase, eti_storage.ViewMixin):
    _tables: tuple[str, ...] = tuple(
        core_tables.collect_table_names(*core_tables.REPEAT_ATTR),
    )

    @property
    def conn(self) -> duckdb.DuckDBPyConnection:
        if self._conn is None:
            super().conn
            sql = """CREATE VIEW IF NOT EXISTS repeat_view AS
                    SELECT 
                        rc.repeat_type AS repeat_type,
                        sr.name AS seqid,
                        rc.repeat_class AS repeat_class,
                        rc.repeat_name AS repeat_name,
                        rf.seq_region_start AS start,
                        rf.seq_region_end AS stop,
                        rf.seq_region_strand AS strand
                    FROM repeat_consensus rc
                    JOIN repeat_feature rf ON rc.repeat_consensus_id = rf.repeat_consensus_id
                    JOIN seq_region sr ON rf.seq_region_id = sr.seq_region_id
                    """
            self._conn.sql(sql)
        return self._conn

    def num_records(self) -> int:
        """returns the number of rows in repeat_feature"""
        sql = "SELECT COUNT(*) FROM repeat_feature"
        return self.conn.sql(sql).fetchone()[0]

    def get_features_matching(
        self,
        *,
        seqid: OptStr = None,
        biotype: OptStr = None,
        name: OptStr = None,
        start: OptInt = None,
        stop: OptInt = None,
        strand: OptStr = None,
        repeat_type: OptStr = None,
        repeat_class: OptStr = None,
        **kwargs,  # noqa: ANN003
    ) -> typing.Iterator[FeatureDataType]:
        limit = kwargs.pop("limit", None)
        local_vars = locals()
        if kwargs := {
            k: v
            for k, v in local_vars.items()
            if k not in ("self", "kwargs", "limit") and v is not None
        }:
            like_cols = "repeat_type", "repeat_class", "repeat_name"
            like_conds = {k: v for k, v in kwargs.items() if k in like_cols}
            equals_conds = {k: v for k, v in kwargs.items() if k not in like_cols}
            columns = "seqid", "start", "stop", "strand"
            sql = _select_records_sql(
                table_name="repeat_view",
                equals_conds=equals_conds,
                like_conds=like_conds,
                columns=columns + like_cols,
            )
            sql += f" LIMIT {limit}" if limit else ""
        else:
            columns = (
                "seqid",
                "start",
                "stop",
                "strand",
                "repeat_type",
                "repeat_class",
                "repeat_name",
            )
            sql = "SELECT * FROM repeat_view LIMIT 10"

        for record in self.conn.sql(sql).fetchall():
            data = dict(zip(columns, record, strict=False))
            spans = numpy.array(
                [(data.pop("start"), data.pop("stop"))],
                dtype=numpy.int32,
            )
            data["spans"] = spans
            data["biotype"] = "repeat"
            yield FeatureDataType(**data)

    def get_children_matching(self, **kwargs):
        return ()

    def count_distinct(
        self,
        seqid: StrOrBool = False,
        repeat_type: OptBool = False,
        repeat_class: OptBool = False,
    ) -> "Table | None":
        if not any((seqid, repeat_type, repeat_class)):
            return None

        local_vars = locals()
        if constraints := {k: v for k, v in local_vars.items() if isinstance(v, str)}:
            where_clause = f"WHERE {_matching_conditions(equals_conds=constraints)}"
        else:
            where_clause = ""

        header = [c for c in ("seqid", "repeat_type", "repeat_class") if local_vars[c]]
        sql = (
            f"SELECT {', '.join(header)}, COUNT(*) as count FROM repeat_view"
            f" {where_clause} GROUP BY {', '.join(header)};"
        )
        return cogent3.make_table(
            header=[*header, "count"],
            data=self.conn.sql(sql).fetchall(),
            column_templates={"count": lambda x: f"{x:,}"},
        )


@dataclasses.dataclass
class Annotations(AnnotationDbABC, eti_storage.ViewMixin):
    """virtual genome annotation database that provides access to gene and repeat features"""

    source: dataclasses.InitVar[pathlib.Path | str]
    _source: pathlib.Path = dataclasses.field(init=False)
    biotypes: BiotypeView | None = dataclasses.field(init=True, default=None)
    genes: GeneView | None = dataclasses.field(init=True, default=None)
    repeats: RepeatView | None = dataclasses.field(init=True, default=None)

    def __post_init__(
        self,
        source: pathlib.Path,
    ) -> None:
        source = pathlib.Path(source)
        self._source = source
        if source.is_dir():
            self.biotypes = BiotypeView(source=source)
            self.genes = GeneView(source=source)
            self.repeats = RepeatView(source=source)

    @property
    def source(self) -> pathlib.Path:
        return self._source

    def get_features_matching(
        self,
        *,
        biotype: str,
        **kwargs,
    ) -> typing.Iterator[FeatureDataType]:
        view = self.repeats if biotype == "repeat" else self.genes
        if biotype != "repeat":
            kwargs["biotype"] = biotype
        yield from view.get_features_matching(**kwargs)

    def __len__(self) -> int:
        return self.num_records()

    def add_records(self, **kwargs):
        raise NotImplementedError

    def add_feature(self, **kwargs):
        raise NotImplementedError

    def get_feature_children(self, **kwargs):
        raise NotImplementedError

    def get_feature_parent(self, **kwargs):
        raise NotImplementedError

    def num_matches(self, **kwargs):
        raise NotImplementedError

    def subset(self, **kwargs):
        # this should return a cogent3 BasicAnnotationDb
        # or we need a AnnotationsUnion class that
        # encloses the individual Annotations instances
        # to avoid copying
        raise NotImplementedError

    def union(self, **kwargs):
        # this should return a cogent3 BasicAnnotationDb
        # or we need a AnnotationsUnion class that
        # encloses the individual Annotations instances
        # to avoid copying
        raise NotImplementedError

    def update(self, **kwargs):
        # this should return a cogent3 BasicAnnotationDb
        # or we need a AnnotationsUnion class that
        # encloses the individual Annotations instances
        # to avoid copying
        raise NotImplementedError

    def to_json(self) -> str:
        raise NotImplementedError

    def to_rich_dict(self) -> dict:
        raise NotImplementedError

    def from_dict(self, data: dict[str, typing.Any]) -> None:
        raise NotImplementedError

    def get_cds(self, **kwargs) -> CdsData:  # noqa: ANN003
        return self.genes.get_cds(**kwargs)

    def get_ids_for_biotype(self, biotype: str, limit: int | None = None) -> list[str]:
        return self.genes.get_ids_for_biotype(biotype=biotype, limit=limit)

    def count_distinct(self, **kwargs) -> "Table | None":
        return None if self.genes is None else self.genes.count_distinct(**kwargs)

    def num_records(self) -> int:
        """returns the total number of genes and repeat features"""
        num_genes = 0 if self.genes is None else self.genes.num_records()
        num_repeats = 0 if self.repeats is None else self.repeats.num_records()
        return num_genes + num_repeats

    def close(self) -> None:
        self.biotypes.close()
        self.genes.close()
        self.repeats.close()
