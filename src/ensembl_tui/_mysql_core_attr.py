import collections
import dataclasses
import functools
import io
import pathlib
import typing

import duckdb
import numpy

from ensembl_tui import _util as eti_util

TRANSCRIPT_ATTR_SCHEMA = (
    "gene_id INTEGER",
    "transcript_id INTEGER",
    "seqid TEXT",
    "start INTEGER",
    "stop INTEGER",
    "strand TINYINT",
    "transcript_spans BLOB",
    "cds_spans BLOB",
    "transcript_stable_id TEXT",
    "cds_stable_id TEXT",
)
TRANSCRIPT_ATTR_COLS = eti_util.make_column_constant(TRANSCRIPT_ATTR_SCHEMA)
GENE_ATTR_SCHEMA = (
    "gene_id INTEGER",
    "stable_id TEXT",
    "biotype TEXT",
    "canonical_transcript_id INTEGER",
    "seqid TEXT",
    "start INTEGER",
    "stop INTEGER",
    "strand TINYINT",
    "symbol TEXT",
)
GENE_ATTR_COLS = eti_util.make_column_constant(GENE_ATTR_SCHEMA)


def array_to_blob(data: numpy.ndarray) -> bytes:
    """writes a numpy array to bytes using numpy.save

    Notes
    -----
    This is robust across platforms as it saves both the type and shape
    """
    with io.BytesIO() as out:
        numpy.save(out, data)
        out.seek(0)
        return out.read()


def blob_to_array(data: bytes) -> numpy.ndarray:
    """reverses the process of array_to_blob"""
    with io.BytesIO(data) as out:
        out.seek(0)
        return numpy.load(out)


# https://asia.ensembl.org/info/docs/api/core/core_schema.html


def collect_table_names(*args: dict[str, str]) -> set[str]:
    """
    returns the set of table names from the provided dictionaries
    """
    tables = set()

    for group in args:
        tables |= set(group.values())

    return tables


@functools.cache
def get_all_tables() -> set[str]:
    return collect_table_names(
        *[v for k, v in globals().items() if k.endswith("_attrs")],
    )


def load_db(db_name: pathlib.Path, table_names: set[str]) -> duckdb.DuckDBPyConnection:
    con = duckdb.connect(":memory:")

    # Import each Parquet file into a table of the same name
    for table_name in table_names:
        parquet_file = db_name / f"{table_name}.parquet"
        con.execute(
            f"CREATE TABLE {table_name} AS SELECT * FROM read_parquet('{parquet_file}')",
        )

    return con


location_attrs = {
    # "coord_names": "coord_system",
    "location": "seq_region",
    # "meta": "meta",
}
gene_attrs = {
    "stableid": "gene",
    "symbol": "xref",
    "description": "gene",
    "biotype": "gene",
    "location": "gene",
    "canonical_transcript": "gene",
}
transcript_attrs = {
    "stableid": "transcript",
    "location": "transcript",
    "status": "transcript",
    "translated_exons": "translation",
}
exon_attrs = {"stableid": "exon", "location": "exon", "transcript": "exon_transcript"}
repeat_attrs = {
    "symbol": "repeat_consensus",
    "repeat_type": "repeat_consensus",
    "repeat_class": "repeat_consensus",
    "consensus": "repeat_consensus",
    "location": "repeat_feature",
}

GENE_ATTR = location_attrs, gene_attrs
TRANSCRIPT_ATTR = location_attrs, transcript_attrs
FULL_GENES_ATTR = location_attrs, gene_attrs, transcript_attrs, exon_attrs
REPEAT_ATTR = location_attrs, repeat_attrs


@functools.cache
def make_mysqldump_names() -> list[str]:
    """makes expected names for the mysqldump files"""
    all_tables = [f"{table}.txt.gz" for table in get_all_tables()]
    all_tables.insert(0, "CHECKSUMS")
    return all_tables


# we integrate transcript groups of tables into the following:
# translation, transcription, exon, exon_transcript and translation tables
# into a single table such that the coordinates of multiple exons are stored as a
# binary blob in a single field -- transcript_spans -- and the coordinates for
# the translated exons in a single field -- cds_spans. This will be stored as
# transcribe_translate.parquet


@dataclasses.dataclass(slots=True)
class LimitExons:
    """stores attr concerning the first and last translated exons"""

    start_rank: int
    stop_rank: int
    rel_start: int
    rel_stop: int
    strand: int
    transcript_id: int


def get_all_limit_exons(
    conn: duckdb.DuckDBPyConnection,
) -> dict[int, list[tuple[int, ...]]]:
    """returns the limiting exons for transcript_id"""
    sql = """SELECT
        t.transcript_id,
        ev.exon_id,
        t.start_exon_id,
        t.end_exon_id,
        ev.rank,
        t.seq_start,
        t.seq_end,
        ev.strand
    FROM exon_view ev
    JOIN translation t ON ev.transcript_id = t.transcript_id
    AND (ev.exon_id = t.start_exon_id OR ev.exon_id = t.end_exon_id);
    """
    limit_exons = collections.defaultdict(list)

    for record in conn.sql(sql).fetchall():
        limit_exons[record[0]].append(record)
    return limit_exons


def get_limit_exons(records: list[tuple[int, ...]]) -> LimitExons:
    if not records:
        msg = "no records"
        raise ValueError(msg)

    transcript_id = records[0][0]
    if len(records) > 2:
        msg = f"too many entries in translation table for {transcript_id}\n{records!r}"
        raise ValueError(msg)

    # the folowing handles cases where a single exon is both start and end
    if records[0][1] == records[0][2]:
        start_exon, end_exon = records[0], records[-1]
    else:
        start_exon, end_exon = records[-1], records[0]

    start_rank = start_exon[4]
    end_rank = end_exon[4]
    rel_start = start_exon[5]
    rel_end = end_exon[6]
    strand = start_exon[7]
    return LimitExons(
        start_rank=start_rank,
        stop_rank=end_rank,
        # subtract 1 to convert to 0-based
        rel_start=rel_start - 1,
        rel_stop=rel_end,
        strand=strand,
        transcript_id=transcript_id,
    )


@dataclasses.dataclass(slots=True)
class TranscriptAttrRecord:
    """one data record for the transcribe table"""

    transcript_id: int
    gene_id: int
    seqid: str
    strand: int
    transcript_spans: numpy.ndarray
    cds_spans: numpy.ndarray | None
    transcript_stable_id: str
    cds_stable_id: str

    @property
    def start(self) -> int:
        return self.transcript_spans.min()

    @property
    def stop(self) -> int:
        return self.transcript_spans.max()

    def to_record(self, columns: tuple[str]) -> tuple:
        cds_blob = array_to_blob(self.cds_spans) if self.cds_spans is not None else None
        mapping = {
            "transcript_id": self.transcript_id,
            "gene_id": self.gene_id,
            "seqid": self.seqid,
            "start": int(self.start),
            "stop": int(self.stop),
            "strand": int(self.strand),
            "transcript_spans": array_to_blob(self.transcript_spans),
            "cds_spans": cds_blob,
            "transcript_stable_id": self.transcript_stable_id,
            "cds_stable_id": self.cds_stable_id,
        }
        return tuple(mapping[c] for c in columns)


def get_transcript_attr_records(
    conn: duckdb.DuckDBPyConnection,
) -> typing.Iterator[TranscriptAttrRecord]:
    """returns a generator of TranscriptAttrRecord"""
    # we use SQL aggregate functions followed by numpy.fromstring to
    # greatly speedup extraction of all exon coords
    sql = """SELECT
    transcript_id, gene_id, strand, seqid, transcript_stable_id, cds_stable_id,
    STRING_AGG(CAST(start AS VARCHAR), ' ') AS agg_start,
    STRING_AGG(CAST(stop AS VARCHAR), ' ') AS agg_stop,
    STRING_AGG(CAST(rank AS VARCHAR), ' ') AS agg_rank
    FROM exon_view
    GROUP BY transcript_id, gene_id, strand, seqid, transcript_stable_id, cds_stable_id
    """
    limit_exons = get_all_limit_exons(conn)
    for (
        transcript_id,
        gene_id,
        strand,
        seqid,
        transcript_stable_id,
        cds_stable_id,
        agg_start,
        agg_stop,
        agg_rank,
    ) in conn.sql(sql).fetchall():
        # Note that the adjustment of start to be 0-based has already been
        # done during the mysqldump import
        starts = numpy.fromstring(agg_start, sep=" ", dtype=numpy.int32)
        stops = numpy.fromstring(agg_stop, sep=" ", dtype=numpy.int32)
        ranks = numpy.fromstring(agg_rank, sep=" ", dtype=numpy.int32)
        # make the transcript exon spans, in rank order
        # to facilitate getting the cds spans
        transcript_spans = numpy.empty((starts.size, 2), dtype=numpy.int32)
        for i, rank in enumerate(ranks):
            transcript_spans[rank - 1] = (starts[i], stops[i])

        if transcript_id not in limit_exons:
            # no translated exons
            yield TranscriptAttrRecord(
                seqid=seqid,
                transcript_id=transcript_id,
                gene_id=gene_id,
                strand=strand,
                transcript_spans=transcript_spans,
                cds_spans=None,
                transcript_stable_id=transcript_stable_id,
                cds_stable_id=cds_stable_id,
            )
            continue

        lex = get_limit_exons(limit_exons[transcript_id])
        start_index = lex.start_rank - 1
        stop_index = lex.stop_rank - 1
        # adjust the start and end using the limiting exons
        # the rel_start and rel_stop are BOTH relative to the
        # 5' end of an exon
        # so the start_exon coords become (exon_start + rel_start, exon_end)
        # the end_exon coords become (exon_start, exon_start + rel_stop)
        cds_spans = transcript_spans.copy()
        start_exon_coords = cds_spans[start_index]
        stop_exon_coords = cds_spans[stop_index]
        if lex.strand == 1:
            start_exon_coords = (
                start_exon_coords[0] + lex.rel_start,
                start_exon_coords[1],
            )
            stop_exon_coords = (
                stop_exon_coords[0],
                stop_exon_coords[0] + lex.rel_stop,
            )
        else:
            start_exon_coords = (
                start_exon_coords[0],
                start_exon_coords[1] - lex.rel_start,
            )
            stop_exon_coords = (
                stop_exon_coords[1] - lex.rel_stop,
                stop_exon_coords[1],
            )

        # update the limit exon coordinates
        cds_spans[start_index] = start_exon_coords
        cds_spans[stop_index] = stop_exon_coords
        cds_spans = cds_spans[start_index : stop_index + 1]

        # sort all spans in ascending numerical order
        # note that the lexsort returns the sorted indices
        cds_spans = cds_spans[numpy.lexsort(cds_spans.T), :]
        transcript_spans = transcript_spans[numpy.lexsort(transcript_spans.T), :]

        yield TranscriptAttrRecord(
            seqid=seqid,
            transcript_id=transcript_id,
            gene_id=gene_id,
            strand=strand,
            transcript_spans=transcript_spans,
            cds_spans=cds_spans,
            transcript_stable_id=transcript_stable_id,
            cds_stable_id=cds_stable_id,
        )

    return


def make_transcript_attr(con: duckdb.DuckDBPyConnection) -> duckdb.DuckDBPyConnection:
    """creates a transcript_attr table from several other tables"""
    sql = """CREATE VIEW IF NOT EXISTS exon_view AS
        SELECT
            et.transcript_id AS transcript_id,
            ex.exon_id AS exon_id,
            sr.name AS seqid,
            ex.seq_region_start AS start,
            ex.seq_region_end AS stop,
            ex.seq_region_strand AS strand,
            et.rank AS rank,
            tr.gene_id as gene_id,
            tr.stable_id as transcript_stable_id,
            tl.stable_id as cds_stable_id,
        FROM exon ex
        JOIN seq_region sr ON ex.seq_region_id = sr.seq_region_id
        JOIN exon_transcript et ON ex.exon_id = et.exon_id
        JOIN transcript tr ON et.transcript_id = tr.transcript_id
        LEFT JOIN translation tl ON tr.transcript_id = tl.transcript_id
        """
    # it's a left join between tr and tl to handle the case where
    # there is no translation and thus a transcript_id is missing from tl
    con.sql(sql)
    # the transcript_attr schema
    value_placeholder = "?, " * len(TRANSCRIPT_ATTR_COLS)
    sql = f"""CREATE TABLE transcript_attr ({','.join(TRANSCRIPT_ATTR_SCHEMA)})"""
    con.sql(sql)
    values = [
        r.to_record(TRANSCRIPT_ATTR_COLS) for r in get_transcript_attr_records(con)
    ]

    sql = f"""INSERT INTO transcript_attr ({','.join(TRANSCRIPT_ATTR_COLS)}) VALUES ({value_placeholder})"""
    con.executemany(sql, values)
    return con


def make_gene_attr(con: duckdb.DuckDBPyConnection) -> duckdb.DuckDBPyConnection:
    """creates a gene_attr 'table' from several other tables"""
    sql = """CREATE VIEW IF NOT EXISTS gene_attr AS
        SELECT 
            g.gene_id AS gene_id,
            g.stable_id AS stable_id,
            g.biotype AS biotype,
            g.canonical_transcript_id AS canonical_transcript_id,
            sr.name AS seqid,
            g.seq_region_start AS start,
            g.seq_region_end AS stop,
            g.seq_region_strand AS strand,
            x.display_label AS symbol,
            g.description as description,
        FROM gene g
        JOIN seq_region sr ON g.seq_region_id = sr.seq_region_id
        LEFT JOIN xref x ON g.display_xref_id = x.xref_id
        """
    con.sql(sql)
    return con


# the repeat and location related tables are left as is
