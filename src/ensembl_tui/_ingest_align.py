import pathlib
import typing

import duckdb
import numpy
import rich.progress as rich_progress
from cogent3 import make_seq
from cogent3.app.composable import LOADER, define_app
from cogent3.app.typing import IdentifierType

from ensembl_tui import _align as eti_align
from ensembl_tui import _config as eti_config
from ensembl_tui import _ingest_annotation as eti_db_ingest
from ensembl_tui import _util as eti_util
from ensembl_tui._maf import parse

_no_gaps = numpy.array([], dtype=numpy.int32)


def seq2gaps(record: dict) -> eti_align.AlignRecord:
    seq = make_seq(record.pop("seq"))
    indel_map, _ = seq.parse_out_gaps()
    if indel_map.num_gaps:
        record["gap_spans"] = numpy.array(
            [indel_map.gap_pos, indel_map.get_gap_lengths()],
            dtype=numpy.int32,
        ).T
    else:
        record["gap_spans"] = _no_gaps
    return eti_align.AlignRecord(**record)


@define_app(app_type=LOADER)
class load_align_records:  # noqa: N801
    def __init__(self, species: set[str] | None = None) -> None:
        self.species = species or {}

    def main(self, path: IdentifierType) -> list[eti_align.AlignRecord]:
        records = []
        for block_id, align in parse(path):
            converted = []
            for maf_name, seq in align.items():
                if self.species and maf_name.species not in self.species:
                    continue
                record = maf_name.to_dict()
                record["block_id"] = block_id
                record["source"] = path.name
                record["seq"] = seq
                converted.append(seq2gaps(record))
            records.extend(converted)
        return records


def make_alignment_aggregator_db() -> duckdb.DuckDBPyConnection:
    conn = duckdb.connect(":memory:")
    conn.sql("CREATE SEQUENCE align_id_seq")
    columns = ", ".join(eti_align.ALIGN_ATTR_SCHEMA)
    conn.sql(f"CREATE TABLE align_blocks ({columns})")
    return conn


def add_records(
    conn: duckdb.DuckDBPyConnection,
    records: typing.Sequence[eti_align.AlignRecord],
    gap_store: eti_align.GapStore,
    progress: rich_progress.Progress | None = None,
) -> None:
    if not records:
        return
    # we need to identify block_id's that have already been used
    col_order = [c for c in eti_align.ALIGN_ATTR_COLS if c != "align_id"]
    block_ids = tuple({r.block_id for r in records})
    val_placeholder = ", ".join("?" * len(block_ids))
    sql = f"SELECT DISTINCT(block_id) from align_blocks WHERE block_id IN ({val_placeholder})"
    used = {r[0] for r in conn.sql(sql, params=block_ids).fetchall()}

    val_placeholder = ", ".join("?" * len(col_order))
    sql = (
        f"INSERT INTO align_blocks ({', '.join(col_order)}) VALUES ({val_placeholder})"
    )
    rowid_sql = "SELECT currval('align_id_seq')"
    if progress is not None:
        msg = "Writings aligns ✍️"
        writing = progress.add_task(total=len(records), description=msg, advance=0)

    for record in records:
        if record.block_id in used:
            continue

        conn.sql(sql, params=[record[c] for c in col_order])
        index = conn.sql(rowid_sql).fetchone()[0]
        gap_store.add_record(index=index, gaps=record.gap_spans)
        if progress is not None:
            progress.update(writing, description=msg, advance=1)

    if progress is not None:
        progress.remove_task(writing)


def install_alignment(
    config: eti_config.Config,
    align_name: str,
    progress: rich_progress.Progress | None = None,
    max_workers: int | None = None,
) -> pathlib.Path:
    src_dir = config.staging_aligns / align_name
    dest_dir = config.install_aligns / align_name

    dest_dir.mkdir(parents=True, exist_ok=True)
    paths = list(src_dir.glob(f"{align_name}*maf*"))
    aln_loader = load_align_records(set(config.db_names))
    agg = make_alignment_aggregator_db()
    records = []
    series = eti_util.get_iterable_tasks(
        func=aln_loader,
        series=paths,
        max_workers=max_workers,
    )
    gap_path = dest_dir / f"{dest_dir.stem}.{eti_align.GAP_STORE_SUFFIX}"

    gap_store = eti_align.GapStore(
        source=gap_path,
        align_name=align_name,
        mode="w",
        in_memory=False,
    )

    if progress is not None:
        msg = "Reading aligns 📖"
        reading = progress.add_task(total=len(paths), description=msg, advance=0)

    for result in series:
        if not result:
            msg = f"{result=}"
            raise RuntimeError(msg)

        records.extend(result)

        if progress is not None:
            progress.update(reading, description=msg, advance=1)

    if progress is not None:
        progress.remove_task(reading)

    add_records(conn=agg, records=records, gap_store=gap_store, progress=progress)
    gap_store.close()
    # write the parquet file, returns patrh to that file
    return eti_db_ingest.export_parquet(
        con=agg,
        table_name="align_blocks",
        dest_dir=dest_dir,
    )
