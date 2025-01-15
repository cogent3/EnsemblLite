import contextlib
import pathlib
import shutil
from collections.abc import Generator

import duckdb
import rich.progress as rich_progress
from cogent3.app.composable import LOADER, define_app

from ensembl_tui import _config as eti_config
from ensembl_tui import _util as eti_util

PARQUET_FORMAT = "(FORMAT PARQUET, COMPRESSION 'zstd', ROW_GROUP_SIZE 100_000)"


def get_start_column(con: duckdb.DuckDBPyConnection, table_name: str) -> set[str]:
    """return list of column names ending with 'start'

    Notes
    -----
    We will be using this to adjust from a 1-based to 0-based counting system
    """
    sql = f"""
    SELECT column_name
    FROM information_schema.columns
    WHERE table_name = '{table_name}' AND column_name LIKE '%_start';
    """
    return {r[0] for r in con.sql(sql).fetchall()}


def show_some_data(
    con: duckdb.DuckDBPyConnection,
    table_name: str,
    limit: int = 5,
) -> None:
    # this is a useful function for debugging
    sql = f"SELECT * FROM {table_name} LIMIT {limit}"
    print(table_name, con.sql(sql), sep="\n")  # noqa: T201
    sql = f"SELECT COUNT(*) as num_rows FROM {table_name}"
    print(con.sql(sql))  # noqa: T201


@contextlib.contextmanager
def tempdb(
    source_db_path: str | pathlib.Path,
) -> Generator[duckdb.DuckDBPyConnection, None, None]:
    """context manager returns a duckdb connection to a source db copy"""
    source_db_path = pathlib.Path(source_db_path)
    if not source_db_path.exists():
        raise FileNotFoundError(source_db_path)

    with eti_util.tempdir() as temp_dir:
        temp_db_path = pathlib.Path(temp_dir) / source_db_path.name
        shutil.copy2(source_db_path, temp_db_path)

        try:
            conn = duckdb.connect(str(temp_db_path))
            conn.sql("PRAGMA disable_progress_bar;")
            yield conn
        finally:
            conn.close()


def migrate_schema(con: duckdb.DuckDBPyConnection, table_name: str) -> None:
    """correct duckdb schema inference to match mysql schema

    Notes
    -----
    duckdb always converts MySQL tinyint to bool. But Ensembl strand field has
    3 possible values (-1, 0, 1). So we explicitly set types for all columns
    whose name ends in "strand".
    """
    sql = f"CREATE TABLE {table_name} AS SELECT * FROM mysqldb.{table_name} LIMIT 0"
    con.sql(sql)
    # assume any column name ending with "strand" is a smallint (-1, 0, 1) for
    # minus, unknown, plus strand
    names = con.sql(f"DESCRIBE {table_name}").to_df()
    for n in names["column_name"].to_list():
        if n.endswith("strand"):
            sql = f"ALTER TABLE {table_name} ALTER COLUMN {n} SET DATA TYPE TINYINT;"
            con.sql(sql)


def make_mysql_connection(
    *,
    db_name: str,
    db_host: str = "ensembldb.ensembl.org",
    db_port: int = 3306,
    db_user: str = "anonymous",
    db_path: pathlib.Path,
) -> duckdb.DuckDBPyConnection:
    # TODO: modify eti config object to support mysqlserver and port
    uk_connect = f"ATTACH 'host={db_host} port={db_port} user={db_user} database={db_name}' AS mysqldb (TYPE mysql)"
    con = duckdb.connect(str(db_path))
    con.sql("PRAGMA disable_progress_bar;")
    # install mysql extension
    con.sql("INSTALL mysql")
    con.sql("LOAD mysql")
    con.sql(uk_connect)
    return con


# convert following to using config object
def make_table_template(
    *,
    dest_dir: pathlib.Path,
    db_name: str,
    table_name: str,
    db_host: str = "ensembldb.ensembl.org",
    db_port: int = 3306,
    db_user: str = "anonymous",
) -> pathlib.Path:
    """make a template db file for a given table

    Notes
    -----
    This is done only ONCE for a given config as we assume all genome db's
    share the same schema.
    """
    outname = dest_dir / f"{table_name}.duckdb"
    if outname.exists():
        return outname

    con = make_mysql_connection(
        db_name=db_name,
        db_host=db_host,
        db_port=db_port,
        db_user=db_user,
        db_path=outname,
    )
    migrate_schema(con, table_name)
    con.close()
    return outname


def import_mysql_table(
    *,
    con: duckdb.DuckDBPyConnection,
    mysql_dump_path: pathlib.Path,
    table_name: str,
    fix_start: bool = True,
) -> None:
    """imports a mysql dump file into a duckdb table

    Parameters
    ----------
    con
        duckdb connection
    mysql_dump_path
        the full path to the dumpy file
    table_name
        the name of the table
    fix_start
        if True, columns ending in "_start" are adjusted from 1-based to 0-based
    """
    sql = f"INSERT INTO {table_name} SELECT * FROM read_csv_auto('{mysql_dump_path}', nullstr='\\N', header=false, ignore_errors=false, delim='\\t')"
    con.sql(sql)
    if fix_start:
        for start_column in get_start_column(con, table_name):
            con.execute(
                f"UPDATE {table_name} SET {start_column} = {start_column} - 1 WHERE {start_column} IS NOT NULL",
            )


def export_parquet(
    *,
    con: duckdb.DuckDBPyConnection,
    table_name: str,
    dest_dir: pathlib.Path,
) -> pathlib.Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    parquet_file = dest_dir / f"{table_name}.parquet"
    con.sql(
        f"COPY {table_name} TO '{parquet_file}' {PARQUET_FORMAT};",
    )
    return parquet_file


def write_parquet(
    *,
    db_templates: pathlib.Path,
    dump_path: pathlib.Path,
    table_name: str,
    dest_dir: pathlib.Path,
    fix_start: bool = True,
) -> pathlib.Path:
    """create a parquet file at `dest_dir/{table_name}.parquet`

    Parameters
    ----------
    db_templates
        directory to the database template files
    dump_path
        exact path to the dump file
    table_name
        name of the table to be created
    dest_dir
        location to write the parquet file
    fix_start
        if True, columns ending in "_start" are adjusted from 1-based to 0-based

    Returns
    -------
    path to the generated parquet file
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    with tempdb(db_templates / f"{table_name}.duckdb") as con:
        import_mysql_table(
            con=con,
            mysql_dump_path=dump_path,
            table_name=table_name,
            fix_start=fix_start,
        )
        return export_parquet(con=con, table_name=table_name, dest_dir=dest_dir)


# we merge tables related to transcription, translation, exons and exons_transcripts
# into a single table. The coordinates (spans) for multiple exons are stored as a binary blob
# for transcripts and the translated exons in a single row. This will improve querying performance
# at the cost of install time and storage.


def _make_db(
    config: eti_config.Config,
    db_name: str,
    table_names: tuple[str, ...],
) -> duckdb.DuckDBPyConnection:
    # utility function for making a virtual db from a series of paquet files
    files = []
    for table_name in table_names:
        parquet_file = config.install_genomes / db_name / f"{table_name}.parquet"
        if not parquet_file.exists():
            msg = "use mysql_dump_to_parquet() app first"
            raise FileNotFoundError(msg)
        files.append(parquet_file)

    conn = duckdb.connect(":memory:")

    for file_path in files:
        conn.execute(
            f"CREATE TABLE {file_path.stem} AS SELECT * FROM read_parquet('{file_path}')",
        )

    return conn


def make_combined_tables(
    *,
    config: eti_config.Config,
    db_name: str,
    cleanup: bool = True,
) -> None:
    """makes combined tables for transcripts and genes and writes as parquet files

    Parameters
    ----------
    config
        an downloaed cfg instance
    db_name
        directory name reprsenting a species
    cleanup
        if provided, the parquet files from which merged tables are built
        are deleted on completion

    Notes
    -----
    Creates transcript_attr.parquet and gene_attr.parquet files.
    """
    from ensembl_tui._mysql_core_attr import make_gene_attr, make_transcript_attr

    # make the transcribed_attr table
    transcribed_tables = (
        "exon",
        "exon_transcript",
        "transcript",
        "translation",
    )
    preserve = (
        "seq_region",
    )  # keep this one separate as we need to keep it for repeats
    # checks these tables already exist in parquet format, fails otherwise
    conn = _make_db(config, db_name, transcribed_tables + preserve)
    _ = make_transcript_attr(con=conn)
    export_parquet(
        con=conn,
        table_name="transcript_attr",
        dest_dir=config.install_genomes / db_name,
    )
    conn.close()
    # make the gene_attr table
    gene_tables = "gene", "xref"
    conn = _make_db(config, db_name, gene_tables + preserve)
    _ = make_gene_attr(con=conn)
    export_parquet(
        con=conn,
        table_name="gene_attr",
        dest_dir=config.install_genomes / db_name,
    )
    conn.close()
    if cleanup:
        for table_name in transcribed_tables + gene_tables:
            (config.install_genomes / db_name / f"{table_name}.parquet").unlink(
                missing_ok=True,
            )


def install_parquet_tables(
    config: eti_config.Config,
    progress: rich_progress.Progress | None = None,
    make_combined: bool = True,
) -> pathlib.Path:
    """installs annotation data as parquet files

    Parameters
    ----------
    config
        a config instance from a downloaded.cfg file
    progress
        rich Progress instance
    make_combined
        makes the combined attr tables, for testing purposes only!

    Returns
    -------
    The path to the installed genomes directory.
    """
    template_dir = config.staging_template_path
    if not template_dir.exists():
        msg = f"no mysql dump dir for {template_dir}"
        raise FileNotFoundError(msg)

    dump_root = config.staging_genomes
    genome_root = config.install_genomes
    genome_root.mkdir(parents=True, exist_ok=True)
    table_names = [fn.stem for fn in template_dir.glob("*.duckdb")]
    if progress is not None:
        msg = "Installing features 📚"
        num_genomes = len(config.species_dbs)
        num_tables = len(table_names)
        writing = progress.add_task(total=num_genomes * num_tables, description=msg)

    for db_name in config.db_names:
        dump_dir = dump_root / db_name / "mysql"
        if not dump_dir.exists():
            msg = f"no mysql dump dir for {db_name}"
            raise FileNotFoundError(msg)

        dest_dir = genome_root / db_name
        dest_dir.mkdir(parents=True, exist_ok=True)
        for table_name in table_names:
            dump_path = dump_dir / f"{table_name}.txt.gz"
            if not dump_path.exists():
                msg = f"no mysqldump file for {table_name}"
                raise FileNotFoundError(msg)

            write_parquet(
                db_templates=template_dir,
                dump_path=dump_path,
                table_name=table_name,
                dest_dir=dest_dir,
                # the translation table has 2 columns that need the
                # correction to 0-based, that coercion will be done in the
                # merging tables function
                fix_start=table_name != "translation",
            )
            if progress is not None:
                progress.update(writing, advance=1)

        # and now we construct the combined attr tables
        if make_combined:
            make_combined_tables(config=config, db_name=db_name, cleanup=True)

    return genome_root


@define_app(app_type=LOADER)
class mysql_dump_to_parquet:  # noqa: N801
    def __init__(
        self,
        config: eti_config.Config,
        verbose: bool = False,
        make_combined: bool = True,
    ) -> None:
        if not config.staging_template_path.exists():
            msg = f"no mysql dump dir for {config.staging_template_path=}"
            raise FileNotFoundError(msg)

        self._config = config
        self._template_dir = config.staging_template_path
        self._staging_dir = config.staging_genomes
        self._install_dir = config.install_genomes
        self._install_dir.mkdir(parents=True, exist_ok=True)
        self._table_names = [fn.stem for fn in self._template_dir.glob("*.duckdb")]
        self._verbose = verbose
        self._make_combined = make_combined

    def main(self, db_name: str) -> pathlib.Path:
        dump_dir = self._staging_dir / db_name / "mysql"
        if not dump_dir.exists():
            msg = f"no mysql dump dir for {db_name}"
            raise FileNotFoundError(msg)

        dest_dir = self._install_dir / db_name
        dest_dir.mkdir(parents=True, exist_ok=True)
        for table_name in self._table_names:
            dump_path = dump_dir / f"{table_name}.txt.gz"
            if not dump_path.exists():
                msg = f"no mysqldump file for {table_name}"
                raise FileNotFoundError(msg)

            write_parquet(
                db_templates=self._template_dir,
                dump_path=dump_path,
                table_name=table_name,
                dest_dir=dest_dir,
            )

        # and now we construct the combined attr tables
        if self._make_combined:
            make_combined_tables(config=self._config, db_name=db_name, cleanup=True)

        return dest_dir
