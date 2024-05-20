import pytest

from cogent3 import load_table

from ensembl_lite._homologydb import (
    _HOMOLOGYDB_NAME,
    HomologyDb,
    HomologyRecord,
    grouped_related,
    homolog_group,
    species_genes,
)
from ensembl_lite._install import LoadHomologies


def _make_expected_o2o(table):
    """return dict with keys stable ID's values list[tuple[str,str]]"""
    data = table.to_list(["species_1", "gene_id_1", "species_2", "gene_id_2"])

    result = {}
    for sp1, g1, sp2, g2 in data:
        value = {(sp1, g1), (sp2, g2)}
        result[g1] = result.get(g1, set()) | value
        result[g2] = result.get(g2, set()) | value

    return result


@pytest.fixture
def o2o_db(DATA_DIR, tmp_dir):
    raw = DATA_DIR / "one2one_homologies.tsv"

    species = {
        "gorilla_gorilla",
        "macaca_mulatta",
        "microcebus_murinus",
        "homo_sapiens",
        "pongo_abelii",
        "pan_troglodytes",
        "macaca_fascicularis",
        "chlorocebus_sabaeus",
        "pan_paniscus",
    }
    loader = LoadHomologies(species)

    table = load_table(raw).get_columns(loader.src_cols)

    table = table.with_new_header(loader.src_cols, loader.dest_col[:-1])

    data = loader([raw])
    homdb = HomologyDb(tmp_dir / _HOMOLOGYDB_NAME)
    homdb.add_records(records=data, col_order=loader.dest_col)
    return homdb, table


@pytest.mark.parametrize(
    "gene_id",
    (
        "ENSGGOG00000026757",
        "ENSGGOG00000025053",
        "ENSGGOG00000022688",
        "ENSGGOG00000026221",
        "ENSGGOG00000024015",
    ),
)
def test_hdb(o2o_db, gene_id):
    homdb, table = o2o_db
    expect = _make_expected_o2o(table)

    got = homdb.get_related_to(gene_id=gene_id, relationship_type="ortholog_one2one")
    assert got == expect[gene_id]


@pytest.fixture
def orth_records():
    common = dict(relationship="ortholog_one2one")
    return [
        HomologyRecord(
            species_1="a", gene_id_1="1", species_2="b", gene_id_2="2", **common
        ),  # grp 1
        HomologyRecord(
            species_1="b", gene_id_1="2", species_2="c", gene_id_2="3", **common
        ),  # grp 1
        HomologyRecord(
            species_1="a", gene_id_1="4", species_2="b", gene_id_2="5", **common
        ),
    ]


@pytest.fixture
def hom_records(orth_records):
    return orth_records + [
        HomologyRecord(
            species_1="a",
            gene_id_1="6",
            species_2="e",
            gene_id_2="7",
            relationship="ortholog_one2many",
        ),
    ]


def _reduced_to_ids(groups):
    result = []
    for group in groups:
        result.append(group.all_gene_ids())
    return result


def test_hdb_get_related_groups(o2o_db):
    homdb, _ = o2o_db
    got = homdb.get_related_groups(relationship_type="ortholog_one2one")
    groups = _reduced_to_ids(got)
    assert len(groups) == 5


@pytest.fixture
def hom_hdb(hom_records):
    col_order = (
        "relationship",
        "species_1",
        "gene_id_1",
        "species_2",
        "gene_id_2",
    )
    records = [[r[c] for c in col_order] for r in hom_records]
    hdb = HomologyDb(source=":memory:")
    hdb.add_records(records=records, col_order=col_order)
    return hdb


def test_group_related(hom_records):
    orths = [r for r in hom_records if r["relationship"] == "ortholog_one2one"]
    got = sorted(grouped_related(orths), key=lambda x: len(x), reverse=True)
    expect = [
        homolog_group(
            relationship="ortholog_one2one",
            data=(
                species_genes(species="a", gene_ids=["1"]),
                species_genes(species="b", gene_ids=["2"]),
                species_genes(species="c", gene_ids=["3"]),
            ),
        ),
        homolog_group(
            relationship="ortholog_one2one",
            data=(
                species_genes(species="a", gene_ids=["4"]),
                species_genes(species="b", gene_ids=["5"]),
            ),
        ),
    ]
    assert got == expect


def test_homology_db(hom_hdb):
    got = sorted(
        hom_hdb.get_related_groups("ortholog_one2one"),
        key=lambda x: len(x),
        reverse=True,
    )

    expect = [
        homolog_group(
            relationship="ortholog_one2one",
            data=(
                species_genes(species="a", gene_ids=["1"]),
                species_genes(species="b", gene_ids=["2"]),
                species_genes(species="c", gene_ids=["3"]),
            ),
        ),
        homolog_group(
            relationship="ortholog_one2one",
            data=(
                species_genes(species="a", gene_ids=["4"]),
                species_genes(species="b", gene_ids=["5"]),
            ),
        ),
    ]
    assert got == expect


def test_species_genes_eq():
    a = species_genes(species="a", gene_ids=["1"])
    b = species_genes(species="a", gene_ids=["1"])
    assert a == b
    c = species_genes(species="a", gene_ids=["2"])
    assert a != c
    d = species_genes(species="b", gene_ids=["1"])
    assert a != d
    e = species_genes(species="b", gene_ids=["2"])
    assert a != e


@pytest.mark.parametrize("col", tuple(HomologyDb._index_columns))
def test_indexing(o2o_db, col):
    db, _ = o2o_db
    expect = ("index", col, db.table_name)
    db.make_indexes()
    sql_template = (
        f"SELECT * FROM sqlite_master WHERE type = 'index' AND "  # nosec B608
        f"tbl_name = {db.table_name!r} and name = {col!r}"  # nosec B608
    )

    result = db._execute_sql(sql_template).fetchone()
    got = tuple(result)[:3]
    assert got == expect
def test_homolog_group_union():
    a = homolog_group(relationship="one2one", gene_ids={"1", "2", "3"})
    b = homolog_group(relationship="one2one", gene_ids={"3", "4"})
    c = a | b
    assert c.gene_ids == {"1", "2", "3", "4"}


def test_homolog_group_union_invalid():
    a = homolog_group(relationship="one2one", gene_ids={"1", "2", "3"})
    b = homolog_group(relationship="one2many", gene_ids={"3", "4"})
    with pytest.raises(ValueError):
        _ = a | b


