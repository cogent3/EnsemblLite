import pytest

from cogent3 import load_table

from ensembl_lite._homologydb import _HOMOLOGYDB_NAME, HomologyDb
from ensembl_lite.install import LoadHomologies


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
    expect = _make_expected_o2o(table)

    data = loader([raw])
    homdb = HomologyDb(tmp_dir / _HOMOLOGYDB_NAME)
    homdb.add_records(records=data, col_order=loader.dest_col)
    return homdb, expect


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
    homdb, expect = o2o_db

    got = homdb.get_related_to(gene_id=gene_id, relationship_type="ortholog_one2one")
    assert got == expect[gene_id]
