import pytest

from ensembl_cli._genomedb import (
    CompressedGenomeDb,
    GenomeDb,
    compress_it,
    decompress_it,
)


@pytest.fixture(scope="function")
def small_data():
    return {"s1": "AAACCCCAA", "s2": "TTGGTT"}


@pytest.fixture(scope="function")
def small_genome(small_data):
    db = GenomeDb(source=":memory:", species="dodo")
    db.add_records(records=small_data.items())
    return db, small_data


@pytest.fixture(scope="function")
def compressed_small_genome(small_data):
    db = CompressedGenomeDb(source=":memory:", species="dodo")
    db.add_records(records=small_data.items())
    return db, small_data


@pytest.mark.parametrize("genome", ("small_genome", "compressed_small_genome"))
@pytest.mark.parametrize(
    "name,start,end", (("s1", 3, 7), ("s1", 3, None), ("s1", None, 7), ("s2", 2, 4))
)
def test_get_seq(genome, request, name, start, end):
    genome, seqs = request.getfixturevalue(genome)
    expect = seqs[name][start:end]
    assert genome.get_seq(coord_name=name, start=start, end=end) == expect


@pytest.mark.parametrize("genome", ("small_genome", "compressed_small_genome"))
@pytest.mark.parametrize("name", ("s1", "s2"))
def test_get_fullseq(genome, request, name):
    genome, seqs = request.getfixturevalue(genome)
    expect = seqs[name]
    assert genome.get_seq(coord_name=name) == expect


def test_add_compressed(small_data):
    db = CompressedGenomeDb(source=":memory:", species="dodo")
    data = {k: compress_it(s) for k, s in small_data.items()}
    db.add_compressed_records(records=data.items())
    assert db.get_seq(coord_name="s1") == small_data["s1"]
    assert db.get_seq(coord_name="s2") == small_data["s2"]
