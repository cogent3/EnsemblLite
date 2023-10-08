import pytest

from ensembl_cli._genome import CompressedGenome, Genome


@pytest.fixture(scope="function")
def small_genome():
    seqs = {"s1": "AAACCCCAA", "s2": "TTGGTT"}
    db = Genome(source=":memory:", species="dodo")
    db.add_records(records=seqs.items())
    return db, seqs


@pytest.fixture(scope="function")
def compressed_small_genome():
    seqs = {"s1": "AAACCCCAA", "s2": "TTGGTT"}
    db = CompressedGenome(source=":memory:", species="dodo")
    db.add_records(records=seqs.items())
    return db, seqs


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
