import pytest

from ensembl_cli._genome import Genome


@pytest.fixture(scope="function")
def small_genome():
    seqs = {"s1": "AAACCCCAA", "s2": "TTGGTT"}
    db = Genome(source=":memory:", species="dodo")
    db.add_records(seqs.items())
    return db, seqs


@pytest.mark.parametrize(
    "name,start,end", (("s1", 3, 7), ("s1", 3, None), ("s1", None, 7), ("s2", 2, 4))
)
def test_get_seq(small_genome, name, start, end):
    genome, seqs = small_genome
    expect = seqs[name][start:end]
    assert genome.get_seq(name, start, end) == expect


@pytest.mark.parametrize("name", ("s1", "s2"))
def test_get_fullseq(small_genome, name):
    genome, seqs = small_genome
    expect = seqs[name]
    assert genome.get_seq(name) == expect
