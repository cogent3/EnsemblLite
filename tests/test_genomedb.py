import pytest

from cogent3.core.annotation_db import GffAnnotationDb

from ensembl_lite._genomedb import CompressedGenomeSeqsDb, Genome, GenomeSeqsDb
from ensembl_lite.util import elt_compress_it


@pytest.fixture(scope="function")
def small_data():
    return {"s1": "TAACCCCAAG", "s2": "TTGGTTGG"}


@pytest.fixture(scope="function")
def small_annots():
    return [
        dict(
            seqid="s1",
            name="gene-01",
            biotype="gene",
            spans=[(1, 3), (7, 9)],
            strand="+",
        ),
        dict(
            seqid="s1",
            name="exon-01",
            biotype="exon",
            spans=[(1, 3)],
            strand="+",
            parent_id="gene-01",
        ),
        dict(
            seqid="s1",
            name="exon-02",
            biotype="exon",
            spans=[(7, 9)],
            strand="+",
            parent_id="gene-01",
        ),
        dict(
            seqid="s2",
            name="gene-02",
            biotype="gene",
            spans=[(2, 4), (6, 8)],
            strand="+",
        ),
    ]


@pytest.fixture(scope="function")
def small_annotdb(small_annots):
    db = GffAnnotationDb(source=":memory:")
    for record in small_annots:
        db.add_feature(**record)
    return db


@pytest.fixture(scope="function")
def small_genome(small_data):
    db = GenomeSeqsDb(source=":memory:", species="dodo")
    db.add_records(records=small_data.items())
    return db, small_data


@pytest.fixture(scope="function")
def compressed_small_genome(small_data):
    db = CompressedGenomeSeqsDb(source=":memory:", species="dodo")
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
    db = CompressedGenomeSeqsDb(source=":memory:", species="dodo")
    data = {k: elt_compress_it(s) for k, s in small_data.items()}
    db.add_compressed_records(records=data.items())
    assert db.get_seq(coord_name="s1") == small_data["s1"]
    assert db.get_seq(coord_name="s2") == small_data["s2"]


def test_annodb(small_annotdb):
    list(small_annotdb.get_features_matching(seqid="s1", biotype="gene"))


def test_selected_seq_is_annotated(small_genome, small_annotdb):
    gen_seqs_db, _ = small_genome
    genome = Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid="s1")
    assert seq.annotation_db == small_annotdb
    gene = list(genome.get_features(seqid="s1", biotype="gene"))[0]
    gene_seq = gene.get_slice()
    assert str(gene_seq) == "AAAA"
    assert gene.name == "gene-01"


@pytest.mark.parametrize("cls", (GenomeSeqsDb, CompressedGenomeSeqsDb))
def test_hashable_genome(cls):
    genome = cls(species="dodo", source=":memory:")
    assert hash(genome) == id(genome)


def test_genome_close(small_genome, small_annotdb):
    import sqlite3

    gen_seqs_db, _ = small_genome
    genome = Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid="s1")
    assert seq
    genome.close()
    with pytest.raises(sqlite3.ProgrammingError):
        genome.get_seq(seqid="s1")
