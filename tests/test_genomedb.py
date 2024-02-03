import pytest

from cogent3 import make_unaligned_seqs
from cogent3.core.annotation_db import GffAnnotationDb

from ensembl_lite._genomedb import (
    CompressedGenomeSeqsDb,
    Genome,
    GenomeSeqsDb,
    get_gene_table_for_species,
)
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
            strand="-",
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
def small_coll(small_data, small_annotdb):
    seqs = make_unaligned_seqs(data=small_data, moltype="dna")
    seqs.annotation_db = small_annotdb
    return seqs


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
    assert genome.get_seq(seqid=name, start=start, end=end) == expect


@pytest.mark.parametrize("genome", ("small_genome", "compressed_small_genome"))
@pytest.mark.parametrize("name", ("s1", "s2"))
def test_get_fullseq(genome, request, name):
    genome, seqs = request.getfixturevalue(genome)
    expect = seqs[name]
    assert genome.get_seq(seqid=name) == expect


def test_add_compressed(small_data):
    db = CompressedGenomeSeqsDb(source=":memory:", species="dodo")
    data = {k: elt_compress_it(s) for k, s in small_data.items()}
    db.add_compressed_records(records=data.items())
    assert db.get_seq(seqid="s1") == small_data["s1"]
    assert db.get_seq(seqid="s2") == small_data["s2"]


def test_annodb(small_annotdb):
    list(small_annotdb.get_features_matching(seqid="s1", biotype="gene"))


def test_selected_seq_is_annotated(small_genome, small_annotdb, namer):
    gen_seqs_db, _ = small_genome
    genome = Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid="s1", namer=namer)
    assert len(seq.annotation_db) == 3
    gene = list(genome.get_features(seqid="s1", biotype="gene"))[0]
    gene_seq = gene.get_slice()
    assert str(gene_seq) == "AAAA"
    assert gene.name == "gene-01"


@pytest.mark.parametrize("cls", (GenomeSeqsDb, CompressedGenomeSeqsDb))
def test_hashable_genome(cls):
    genome = cls(species="dodo", source=":memory:")
    assert hash(genome) == id(genome)


def test_genome_close(small_genome, small_annotdb, namer):
    import sqlite3

    gen_seqs_db, _ = small_genome
    genome = Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid="s1", namer=namer)
    assert seq
    genome.close()
    with pytest.raises(sqlite3.ProgrammingError):
        genome.get_seq(seqid="s1")


@pytest.mark.parametrize("seqid", ("s1", "s2"))
def test_get_seq_num_annotations_correct(
    small_genome, small_annotdb, small_coll, seqid, namer
):
    gen_seqs_db, small_data = small_genome
    genome = Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid=seqid, namer=namer)
    expect = list(small_coll.get_features(seqid=seqid))
    assert len(seq.annotation_db) == len(expect)


@pytest.mark.parametrize(
    "seqid,feature_name,start,end",
    (
        ("s1", None, None, None),
        ("s1", "gene-01", 2, 8),
        ("s2", "gene-02", 1, 8),
    ),
)
def test_get_seq_feature_seq_correct(
    small_genome, small_annotdb, small_coll, seqid, feature_name, start, end, namer
):
    gen_seqs_db, small_data = small_genome
    genome = Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid=seqid, start=start, end=end, namer=namer)
    coll_seq = small_coll.get_seq(seqid)
    assert seq == coll_seq[start:end]
    expect = list(coll_seq[start:end].get_features(allow_partial=True))[0]
    got = list(seq.get_features(allow_partial=True))[0]
    # should also get the same slice
    assert got.get_slice() == expect.get_slice()


def test_get_gene_table_for_species(small_annotdb):
    from cogent3.util.table import Table

    # we do not check values here, only the Type and that we have > 0 records
    got = get_gene_table_for_species(annot_db=small_annotdb, limit=None)
    assert isinstance(got, Table)
    assert len(got) > 0
