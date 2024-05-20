import pytest

from cogent3 import make_unaligned_seqs
from cogent3.core.annotation_db import GffAnnotationDb

from ensembl_lite._genomedb import (
    _SEQDB_NAME,
    Genome,
    SeqsDataHdf5,
    get_gene_table_for_species,
    get_species_summary,
    str2arr,
)


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
def small_coll(small_data, small_annotdb):
    seqs = make_unaligned_seqs(data=small_data, moltype="dna")
    seqs.annotation_db = small_annotdb
    return seqs


@pytest.fixture(scope="function")
def h5_genome(tmp_path):
    # in memory db
    return SeqsDataHdf5(
        source=tmp_path / "small-hd5f.genome-h5",
        mode="w",
        species="Human",
        in_memory=True,
    )


@pytest.fixture
def small_h5_genome(small_data, h5_genome):
    # in memory db
    h5_genome.add_records(records=small_data.items())
    return h5_genome, small_data


@pytest.mark.parametrize(
    "name,start,stop", (("s1", 3, 7), ("s1", 3, None), ("s1", None, 7), ("s2", 2, 4))
)
def test_get_seq(small_h5_genome, name, start, stop):
    genome, seqs = small_h5_genome
    expect = seqs[name][start:stop]
    assert genome.get_seq_str(seqid=name, start=start, stop=stop) == expect


@pytest.mark.parametrize("name", ("s1", "s2"))
def test_get_fullseq(small_h5_genome, name):
    genome, seqs = small_h5_genome
    expect = seqs[name]
    assert genome.get_seq_str(seqid=name) == expect


def test_annodb(small_annotdb):
    list(small_annotdb.get_features_matching(seqid="s1", biotype="gene"))


def test_selected_seq_is_annotated(small_h5_genome, small_annotdb, namer):
    gen_seqs_db, _ = small_h5_genome
    genome = Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid="s1", namer=namer)
    assert len(seq.annotation_db) == 3
    gene = list(genome.get_features(seqid="s1", biotype="gene"))[0]
    gene_seq = gene.get_slice()
    assert str(gene_seq) == "AAAA"
    assert gene.name == "gene-01"


def test_hashable_genome_seqs(h5_genome):
    assert hash(h5_genome) == id(h5_genome)


def test_genome_close(small_h5_genome, small_annotdb, namer):
    gen_seqs_db, _ = small_h5_genome
    genome = Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid="s1", namer=namer)
    assert seq
    genome.close()
    with pytest.raises(OSError):
        genome.get_seq(seqid="s1")


@pytest.mark.parametrize("seqid", ("s1", "s2"))
def test_get_seq_num_annotations_correct(
    small_h5_genome, small_annotdb, small_coll, seqid, namer
):
    gen_seqs_db, small_data = small_h5_genome
    genome = Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid=seqid, namer=namer)
    expect = list(small_coll.get_features(seqid=seqid))
    assert len(seq.annotation_db) == len(expect)


@pytest.mark.parametrize(
    "seqid,feature_name,start,stop",
    (
        ("s1", None, None, None),
        ("s1", "gene-01", 2, 8),
        ("s2", "gene-02", 1, 8),
    ),
)
def test_get_seq_feature_seq_correct(
    small_h5_genome, small_annotdb, small_coll, seqid, feature_name, start, stop, namer
):
    gen_seqs_db, small_data = small_h5_genome
    genome = Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid=seqid, start=start, stop=stop, namer=namer)
    coll_seq = small_coll.get_seq(seqid)
    assert seq == coll_seq[start:stop]
    expect = list(coll_seq[start:stop].get_features(allow_partial=True))[0]
    got = list(seq.get_features(allow_partial=True))[0]
    # should also get the same slice
    assert got.get_slice() == expect.get_slice()


def test_get_gene_table_for_species(small_annotdb):
    from cogent3.util.table import Table

    # we do not check values here, only the Type and that we have > 0 records
    got = get_gene_table_for_species(annot_db=small_annotdb, limit=None, species="none")
    assert isinstance(got, Table)
    assert len(got) > 0


def test_get_species_summary(small_annotdb):
    from cogent3.util.table import Table

    got = get_species_summary(annot_db=small_annotdb, species="none")
    # we do not check values here, only the Type and that we have > 0 records
    assert isinstance(got, Table)
    assert len(got) > 0


def test_hdf5_genome_skip_duplicates(small_h5_genome):
    genome, data = small_h5_genome
    # should not fail
    genome.add_records(records=data.items())


def test_hdf5_genome_errors_sameid_diff_seq(small_h5_genome):
    genome, data = small_h5_genome
    # same eqid but diff seq should fail
    data = {"s1": "AAA"}
    with pytest.raises(ValueError):
        genome.add_records(records=data.items())


def test_hdf5_genome_error_duplicate_names(small_h5_genome):
    genome, data = small_h5_genome
    with pytest.raises(ValueError):
        # duplicate name, but seq is different
        genome.add_record(seqid="s1", seq=data["s1"][:-2])


def test_hdf5_genome_coord_names(small_h5_genome):
    genome, data = small_h5_genome
    assert genome.get_coord_names() == tuple(data)


def test_empty_hdf5_genome_coord_names(h5_genome):
    assert h5_genome.get_coord_names() == ()


@pytest.mark.parametrize(
    "name,start,stop",
    (
        ("s1", 3, 7),
        ("s1", 3, None),
        ("s1", None, 7),
        ("s2", 2, 4),
        ("s1", None, None),
        ("s2", None, None),
    ),
)
def test_h5_get_seq(small_h5_genome, name, start, stop):
    genome, seqs = small_h5_genome
    expect = seqs[name][start:stop]
    assert genome.get_seq_str(seqid=name, start=start, stop=stop) == expect
    convert = str2arr(moltype="dna")
    assert (
        genome.get_seq_arr(seqid=name, start=start, stop=stop) == convert(expect)
    ).all()


def test_pickling_round_trip(small_data, tmp_path):
    import pickle  # nosec B403

    path = tmp_path / f"small.{_SEQDB_NAME}"
    kwargs = dict(source=path, species="human")
    genome = SeqsDataHdf5(mode="w", **kwargs)
    genome.add_records(records=small_data.items())
    with pytest.raises(NotImplementedError):
        pickle.dumps(genome)  # nosec B301

    ro = SeqsDataHdf5(mode="r", **kwargs)
    assert ro.get_seq_str(seqid="s1") == small_data["s1"]
    unpkl = pickle.loads(pickle.dumps(ro))  # nosec B301
    got = unpkl.get_seq_str(seqid="s1")
    assert got == small_data["s1"]


def test_species_setting(small_data, tmp_path):
    path = tmp_path / f"small.{_SEQDB_NAME}"
    kwargs = dict(source=path, species="human")
    genome = SeqsDataHdf5(mode="w", **kwargs)
    genome.add_records(records=small_data.items())
    genome.close()

    genome = SeqsDataHdf5(mode="r", source=path)
    # note that species are converted into the Ensembl db prefix
    assert genome.species == "homo_sapiens"
    with pytest.raises(ValueError):
        _ = SeqsDataHdf5(mode="r", source=path, species="cat")


def test_has_of_seqsdata(h5_genome):
    assert hash(h5_genome) == id(h5_genome)
