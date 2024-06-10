import pytest

from cogent3 import make_unaligned_seqs
from cogent3.core.annotation_db import GffAnnotationDb

from ensembl_lite._genomedb import (
    _ANNOTDB_NAME,
    _SEQDB_NAME,
    EnsemblGffDb,
    EnsemblGffRecord,
    Genome,
    SeqsDataHdf5,
    custom_gff_parser,
    get_gene_table_for_species,
    get_species_summary,
    make_annotation_db,
    make_gene_relationships,
    str2arr,
    tidy_gff3_stableids,
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
            start=1,
            stop=9,
            strand="+",
        ),
        dict(
            seqid="s1",
            name="exon-01",
            biotype="exon",
            spans=[(1, 3)],
            start=1,
            stop=3,
            strand="+",
            parent_id="gene-01",
        ),
        dict(
            seqid="s1",
            name="exon-02",
            biotype="exon",
            spans=[(7, 9)],
            start=7,
            stop=9,
            strand="+",
            parent_id="gene-01",
        ),
        dict(
            seqid="s2",
            name="gene-02",
            biotype="gene",
            spans=[(2, 4), (6, 8)],
            start=2,
            stop=8,
            strand="-",
        ),
    ]


@pytest.fixture(scope="function")
def small_annotdb(small_annots):
    db = EnsemblGffDb(source=":memory:")
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
    assert len(seq.annotation_db) == 4
    genes = list(genome.get_features(seqid="s1", biotype="gene"))
    gene = genes[0]
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
    assert len(list(seq.get_features())) == len(expect)


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
    got = list(seq.get_features(allow_partial=True))
    got = got[0]
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


def test_tidying_stableids_in_gff3():

    orig = (
        "ID=Transcript:ENST00000461467;Parent=Gene:ENSG00000237613;Name=FAM138A-202;bio"
    )
    expect = (
        "ID=transcript:ENST00000461467;Parent=gene:ENSG00000237613;Name=FAM138A-202;bio"
    )
    assert tidy_gff3_stableids(orig) == expect


def test_custom_gff3_parser(DATA_DIR):
    path = DATA_DIR / "c_elegans_WS199_shortened.gff3"
    records, _ = custom_gff_parser(path, 0)

    rel = make_gene_relationships(records.values())
    children = rel["gene:WBGene00000138"]
    # as the records are hashable by their .name attribute, we can just
    # check returned value against their names
    assert children == {"cds:B0019.1", "transcript:B0019.1"}
    # check that multi row records have the correct spans, start, stop and strand
    assert_allclose(
        records["cds:B0019.1"].spans, numpy.array([(9, 20), (29, 45), (59, 70)])
    )
    assert records["cds:B0019.1"].start == 9
    assert records["cds:B0019.1"].stop == 70
    assert records["cds:B0019.1"].strand == "-"


def test_gff_record_size(DATA_DIR):
    merged, _ = custom_gff_parser(DATA_DIR / "c_elegans_WS199_shortened.gff3", 0)
    # record CDS:B0019.1 has spans [(9, 20), (29, 45), (59, 70)] which sum to 38
    starts, stops = numpy.array([(9, 20), (29, 45), (59, 70)]).T
    expect = (stops - starts).sum()
    assert merged["cds:B0019.1"].size == expect


@pytest.mark.parametrize("val", ((), [], numpy.array([]), None))
def tess_gff_record_size_zero(val):
    record = EnsemblGffRecord(spans=val)
    assert record.size == 0


@pytest.mark.parametrize(
    "attrs", ("Ensembl_canonical", "text;other;Ensembl_canonical;than")
)
def test_is_canonical(attrs):
    f = EnsemblGffRecord(attrs=attrs)
    assert f.is_canonical


@pytest.mark.parametrize("attrs", (None, "text;other;than"))
def test_not_is_canonical(attrs):
    f = EnsemblGffRecord(attrs=attrs)
    assert not f.is_canonical


@pytest.mark.parametrize(
    "val",
    (
        [(10, 48)],
        [(9, 20), (29, 45), (59, 70)],
        numpy.array([(9, 20), (29, 45), (59, 70)]),
    ),
)
def tess_gff_record_size_nonzero(val):
    record = EnsemblGffRecord(spans=val)
    assert record.size == 38


def test_gff_record_hashing():
    name = "abcd"
    record = EnsemblGffRecord(name=name)
    assert hash(record) == hash(name)
    v = {record: 21}
    assert v[name] == 21
    n = {name: 21}
    assert v == n


@pytest.fixture
def ensembl_gff_records(DATA_DIR):
    records, _ = custom_gff_parser(DATA_DIR / "c_elegans_WS199_shortened.gff3", 0)
    return records


@pytest.fixture
def non_canonical_related(ensembl_gff_records):
    return make_gene_relationships(ensembl_gff_records.values())


@pytest.fixture
def canonical_related(ensembl_gff_records):
    transcript = ensembl_gff_records["transcript:B0019.1"]
    transcript.attrs = f"Ensembl_canonical;{transcript.attrs}"
    return ensembl_gff_records, make_gene_relationships(ensembl_gff_records.values())


def test_make_gene_relationships(ensembl_gff_records):
    # make the mRNA is_canonical
    transcript = ensembl_gff_records["transcript:B0019.1"]
    transcript.attrs = f"Ensembl_canonical;{transcript.attrs}"
    # at this point the related CDS is not canonical
    assert not ensembl_gff_records["cds:B0019.1"].is_canonical
    related = make_gene_relationships(ensembl_gff_records.values())
    got = {c.is_canonical for c in related["gene:WBGene00000138"]}
    assert got == {True}
    # the related CDS is now canonical
    assert ensembl_gff_records["cds:B0019.1"].is_canonical


def test_featuredb(canonical_related):
    records, related = canonical_related
    db = EnsemblGffDb(source=":memory:")
    db.add_records(records=records.values(), gene_relations=related)
    cds = list(
        db.get_feature_children(name="WBGene00000138", biotype="cds", is_canonical=True)
    )[0]
    assert cds["name"] == "B0019.1"


def test_featuredb_num_records(canonical_related):
    records, related = canonical_related
    db = EnsemblGffDb(source=":memory:")
    assert db.num_records() == 0
    db.add_records(records=records.values(), gene_relations=related)
    assert db.num_records() == 11


def test_make_annotation_db(DATA_DIR, tmp_path):
    src = DATA_DIR / "c_elegans_WS199_shortened.gff3"
    dest = tmp_path / _ANNOTDB_NAME
    make_annotation_db((src, dest))
    got = EnsemblGffDb(source=dest)
    assert got.num_records() == 11


@pytest.mark.parametrize("table_name", tuple(EnsemblGffDb._index_columns))
def test_indexing(canonical_related, table_name):
    records, related = canonical_related
    db = EnsemblGffDb(source=":memory:")
    db.add_records(records=records.values(), gene_relations=related)
    col = EnsemblGffDb._index_columns[table_name][0]
    expect = ("index", f"{col}_index", table_name)
    db.make_indexes()
    sql_template = (
        f"SELECT * FROM sqlite_master WHERE type = 'index' AND "  # nosec B608
        f"tbl_name = {table_name!r} and name = '{col}_index'"  # nosec B608
    )

    result = db._execute_sql(sql_template).fetchone()
    got = tuple(result)[:3]
    assert got == expect
