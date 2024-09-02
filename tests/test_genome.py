import pathlib

import numpy
import pytest
from cogent3 import make_unaligned_seqs
from numpy.testing import assert_allclose

from ensembl_tui import _genome as elt_genome
from ensembl_tui import _storage_mixin as elt_mixin


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
    db = elt_genome.EnsemblGffDb(source=":memory:")
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
    return elt_genome.SeqsDataHdf5(
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
    "name,start,stop",
    (("s1", 3, 7), ("s1", 3, None), ("s1", None, 7), ("s2", 2, 4)),
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
    genome = elt_genome.Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
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
    genome = elt_genome.Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
    seq = genome.get_seq(seqid="s1", namer=namer)
    assert seq
    genome.close()
    with pytest.raises(OSError):
        genome.get_seq(seqid="s1")


@pytest.mark.parametrize("seqid", ("s1", "s2"))
def test_get_seq_num_annotations_correct(
    small_h5_genome,
    small_annotdb,
    small_coll,
    seqid,
    namer,
):
    gen_seqs_db, small_data = small_h5_genome
    genome = elt_genome.Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
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
    small_h5_genome,
    small_annotdb,
    small_coll,
    seqid,
    feature_name,
    start,
    stop,
    namer,
):
    gen_seqs_db, small_data = small_h5_genome
    genome = elt_genome.Genome(species="dodo", seqs=gen_seqs_db, annots=small_annotdb)
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
    got = elt_genome.get_gene_table_for_species(
        annot_db=small_annotdb,
        limit=None,
        species="none",
    )
    assert isinstance(got, Table)
    assert len(got) > 0


def test_get_species_summary(small_annotdb):
    from cogent3.util.table import Table

    got = elt_genome.get_species_summary(annot_db=small_annotdb, species="none")
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
        genome.add_record(data["s1"][:-2], "s1")


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
    convert = elt_genome.str2arr(moltype="dna")
    assert (
        genome.get_seq_arr(seqid=name, start=start, stop=stop) == convert(expect)
    ).all()


def test_pickling_round_trip(small_data, tmp_path):
    import pickle  # nosec B403

    path = tmp_path / f"small.{elt_genome.SEQ_STORE_NAME}"
    kwargs = dict(source=path, species="human")
    genome = elt_genome.SeqsDataHdf5(mode="w", **kwargs)
    genome.add_records(records=small_data.items())
    with pytest.raises(NotImplementedError):
        pickle.dumps(genome)  # nosec B301

    ro = elt_genome.SeqsDataHdf5(mode="r", **kwargs)
    assert ro.get_seq_str(seqid="s1") == small_data["s1"]
    unpkl = pickle.loads(pickle.dumps(ro))  # nosec B301
    got = unpkl.get_seq_str(seqid="s1")
    assert got == small_data["s1"]


def test_species_setting(small_data, tmp_path):
    path = tmp_path / f"small.{elt_genome.SEQ_STORE_NAME}"
    kwargs = dict(source=path, species="human")
    genome = elt_genome.SeqsDataHdf5(mode="w", **kwargs)
    genome.add_records(records=small_data.items())
    genome.close()

    genome = elt_genome.SeqsDataHdf5(mode="r", source=path)
    # note that species are converted into the Ensembl db prefix
    assert genome.species == "homo_sapiens"
    with pytest.raises(ValueError):
        _ = elt_genome.SeqsDataHdf5(mode="r", source=path, species="cat")


def test_has_of_seqsdata(h5_genome):
    assert hash(h5_genome) == id(h5_genome)


def test_tidying_stableids_in_gff3():
    orig = (
        "ID=Transcript:ENST00000461467;Parent=Gene:ENSG00000237613;Name=FAM138A-202;bio"
    )
    expect = (
        "ID=transcript:ENST00000461467;Parent=gene:ENSG00000237613;Name=FAM138A-202;bio"
    )
    assert elt_genome.tidy_gff3_stableids(orig) == expect


def test_custom_gff3_parser(DATA_DIR):
    path = DATA_DIR / "c_elegans_WS199_shortened.gff3"
    records, _ = elt_genome.custom_gff_parser(path, 0)

    rel = elt_genome.make_gene_relationships(records.values())
    children = rel["gene:WBGene00000138"]
    # as the records are hashable by their .name attribute, we can just
    # check returned value against their names
    assert children == {"cds:B0019.1", "transcript:B0019.1"}
    # check that multi row records have the correct spans, start, stop and strand
    assert_allclose(
        records["cds:B0019.1"].spans,
        numpy.array([(9, 20), (29, 45), (59, 70)]),
    )
    assert records["cds:B0019.1"].start == 9
    assert records["cds:B0019.1"].stop == 70
    assert records["cds:B0019.1"].strand == "-"


def test_gff_record_size(DATA_DIR):
    merged, _ = elt_genome.custom_gff_parser(
        DATA_DIR / "c_elegans_WS199_shortened.gff3",
        0,
    )
    # record CDS:B0019.1 has spans [(9, 20), (29, 45), (59, 70)] which sum to 38
    starts, stops = numpy.array([(9, 20), (29, 45), (59, 70)]).T
    expect = (stops - starts).sum()
    assert merged["cds:B0019.1"].size == expect


@pytest.mark.parametrize("val", ((), [], numpy.array([]), None))
def tess_gff_record_size_zero(val):
    record = elt_genome.EnsemblGffRecord(spans=val)
    assert record.size == 0


@pytest.mark.parametrize(
    "attrs",
    ("Ensembl_canonical", "text;other;Ensembl_canonical;than"),
)
def test_is_canonical(attrs):
    f = elt_genome.EnsemblGffRecord(attrs=attrs)
    assert f.is_canonical


@pytest.mark.parametrize("attrs", (None, "text;other;than"))
def test_not_is_canonical(attrs):
    f = elt_genome.EnsemblGffRecord(attrs=attrs)
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
    record = elt_genome.EnsemblGffRecord(spans=val)
    assert record.size == 38


@pytest.mark.parametrize("symbol,expect", (("Name=ATAD3B;", "ATAD3B"), ("", None)))
def test_gff_record_symbol(symbol, expect):
    data = {"attrs": f"ID=gene:ENSG00000160072;{symbol}biotype=protein_coding;"}
    record = elt_genome.EnsemblGffRecord(**data)
    assert record.symbol == expect


@pytest.mark.parametrize(
    "descr,expect",
    (
        (
            "description=...domain containing 3B [Source:HGNC Symbol%3BAcc:HGNC:24007];",
            "...domain containing 3B [Source:HGNC Symbol%3BAcc:HGNC:24007]",
        ),
        ("", None),
    ),
)
def test_gff_record_description(descr, expect):
    data = {"attrs": f"ID=gene:ENSG00000160072;{descr}biotype=protein_coding;"}
    record = elt_genome.EnsemblGffRecord(**data)
    assert record.description == expect


def test_gff_record_hashing():
    name = "abcd"
    record = elt_genome.EnsemblGffRecord(name=name)
    assert hash(record) == hash(name)
    v = {record: 21}
    assert v[name] == 21
    n = {name: 21}
    assert v == n


@pytest.mark.parametrize("exclude_null", (True, False))
def test_gff_record_to_record(exclude_null):
    data = {
        "seqid": "s1",
        "name": "gene-01",
        "biotype": "gene",
        "spans": [(1, 3), (7, 9)],
        "start": 1,
        "stop": 9,
        "strand": "+",
    }
    all_fields = (
        {} if exclude_null else {s: None for s in elt_genome.EnsemblGffRecord.__slots__}
    )
    all_fields.pop("_is_updated", None)
    record = elt_genome.EnsemblGffRecord(**data)
    got = record.to_record(exclude_null=exclude_null)
    expect = all_fields | data
    expect.pop("spans")
    got_spans = got.pop("spans")
    assert got == expect
    assert numpy.array_equal(elt_mixin.blob_to_array(got_spans), data["spans"])


@pytest.mark.parametrize("exclude_null", (True, False))
def test_gff_record_to_record_selected_fields(exclude_null):
    data = {
        "seqid": "s1",
        "name": "gene-01",
        "start": None,
        "stop": None,
    }
    fields = list(data)
    record = elt_genome.EnsemblGffRecord(**data)
    got = record.to_record(fields=fields, exclude_null=exclude_null)
    expect = {f: data[f] for f in fields if data[f] is not None or not exclude_null}
    assert got == expect


@pytest.fixture
def ensembl_gff_records(DATA_DIR):
    records, _ = elt_genome.custom_gff_parser(
        DATA_DIR / "c_elegans_WS199_shortened.gff3",
        0,
    )
    return records


@pytest.fixture
def non_canonical_related(ensembl_gff_records):
    return elt_genome.make_gene_relationships(ensembl_gff_records.values())


@pytest.fixture
def canonical_related(ensembl_gff_records):
    transcript = ensembl_gff_records["transcript:B0019.1"]
    transcript.attrs = f"Ensembl_canonical;{transcript.attrs}"
    return ensembl_gff_records, elt_genome.make_gene_relationships(
        ensembl_gff_records.values(),
    )


def test_make_gene_relationships(ensembl_gff_records):
    # make the mRNA is_canonical
    transcript = ensembl_gff_records["transcript:B0019.1"]
    transcript.attrs = f"Ensembl_canonical;{transcript.attrs}"
    # at this point the related CDS is not canonical
    assert not ensembl_gff_records["cds:B0019.1"].is_canonical
    related = elt_genome.make_gene_relationships(ensembl_gff_records.values())
    got = {c.is_canonical for c in related["gene:WBGene00000138"]}
    assert got == {True}
    # the related CDS is now canonical
    assert ensembl_gff_records["cds:B0019.1"].is_canonical


def test_featuredb(canonical_related):
    records, related = canonical_related
    db = elt_genome.EnsemblGffDb(source=":memory:")
    db.add_records(records=records.values(), gene_relations=related)
    cds = list(
        db.get_feature_children(
            name="WBGene00000138",
            biotype="cds",
            is_canonical=True,
        ),
    )[0]
    assert cds["name"] == "B0019.1"


def test_featuredb_num_records(canonical_related):
    records, related = canonical_related
    db = elt_genome.EnsemblGffDb(source=":memory:")
    assert db.num_records() == 0
    db.add_records(records=records.values(), gene_relations=related)
    assert db.num_records() == 11


def test_make_annotation_db(DATA_DIR, tmp_path):
    src = DATA_DIR / "c_elegans_WS199_shortened.gff3"
    dest = tmp_path / elt_genome.ANNOT_STORE_NAME
    elt_genome.make_annotation_db((src, dest))
    got = elt_genome.EnsemblGffDb(source=dest)
    assert got.num_records() == 11


def test_get_features_matching(canonical_related):
    records, related = canonical_related
    db = elt_genome.EnsemblGffDb(source=":memory:")
    db.add_records(records=records.values(), gene_relations=related)
    got = list(db.get_features_matching(biotype="cds"))
    assert got[0]["name"] == "B0019.1"
    assert got[0]["biotype"] == "cds"


@pytest.mark.parametrize("table_name", tuple(elt_genome.EnsemblGffDb._index_columns))
def test_indexing(canonical_related, table_name):
    records, related = canonical_related
    db = elt_genome.EnsemblGffDb(source=":memory:")
    db.add_records(records=records.values(), gene_relations=related)
    col = elt_genome.EnsemblGffDb._index_columns[table_name][0]
    expect = ("index", f"{col}_index", table_name)
    db.make_indexes()
    sql_template = (
        f"SELECT * FROM sqlite_master WHERE type = 'index' AND "  # nosec B608
        f"tbl_name = {table_name!r} and name = '{col}_index'"  # nosec B608
    )

    result = db._execute_sql(sql_template).fetchone()
    got = tuple(result)[:3]
    assert got == expect


def test_get_feature_parent(canonical_related):
    records, related = canonical_related
    db = elt_genome.EnsemblGffDb(source=":memory:")
    db.add_records(records=records.values(), gene_relations=related)
    got = list(db.get_feature_parent(name="B0019.1"))[0]
    assert got["name"] == "WBGene00000138"


def test_get_feature_children(canonical_related):
    records, related = canonical_related
    db = elt_genome.EnsemblGffDb(source=":memory:")
    db.add_records(records=records.values(), gene_relations=related)
    got = list(
        db.get_feature_children(
            name="WBGene00000138",
            biotype="cds",
            is_canonical=True,
        ),
    )[0]
    assert got["name"] == "B0019.1"


def test_add_feature():
    db = elt_genome.EnsemblGffDb(source=":memory:")
    feature = elt_genome.EnsemblGffRecord(
        start=2,
        stop=3,
        seqid="s0",
        name="demo",
        spans=[(2, 3)],
        biotype="gene",
    )
    db.add_feature(feature=feature)
    got = list(db.get_features_matching(seqid="s0"))
    assert len(got) == 1
    assert got[0]["name"] == "demo"


@pytest.fixture(params=("\n", "\r\n"))
def fasta_data(DATA_DIR, tmp_path, request):
    data = pathlib.Path(DATA_DIR / "c_elegans_WS199_shortened.fasta").read_text()
    outpath = tmp_path / "demo.fasta"
    outpath.write_text(data, newline=request.param)
    return outpath


def test_faster_fasta(fasta_data):
    from cogent3.parse.fasta import MinimalFastaParser

    from ensembl_tui._faster_fasta import bytes_to_array, quicka_parser

    expect = {
        n: bytes_to_array(s.encode("utf8")) for n, s in MinimalFastaParser(fasta_data)
    }
    got = dict(quicka_parser(fasta_data))
    assert (got["I"] == expect["I"]).all()


def test_gff_parse_merge(DATA_DIR):
    records, _ = elt_genome.custom_gff_parser(
        DATA_DIR / "gene-multi-transcript.gff3",
        0,
    )
    related = elt_genome.make_gene_relationships(list(records.values()))
    homologs = related["gene:ENSG00000160072"]
    assert len(homologs) == 4
    # we know that cds:ENSP00000500094 and transcript:ENST00000673477 are canonical
    assert records["cds:ENSP00000500094"].is_canonical
    assert records["transcript:ENST00000673477"].is_canonical
    # and cds:ENSP00000311766 and transcript:ENST00000308647 are not
    assert not records["cds:ENSP00000311766"].is_canonical
    assert not records["transcript:ENST00000308647"].is_canonical
    assert {
        "transcript:ENST00000308647",
        "cds:ENSP00000311766",
        "transcript:ENST00000673477",
        "cds:ENSP00000500094",
    } == homologs


def test_genome_segment():
    segment = elt_genome.genome_segment(
        species="abcd_efg",
        seqid="1",
        start=20,
        stop=40,
        strand="+",
    )
    assert segment.unique_id == "abcd_efg-1-20-40"
    segment = elt_genome.genome_segment(
        species="abcd_efg",
        seqid="1",
        start=20,
        stop=40,
        strand="+",
        unique_id="gene:NICE",
    )
    assert segment.unique_id == "NICE"


def test_get_gene_segments(small_annotdb):
    segments = elt_genome.get_gene_segments(annot_db=small_annotdb, species="dodo")
    assert len(segments) == len(
        list(small_annotdb.get_features_matching(biotype="gene")),
    )
    assert {s.unique_id for s in segments} == {"gene-01", "gene-02"}


def test_get_gene_segments_stableids(small_annotdb):
    segments = elt_genome.get_gene_segments(
        annot_db=small_annotdb,
        species="dodo",
        stableids=["gene-02"],
    )
    assert len(segments) == 1
    segment = segments[0]
    assert segment.unique_id == "gene-02"
    assert segment.source == "gene-02"
