import pathlib

import pytest

from ensembl_tui import _genome as eti_genome


def test_get_gene_segments(yeast_db):
    segs = eti_genome.get_gene_segments(
        annot_db=yeast_db,
        species="saccharomyces_cerevisiae",
        limit=10,
    )
    assert len(segs) == 10


def test_get_seq_feature_seq_correct(yeast):
    seq = yeast.get_seq(seqid="III", start=309069, stop=310155)
    raw_seq = str(seq)
    assert raw_seq.startswith("ATGCTTTACCCAGAAAAATTTCA")  # expected from ensembl.org
    assert raw_seq.endswith("ATAAGAAATTCCATAAATAG")
    got = next(iter(seq.get_features()))
    feat_seq = got.get_slice(apply_name=True)
    assert feat_seq.name == "YCR105W"


def test_get_gene_table_for_species(yeast_db):
    from cogent3.util.table import Table

    # we do not check values here, only the Type and that we have > 0 records
    got = eti_genome.get_gene_table_for_species(
        annot_db=yeast_db,
        limit=None,
        species="saccharomyces_cerevisiae",
    )
    assert isinstance(got, Table)
    assert len(got) > 0


def test_get_species_summary(yeast_db):
    from cogent3.util.table import Table

    got = eti_genome.get_species_summary(annot_db=yeast_db)
    # we do not check values here, only the Type and that we have > 0 records
    assert isinstance(got, Table)
    assert len(got) > 0


def test_genome_coord_names(yeast_db):
    counts = yeast_db.count_distinct(seqid=True)
    assert counts.shape[0] == 17


def test_empty_hdf5_genome_coord_names(yeast_seqs):
    assert yeast_seqs.get_coord_names()
    assert "III" in yeast_seqs.get_coord_names()


def test_get_seq(yeast_seqs):
    mt = yeast_seqs.get_seq_str(seqid="Mito")
    assert len(mt) == 85779  # number for ensembl.org


def test_pickling_round_trip(yeast_seqs):
    import pickle  # nosec B403

    ro = yeast_seqs
    kwargs = {"seqid": "Mito", "start": 200, "stop": 220}
    small_seq = yeast_seqs.get_seq_str(**kwargs)
    # expected value from one run
    assert small_seq == "AAAGATAAAAAAAATAATGT"

    unpkl = pickle.loads(pickle.dumps(ro))  # nosec B301  # noqa: S301
    got = unpkl.get_seq_str(**kwargs)
    assert got == small_seq


def test_species_setting(yeast_seqs):
    # note that species are converted into the Ensembl db prefix
    assert yeast_seqs.species == "saccharomyces_cerevisiae"
    with pytest.raises(ValueError):  # noqa: PT011
        _ = eti_genome.SeqsDataHdf5(mode="r", source=yeast_seqs.source, species="cat")


def test_hash_of_seqs_data(yeast_seqs):
    assert hash(yeast_seqs) == id(yeast_seqs)


def test_gene_description(worm_db):
    gene = next(iter(worm_db.genes.get_by_stable_id(stable_id="WBGene00010790")))
    assert "alcohol dehydrogenase" in gene.description.lower()


def test_featuredb_num_records(worm_db):
    db = worm_db
    num_repeats = db.repeats.num_records()
    num_genes = db.genes.num_records()
    assert db.num_records() == num_repeats + num_genes
    assert num_repeats > 1000
    assert num_genes > 1000
    num_biotypes = db.biotypes.num_records()
    assert num_biotypes > 10


def test_get_features_matching(yeast_db):
    got = list(yeast_db.get_features_matching(biotype="rRNA"))
    assert len(got) > 20  # number from ensembl.org


def test_get_feature_child_parent(worm_db):
    db = worm_db.genes
    gene = next(iter(db.get_by_stable_id(stable_id="WBGene00010790")))
    transcript = next(iter(db.get_feature_children(gene)))
    cds = next(iter(db.get_feature_children(transcript)))
    assert gene.canonical_transcript_id == transcript.transcript_id == cds.transcript_id
    got = db.get_feature_parent(cds)
    assert isinstance(got, type(transcript))
    assert got.stable_id == transcript.stable_id
    got = db.get_feature_parent(transcript)
    assert isinstance(got, type(gene))
    assert got.stable_id == gene.stable_id


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


def test_genome_segment():
    segment = eti_genome.genome_segment(
        species="abcd_efg",
        seqid="1",
        start=20,
        stop=40,
        strand="+",
    )
    assert segment.unique_id == "abcd_efg-1-20-40"
    segment = eti_genome.genome_segment(
        species="abcd_efg",
        seqid="1",
        start=20,
        stop=40,
        strand="+",
        unique_id="gene:NICE",
    )
    assert segment.unique_id == "NICE"


def test_get_gene_segments(worm_db):
    segments = eti_genome.get_gene_segments(
        annot_db=worm_db,
        species="caenorhabditis_elegans",
    )
    assert len(segments) == len(
        list(worm_db.get_features_matching(biotype="protein_coding")),
    )
    assert "WBGene00010790" in {s.unique_id for s in segments}


def test_get_gene_segments_stableids(worm_db):
    segments = eti_genome.get_gene_segments(
        annot_db=worm_db,
        species="caenorhabditis_elegans",
        stableids=["WBGene00010790"],
    )
    assert len(segments) == 1
    segment = segments[0]
    assert segment.unique_id == "WBGene00010790"
    assert segment.source == "WBGene00010790"


def test_get_features(yeast):
    features = list(yeast.get_features(biotype="rRNA", limit=10))
    assert len(features) == 10


def test_get_ids_for_biotype(yeast):
    features = list(yeast.get_ids_for_biotype(biotype="rRNA", limit=10))
    assert len(features) == 10
