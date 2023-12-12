import pytest

from cogent3 import make_seq

from ensembl_lite._aligndb import (
    AlignDb,
    AlignRecordType,
    GapPositions,
    get_alignment,
)
from ensembl_lite.convert import seq_to_gap_coords


def small_seqs():
    from cogent3 import make_aligned_seqs

    seqs = {
        "s1": "GTTGAAGTAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGTTCCAAATAATGAA",
        "s3": "GCTGAAGTAGTGGAAGTTGCAAAT---GAA",
    }
    return make_aligned_seqs(
        data=seqs,
        moltype="dna",
        array_align=False,
        info=dict(species=dict(s1="human", s2="mouse", s3="dog")),
    )


@pytest.fixture
def small_records():
    start, end = 1, 5
    aln = small_seqs()[start:end]
    records = []
    species = aln.info.species
    for seq in aln.seqs:
        gs = seq.get_gapped_seq()
        c, s = seq_to_gap_coords(gs)
        record = AlignRecordType(
            source="blah",
            species=species[seq.name],
            block_id=0,
            coord_name=seq.name,
            start=seq.map.start,
            end=seq.map.end,
            strand="+",
            gap_spans=c,
        )
        records.append(record)
    return records


def test_aligndb_records_match_input(small_records):
    import copy

    orig_records = copy.deepcopy(small_records)
    db = AlignDb(source=":memory:")
    db.add_records(records=small_records)
    got = list(db.get_records_matching(species="human", coord_name="s1"))[0]
    for g, o in zip(got, orig_records):
        g_spans = g.pop("gap_spans")
        o_spans = o.pop("gap_spans")
        assert g == o
        assert (g_spans == o_spans).all()


@pytest.mark.parametrize("data", ("AC---GT--TT", "---GT--TT", "AC---GT--"))
def test_gapped_convert_seq2aln(data):
    seq = make_seq(data, moltype="dna")
    g, s = seq_to_gap_coords(seq)
    gaps = GapPositions(g, len(seq))
    idx = gaps.from_seq_to_align_index(3)
    assert seq[idx] == data[idx]


def test_gapped_convert_aln2seq():
    seq = make_seq("AC---GT--TT", moltype="dna")
    g, s = seq_to_gap_coords(seq)
    gaps = GapPositions(g, len(seq))
    idx = gaps.from_align_to_seq_index(6)
    assert idx == 3


# fixture to make synthetic GenomeSeqsDb and alignment db
# based on a given alignment
@pytest.fixture
def genomedbs_aligndb(small_records):
    from ensembl_lite._genomedb import CompressedGenomeSeqsDb

    align_db = AlignDb(source=":memory:")
    align_db.add_records(records=small_records)
    seqs = small_seqs().degap()
    species = seqs.info.species
    data = seqs.to_dict()
    genomes = {}
    for name, seq in data.items():
        genome = CompressedGenomeSeqsDb(source=":memory:", species=name)
        genome.add_records(records=[(name, seq)])
        genomes[species[name]] = genome

    return genomes, align_db


def test_building_alignment(genomedbs_aligndb):
    genomes, align_db = genomedbs_aligndb
    # todo add test for incorrect species name, incorrect coord name
    got = get_alignment(align_db, genomes, species="mouse", coord_name="s2")
    orig = small_seqs()[1:5]
    assert got.to_dict() == orig.to_dict()
