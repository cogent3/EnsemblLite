import numpy
import pytest

from cogent3 import make_seq

from ensembl_lite._aligndb import (
    AlignDb,
    AlignRecordType,
    GapPositions,
    get_alignment,
)
from ensembl_lite._genomedb import CompressedGenomeSeqsDb
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


def make_records(start, end):
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


@pytest.fixture
def small_records():
    records = make_records(1, 5)
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


@pytest.mark.parametrize("data", ("AB---CD--EF", "---ABCD--EF", "ABCD---EF--"))
@pytest.mark.parametrize("index", range(6))  # the ungapped sequence is 6 long
def test_gapped_convert_seq2aln(data, index):
    # converting a sequence index to alignment index
    ungapped = data.replace("-", "")
    seq = make_seq(data, moltype="text")
    g, s = seq_to_gap_coords(seq)
    gaps = GapPositions(g, len(seq))
    idx = gaps.from_seq_to_align_index(index)
    assert data[idx] == ungapped[index]


@pytest.mark.parametrize("data", ("AC--GTA-TG", "--GTA-TGAA", "AC--GTA---"))
@pytest.mark.parametrize("index", range(10))
def test_gapped_convert_aln2seq(data, index):
    seq = make_seq(data, moltype="dna")
    g, s = seq_to_gap_coords(seq)
    gaps = GapPositions(g, len(seq))
    expect = data[:index].replace("-", "")
    idx = gaps.from_align_to_seq_index(index)
    assert idx == len(expect)


def test_gapped_convert_aln2seq_invalid():
    seq = make_seq("AC--GTA-TG", moltype="dna")
    g, s = seq_to_gap_coords(seq)
    gaps = GapPositions(g, len(seq))
    with pytest.raises(NotImplementedError):
        gaps.from_align_to_seq_index(-1)


# fixture to make synthetic GenomeSeqsDb and alignment db
# based on a given alignment
@pytest.fixture
def genomedbs_aligndb(small_records):
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
    got = get_alignment(align_db, genomes, species="mouse", coord_name="s2")
    orig = small_seqs()[1:5]
    assert got.to_dict() == orig.to_dict()


@pytest.mark.parametrize(
    "kwargs",
    (dict(species="dodo", coord_name="s2"), dict(species="mouse", coord_name="s222")),
)
def test_building_alignment_invalid_details(genomedbs_aligndb, kwargs):
    genomes, align_db = genomedbs_aligndb
    with pytest.raises(ValueError):
        get_alignment(align_db, genomes, **kwargs)


@pytest.mark.parametrize(
    "invalid_slice",
    (slice(None, None, -1), slice(None, -1, None), slice(-1, None, None)),
)
def test_gap_pos_invalid_slice(invalid_slice):
    gp = GapPositions(numpy.array([[1, 3]], dtype=numpy.int32), 20)
    with pytest.raises(NotImplementedError):
        gp[invalid_slice]


@pytest.mark.parametrize(
    "slice",
    (
        slice(3, 7),
        slice(20, None),
    ),
)
def test_no_gaps_in_slice(slice):
    # aligned length is 25
    seq_length = 20
    gap_length = 5
    gp = GapPositions(
        gaps=numpy.array([[10, gap_length]], dtype=numpy.int32), seq_length=seq_length
    )
    got = gp[slice]
    assert not len(got.gaps)
    start = slice.start or 0
    stop = slice.stop or (seq_length + gap_length)
    assert len(got) == stop - start


def test_len_gapped():
    seq_length = 20
    gap_length = 5

    gp = GapPositions(
        gaps=numpy.array([[10, gap_length]], dtype=numpy.int32), seq_length=seq_length
    )
    assert len(gp) == (seq_length + gap_length)


def test_all_gaps_in_slice():
    # slicing GapPositions
    # sample seq 1
    data = "AC--GTA-TG"
    seq = make_seq(data, moltype="dna")
    g, s = seq_to_gap_coords(seq)
    gp = GapPositions(g, len(data.replace("-", "")))
    sl = slice(1, 9)

    got = gp[sl]
    expect_gaps, expect_seq = seq_to_gap_coords(make_seq(data[sl], moltype="dna"))
    assert (got.gaps == expect_gaps).all()
    assert got.seq_length == 5


@pytest.mark.parametrize(
    "data",
    (
        "----GTA-TG",
        "AC--GTA---",
        "AC--GTA-TG",
        "A-C-G-T-A-",
        "-A-C-G-T-A",
        "ACGTAACGTA",
        "----------",
    ),
)
@pytest.mark.parametrize(
    "slice",
    (
        slice(0, 2, None),
        slice(0, 5, None),
        slice(0, 7, None),
        slice(0, 8, None),
        slice(1, 9, None),
        slice(2, 3, None),
        slice(2, 4, None),
        slice(2, 9, None),
        slice(3, 8, None),
        slice(3, 9, None),
        slice(4, 6, None),
        slice(4, 9, None),
        slice(6, 9, None),
        slice(8, 10, None),
    ),
)
def test_variant_slices(data, slice):
    seq = make_seq(data, moltype="dna")
    g, s = seq_to_gap_coords(seq)
    gaps = GapPositions(g, len(s))
    orig = gaps.gaps.copy()
    got = gaps[slice]

    expect_gaps, expect_seq = seq_to_gap_coords(make_seq(data[slice], moltype="dna"))
    assert got.seq_length == len(expect_seq)
    assert (got.gaps == expect_gaps).all()
    # make sure original data unmodified
    assert (orig == gaps.gaps).all()


def make_sample():
    seqs = small_seqs()
    # we will reverse complement the s2 genome compared to the original
    # this means our coordinates for alignment records from that genome
    # also need to be rc'ed
    species = seqs.info.species
    genomes = {}
    for seq in seqs.seqs:
        name = seq.name
        seq = seq.data.degap()
        if seq.name == "s2":
            seq = seq.rc()
            s2_genome = str(seq)
        genome = CompressedGenomeSeqsDb(source=":memory:", species=species[seq.name])
        genome.add_records(records=[(name, str(seq))])
        genomes[species[name]] = genome

    # define a slice
    start, end = 1, 12
    align_records = make_records(start, end)
    # identify the rc segment for s2, which will be reverse complemented
    # relative to "genome"
    s2 = seqs.get_gapped_seq("s2")
    selected = s2[start:end].rc().degap()
    start = s2_genome.find(str(selected))
    end = start + len(selected)
    for record in align_records:
        if record["coord_name"] == "s2":
            record["start"] = start
            record["end"] = end
            record["strand"] = "-"
            break
    align_db = AlignDb(source=":memory:")
    align_db.add_records(records=align_records)

    return genomes, align_db


def test_select_alignment_rc():
    expect = small_seqs()[1:12]
    # one sequence is stored in reverse complement
    genomes, align_db = make_sample()
    got = get_alignment(
        align_db=align_db, genomes=genomes, species="human", coord_name="s1"
    )
    assert got.to_dict() == expect.to_dict()


@pytest.mark.parametrize(
    "coord",
    (
        ("human", "s1", None, 11),  # finish within
        ("human", "s1", 3, None),  # start within
        ("human", "s1", 3, 9),  # within
        ("human", "s1", 3, 13),  # extends past
    ),
)
def test_align_db_get_records(coord):
    kwargs = dict(zip(("species", "coord_name", "start", "end"), coord))
    # records are, we should get a single hit from each query
    # [('blah', 0, 'human', 's1', 1, 12, '+', array([], dtype=int32)),
    _, align_db = make_sample(two_aligns=True)
    got = list(align_db.get_records_matching(**kwargs))
    assert len(got) == 1


@pytest.mark.parametrize(
    "coord",
    (
        ("human", "s1"),
        ("mouse", "s2"),
        ("dog", "s3"),
    ),
)
def test_align_db_get_records_required_only(coord):
    kwargs = dict(zip(("species", "coord_name"), coord))
    # two hits for each species
    _, align_db = make_sample(two_aligns=True)
    got = list(align_db.get_records_matching(**kwargs))
    assert len(got) == 2


@pytest.mark.parametrize(
    "coord",
    (
        ("human", "s2"),
        ("mouse", "xx"),
        ("blah", "s3"),
    ),
)
def test_align_db_get_records_no_matches(coord):
    kwargs = dict(zip(("species", "coord_name"), coord))
    # no hits at all
    _, align_db = make_sample()
    got = list(align_db.get_records_matching(**kwargs))
    assert not len(got)
