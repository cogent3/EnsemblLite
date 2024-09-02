import numpy
import pytest

from ensembl_tui import _align as elt_align
from ensembl_tui import _genome as elt_genome


def small_seqs():
    from cogent3 import make_aligned_seqs

    seqs = {
        "s1": "GTTGAAGTAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGTTCCAAATAATGAA",
        "s3": "GCTGAAGTAGTGGAAGTTGCAAAT---GAA",
    }
    seqs = make_aligned_seqs(
        data=seqs,
        moltype="dna",
        array_align=False,
        info=dict(species=dict(s1="human", s2="mouse", s3="dog")),
    )
    annot_db = elt_genome.EnsemblGffDb(source=":memory:")
    annot_db.add_feature(seqid="s1", biotype="gene", name="not-on-s2", spans=[(4, 7)])
    annot_db.add_feature(
        seqid="s2",
        biotype="gene",
        name="includes-s2-gap",
        spans=[(2, 6)],
    )
    annot_db.add_feature(
        seqid="s3",
        biotype="gene",
        name="includes-s3-gap",
        spans=[(22, 27)],
    )
    seqs.annotation_db = annot_db
    return seqs


def make_records(start, end, block_id):
    aln = small_seqs()[start:end]
    records = []
    species = aln.info.species
    for seq in aln.seqs:
        seqid, seq_start, seq_end, seq_strand = seq.data.parent_coordinates()
        gs = seq.get_gapped_seq()
        imap, s = gs.parse_out_gaps()
        if imap.num_gaps:
            gap_spans = numpy.array(
                [imap.gap_pos, imap.get_gap_lengths()],
                dtype=numpy.int32,
            ).T
        else:
            gap_spans = numpy.array([], dtype=numpy.int32)
        record = elt_align.AlignRecord(
            source="blah",
            species=species[seq.name],
            block_id=block_id,
            seqid=seqid,
            start=seq_start,
            stop=seq_end,
            strand="-" if seq_strand == -1 else "+",
            gap_spans=gap_spans,
        )
        records.append(record)
    return records


@pytest.fixture
def small_records():
    records = make_records(1, 5, 0)
    return records


def test_aligndb_records_match_input(small_records):
    import copy

    orig_records = copy.deepcopy(small_records)

    db = elt_align.AlignDb(source=":memory:")
    db.add_records(records=small_records)
    got = list(db.get_records_matching(species="human", seqid="s1"))[0]
    assert got == set(orig_records)


def test_aligndb_records_skip_duplicated_block_ids(small_records):
    db = elt_align.AlignDb(source=":memory:")
    db.add_records(records=small_records)
    orig = list(db.get_records_matching(species="human", seqid="s1"))
    db.add_records(records=small_records)
    got = list(db.get_records_matching(species="human", seqid="s1"))
    assert len(got) == len(orig)


def _find_nth_gap_index(data: str, n: int) -> int:
    num = -1
    for i, c in enumerate(data):
        if c == "-":
            num += 1
        if num == n:
            return i
    raise ValueError(f"{data=}, {n=}")


def _get_expected_seqindex(data: str, align_index: int) -> int:
    # compute the expected seqindex
    refseq = data.replace("-", "")
    got = data[align_index:].lstrip("-")
    return refseq.find(got[0]) if got else len(refseq)


# fixture to make synthetic GenomeSeqsDb and alignment db
# based on a given alignment
@pytest.fixture
def genomedbs_aligndb(small_records):
    align_db = elt_align.AlignDb(source=":memory:")
    align_db.add_records(records=small_records)
    seqs = small_seqs().degap()
    species = seqs.info.species
    data = seqs.to_dict()
    genomes = {}
    for name, seq in data.items():
        genome = elt_genome.SeqsDataHdf5(
            source=f"{name}",
            species=species[name],
            mode="w",
            in_memory=True,
        )
        genome.add_records(records=[(name, seq)])
        genomes[species[name]] = elt_genome.Genome(
            seqs=genome,
            annots=None,
            species=species[name],
        )

    return genomes, align_db


def test_building_alignment(genomedbs_aligndb, namer):
    genomes, align_db = genomedbs_aligndb
    got = list(
        elt_align.get_alignment(
            align_db,
            genomes,
            ref_species="mouse",
            seqid="s2",
            namer=namer,
        ),
    )[0]
    orig = small_seqs()[1:5]
    assert got.to_dict() == orig.to_dict()


@pytest.mark.parametrize(
    "kwargs",
    (dict(ref_species="dodo", seqid="s2"),),
)
def test_building_alignment_invalid_details(genomedbs_aligndb, kwargs):
    genomes, align_db = genomedbs_aligndb
    with pytest.raises(ValueError):
        list(elt_align.get_alignment(align_db, genomes, **kwargs))


def make_sample(two_aligns=False):
    aln = small_seqs()
    species = aln.info.species
    # make annotation db's
    annot_dbs = {}
    for name in aln.names:
        feature_db = aln.annotation_db.subset(seqid=name)
        annot_dbs[name] = feature_db

    # we will reverse complement the s2 genome compared to the original
    # this means our coordinates for alignment records from that genome
    # also need to be rc'ed
    genomes = {}
    for seq in aln.seqs:
        name = seq.name
        seq = seq.data.degap()
        if seq.name == "s2":
            seq = seq.rc()
            s2_genome = str(seq)
        genome = elt_genome.SeqsDataHdf5(
            source=f"{name}",
            mode="w",
            in_memory=True,
            species=species[seq.name],
        )
        genome.add_records(records=[(name, str(seq))])
        genomes[species[name]] = elt_genome.Genome(
            seqs=genome,
            annots=annot_dbs[name],
            species=species[name],
        )

    # define two alignment blocks that incorporate features
    align_records = _update_records(s2_genome, aln, "0", 1, 12)
    if two_aligns:
        align_records += _update_records(s2_genome, aln, "1", 22, 30)
    align_db = elt_align.AlignDb(source=":memory:")
    align_db.add_records(records=align_records)

    return genomes, align_db


def _update_records(s2_genome, aln, block_id, start, end):
    # start, stop are the coordinates used to slice the alignment
    align_records = make_records(start, end, block_id)
    # in the alignment, s2 is in reverse complement relative to its genome
    # In order to be sure what "genome" coordinates are for s2, we first slice
    # the alignment
    aln = aln[start:end]
    # then get the ungapped sequence
    s2 = aln.get_seq("s2")
    # and reverse complement it ...
    selected = s2.rc()
    # so we can get the genome coordinates for this segment on the s2 genome
    start = s2_genome.find(str(selected))
    end = start + len(selected)
    for record in align_records:
        if record.seqid == "s2":
            record.start = start
            record.stop = end
            record.strand = "-"
            break
    return align_records


@pytest.mark.parametrize(
    "start_end",
    (
        (None, None),
        (None, 11),
        (3, None),
        (3, 13),
    ),
)
@pytest.mark.parametrize(
    "species_coord",
    (
        ("human", "s1"),
        ("dog", "s3"),
    ),
)
def test_select_alignment_plus_strand(species_coord, start_end, namer):
    species, seqid = species_coord
    start, end = start_end
    aln = small_seqs()
    # the sample alignment db has an alignment block that
    # starts at 1 and ends at 12. The following slice is to
    # get the expected answer
    expect = aln[max(1, start or 1) : min(end or 12, 12)]
    # one sequence is stored in reverse complement
    genomes, align_db = make_sample()
    got = list(
        elt_align.get_alignment(
            align_db=align_db,
            genomes=genomes,
            ref_species=species,
            seqid=seqid,
            ref_start=start,
            ref_end=end,
            namer=namer,
        ),
    )
    assert len(got) == 1
    assert got[0].to_dict() == expect.to_dict()


@pytest.mark.parametrize(
    "start_end",
    (
        (None, None),
        (None, 5),
        (2, None),
        (2, 7),
    ),
)
def test_select_alignment_minus_strand(start_end, namer):
    species, seqid = "mouse", "s2"
    start, end = start_end
    aln = small_seqs()
    ft = aln.add_feature(
        biotype="custom",
        name="selected",
        seqid="s2",
        spans=[(max(1, start or 0), min(end or 12, 12))],
    )
    expect = aln[ft.map.start : min(ft.map.end, 12)]
    # mouse sequence is on minus strand, so need to adjust
    # coordinates for query
    s2 = aln.get_seq("s2")
    s2_ft = list(s2.get_features(name="selected"))[0]
    if not any([start is None, end is None]):
        start = len(s2) - s2_ft.map.end
        end = len(s2) - s2_ft.map.start
    elif start == None != end:  # noqa E711
        start = len(s2) - s2_ft.map.end
        end = None
    elif start != None == end:  # noqa E711
        end = len(s2) - s2_ft.map.start
        start = None

    # mouse sequence is on minus strand, so need to adjust
    # coordinates for query

    genomes, align_db = make_sample(two_aligns=False)
    got = list(
        elt_align.get_alignment(
            align_db=align_db,
            genomes=genomes,
            ref_species=species,
            seqid=seqid,
            ref_start=start,
            ref_end=end,
            namer=namer,
        ),
    )
    # drop the strand info
    assert len(got) == 1, f"{s2_ft=}"
    assert got[0].to_dict() == expect.to_dict()


@pytest.mark.parametrize(
    "coord",
    (
        ("human", "s1", None, 11),  # finish within
        ("human", "s1", 3, None),  # start within
        ("human", "s1", 3, 9),  # within
        ("human", "s1", 3, 13),  # extends past
    ),
)
def test_get_alignment_features(coord):
    kwargs = dict(
        zip(("ref_species", "seqid", "ref_start", "ref_end"), coord, strict=False),
    )
    genomes, align_db = make_sample(two_aligns=False)
    got = list(elt_align.get_alignment(align_db=align_db, genomes=genomes, **kwargs))[0]
    assert len(got.annotation_db) == 1


@pytest.mark.parametrize(
    "coord",
    (
        ("human", "s1", None, 11),  # finish within
        ("human", "s1", 3, None),  # start within
        ("human", "s1", 3, 9),  # within
        ("human", "s1", 3, 13),  # extends past
    ),
)
def test_get_alignment_masked_features(coord):
    kwargs = dict(
        zip(("ref_species", "seqid", "ref_start", "ref_end"), coord, strict=False),
    )
    kwargs["mask_features"] = ["gene"]
    genomes, align_db = make_sample(two_aligns=False)
    got = list(elt_align.get_alignment(align_db=align_db, genomes=genomes, **kwargs))[0]
    assert len(got.annotation_db) == 1


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
    kwargs = dict(zip(("species", "seqid", "start", "stop"), coord, strict=False))
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
    kwargs = dict(zip(("species", "seqid"), coord, strict=False))
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
    kwargs = dict(zip(("species", "seqid"), coord, strict=False))
    # no hits at all
    _, align_db = make_sample()
    got = list(align_db.get_records_matching(**kwargs))
    assert not len(got)


def test_get_species():
    _, align_db = make_sample()
    assert set(align_db.get_species_names()) == {"dog", "human", "mouse"}


def test_write_alignments(tmp_path):
    genomes, align_db = make_sample(two_aligns=True)
    locations = elt_genome.get_gene_segments(
        annot_db=genomes["human"].annotation_db,
        species="human",
        stableids=["not-on-s2"],
    )
    app = elt_align.construct_alignment(
        align_db=align_db,
        genomes=genomes,
    )
    aln = app(locations[0])  # pylint: disable=not-callable
    assert len(aln[0]) == 3


@pytest.mark.parametrize(
    "gaps",
    [numpy.array([[0, 24], [2, 3]], dtype=int), numpy.array([], dtype=int)],
)
def test_gapstore_add_retrieve(gaps):
    gap_store = elt_align.GapStore(
        source="stuff",
        in_memory=True,
        mode="w",
        align_name="demo",
    )
    gap_store.add_record(index=20, gaps=gaps)
    got = gap_store.get_record(index=20)
    assert got is not gaps
    assert (got == gaps).all()


def test_gapstore_add_duplicate():
    a = numpy.array([[0, 24], [2, 3]], dtype=int)
    gap_store = elt_align.GapStore(
        source="stuff",
        in_memory=True,
        mode="w",
        align_name="demo",
    )
    gap_store.add_record(index=20, gaps=a)
    # adding it again has now effect
    gap_store.add_record(index=20, gaps=a)
    got = gap_store.get_record(index=20)
    assert got is not a
    assert (got == a).all()


def test_gapstore_add_invalid_duplicate():
    a = numpy.array([[0, 24], [2, 3]], dtype=int)
    gap_store = elt_align.GapStore(
        source="stuff",
        in_memory=True,
        mode="w",
        align_name="demo",
    )
    gap_store.add_record(index=20, gaps=a)
    with pytest.raises(ValueError):
        gap_store.add_record(index=20, gaps=a[:1])


@pytest.fixture
def small_db(small_records):
    import copy

    db = elt_align.AlignDb(source=":memory:")
    db.add_records(records=copy.deepcopy(small_records))
    return db


@pytest.mark.parametrize("col", tuple(elt_align.AlignDb._index_columns["align"]))
def test_indexing(small_db, col):
    table_name = "align"
    expect = ("index", f"{col}_index", table_name)
    small_db.make_indexes()
    sql_template = (
        f"SELECT * FROM sqlite_master WHERE type = 'index' AND "  # nosec B608
        f"tbl_name = {table_name!r} and name = '{col}_index'"  # nosec B608
    )

    result = small_db._execute_sql(sql_template).fetchone()
    got = tuple(result)[:3]
    assert got == expect
