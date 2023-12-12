import pytest

from cogent3 import make_seq

from ensembl_lite.convert import gap_coords_to_seq, seq_to_gap_coords


@pytest.mark.parametrize(
    "seq", ("----", "---AGC--TGC--", "AGC--TGC--", "---AGC--TGC", "AGCTGC")
)
def test_roundtrip_gapped_seqs(seq):
    seq = make_seq(seq, moltype="dna")
    c, ug = seq_to_gap_coords(seq)
    aligned = gap_coords_to_seq(c, ug)
    assert str(aligned) == str(seq)


def test_roundtrip_sliced_gapped_seqs():
    aligned = make_seq("---AGC--TGC", moltype="dna")[1:8]
    c, s = seq_to_gap_coords(aligned)
    got = gap_coords_to_seq(c, s)
    assert str(got) == str(aligned)
