import pytest

from ensembl_tui import _emf as elt_emf
from ensembl_tui import _name as elt_name


def test_load(DATA_DIR):
    path = DATA_DIR / "sample.emf"
    got = list(elt_emf.parse_emf(path))[0]
    expect = {
        elt_name.EmfName(
            "human",
            "4",
            "450000",
            "560000",
            "1",
            "(chr_length=201709)",
        ): "-TCGC",
        elt_name.EmfName(
            "mouse",
            "17",
            "780000",
            "790000",
            "-1",
            "(chr_length=201709)",
        ): "AT--G",
        elt_name.EmfName(
            "rat",
            "12",
            "879999",
            "889998",
            "1",
            "(chr_length=201709)",
        ): "AAA--",
    }
    assert got == expect


def test_unsupported_format(tmp_path, DATA_DIR):
    data = (DATA_DIR / "sample.emf").read_text().splitlines(keepends=True)
    data[0] = data[0].replace("compara", "resequencing")
    outpath = tmp_path / "sample.emf"
    outpath.write_text("".join(data))
    with pytest.raises(NotImplementedError):
        list(elt_emf.parse_emf(outpath))
