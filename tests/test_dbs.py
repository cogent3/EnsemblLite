import numpy
import pytest

from ensembl_tui import _storage_mixin as eti_storage


@pytest.mark.parametrize(
    "data",
    [numpy.array([], dtype=numpy.int32), numpy.array([0, 3], dtype=numpy.uint8)],
)
def test_array_blob_roundtrip(data):
    blob = eti_storage.array_to_blob(data)
    assert isinstance(blob, bytes)
    inflated = eti_storage.blob_to_array(blob)
    assert numpy.array_equal(inflated, data)
    assert inflated.dtype is data.dtype


@pytest.mark.parametrize(
    "data",
    [
        numpy.array([0, 3], dtype=numpy.uint8),
        eti_storage.array_to_blob(numpy.array([0, 3], dtype=numpy.uint8)),
    ],
)
def test_blob_array(data):
    # handles array or bytes as input
    inflated = eti_storage.blob_to_array(data)
    assert numpy.array_equal(inflated, numpy.array([0, 3], dtype=numpy.uint8))
