import pathlib
import typing

import numpy
from cogent3 import get_moltype, open_

# old style moltype
alphabet = get_moltype("dna").alphabets.degen_gapped


class converter:
    """Defines a linear mapping from provided characters to uint8.
    The resulting object is callable, taking a bytes object and returning a
    numpy array."""

    def __init__(self, dtype=numpy.uint8):
        self._tr = b"".maketrans(
            "".join(alphabet).encode("utf8"),
            bytes(bytearray(range(len(alphabet)))),
        )
        self.dtype = dtype

    def __call__(self, seq: bytes) -> numpy.ndarray:
        b = seq.translate(self._tr, delete=b" \n\r")
        return numpy.array(memoryview(b), dtype=self.dtype)


bytes_to_array = converter()


def quicka_parser(
    path: pathlib.Path,
    converter: typing.Callable[[bytes], bytes] = bytes_to_array,
):
    """generator returning sequence labels and sequences converted bytes from a fasta file

    Parameters
    ----------
    path
        location of the fasta file
    converter
        a callable that uses converts sequence characters into nominated bytes,
        deleting unwanted characters. Must handle newlines. Whatever type this
        callable returns will be the type of the sequence returned.

    Returns
    -------
    the sequence label as a string and the sequence as transformed by converter
    """
    with open_(path, mode="rb") as infile:
        data: bytes = infile.read()

    records = data.split(b">")
    for record in records:
        eol = record.find(b"\n")
        if eol == -1:
            continue
        label = record[:eol].strip().decode("utf8")
        seq = converter(record[eol + 1 :])
        yield label, seq
