from __future__ import annotations

import functools
import os
import pathlib
import re
import shutil
import subprocess
import sys
import uuid

from hashlib import md5
from tempfile import mkdtemp
from typing import IO, Callable, Union

import blosc2
import numba
import numpy

from cogent3.app.composable import define_app


def md5sum(data: bytes, *args) -> str:
    """computes MD5SUM

    Notes
    -----
    *args is for signature compatability with checksum
    """
    return md5(data).hexdigest()


# based on https://www.reddit.com/r/learnpython/comments/9bpgjl/implementing_bsd_16bit_checksum/
# and https://www.gnu.org/software/coreutils/manual/html_node/sum-invocation.html#sum-invocation
@numba.jit(nopython=True)
def checksum(data: bytes, size: int):
    """computes BSD style checksum"""
    # equivalent to command line BSD sum
    nb = numpy.ceil(size / 1024)
    cksum = 0
    for c in data:
        cksum = (cksum >> 1) + ((cksum & 1) << 15)
        cksum += c
        cksum &= 0xFFFF
    return cksum, int(nb)


def _get_resource_dir() -> os.PathLike:
    """returns path to resource directory"""
    if "ENSEMBLDBRC" in os.environ:
        path = os.environ["ENSEMBLDBRC"]
    else:
        from ensembl_lite import data

        path = pathlib.Path(data.__file__).parent

    path = pathlib.Path(path).expanduser().absolute()
    if not path.exists():
        raise ValueError("ENSEMBLDBRC directory '%s' does not exist")

    return pathlib.Path(path)


def get_resource_path(resource: Union[str, os.PathLike]) -> os.PathLike:
    path = ENSEMBLDBRC / resource
    assert path.exists()
    return path


# the following is where essential files live, such as
# the species/common name map and sample download.cfg
ENSEMBLDBRC = _get_resource_dir()


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0

    Parameters
    ----------

    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""
    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        msg = err
        sys.stderr.writelines(f"FAILED: {cmnd}\n{msg}")
        sys.exit(proc.returncode)
    return out.decode("utf8") if out is not None else None


class CaseInsensitiveString(str):
    """A case-insensitive string class. Comparisons are also case-insensitive."""

    def __new__(cls, arg, h=None):
        n = str.__new__(cls, str(arg))
        n._lower = "".join(list(n)).lower()
        n._hash = hash(n._lower)
        return n

    def __eq__(self, other):
        return self._lower == "".join(list(other)).lower()

    def __hash__(self):
        # dict hashing done via lower case
        return self._hash

    def __str__(self):
        return "".join(list(self))


def load_ensembl_checksum(path: os.PathLike) -> dict:
    """loads the BSD checksums from Ensembl CHECKSUMS file"""
    result = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        s, b, p = line.split()
        result[p] = int(s), int(b)
    result.pop("README", None)
    return result


def load_ensembl_md5sum(path: os.PathLike) -> dict:
    """loads the md5 sum from Ensembl MD5SUM file"""
    result = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        s, p = line.split()
        result[p] = s
    result.pop("README", None)
    return result


class atomic_write:
    """performs atomic write operations, cleans up if fails"""

    def __init__(self, path: os.PathLike, tmpdir=None, mode="wb", encoding=None):
        """

        Parameters
        ----------
        path
            path to file
        tmpdir
            directory where temporary file will be created
        mode
            file writing mode
        encoding
            text encoding
        """
        path = pathlib.Path(path).expanduser()

        self._path = path
        self._mode = mode
        self._file = None
        self._encoding = encoding
        self._tmppath = self._make_tmppath(tmpdir)

        self.succeeded = None
        self._close_func = self._close_rename_standard

    def _make_tmppath(self, tmpdir):
        """returns path of temporary file

        Parameters
        ----------
        tmpdir: Path
            to directory

        Returns
        -------
        full path to a temporary file

        Notes
        -----
        Uses a random uuid as the file name, adds suffixes from path
        """
        suffixes = "".join(self._path.suffixes)
        parent = self._path.parent
        name = f"{uuid.uuid4()}{suffixes}"
        tmpdir = (
            pathlib.Path(mkdtemp(dir=parent))
            if tmpdir is None
            else pathlib.Path(tmpdir)
        )

        if not tmpdir.exists():
            raise FileNotFoundError(f"{tmpdir} directory does not exist")

        return tmpdir / name

    def _get_fileobj(self):
        """returns file to be written to"""
        if self._file is None:
            self._file = open(self._tmppath, self._mode)

        return self._file

    def __enter__(self) -> IO:
        return self._get_fileobj()

    def _close_rename_standard(self, src):
        dest = pathlib.Path(self._path)
        try:
            dest.unlink()
        except FileNotFoundError:
            pass
        finally:
            src.rename(dest)

        shutil.rmtree(src.parent)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file.close()
        if exc_type is None:
            self._close_func(self._tmppath)
            self.succeeded = True
        else:
            self.succeeded = False

        shutil.rmtree(self._tmppath.parent, ignore_errors=True)

    def write(self, text):
        """writes text to file"""
        fileobj = self._get_fileobj()
        fileobj.write(text)

    def close(self):
        """closes file"""
        self.__exit__(None, None, None)


_sig_load_funcs = dict(CHECKSUMS=load_ensembl_checksum, MD5SUM=load_ensembl_md5sum)
_sig_calc_funcs = dict(CHECKSUMS=checksum, MD5SUM=md5sum)
_dont_checksum = re.compile("(CHECKSUMS|MD5SUM|README)")
_sig_file = re.compile("(CHECKSUMS|MD5SUM)")


def dont_checksum(path: os.PathLike) -> bool:
    return _dont_checksum.search(str(path)) is not None


@functools.singledispatch
def is_signature(path: os.PathLike) -> bool:
    return _sig_file.search(path.name) is not None


@is_signature.register
def _(path: str) -> bool:
    return _sig_file.search(path) is not None


@functools.singledispatch
def get_sig_calc_func(sig_path: os.PathLike) -> Callable:
    return _sig_calc_funcs[sig_path.name]


@get_sig_calc_func.register
def _(sig_path: str) -> Callable:
    return _sig_calc_funcs[sig_path]


def get_signature_data(path: os.PathLike) -> Callable:
    return _sig_load_funcs[path.name](path)


def rich_display(c3t, title_justify="left"):
    """converts a cogent3 Table to a Rich Table and displays it"""
    from cogent3.format.table import formatted_array
    from rich.console import Console
    from rich.table import Table

    cols = c3t.columns
    columns = [formatted_array(cols[c], pad=False)[0] for c in c3t.header]
    rich_table = Table(
        title=c3t.title,
        highlight=True,
        title_justify=title_justify,
        title_style="bold blue",
    )
    for col in c3t.header:
        numeric_type = any(v in cols[col].dtype.name for v in ("int", "float"))
        j = "right" if numeric_type else "left"
        rich_table.add_column(col, justify=j, no_wrap=numeric_type)

    for row in zip(*columns):
        rich_table.add_row(*row)

    console = Console()
    console.print(rich_table)


_seps = re.compile(r"[-._\s]")


def _name_parts(path: str) -> list[str]:
    return _seps.split(pathlib.Path(path).name.lower())


def _simple_check(align_parts: str, tree_parts: str) -> int:
    """evaluates whether the start of the two paths match"""
    matches = 0
    for a, b in zip(align_parts, tree_parts):
        if a != b:
            break
        matches += 1

    return matches


def trees_for_aligns(aligns, trees) -> dict[str, str]:
    from cogent3.maths.distance_transform import jaccard

    aligns = {p: _name_parts(p) for p in aligns}
    trees = {p: _name_parts(p) for p in trees}
    result = {}
    for align, align_parts in aligns.items():
        dists = [
            (_simple_check(align_parts, tree_parts), tree)
            for tree, tree_parts in trees.items()
        ]
        v, p = max(dists)
        if v == 0:
            raise ValueError(f"no tree for {align}")

        result[align] = p

    return result


@define_app
def _str_to_bytes(data: str) -> bytes:
    """converts string to bytes"""
    return data.encode("utf8")


@define_app
def _bytes_to_str(data: bytes) -> str:
    """converts bytes into string"""
    return data.decode("utf8")


@define_app
def blosc_compress_it(data: bytes) -> bytes:
    return blosc2.compress(data, clevel=9, filter=blosc2.Filter.SHUFFLE)


@define_app
def blosc_decompress_it(data: bytes, as_bytearray=True) -> bytes:
    return bytes(blosc2.decompress(data, as_bytearray=as_bytearray))


elt_compress_it = _str_to_bytes() + blosc_compress_it()
elt_decompress_it = blosc_decompress_it() + _bytes_to_str()
