import contextlib
import functools
import inspect
import os
import pathlib
import re
import shutil
import subprocess
import sys
import typing
import uuid
from collections.abc import Callable
from hashlib import md5
from tempfile import mkdtemp
from typing import IO, Union

import blosc2
import hdf5plugin
import numba
import numpy
from cogent3.app.composable import define_app
from cogent3.util.parallel import as_completed

PathType = Union[str, pathlib.Path, os.PathLike]

_HDF5_BLOSC2_KWARGS = hdf5plugin.Blosc2(
    cname="blosclz",
    clevel=9,
    filters=hdf5plugin.Blosc2.BITSHUFFLE,
)


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
def checksum(data: bytes, size: int):  # pragma: no cover
    """computes BSD style checksum"""
    # equivalent to command line BSD sum
    nb = numpy.ceil(size / 1024)
    cksum = 0
    for c in data:
        cksum = (cksum >> 1) + ((cksum & 1) << 15)
        cksum += c
        cksum &= 0xFFFF
    return cksum, int(nb)


def _get_resource_dir() -> PathType:
    """returns path to resource directory"""
    if "ENSEMBLDBRC" in os.environ:
        path = os.environ["ENSEMBLDBRC"]
    else:
        from ensembl_tui import data

        path = pathlib.Path(data.__file__).parent

    path = pathlib.Path(path).expanduser().absolute()
    if not path.exists():
        raise ValueError(f"ENSEMBLDBRC directory {str(path)!r} does not exist")

    return pathlib.Path(path)


def get_resource_path(resource: PathType) -> PathType:
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


def load_ensembl_checksum(path: PathType) -> dict:
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


def load_ensembl_md5sum(path: PathType) -> dict:
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

    def __init__(self, path: PathType, tmpdir=None, mode="wb", encoding=None):
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


def dont_checksum(path: PathType) -> bool:
    return _dont_checksum.search(str(path)) is not None


@functools.singledispatch
def is_signature(path: PathType) -> bool:
    return _sig_file.search(path.name) is not None


@is_signature.register
def _(path: str) -> bool:
    return _sig_file.search(path) is not None


@functools.singledispatch
def get_sig_calc_func(sig_path) -> Callable:
    """returns signature calculating function based on Ensembl path name"""
    raise NotImplementedError(f"{type(sig_path)} not supported")


@get_sig_calc_func.register
def _(sig_path: str) -> Callable:
    return _sig_calc_funcs[sig_path]


def get_signature_data(path: PathType) -> Callable:
    return _sig_load_funcs[path.name](path)


def rich_display(c3t, title_justify="left"):
    """converts a cogent3 Table to a Rich Table and displays it"""
    from rich.console import Console
    from rich.table import Table

    cols = c3t.columns
    columns = []
    for c in c3t.header:
        if tmplt := c3t._column_templates.get(c, None):
            col = [tmplt(v) for v in cols[c]]
        else:
            col = cols[c]
        columns.append(col)

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

    for row in zip(*columns, strict=False):
        rich_table.add_row(*row)

    console = Console()
    console.print(rich_table)


_seps = re.compile(r"[-._\s]")


def _name_parts(path: str) -> list[str]:
    return _seps.split(pathlib.Path(path).name.lower())


def _simple_check(align_parts: str, tree_parts: str) -> int:
    """evaluates whether the start of the two paths match"""
    matches = 0
    for a, b in zip(align_parts, tree_parts, strict=False):
        if a != b:
            break
        matches += 1

    return matches


def trees_for_aligns(aligns, trees) -> dict[str, str]:
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

_biotypes = re.compile(r"(gene|transcript|exon|mRNA|rRNA|protein):")


def sanitise_stableid(stableid: str) -> str:
    """remove <biotype>:E.. from Ensembl stable ID

    Notes
    -----
    The GFF3 files from Ensembl store identifiers as <biotype>:<identifier>,
    this function removes redundant biotype component.
    """
    return _biotypes.sub("", stableid)


@contextlib.contextmanager
def fake_wake(*args, **kwargs):
    yield


class SerialisableMixin:
    """mixin class, adds a self._init_vals dict attribute which
    contains the keyword/arg mapping of arguments provided to the
    constructor"""

    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        init_sig = inspect.signature(cls.__init__)
        bargs = init_sig.bind_partial(cls, *args, **kwargs)
        bargs.apply_defaults()
        init_vals = bargs.arguments
        init_vals.pop("self", None)
        obj._init_vals = init_vals
        return obj


_quotes = re.compile(r"^[\'\"]|[\'\"]$")


def strip_quotes(text: str):
    return _quotes.sub("", text)


def get_iterable_tasks(
    *,
    func: typing.Callable,
    series: typing.Sequence,
    max_workers: int | None,
    **kwargs,
) -> typing.Iterator:
    if max_workers == 1:
        return map(func, series)
    return as_completed(func, series, max_workers=max_workers, **kwargs)


# From http://mart.ensembl.org/info/genome/stable_ids/prefixes.html
# The Ensembl stable id structure is
# [species prefix][feature type prefix][a unique eleven digit number]
# feature type prefixes are
# E exon
# FM Ensembl protein family
# G gene
# GT gene tree
# P protein
# R regulatory feature
# T transcript
_feature_type_1 = {"E", "G", "P", "R", "T"}
_feature_type_2 = {"FM", "GT"}


def get_stableid_prefix(stableid: str) -> str:
    """returns the prefix component of a stableid"""
    if len(stableid) < 15:
        raise ValueError(f"{stableid!r} too short")

    if stableid[-13:-11] in _feature_type_2:
        return stableid[:-13]
    if stableid[-12] not in _feature_type_1:
        raise ValueError(f"{stableid!r} has unknown feature type {stableid[-13]!r}")
    return stableid[:-12]
