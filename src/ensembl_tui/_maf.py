# parser for MAF, defined at
# https://genome.ucsc.edu/FAQ/FAQformat.html#format5

import re
import typing

from cogent3 import open_

from ensembl_tui import _name as eti_name
from ensembl_tui import _util as eti_util

_id_pattern = re.compile(r"(?<=id[:])\s*\d+")


def _get_alignment_block_indices(data: list[str]) -> list[tuple[int, int]]:
    blocks = []
    start = None
    for i, line in enumerate(data):
        if _id_pattern.search(line):
            if start is not None:
                blocks.append((start, i))
            start = i

    if start is None:
        return []

    blocks.append((start, i))
    return blocks


def process_id_line(line: str) -> int:
    if match := _id_pattern.search(line):
        return int(match.group())

    msg = f"{line=} is not a tree id line"
    raise ValueError(msg)


def process_maf_line(line: str) -> tuple[eti_name.MafName, str]:
    # after the s token we have src.seqid, start, size, strand, src_size, seq
    _, src_coord, start, size, strand, coord_length, seq = line.strip().split()
    species, coord = src_coord.split(".", maxsplit=1)
    start, size, coord_length = int(start), int(size), int(coord_length)
    if strand == "-":
        start = coord_length - (start + size)

    stop = start + size
    n = eti_name.MafName(
        species=species,
        seqid=coord,
        start=start,
        stop=stop,
        strand=strand,
        coord_length=coord_length,
    )
    return n, seq


def _get_seqs(lines: list[str]) -> dict[eti_name.MafName, str]:
    alignment = {}
    for line in lines:
        if not line.startswith("s") or "ancestral" in line[:100]:
            continue
        n, seq = process_maf_line(line)
        alignment[n] = seq
    return alignment


def parse(
    path: eti_util.PathType,
) -> typing.Iterable[tuple[int, dict[eti_name.MafName, str]]]:
    with open_(path) as infile:
        data = infile.readlines()

    blocks = _get_alignment_block_indices(data)
    for block_start, block_end in blocks:
        block_id = process_id_line(data[block_start])
        yield block_id, _get_seqs(data[block_start + 1 : block_end])
