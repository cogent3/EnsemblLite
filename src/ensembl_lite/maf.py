# parser for MAF, defined at
# https://genome.ucsc.edu/FAQ/FAQformat.html#format5

import os

from cogent3 import open_

from ensembl_lite.name import MafName


def _get_alignment_block_indices(data: list[str]) -> list[tuple[int]]:
    blocks = []
    start = None
    for i, line in enumerate(data):
        if line.startswith("a"):
            if start is not None:
                blocks.append((start, i))
            start = i

    if start is None:
        return []

    blocks.append((start, i))
    return blocks


def _get_seqs(lines: list[str]) -> dict[MafName, str]:
    alignment = {}
    for line in lines:
        if not line.startswith("s") or "ancestral" in line[:100]:
            continue
        # after the s token we have src.seqid, start, size, strand, src_size, seq
        _, src_coord, start, size, strand, coord_length, seq = line.strip().split()
        species, coord = src_coord.split(".", maxsplit=1)
        start, size, coord_length = int(start), int(size), int(coord_length)
        n = MafName(
            species=species,
            seqid=coord,
            start=start,
            end=start + start,
            strand=strand,
            coord_length=coord_length,
        )
        alignment[n] = seq
    return alignment


def parse(path: os.PathLike) -> dict[MafName, str]:
    with open_(path) as infile:
        data = infile.readlines()

    blocks = _get_alignment_block_indices(data)
    for block_start, block_end in blocks:
        yield _get_seqs(data[block_start:block_end])
