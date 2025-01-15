# parser for MAF, defined at
# https://genome.ucsc.edu/FAQ/FAQformat.html#format5

import re
import typing

import numpy
from cogent3 import make_seq, open_
from cogent3.app.composable import LOADER, define_app
from cogent3.app.typing import IdentifierType

from ensembl_tui import _align as eti_align
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


def seq2gaps(record: dict) -> eti_align.AlignRecord:
    seq = make_seq(record.pop("seq"))
    indel_map, _ = seq.parse_out_gaps()
    if indel_map.num_gaps:
        record["gap_spans"] = numpy.array(
            [indel_map.gap_pos, indel_map.get_gap_lengths()],
            dtype=numpy.int32,
        ).T
    else:
        record["gap_spans"] = numpy.array([], dtype=numpy.int32)
    return eti_align.AlignRecord(**record)


@define_app(app_type=LOADER)
class load_align_records:
    def __init__(self, species: set[str] | None = None):
        self.species = species or {}

    def main(self, path: IdentifierType) -> list[eti_align.AlignRecord]:
        records = []
        for block_id, align in parse(path):
            converted = []
            for maf_name, seq in align.items():
                if self.species and maf_name.species not in self.species:
                    continue
                record = maf_name.to_dict()
                record["block_id"] = block_id
                record["source"] = path.name
                record["seq"] = seq
                converted.append(seq2gaps(record))
            records.extend(converted)
        return records
