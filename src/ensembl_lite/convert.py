import numpy

from cogent3 import Sequence
from cogent3.core.alignment import Aligned
from cogent3.core.location import LostSpan, Map, Span


O = tuple[numpy.ndarray, Sequence]


def seq_to_gap_coords(seq: Sequence) -> O:
    """returns coordinates of sequence gaps"""
    m, s = seq.parse_out_gaps()
    # Assuming the maximum integer is < 2^31
    return numpy.array(m.get_gap_coordinates(), dtype=numpy.int32), s


def gap_coords_to_seq(coords: numpy.ndarray, ungapped: Sequence) -> Aligned:
    """returns Aligned instance

    Parameters
    ----------
    coords
        2D array with first column being gap insertion point and
        second column being gap length
    ungapped
        the ungapped sequence instance the coordinates correspond to
    """
    segment_start = None
    spans = []
    insert = 0
    for insert, length in coords:
        gap = LostSpan(length)
        if segment_start is None:
            # this is the first gap
            if insert:
                # gap is within seq, so we include segment span first
                spans.append(Span(start=0, end=insert))
                # followed by the lost span

            # because alternate (insert == 0) means
            # gap is before seq, we just append gap for both cases
            spans.append(gap)
            # next segment_start is current insert point
            segment_start = insert
            continue

        spans.extend((Span(start=segment_start, end=insert), gap))
        segment_start = insert

    segment_total = sum(len(s) for s in spans if isinstance(s, Span))
    if len(ungapped) - segment_total:
        # last gap was also internal, so we add span for remainder of
        # sequence
        spans.append(Span(start=insert, end=len(ungapped)))

    m = Map(spans=spans, parent_length=sum(len(s) for s in spans))
    return Aligned(m, ungapped)
