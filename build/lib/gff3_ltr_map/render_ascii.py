"""Plain-text postcard renderer for RepeatRegion objects."""

from __future__ import annotations

from typing import List

from .model import RepeatRegion


def _scale(pos: int, start: int, end: int, width: int) -> int:
    span = max(1, end - start)
    rel = (pos - start) / span
    return max(0, min(width - 1, round(rel * (width - 1))))


def _fill(line: List[str], start: int, end: int, char: str, elem: RepeatRegion, width: int) -> None:
    ia = _scale(start, elem.start, elem.end, width)
    ib = _scale(end, elem.start, elem.end, width)
    if ia > ib:
        ia, ib = ib, ia
    for idx in range(ia, ib + 1):
        line[idx] = char


def ascii_map(elem: RepeatRegion, width: int = 100, ruler: bool = False) -> str:
    width = max(10, width)
    line = [" "] * width
    line[0] = "|"
    line[-1] = "|"

    if elem.ltr5:
        _fill(line, elem.ltr5.start, elem.ltr5.end, "=", elem, width)
    if elem.ltr3:
        _fill(line, elem.ltr3.start, elem.ltr3.end, "=", elem, width)
    if elem.ltr5 and elem.ltr3 and elem.ltr3.start - elem.ltr5.end > 1:
        _fill(line, elem.ltr5.end + 1, elem.ltr3.start - 1, "-", elem, width)

    if elem.tsd5 is not None:
        line[_scale(elem.tsd5, elem.start, elem.end, width)] = "T"
    if elem.tsd3 is not None:
        line[_scale(elem.tsd3, elem.start, elem.end, width)] = "D"

    meta = (
        f"> {elem.scaffold}:{elem.id}  fam:{elem.superfamily or 'NA'}  strand:{elem.strand}  "
        f"span:{elem.start}-{elem.end}  ltr_id:{elem.ltr_identity if elem.ltr_identity is not None else 'NA'}  "
        f"motif:{elem.motif or 'NA'}  tsd:{elem.tsd or 'NA'}"
    )

    bar = "".join(line)
    if not ruler:
        return f"{meta}\n{bar}\n"

    start_label = str(elem.start)
    end_label = str(elem.end)
    pad = max(0, width - len(start_label) - len(end_label))
    ruler_line = start_label + (" " * pad) + end_label
    return f"{meta}\n{bar}\n{ruler_line}\n"
