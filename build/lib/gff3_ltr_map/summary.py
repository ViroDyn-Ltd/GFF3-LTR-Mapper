"""Summary TSV writer."""

from __future__ import annotations

from typing import Iterable, TextIO

from .model import RepeatRegion

HEADER = [
    "element_id",
    "scaffold",
    "start",
    "end",
    "strand",
    "superfamily",
    "ltr_identity",
    "motif",
    "tsd",
    "length_bp",
    "ltr5_len",
    "ltr3_len",
    "internal_len",
    "has_both_tsd",
    "n_children_warn",
]


def write_summary(rows: Iterable[RepeatRegion], handle: TextIO) -> None:
    handle.write("\t".join(HEADER) + "\n")
    for row in rows:
        handle.write(
            "\t".join(
                [
                    row.id,
                    row.scaffold,
                    str(row.start),
                    str(row.end),
                    row.strand,
                    row.superfamily or "",
                    "" if row.ltr_identity is None else str(row.ltr_identity),
                    row.motif or "",
                    row.tsd or "",
                    str(row.length_bp),
                    "" if row.ltr5_len is None else str(row.ltr5_len),
                    "" if row.ltr3_len is None else str(row.ltr3_len),
                    "" if row.internal_len is None else str(row.internal_len),
                    "true" if row.has_both_tsd else "false",
                    str(row.n_children_warn),
                ]
            )
            + "\n"
        )
