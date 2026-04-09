"""Summary and validation TSV writers."""

from __future__ import annotations

from typing import Iterable, TextIO

from .model import RepeatRegion

SUMMARY_HEADER = [
    "element_id",
    "scaffold",
    "start",
    "end",
    "strand",
    "source",
    "name",
    "classification",
    "superfamily",
    "method",
    "retrotransposon_type",
    "ltr_identity",
    "motif",
    "tsd",
    "tsd_len",
    "length_bp",
    "ltr5_len",
    "ltr3_len",
    "internal_len",
    "has_both_tsd",
    "is_intact",
    "qc_status",
    "validation_warnings",
    "validation_errors",
    "n_children",
    "child_signature",
]

VALIDATION_HEADER = [
    "element_id",
    "scaffold",
    "classification",
    "superfamily",
    "is_intact",
    "qc_status",
    "warning_count",
    "error_count",
    "validation_warnings",
    "validation_errors",
    "child_signature",
]


def write_summary(rows: Iterable[RepeatRegion], handle: TextIO) -> None:
    handle.write("\t".join(SUMMARY_HEADER) + "\n")
    for row in rows:
        handle.write(
            "\t".join(
                [
                    row.id,
                    row.scaffold,
                    str(row.start),
                    str(row.end),
                    row.strand,
                    row.source,
                    row.name or "",
                    row.classification or "",
                    row.superfamily or "",
                    row.method or "",
                    row.retrotransposon_type or "",
                    "" if row.ltr_identity is None else str(row.ltr_identity),
                    row.motif or "",
                    row.tsd or "",
                    "" if row.tsd_len is None else str(row.tsd_len),
                    str(row.length_bp),
                    "" if row.ltr5_len is None else str(row.ltr5_len),
                    "" if row.ltr3_len is None else str(row.ltr3_len),
                    "" if row.internal_len is None else str(row.internal_len),
                    "true" if row.has_both_tsd else "false",
                    "true" if row.is_intact else "false",
                    row.qc_status,
                    row.validation_warning_text,
                    row.validation_error_text,
                    str(row.n_children_warn),
                    row.child_signature,
                ]
            )
            + "\n"
        )


def write_validation_report(rows: Iterable[RepeatRegion], handle: TextIO) -> None:
    handle.write("\t".join(VALIDATION_HEADER) + "\n")
    for row in rows:
        handle.write(
            "\t".join(
                [
                    row.id,
                    row.scaffold,
                    row.classification or "",
                    row.superfamily or "",
                    "true" if row.is_intact else "false",
                    row.qc_status,
                    str(len(row.validation_warnings)),
                    str(len(row.validation_errors)),
                    row.validation_warning_text,
                    row.validation_error_text,
                    row.child_signature,
                ]
            )
            + "\n"
        )
