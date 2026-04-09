"""Batch summaries for collections of EDTA intact LTR reports."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from statistics import median
from typing import Iterable, List, Sequence

from .aggregates import format_top_counts
from .model import RepeatRegion

BATCH_SAMPLE_COLUMNS = [
    "sample",
    "input_path",
    "n_repeat_regions",
    "n_intact",
    "pct_intact",
    "n_pass",
    "n_warn",
    "n_fail",
    "n_scaffolds",
    "n_superfamilies",
    "median_ltr_identity",
    "median_length_bp",
    "median_internal_len",
    "pct_with_both_tsd",
    "top_motifs",
    "top_tsd",
    "notes",
]

BATCH_SUPERFAMILY_COLUMNS = [
    "sample",
    "superfamily",
    "n_elements",
    "pct_sample_elements",
    "median_ltr_identity",
    "median_length_bp",
    "median_internal_len",
    "top_motifs",
    "top_tsd",
]


@dataclass
class BatchSampleRow:
    sample: str
    input_path: str
    n_repeat_regions: int
    n_intact: int
    pct_intact: float
    n_pass: int
    n_warn: int
    n_fail: int
    n_scaffolds: int
    n_superfamilies: int
    median_ltr_identity: float | None
    median_length_bp: float | None
    median_internal_len: float | None
    pct_with_both_tsd: float
    top_motifs: str
    top_tsd: str
    notes: str


@dataclass
class BatchSuperfamilyRow:
    sample: str
    superfamily: str
    n_elements: int
    pct_sample_elements: float
    median_ltr_identity: float | None
    median_length_bp: float | None
    median_internal_len: float | None
    top_motifs: str
    top_tsd: str


def summarize_sample(
    sample: str,
    input_path: str,
    all_elements: Sequence[RepeatRegion],
    intact_elements: Sequence[RepeatRegion],
    *,
    top_k: int,
) -> BatchSampleRow:
    n_total = len(all_elements)
    n_intact = len(intact_elements)
    identities = [elem.ltr_identity for elem in intact_elements if elem.ltr_identity is not None]
    lengths = [elem.length_bp for elem in intact_elements]
    internals = [elem.internal_len for elem in intact_elements if elem.internal_len is not None]
    motifs_text, _ = format_top_counts([elem.motif or "" for elem in intact_elements], top_k, n_intact)
    tsd_text, _ = format_top_counts([elem.tsd or "" for elem in intact_elements], top_k, n_intact)
    n_pass = sum(1 for elem in all_elements if elem.qc_status == "PASS")
    n_warn = sum(1 for elem in all_elements if elem.qc_status == "WARN")
    n_fail = sum(1 for elem in all_elements if elem.qc_status == "FAIL")
    notes: List[str] = []
    if n_fail:
        notes.append(f"{n_fail} QC fail")
    if n_warn:
        notes.append(f"{n_warn} QC warn")
    return BatchSampleRow(
        sample=sample,
        input_path=input_path,
        n_repeat_regions=n_total,
        n_intact=n_intact,
        pct_intact=(n_intact / n_total * 100) if n_total else 0.0,
        n_pass=n_pass,
        n_warn=n_warn,
        n_fail=n_fail,
        n_scaffolds=len({elem.scaffold for elem in intact_elements}),
        n_superfamilies=len({elem.superfamily for elem in intact_elements if elem.superfamily}),
        median_ltr_identity=median(identities) if identities else None,
        median_length_bp=median(lengths) if lengths else None,
        median_internal_len=median(internals) if internals else None,
        pct_with_both_tsd=(sum(1 for elem in intact_elements if elem.has_both_tsd) / n_intact * 100)
        if n_intact
        else 0.0,
        top_motifs=motifs_text,
        top_tsd=tsd_text,
        notes="; ".join(notes),
    )


def summarize_sample_superfamilies(
    sample: str,
    intact_elements: Sequence[RepeatRegion],
    *,
    top_k: int,
) -> List[BatchSuperfamilyRow]:
    groups = {}
    for elem in intact_elements:
        key = elem.superfamily or "NA"
        groups.setdefault(key, []).append(elem)

    rows: List[BatchSuperfamilyRow] = []
    total = len(intact_elements)
    for superfamily, elems in sorted(groups.items()):
        identities = [elem.ltr_identity for elem in elems if elem.ltr_identity is not None]
        lengths = [elem.length_bp for elem in elems]
        internals = [elem.internal_len for elem in elems if elem.internal_len is not None]
        motifs_text, _ = format_top_counts([elem.motif or "" for elem in elems], top_k, len(elems))
        tsd_text, _ = format_top_counts([elem.tsd or "" for elem in elems], top_k, len(elems))
        rows.append(
            BatchSuperfamilyRow(
                sample=sample,
                superfamily=superfamily,
                n_elements=len(elems),
                pct_sample_elements=(len(elems) / total * 100) if total else 0.0,
                median_ltr_identity=median(identities) if identities else None,
                median_length_bp=median(lengths) if lengths else None,
                median_internal_len=median(internals) if internals else None,
                top_motifs=motifs_text,
                top_tsd=tsd_text,
            )
        )
    return rows


def write_batch_sample_tsv(rows: Iterable[BatchSampleRow], handle) -> None:
    handle.write("\t".join(BATCH_SAMPLE_COLUMNS) + "\n")
    for row in rows:
        payload = asdict(row)
        handle.write("\t".join(_format_value(payload[column]) for column in BATCH_SAMPLE_COLUMNS) + "\n")


def write_batch_superfamily_tsv(rows: Iterable[BatchSuperfamilyRow], handle) -> None:
    handle.write("\t".join(BATCH_SUPERFAMILY_COLUMNS) + "\n")
    for row in rows:
        payload = asdict(row)
        handle.write(
            "\t".join(_format_value(payload[column]) for column in BATCH_SUPERFAMILY_COLUMNS) + "\n"
        )


def _format_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)
