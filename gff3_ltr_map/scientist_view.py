"""Scientist-facing TSV rows and CLI printers."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from statistics import median
from typing import Iterable, List, Sequence

from .model import RepeatRegion

SCIENTIST_REGION_COLUMNS = [
    "sample",
    "region",
    "window_bp",
    "intact_elements",
    "intact_per_mbp",
    "qc_warn",
    "median_ltr_identity",
    "identity_band",
    "median_length_bp",
    "coverage_pct",
    "pct_with_both_tsd",
    "dominant_superfamily",
    "dominant_motif",
    "note",
]

SCIENTIST_ELEMENT_COLUMNS = [
    "scaffold",
    "start",
    "end",
    "superfamily",
    "ltr_identity",
    "identity_band",
    "length_bp",
    "internal_len",
    "motif",
    "tsd",
    "qc_status",
    "qc_note",
    "edta_id",
]


@dataclass
class ScientistRegionRow:
    sample: str
    region: str
    window_bp: int | None
    intact_elements: int
    intact_per_mbp: float | None
    qc_warn: int
    median_ltr_identity: float | None
    identity_band: str
    median_length_bp: float | None
    coverage_pct: float | None
    pct_with_both_tsd: float
    dominant_superfamily: str
    dominant_motif: str
    note: str


@dataclass
class ScientistElementRow:
    scaffold: str
    start: int
    end: int
    superfamily: str
    ltr_identity: float | None
    identity_band: str
    length_bp: int
    internal_len: int | None
    motif: str
    tsd: str
    qc_status: str
    qc_note: str
    edta_id: str


def build_scientist_region_row(
    sample: str,
    region: str,
    all_elements: Sequence[RepeatRegion],
    intact_elements: Sequence[RepeatRegion],
    *,
    window_bp: int | None = None,
    coverage_pct: float | None = None,
    extra_note: str | None = None,
) -> ScientistRegionRow:
    identities = [elem.ltr_identity for elem in intact_elements if elem.ltr_identity is not None]
    lengths = [elem.length_bp for elem in intact_elements]
    superfamilies = Counter(elem.superfamily or "NA" for elem in intact_elements)
    motifs = Counter(elem.motif or "NA" for elem in intact_elements)
    dominant_superfamily = superfamilies.most_common(1)[0][0] if superfamilies else "NA"
    dominant_motif = motifs.most_common(1)[0][0] if motifs else "NA"
    warn_count = sum(1 for elem in all_elements if elem.qc_status == "WARN")
    note_parts: List[str] = []
    if extra_note:
        note_parts.append(extra_note)
    if warn_count:
        note_parts.append(f"{warn_count} warned elements")
    if not intact_elements:
        note_parts.append("no intact elements")
    median_identity = median(identities) if identities else None
    return ScientistRegionRow(
        sample=sample,
        region=region,
        window_bp=window_bp,
        intact_elements=len(intact_elements),
        intact_per_mbp=(len(intact_elements) / (window_bp / 1_000_000.0)) if window_bp else None,
        qc_warn=warn_count,
        median_ltr_identity=median_identity,
        identity_band=identity_band(median_identity),
        median_length_bp=median(lengths) if lengths else None,
        coverage_pct=coverage_pct,
        pct_with_both_tsd=(
            sum(1 for elem in intact_elements if elem.has_both_tsd) / len(intact_elements) * 100
            if intact_elements
            else 0.0
        ),
        dominant_superfamily=dominant_superfamily,
        dominant_motif=dominant_motif,
        note="; ".join(note_parts),
    )


def build_scientist_element_rows(elements: Sequence[RepeatRegion]) -> List[ScientistElementRow]:
    rows: List[ScientistElementRow] = []
    for elem in elements:
        rows.append(
            ScientistElementRow(
                scaffold=elem.scaffold,
                start=elem.start,
                end=elem.end,
                superfamily=elem.superfamily or "NA",
                ltr_identity=elem.ltr_identity,
                identity_band=identity_band(elem.ltr_identity),
                length_bp=elem.length_bp,
                internal_len=elem.internal_len,
                motif=elem.motif or "NA",
                tsd=elem.tsd or "NA",
                qc_status=elem.qc_status,
                qc_note=_qc_note(elem),
                edta_id=elem.id,
            )
        )
    return rows


def write_scientist_region_tsv(rows: Iterable[ScientistRegionRow], handle) -> None:
    handle.write("\t".join(SCIENTIST_REGION_COLUMNS) + "\n")
    for row in rows:
        payload = asdict(row)
        handle.write("\t".join(_format_value(payload[column]) for column in SCIENTIST_REGION_COLUMNS) + "\n")


def write_scientist_element_tsv(rows: Iterable[ScientistElementRow], handle) -> None:
    handle.write("\t".join(SCIENTIST_ELEMENT_COLUMNS) + "\n")
    for row in rows:
        payload = asdict(row)
        handle.write("\t".join(_format_value(payload[column]) for column in SCIENTIST_ELEMENT_COLUMNS) + "\n")


def print_scientist_region_table(rows: Sequence[ScientistRegionRow], *, max_rows: int = 200) -> None:
    ordered_rows = sorted(rows, key=lambda row: (-row.intact_elements, row.sample))
    _print_table(
        title="[gff3-ltr-map] Simplified region table:",
        columns=[
            ("sample", lambda row: row.sample),
            ("region", lambda row: row.region),
            ("window_bp", lambda row: _fmt_int(row.window_bp)),
            ("intact", lambda row: str(row.intact_elements)),
            ("intact_per_mbp", lambda row: _fmt_float(row.intact_per_mbp)),
            ("qc_warn", lambda row: str(row.qc_warn)),
            ("median_id", lambda row: _fmt_float(row.median_ltr_identity)),
            ("band", lambda row: row.identity_band),
            ("median_len", lambda row: _fmt_int(row.median_length_bp)),
            ("coverage_pct", lambda row: _fmt_float(row.coverage_pct)),
            ("dominant_sf", lambda row: row.dominant_superfamily),
            ("motif", lambda row: row.dominant_motif),
            ("note", lambda row: row.note),
        ],
        rows=ordered_rows,
        max_rows=max_rows,
    )


def print_scientist_element_table(rows: Sequence[ScientistElementRow], *, max_rows: int = 25) -> None:
    _print_table(
        title="[gff3-ltr-map] Detailed element table:",
        columns=[
            ("scaffold", lambda row: row.scaffold),
            ("start", lambda row: str(row.start)),
            ("end", lambda row: str(row.end)),
            ("superfamily", lambda row: row.superfamily),
            ("ltr_identity", lambda row: _fmt_float(row.ltr_identity)),
            ("identity_band", lambda row: row.identity_band),
            ("length_bp", lambda row: _fmt_int(row.length_bp)),
            ("motif", lambda row: row.motif),
            ("tsd", lambda row: row.tsd),
            ("qc", lambda row: row.qc_status),
            ("edta_id", lambda row: row.edta_id),
        ],
        rows=rows,
        max_rows=max_rows,
    )


def identity_band(identity: float | None) -> str:
    if identity is None:
        return "NA"
    if identity >= 0.98:
        return "very_high"
    if identity >= 0.95:
        return "high"
    if identity >= 0.90:
        return "moderate"
    return "low"


def _qc_note(elem: RepeatRegion) -> str:
    if elem.validation_warnings:
        return "; ".join(elem.validation_warnings)
    if elem.validation_errors:
        return "; ".join(elem.validation_errors)
    return ""


def _print_table(title: str, columns, rows, *, max_rows: int) -> None:
    print(title)
    print("\t".join(name for name, _ in columns))
    shown = rows[:max_rows]
    for row in shown:
        print("\t".join(getter(row) for _, getter in columns))
    if len(rows) > max_rows:
        print(f"[gff3-ltr-map] Showing first {max_rows} rows of {len(rows)}")


def _format_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _fmt_float(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.3f}"


def _fmt_int(value: float | int | None) -> str:
    if value is None:
        return "NA"
    return f"{int(round(value)):,}"
