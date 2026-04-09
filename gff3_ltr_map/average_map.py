"""Utilities to build ASCII average maps from multiple RepeatRegion entries."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from statistics import mean, median
from typing import List, Optional, Sequence, Tuple

from .model import RepeatRegion


@dataclass
class AverageMapProfile:
    """Aggregated metrics describing a cohort of repeats."""

    label: str
    group_label: str
    n_elements: int
    mean_identity: Optional[float]
    median_identity: Optional[float]
    ltr5_len: float
    internal_len: float
    ltr3_len: float
    total_len: float
    identity_range: Optional[Tuple[Optional[float], Optional[float]]]
    low_n: bool
    note: str
    q25_len: Optional[float]
    q75_len: Optional[float]
    strand_summary: str
    pct_with_both_tsd: float
    median_tsd_len: Optional[float]

    @property
    def has_content(self) -> bool:
        return self.n_elements > 0 and self.total_len > 0


def build_average_profile(
    elements: Sequence[RepeatRegion],
    *,
    label: str,
    group_label: str,
    identity_range: Optional[Tuple[Optional[float], Optional[float]]] = None,
    min_n: int = 1,
) -> Optional[AverageMapProfile]:
    """Compute average lengths and metadata for a cohort."""

    if identity_range:
        min_id, max_id = identity_range
        filtered = [elem for elem in elements if _within_identity(elem, min_id, max_id)]
    else:
        filtered = list(elements)

    if not filtered:
        return None

    identities = [elem.ltr_identity for elem in filtered if elem.ltr_identity is not None]
    ltr5 = [elem.ltr5_len for elem in filtered if elem.ltr5_len is not None]
    ltr3 = [elem.ltr3_len for elem in filtered if elem.ltr3_len is not None]
    internal = [elem.internal_len for elem in filtered if elem.internal_len is not None]
    totals = [elem.length_bp for elem in filtered]
    tsd_lengths = [elem.tsd_len for elem in filtered if elem.tsd_len is not None]

    ltr5_len = median(ltr5) if ltr5 else 0.0
    ltr3_len = median(ltr3) if ltr3 else 0.0
    total_len = median(totals) if totals else ltr5_len + ltr3_len
    internal_len = median(internal) if internal else max(total_len - ltr5_len - ltr3_len, 0.0)

    if total_len <= 0:
        return None
    remainder = total_len - (ltr5_len + internal_len + ltr3_len)
    if abs(remainder) > 1e-6:
        internal_len = max(0.0, internal_len + remainder)

    q25 = _quantile(totals, 0.25)
    q75 = _quantile(totals, 0.75)
    low_n = len(filtered) < min_n
    note = ""
    if low_n:
        note = f"LOW N (n={len(filtered)})"

    strand_counts = Counter(elem.strand for elem in filtered if elem.strand in {"+", "-"})
    strand_summary = _format_strand_summary(strand_counts, len(filtered))

    return AverageMapProfile(
        label=label,
        group_label=group_label,
        n_elements=len(filtered),
        mean_identity=mean(identities) if identities else None,
        median_identity=median(identities) if identities else None,
        ltr5_len=ltr5_len,
        internal_len=internal_len,
        ltr3_len=ltr3_len,
        total_len=total_len,
        identity_range=identity_range,
        low_n=low_n,
        note=note,
        q25_len=q25,
        q75_len=q75,
        strand_summary=strand_summary,
        pct_with_both_tsd=(sum(1 for elem in filtered if elem.has_both_tsd) / len(filtered) * 100),
        median_tsd_len=median(tsd_lengths) if tsd_lengths else None,
    )


def average_ascii_map(
    profile: AverageMapProfile,
    *,
    width: int = 100,
    ruler: bool = True,
    show_quantiles: bool = False,
) -> str:
    """Render an aggregate profile as ASCII art."""

    width = max(20, width)
    line = [" "] * width
    line[0] = "|"
    line[-1] = "|"

    def fill(length: float, char: str) -> None:
        start_idx = _scale(fill.position, profile.total_len, width)
        fill.position += length
        end_idx = _scale(fill.position, profile.total_len, width)
        if end_idx < start_idx:
            start_idx, end_idx = end_idx, start_idx
        for idx in range(start_idx, min(width, end_idx + 1)):
            line[idx] = char

    fill.position = 0.0  # type: ignore[attr-defined]
    fill(profile.ltr5_len, "=")
    fill(profile.internal_len, "-")
    fill(profile.ltr3_len, "=")

    identity_span = _format_identity_range(profile.identity_range)
    meta = f"> AVG {profile.group_label}:{profile.label}  range:{identity_span}"
    table_lines = _build_table_rows(profile, identity_span)
    bar = "".join(line)
    if not ruler:
        lines = [meta]
        lines.extend(f"  {row}" for row in table_lines)
        lines.extend(["", bar])
        return "\n".join(lines) + "\n"

    ticks = [
        ("LTR5", profile.ltr5_len),
        ("INT", profile.internal_len),
        ("LTR3", profile.ltr3_len),
    ]
    cumulative = 0.0
    labels: List[str] = []
    for name, length in ticks:
        cumulative += length
        labels.append(f"{name}:{cumulative:.0f}")
    ruler_line = " | ".join(labels)
    lines = [meta]
    lines.extend(f"  {row}" for row in table_lines)
    lines.extend(["", bar, ruler_line])
    if show_quantiles and profile.q25_len is not None and profile.q75_len is not None:
        lines.append(f"Q25:{profile.q25_len:.0f}bp  Q75:{profile.q75_len:.0f}bp")
    return "\n".join(lines) + "\n"


def _scale(value: float, total: float, width: int) -> int:
    if total <= 0:
        return 0
    rel = value / total
    return max(0, min(width - 1, int(round(rel * (width - 1)))))


def _within_identity(
    elem: RepeatRegion,
    min_id: Optional[float],
    max_id: Optional[float],
) -> bool:
    if elem.ltr_identity is None:
        return False
    if min_id is not None and elem.ltr_identity < min_id:
        return False
    if max_id is not None and elem.ltr_identity > max_id:
        return False
    return True


def _format_identity_range(range_tuple: Optional[Tuple[Optional[float], Optional[float]]]) -> str:
    if not range_tuple:
        return "all"
    min_id, max_id = range_tuple
    min_str = "-" if min_id is None else f"{min_id:.3f}"
    max_str = "-" if max_id is None else f"{max_id:.3f}"
    return f"{min_str}–{max_str}"


def _format_float(value: Optional[float]) -> str:
    return "NA" if value is None else f"{value:.3f}"


def _quantile(values: Sequence[float], q: float) -> Optional[float]:
    if not values:
        return None
    if q <= 0:
        return float(min(values))
    if q >= 1:
        return float(max(values))
    ordered = sorted(values)
    idx = (len(ordered) - 1) * q
    lower = int(idx)
    upper = min(len(ordered) - 1, lower + 1)
    weight = idx - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def _format_strand_summary(counts: Counter, total: int) -> str:
    plus = counts.get("+", 0)
    minus = counts.get("-", 0)
    if total <= 0 or (plus == 0 and minus == 0):
        return "strand:NA"
    parts = []
    for label, count in [("+", plus), ("-", minus)]:
        if count:
            perc = (count / total) * 100
            parts.append(f"{label}:{count} ({perc:.0f}%)")
    return "strand " + "/".join(parts)


def _build_table_rows(profile: AverageMapProfile, identity_span: str) -> List[str]:
    rows: List[str] = []
    n_line = f"n: {profile.n_elements}"
    if profile.note:
        n_line += f" ({profile.note})"
    rows.append(n_line)
    rows.append(f"median len: {profile.total_len:.1f} bp")
    if profile.q25_len is not None and profile.q75_len is not None:
        rows.append(f"IQR: {profile.q25_len:.0f}–{profile.q75_len:.0f} bp")
    rows.append(
        f"identity mean:{_format_float(profile.mean_identity)} median:{_format_float(profile.median_identity)} range:{identity_span}"
    )
    rows.append(f"LTR5 median: {profile.ltr5_len:.1f} bp   LTR3 median: {profile.ltr3_len:.1f} bp")
    rows.append(f"Internal median: {profile.internal_len:.1f} bp")
    rows.append(f"TSD pair: {profile.pct_with_both_tsd:.1f}%   median TSD len: {profile.median_tsd_len or 'NA'}")
    rows.append(profile.strand_summary)
    return rows
