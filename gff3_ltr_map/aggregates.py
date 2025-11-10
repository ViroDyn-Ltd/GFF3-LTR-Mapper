"""Aggregate statistics over collections of RepeatRegion objects."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from statistics import mean, median, pstdev
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

from .model import RepeatRegion

GroupType = str

AGGREGATE_COLUMNS = [
    "group_type",
    "group",
    "n_elements",
    "pct_with_both_tsd",
    "ltr_identity_mean",
    "ltr_identity_median",
    "ltr_identity_stdev",
    "ltr5_len_mean",
    "ltr5_len_median",
    "ltr5_len_stdev",
    "ltr3_len_mean",
    "ltr3_len_median",
    "ltr3_len_stdev",
    "internal_len_mean",
    "internal_len_median",
    "internal_len_stdev",
    "length_bp_mean",
    "length_bp_median",
    "length_bp_stdev",
    "ltr_asymmetry_mean",
    "density_per_Mb",
    "coverage_pct",
    "top_motifs",
    "top_tsd",
    "approx_age_median_Myr",
    "notes",
]


@dataclass
class AggregateRow:
    """Single aggregate row ready for serialization."""

    group_type: GroupType
    group: str
    n_elements: int
    pct_with_both_tsd: float
    ltr_identity_mean: Optional[float]
    ltr_identity_median: Optional[float]
    ltr_identity_stdev: Optional[float]
    ltr5_len_mean: Optional[float]
    ltr5_len_median: Optional[float]
    ltr5_len_stdev: Optional[float]
    ltr3_len_mean: Optional[float]
    ltr3_len_median: Optional[float]
    ltr3_len_stdev: Optional[float]
    internal_len_mean: Optional[float]
    internal_len_median: Optional[float]
    internal_len_stdev: Optional[float]
    length_bp_mean: Optional[float]
    length_bp_median: Optional[float]
    length_bp_stdev: Optional[float]
    ltr_asymmetry_mean: Optional[float]
    density_per_Mb: Optional[float]
    coverage_pct: Optional[float]
    top_motifs: str
    top_tsd: str
    approx_age_median_Myr: Optional[float]
    notes: str


def _describe(values: List[float]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    if not values:
        return None, None, None
    mu = mean(values)
    med = median(values)
    sd = pstdev(values) if len(values) > 1 else 0.0
    return mu, med, sd


def _format_top_counts(values: Iterable[str], limit: int, total: int) -> Tuple[str, bool]:
    counter = Counter(v for v in values if v)
    if not counter:
        return "", False
    ordered = sorted(counter.items(), key=lambda item: (-item[1], item[0]))
    top_items = ordered[:limit]
    parts = []
    for value, count in top_items:
        perc = count / total if total else 0.0
        parts.append(f"{value} ({count}, {perc:.0%})")
    consensus_warn = False
    if total > 0:
        most_common_count = ordered[0][1]
        if (most_common_count / total) < 0.4:
            consensus_warn = True
    text = ", ".join(parts)
    if consensus_warn and text:
        text = f"{text} (no single consensus)"
    return text, consensus_warn


def _ltr_asymmetry(elem: RepeatRegion) -> Optional[float]:
    if not elem.ltr5 or not elem.ltr3:
        return None
    l5 = elem.ltr5.length
    l3 = elem.ltr3.length
    denom = (l5 + l3) / 2
    if denom <= 0:
        return None
    return abs(l5 - l3) / denom


def _age_my(elem: RepeatRegion, substitution_rate: float) -> Optional[float]:
    if elem.ltr_identity is None:
        return None
    identity = elem.ltr_identity
    if identity < 0 or identity > 1:
        return None
    if substitution_rate <= 0:
        return None
    return (1 - identity) / (2 * substitution_rate) / 1_000_000


def _merged_coverage_bp(elements: Sequence[RepeatRegion]) -> int:
    intervals = sorted((e.start, e.end) for e in elements)
    if not intervals:
        return 0
    merged = []
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end
    merged.append((cur_start, cur_end))
    return sum(end - start + 1 for start, end in merged)


def compute_aggregates(
    elements: Sequence[RepeatRegion],
    group_types: Sequence[GroupType],
    *,
    scaffold_lengths: Optional[Mapping[str, int]] = None,
    substitution_rate: Optional[float] = None,
    top_k: int = 3,
    min_n: int = 1,
) -> List[AggregateRow]:
    if not elements or not group_types:
        return []

    buckets: MutableMapping[Tuple[GroupType, str], List[RepeatRegion]] = defaultdict(list)
    for elem in elements:
        for group_type in group_types:
            key = _group_key(elem, group_type)
            if key is None:
                continue
            buckets[(group_type, key)].append(elem)

    rows: List[AggregateRow] = []
    for (group_type, group_name), elems in sorted(buckets.items()):
        rows.append(
            _summarize_group(
                group_type,
                group_name,
                elems,
                scaffold_lengths=scaffold_lengths,
                substitution_rate=substitution_rate,
                top_k=top_k,
                min_n=min_n,
            )
        )
    return rows


def _group_key(elem: RepeatRegion, group_type: GroupType) -> Optional[str]:
    if group_type == "genome":
        return "genome"
    if group_type == "superfamily":
        return elem.superfamily or "NA"
    if group_type == "scaffold":
        return elem.scaffold
    return None


def summarize_cohort(
    label: str,
    elems: Sequence[RepeatRegion],
    *,
    substitution_rate: Optional[float],
    top_k: int,
    min_n: int,
) -> AggregateRow:
    return _summarize_group(
        "identity",
        label,
        elems,
        scaffold_lengths=None,
        substitution_rate=substitution_rate,
        top_k=top_k,
        min_n=min_n,
    )


def _summarize_group(
    group_type: GroupType,
    group_name: str,
    elems: Sequence[RepeatRegion],
    *,
    scaffold_lengths: Optional[Mapping[str, int]] = None,
    substitution_rate: Optional[float] = None,
    top_k: int = 3,
    min_n: int = 1,
) -> AggregateRow:
    n = len(elems)
    pct_both = (sum(1 for e in elems if e.has_both_tsd) / n) * 100 if n else 0.0

    ltr_identity = [e.ltr_identity for e in elems if e.ltr_identity is not None]
    ltr5 = [e.ltr5_len for e in elems if e.ltr5_len is not None]
    ltr3 = [e.ltr3_len for e in elems if e.ltr3_len is not None]
    internal = [e.internal_len for e in elems if e.internal_len is not None]
    length_bp = [e.length_bp for e in elems]
    asymmetries = [val for e in elems if (val := _ltr_asymmetry(e)) is not None]
    ages = (
        [val for e in elems if substitution_rate and (val := _age_my(e, substitution_rate)) is not None]
        if substitution_rate
        else []
    )

    identity_stats = _describe(ltr_identity)
    ltr5_stats = _describe(ltr5)
    ltr3_stats = _describe(ltr3)
    internal_stats = _describe(internal)
    length_stats = _describe(length_bp)

    density = None
    coverage_pct = None
    if group_type == "scaffold" and scaffold_lengths:
        scaffold_len = scaffold_lengths.get(group_name)
        if scaffold_len and scaffold_len > 0:
            density = n / (scaffold_len / 1_000_000)
            covered = _merged_coverage_bp(elems)
            coverage_pct = min(100.0, (covered / scaffold_len) * 100)

    notes = []
    if n and n < min_n:
        notes.append(f"LOW N (n={n})")
    elif n == 0:
        notes.append("NO DATA")

    motifs_text, motifs_warn = _format_top_counts(((e.motif or "") for e in elems), top_k, n)
    tsd_text, tsd_warn = _format_top_counts(((e.tsd or "") for e in elems), top_k, n)
    if motifs_warn:
        notes.append("motif consensus <40%")
    if tsd_warn:
        notes.append("TSD consensus <40%")

    return AggregateRow(
        group_type=group_type,
        group=group_name,
        n_elements=n,
        pct_with_both_tsd=pct_both,
        ltr_identity_mean=identity_stats[0],
        ltr_identity_median=identity_stats[1],
        ltr_identity_stdev=identity_stats[2],
        ltr5_len_mean=ltr5_stats[0],
        ltr5_len_median=ltr5_stats[1],
        ltr5_len_stdev=ltr5_stats[2],
        ltr3_len_mean=ltr3_stats[0],
        ltr3_len_median=ltr3_stats[1],
        ltr3_len_stdev=ltr3_stats[2],
        internal_len_mean=internal_stats[0],
        internal_len_median=internal_stats[1],
        internal_len_stdev=internal_stats[2],
        length_bp_mean=length_stats[0],
        length_bp_median=length_stats[1],
        length_bp_stdev=length_stats[2],
        ltr_asymmetry_mean=mean(asymmetries) if asymmetries else None,
        density_per_Mb=density,
        coverage_pct=coverage_pct,
        top_motifs=motifs_text,
        top_tsd=tsd_text,
        approx_age_median_Myr=median(ages) if ages else None,
        notes="; ".join(notes),
    )


def write_aggregate_tsv(rows: Sequence[AggregateRow], handle) -> None:
    handle.write("\t".join(AGGREGATE_COLUMNS) + "\n")
    for row in rows:
        payload = asdict(row)
        handle.write(
            "\t".join(_format_value(payload[column]) for column in AGGREGATE_COLUMNS) + "\n"
        )


def write_aggregate_json(rows: Sequence[AggregateRow], handle) -> None:
    json.dump([asdict(row) for row in rows], handle, indent=2)
    handle.write("\n")


def _format_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)
