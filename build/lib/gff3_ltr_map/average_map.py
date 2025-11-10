"""Utilities to build ASCII/SVG average maps from multiple RepeatRegion entries."""

from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
from statistics import mean, median
from typing import Iterable, List, Optional, Sequence, Tuple

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
        filtered = [
            elem
            for elem in elements
            if _within_identity(elem, min_id, max_id)
        ]
    else:
        filtered = list(elements)

    if not filtered:
        return None

    identities = [elem.ltr_identity for elem in filtered if elem.ltr_identity is not None]
    ltr5 = [elem.ltr5_len for elem in filtered if elem.ltr5_len is not None]
    ltr3 = [elem.ltr3_len for elem in filtered if elem.ltr3_len is not None]
    internal = [elem.internal_len for elem in filtered if elem.internal_len is not None]
    totals = [elem.length_bp for elem in filtered]

    ltr5_len = median(ltr5) if ltr5 else 0.0
    ltr3_len = median(ltr3) if ltr3 else 0.0
    total_len = median(totals) if totals else ltr5_len + ltr3_len
    internal_len = median(internal) if internal else max(total_len - ltr5_len - ltr3_len, 0.0)

    # Ensure components add up sensibly
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


def average_svg_map(
    profile: AverageMapProfile,
    *,
    width: int = 800,
    height: int = 120,
    palette: str = "classic",
    show_quantiles: bool = False,
) -> "svgwrite.Drawing":
    """Render an aggregate profile as SVG."""

    try:
        import svgwrite
    except ImportError as exc:  # pragma: no cover - import-time guard
        raise RuntimeError("svgwrite is required for average SVG maps") from exc

    from .render_svg import PALETTES, _scale as svg_scale

    colors = PALETTES.get(palette, PALETTES["classic"])
    pad = 30
    ltr_height = 30
    internal_height = 18
    identity_span = _format_identity_range(profile.identity_range)
    table_lines = _build_table_rows(profile, identity_span)
    header = f"AVG {profile.group_label}:{profile.label}"

    header_y = 20
    table_line_height = 16
    table_start_y = header_y + 20
    table_block_height = len(table_lines) * table_line_height
    body_top = table_start_y + table_block_height + 20
    required_height = body_top + ltr_height + 50
    canvas_height = max(height, required_height)
    body_width = width - 2 * pad
    mid_y = body_top + ltr_height / 2

    dwg = svgwrite.Drawing(size=(width, canvas_height))
    dwg.add(dwg.rect(insert=(0, 0), size=(width, canvas_height), fill=colors["BG"]))

    dwg.add(dwg.text(header, insert=(pad, header_y), font_size=16, fill=colors["TEXT"]))
    for idx, row in enumerate(table_lines):
        y = table_start_y + idx * table_line_height
        dwg.add(dwg.text(row, insert=(pad, y), font_size=12, fill=colors["TEXT"]))

    def draw(length: float, color_key: str, height_px: float, *, label: Optional[str] = None, label_color: str = "#FFFFFF", label_offset: float = 0.0) -> Tuple[float, float]:
        start = draw.position
        draw.position += length
        x1 = pad + svg_scale(start, 0, profile.total_len, body_width)
        x2 = pad + svg_scale(draw.position, 0, profile.total_len, body_width)
        rect_width = max(1.0, x2 - x1)
        dwg.add(
            dwg.rect(
                insert=(x1, mid_y - height_px / 2),
                size=(rect_width, height_px),
                rx=6,
                ry=6,
                fill=colors[color_key],
            )
        )
        if label:
            center = (x1 + x2) / 2
            text = dwg.text(label, insert=(center, mid_y + label_offset), font_size=10, fill=label_color)
            text.attribs["text-anchor"] = "middle"
            text.attribs["dominant-baseline"] = "middle"
            dwg.add(text)
        return x1, x2

    draw.position = 0.0  # type: ignore[attr-defined]
    draw(
        profile.ltr5_len,
        "LTR",
        ltr_height,
        label=f"{profile.ltr5_len:.0f}bp",
        label_color="#FFFFFF",
    )
    draw(
        profile.internal_len,
        "INTERNAL",
        internal_height,
        label=f"{profile.internal_len:.0f}bp",
        label_color=colors["TEXT"],
        label_offset=-ltr_height,
    )
    draw(
        profile.ltr3_len,
        "LTR",
        ltr_height,
        label=f"{profile.ltr3_len:.0f}bp",
        label_color="#FFFFFF",
    )

    if show_quantiles and profile.q25_len is not None and profile.q75_len is not None:
        for q_label, q_val in (("Q25", profile.q25_len), ("Q75", profile.q75_len)):
            q_val_clamped = max(0.0, min(profile.total_len, q_val))
            x = pad + svg_scale(q_val_clamped, 0, profile.total_len, body_width)
            dwg.add(
                dwg.line(
                    start=(x, mid_y - ltr_height),
                    end=(x, mid_y + ltr_height),
                    stroke=colors["TSD"],
                    stroke_dasharray="4,2",
                    stroke_width=2,
                )
            )
            dwg.add(dwg.text(q_label, insert=(x + 4, mid_y - ltr_height - 4), font_size=10, fill=colors["TEXT"]))

    return dwg


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
    rows.append(profile.strand_summary)
    return rows
