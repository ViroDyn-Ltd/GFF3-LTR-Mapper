"""SVG postcard renderer built on svgwrite."""

from __future__ import annotations

from typing import Dict

import svgwrite

from .model import RepeatRegion

PALETTES: Dict[str, Dict[str, str]] = {
    "classic": {
        "LTR": "#4E79A7",
        "INTERNAL": "#A0CBE8",
        "TSD": "#333333",
        "TEXT": "#111111",
        "BG": "#FFFFFF",
    },
    "mono": {
        "LTR": "#666666",
        "INTERNAL": "#BBBBBB",
        "TSD": "#111111",
        "TEXT": "#000000",
        "BG": "#FFFFFF",
    },
    "protanopia": {
        "LTR": "#0072B2",
        "INTERNAL": "#56B4E9",
        "TSD": "#444444",
        "TEXT": "#111111",
        "BG": "#FFFFFF",
    },
}


def _scale(pos: int, start: int, end: int, width: float) -> float:
    span = max(1, end - start)
    return (pos - start) / span * width


def svg_map(
    elem: RepeatRegion,
    *,
    width: int = 800,
    height: int = 80,
    ruler: bool = False,
    palette: str = "classic",
) -> svgwrite.Drawing:
    colors = PALETTES.get(palette, PALETTES["classic"])
    pad_x = 12
    body_width = width - 2 * pad_x
    mid_y = height / 2
    ltr_height = 24
    internal_height = 14

    dwg = svgwrite.Drawing(size=(width, height))
    dwg.add(dwg.rect(insert=(0, 0), size=(width, height), fill=colors["BG"]))

    arrow = "\u2192" if elem.strand == "+" else "\u2190"
    meta = (
        f"{elem.scaffold}:{elem.id}  fam:{elem.superfamily or 'NA'}  span:{elem.start}-{elem.end}  "
        f"ltr_id:{elem.ltr_identity if elem.ltr_identity is not None else 'NA'}  "
        f"motif:{elem.motif or 'NA'}  tsd:{elem.tsd or 'NA'}"
    )
    dwg.add(
        dwg.text(
            arrow,
            insert=(pad_x, 20),
            font_size=16,
            fill=colors["TEXT"],
        )
    )
    dwg.add(
        dwg.text(
            meta,
            insert=(pad_x + 20, 20),
            font_size=12,
            fill=colors["TEXT"],
        )
    )

    def draw_span(start: int, end: int, *, color: str, height_px: float, y_center: float) -> None:
        x1 = pad_x + _scale(start, elem.start, elem.end, body_width)
        x2 = pad_x + _scale(end, elem.start, elem.end, body_width)
        w = max(1.0, x2 - x1)
        dwg.add(
            dwg.rect(
                insert=(x1, y_center - height_px / 2),
                size=(w, height_px),
                rx=4,
                ry=4,
                fill=color,
            )
        )

    if elem.ltr5:
        draw_span(elem.ltr5.start, elem.ltr5.end, color=colors["LTR"], height_px=ltr_height, y_center=mid_y)
    if elem.ltr3:
        draw_span(elem.ltr3.start, elem.ltr3.end, color=colors["LTR"], height_px=ltr_height, y_center=mid_y)
    if elem.ltr5 and elem.ltr3 and elem.ltr3.start - elem.ltr5.end > 1:
        draw_span(
            elem.ltr5.end + 1,
            elem.ltr3.start - 1,
            color=colors["INTERNAL"],
            height_px=internal_height,
            y_center=mid_y,
        )

    for pos in [elem.tsd5, elem.tsd3]:
        if pos is None:
            continue
        x = pad_x + _scale(pos, elem.start, elem.end, body_width)
        dwg.add(
            dwg.line(
                start=(x, mid_y - ltr_height / 2 - 6),
                end=(x, mid_y + ltr_height / 2 + 6),
                stroke=colors["TSD"],
                stroke_width=2,
            )
        )

    if ruler:
        dwg.add(dwg.text(str(elem.start), insert=(pad_x, height - 6), font_size=10, fill=colors["TEXT"]))
        end_text = dwg.text(str(elem.end), insert=(width - pad_x, height - 6), font_size=10, fill=colors["TEXT"])
        end_text.attribs["text-anchor"] = "end"
        dwg.add(end_text)

    return dwg
