"""GFF3 reader specialized for EDTA "intact" LTR annotations."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, MutableMapping, Optional, Tuple

from .model import LTRChild, RepeatRegion

AttrMap = MutableMapping[str, str]


def parse_attrs(attr_str: str) -> AttrMap:
    """Parse a GFF3 attribute column into a case-insensitive dict."""

    attrs: AttrMap = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
        else:
            key, value = part, ""
        attrs[key.strip().lower()] = value.strip()
    return attrs


@dataclass
class GFFRow:
    """Minimal representation of a GFF3 row."""

    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: AttrMap

    @classmethod
    def from_line(cls, line: str) -> "GFFRow":
        parts = line.rstrip().split("\t")
        if len(parts) != 9:
            raise ValueError(f"Expected 9 columns, got {len(parts)} in line: {line!r}")
        seqid, source, type_, start, end, score, strand, phase, attrs = parts
        return cls(
            seqid=seqid,
            source=source,
            type=type_,
            start=int(start),
            end=int(end),
            score=score,
            strand=strand,
            phase=phase,
            attrs=parse_attrs(attrs),
        )

    def attr(self, key: str) -> Optional[str]:
        return self.attrs.get(key.lower())


Regions = Dict[str, GFFRow]
Children = Dict[str, List[GFFRow]]


def load_gff(gff_path: str) -> Tuple[Regions, Children]:
    """Read GFF3 and partition rows by repeat_region parent."""

    regions: Regions = {}
    by_parent: Children = defaultdict(list)
    with open(gff_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            row = GFFRow.from_line(line)
            rid = row.attr("id")
            if row.type == "repeat_region" and rid:
                regions[rid] = row
                continue
            parent_attr = row.attr("parent")
            if not parent_attr:
                continue
            for parent in parent_attr.split(","):
                parent = parent.strip()
                if parent:
                    by_parent[parent].append(row)
    return regions, by_parent


def classification_of(children: Iterable[GFFRow]) -> Optional[str]:
    for child in children:
        cls_attr = child.attr("classification") or child.attr("class")
        if cls_attr and "LTR/" in cls_attr:
            return cls_attr.split("LTR/")[-1].split(";")[0]
        type_lower = child.type.lower()
        if "gypsy" in type_lower:
            return "Gypsy"
        if "copia" in type_lower:
            return "Copia"
    return None


def ltr_identity_of(children: Iterable[GFFRow]) -> Optional[float]:
    for child in children:
        val = child.attr("ltr_identity") or child.attr("ltrid") or child.attr("identity")
        if not val:
            continue
        try:
            return float(val)
        except ValueError:
            continue
    return None


def motif_of(children: Iterable[GFFRow]) -> Optional[str]:
    for child in children:
        val = child.attr("motif")
        if val:
            return val
    return None


def tsd_of(children: Iterable[GFFRow]) -> Optional[str]:
    for child in children:
        val = child.attr("tsd")
        if val:
            return val
    return None


def pick_ltrs(children: Iterable[GFFRow]) -> Tuple[Optional[LTRChild], Optional[LTRChild]]:
    ltrs = [child for child in children if child.type.lower() == "long_terminal_repeat"]
    if not ltrs:
        return None, None
    ltrs.sort(key=lambda row: (row.start, row.end))
    if len(ltrs) == 1:
        first = ltrs[0]
        return LTRChild(first.start, first.end), None
    left = ltrs[0]
    right = ltrs[-1]
    return LTRChild(left.start, left.end), LTRChild(right.start, right.end)


def pick_tsd_positions(parent: GFFRow, children: Iterable[GFFRow]) -> Tuple[Optional[int], Optional[int]]:
    tsds = [child for child in children if child.type.lower() == "target_site_duplication"]
    if not tsds:
        return None, None
    left = min(tsds, key=lambda row: row.start)
    right = max(tsds, key=lambda row: row.end)
    tsd5 = left.start if abs(left.start - parent.start) < abs(left.end - parent.start) else left.end
    tsd3 = right.end if abs(parent.end - right.end) < abs(parent.end - right.start) else right.start
    return tsd5, tsd3


def to_repeat_region_object(rid: str, parent: GFFRow, children: List[GFFRow]) -> RepeatRegion:
    ltr5, ltr3 = pick_ltrs(children)
    tsd5, tsd3 = pick_tsd_positions(parent, children)
    return RepeatRegion(
        id=rid,
        scaffold=parent.seqid,
        start=parent.start,
        end=parent.end,
        strand=parent.strand,
        superfamily=classification_of(children),
        ltr_identity=ltr_identity_of(children),
        motif=motif_of(children),
        tsd=tsd_of(children),
        ltr5=ltr5,
        ltr3=ltr3,
        tsd5=tsd5,
        tsd3=tsd3,
        n_children_warn=len(children),
    )
