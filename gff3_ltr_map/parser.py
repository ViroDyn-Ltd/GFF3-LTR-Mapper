"""GFF3 reader specialized for EDTA intact LTR annotations."""

from __future__ import annotations

import gzip
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, MutableMapping, Optional, Tuple
from urllib.parse import unquote

from .model import FeatureSpan, LTRChild, RepeatRegion

AttrMap = MutableMapping[str, str]

RETROTRANSPOSON_TYPES = {
    "ltr_retrotransposon",
    "gypsy_ltr_retrotransposon",
    "copia_ltr_retrotransposon",
}
VALID_STRANDS = {"+", "-", ".", "?"}


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
        attrs[key.strip().lower()] = unquote(value.strip())
    return attrs


def _open_text(path: str) -> Iterator[str]:
    gff_path = Path(path)
    if gff_path.suffix == ".gz":
        with gzip.open(gff_path, "rt", encoding="utf-8") as handle:
            yield from handle
        return
    with gff_path.open("r", encoding="utf-8") as handle:
        yield from handle


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
    line_number: int

    @classmethod
    def from_line(cls, line: str, line_number: int) -> "GFFRow":
        parts = line.rstrip().split("\t")
        if len(parts) != 9:
            raise ValueError(
                f"Expected 9 columns at line {line_number}, got {len(parts)} in line: {line!r}"
            )
        seqid, source, type_, start_text, end_text, score, strand, phase, attrs = parts
        try:
            start = int(start_text)
            end = int(end_text)
        except ValueError as exc:
            raise ValueError(
                f"Invalid coordinates at line {line_number}: start={start_text!r}, end={end_text!r}"
            ) from exc
        if start <= 0 or end <= 0 or end < start:
            raise ValueError(
                f"Invalid coordinate span at line {line_number}: start={start}, end={end}"
            )
        return cls(
            seqid=seqid,
            source=source,
            type=type_,
            start=start,
            end=end,
            score=score,
            strand=strand,
            phase=phase,
            attrs=parse_attrs(attrs),
            line_number=line_number,
        )

    def attr(self, key: str) -> Optional[str]:
        return self.attrs.get(key.lower())


Regions = Dict[str, GFFRow]
Children = Dict[str, List[GFFRow]]


def load_gff(gff_path: str) -> Tuple[Regions, Children]:
    """Read GFF3 and partition rows by repeat_region parent."""

    regions: Regions = {}
    by_parent: Children = defaultdict(list)
    for line_number, line in enumerate(_open_text(gff_path), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        row = GFFRow.from_line(line, line_number)
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


def _first_attr(parent: GFFRow, children: Iterable[GFFRow], *keys: str) -> Optional[str]:
    lowered = tuple(key.lower() for key in keys)
    for key in lowered:
        value = parent.attr(key)
        if value:
            return value
    for child in children:
        for key in lowered:
            value = child.attr(key)
            if value:
                return value
    return None


def _first_float(parent: GFFRow, children: Iterable[GFFRow], *keys: str) -> Optional[float]:
    text = _first_attr(parent, children, *keys)
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _normalize_classification(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    parts = [part for part in value.split("/") if part]
    if not parts:
        return value
    prefix = parts[0].upper()
    remainder = "/".join(parts[1:])
    if not remainder:
        return prefix
    if remainder.lower() == "unknown":
        return f"{prefix}/unknown"
    return f"{prefix}/{remainder[0].upper()}{remainder[1:]}"


def _superfamily_from_classification(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    lowered = value.lower()
    if not lowered.startswith("ltr/"):
        return None
    tail = value.split("/", 1)[1]
    if not tail:
        return None
    if tail.lower() == "unknown":
        return "Unknown"
    return f"{tail[0].upper()}{tail[1:]}"


def classification_of(parent: GFFRow, children: Iterable[GFFRow]) -> Optional[str]:
    return _normalize_classification(_first_attr(parent, children, "classification", "class"))


def ltr_identity_of(parent: GFFRow, children: Iterable[GFFRow]) -> Optional[float]:
    return _first_float(parent, children, "ltr_identity", "ltrid", "identity")


def motif_of(parent: GFFRow, children: Iterable[GFFRow]) -> Optional[str]:
    return _first_attr(parent, children, "motif")


def tsd_of(parent: GFFRow, children: Iterable[GFFRow]) -> Optional[str]:
    value = _first_attr(parent, children, "tsd")
    if not value:
        return None
    if "_" in value:
        parts = [chunk for chunk in value.split("_") if chunk]
        if len(parts) >= 2 and parts[0] == parts[1]:
            return parts[0]
    return value


def _feature_span(row: GFFRow) -> FeatureSpan:
    return FeatureSpan(
        feature_id=row.attr("id") or "",
        feature_type=row.type,
        start=row.start,
        end=row.end,
    )


def pick_retrotransposon(children: Iterable[GFFRow]) -> Optional[FeatureSpan]:
    retro_rows = [
        child
        for child in children
        if child.type.lower() in RETROTRANSPOSON_TYPES
    ]
    if len(retro_rows) != 1:
        return None
    return _feature_span(retro_rows[0])


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


def pick_tsd_features(children: Iterable[GFFRow]) -> Tuple[Optional[FeatureSpan], Optional[FeatureSpan]]:
    tsds = [child for child in children if child.type.lower() == "target_site_duplication"]
    if not tsds:
        return None, None
    tsds.sort(key=lambda row: (row.start, row.end))
    if len(tsds) == 1:
        return _feature_span(tsds[0]), None
    return _feature_span(tsds[0]), _feature_span(tsds[-1])


def pick_tsd_positions(children: Iterable[GFFRow]) -> Tuple[Optional[int], Optional[int]]:
    tsd5_span, tsd3_span = pick_tsd_features(children)
    tsd5 = tsd5_span.end if tsd5_span else None
    tsd3 = tsd3_span.start if tsd3_span else None
    return tsd5, tsd3


def _child_signature(children: Iterable[GFFRow]) -> str:
    counts = Counter(child.type for child in children)
    if not counts:
        return ""
    return "; ".join(f"{feature}={count}" for feature, count in sorted(counts.items()))


def _validate_structure(
    parent: GFFRow,
    *,
    classification: Optional[str],
    method: Optional[str],
    ltr_identity: Optional[float],
    retro: Optional[FeatureSpan],
    ltr5: Optional[LTRChild],
    ltr3: Optional[LTRChild],
    tsd5_span: Optional[FeatureSpan],
    tsd3_span: Optional[FeatureSpan],
    children: List[GFFRow],
) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    warnings: List[str] = []
    errors: List[str] = []

    if classification is None:
        warnings.append("MISSING_CLASSIFICATION")
    elif not classification.lower().startswith("ltr/"):
        errors.append("NON_LTR_CLASSIFICATION")

    if method and method.lower() != "structural":
        warnings.append("NON_STRUCTURAL_METHOD")

    retro_count = sum(1 for child in children if child.type.lower() in RETROTRANSPOSON_TYPES)
    if retro_count != 1:
        errors.append("EXPECTED_ONE_RETROTRANSPOSON")

    ltr_rows = [child for child in children if child.type.lower() == "long_terminal_repeat"]
    if len(ltr_rows) != 2:
        errors.append("EXPECTED_TWO_LTRS")

    tsd_rows = [child for child in children if child.type.lower() == "target_site_duplication"]
    if len(tsd_rows) not in (0, 2):
        errors.append("EXPECTED_ZERO_OR_TWO_TSDS")
    if len(tsd_rows) == 0:
        warnings.append("MISSING_TSD_PAIR")

    if parent.strand not in VALID_STRANDS:
        errors.append("INVALID_STRAND")

    if ltr_identity is None:
        warnings.append("MISSING_LTR_IDENTITY")
    elif not (0.0 <= ltr_identity <= 1.0):
        errors.append("LTR_IDENTITY_OUT_OF_RANGE")

    if ltr5 and ltr3 and ltr5.end >= ltr3.start:
        errors.append("LTRS_OVERLAP")

    if retro and ltr5 and retro.start > ltr5.start:
        errors.append("RETRO_START_AFTER_LTR5")
    if retro and ltr3 and retro.end < ltr3.end:
        errors.append("RETRO_END_BEFORE_LTR3")
    if retro and ltr5 and retro.start != ltr5.start:
        warnings.append("RETRO_START_DIFFERS_FROM_LTR5")
    if retro and ltr3 and retro.end != ltr3.end:
        warnings.append("RETRO_END_DIFFERS_FROM_LTR3")

    if ltr5 and (parent.start > ltr5.start or parent.end < ltr5.end):
        errors.append("PARENT_DOES_NOT_COVER_LTR5")
    if ltr3 and (parent.start > ltr3.start or parent.end < ltr3.end):
        errors.append("PARENT_DOES_NOT_COVER_LTR3")

    if tsd5_span and ltr5 and tsd5_span.end > ltr5.start:
        warnings.append("LEFT_TSD_OVERLAPS_LTR")
    if tsd3_span and ltr3 and tsd3_span.start < ltr3.end:
        warnings.append("RIGHT_TSD_OVERLAPS_LTR")

    return tuple(sorted(set(errors))), tuple(sorted(set(warnings)))


def to_repeat_region_object(rid: str, parent: GFFRow, children: List[GFFRow]) -> RepeatRegion:
    classification = classification_of(parent, children)
    ltr5, ltr3 = pick_ltrs(children)
    tsd5_span, tsd3_span = pick_tsd_features(children)
    tsd5, tsd3 = pick_tsd_positions(children)
    retro = pick_retrotransposon(children)
    method = _first_attr(parent, children, "method")
    ltr_identity = ltr_identity_of(parent, children)
    warnings, errors = (), ()
    errors, warnings = _validate_structure(
        parent,
        classification=classification,
        method=method,
        ltr_identity=ltr_identity,
        retro=retro,
        ltr5=ltr5,
        ltr3=ltr3,
        tsd5_span=tsd5_span,
        tsd3_span=tsd3_span,
        children=children,
    )
    return RepeatRegion(
        id=rid,
        scaffold=parent.seqid,
        start=parent.start,
        end=parent.end,
        strand=parent.strand,
        source=parent.source,
        name=parent.attr("name"),
        classification=classification,
        superfamily=_superfamily_from_classification(classification),
        method=method,
        ltr_identity=ltr_identity,
        motif=motif_of(parent, children),
        tsd=tsd_of(parent, children),
        retrotransposon=retro,
        ltr5=ltr5,
        ltr3=ltr3,
        tsd5=tsd5,
        tsd3=tsd3,
        tsd5_span=tsd5_span,
        tsd3_span=tsd3_span,
        child_signature=_child_signature(children),
        validation_errors=errors,
        validation_warnings=warnings,
        n_children_warn=len(children),
    )
