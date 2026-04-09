"""Data structures used across parsing, validation, and rendering."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple


@dataclass(frozen=True)
class LTRChild:
    """Basic coordinate span for an LTR child feature."""

    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class FeatureSpan:
    """Coordinate span with feature metadata."""

    feature_id: str
    feature_type: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class RepeatRegion:
    """Normalized representation of a repeat_region block."""

    id: str
    scaffold: str
    start: int
    end: int
    strand: str
    source: str
    name: Optional[str]
    classification: Optional[str]
    superfamily: Optional[str]
    method: Optional[str]
    ltr_identity: Optional[float]
    motif: Optional[str]
    tsd: Optional[str]
    retrotransposon: Optional[FeatureSpan]
    ltr5: Optional[LTRChild]
    ltr3: Optional[LTRChild]
    tsd5: Optional[int]
    tsd3: Optional[int]
    tsd5_span: Optional[FeatureSpan] = None
    tsd3_span: Optional[FeatureSpan] = None
    child_signature: str = ""
    validation_errors: Tuple[str, ...] = field(default_factory=tuple)
    validation_warnings: Tuple[str, ...] = field(default_factory=tuple)
    n_children_warn: int = 0

    @property
    def length_bp(self) -> int:
        return self.end - self.start + 1

    @property
    def retrotransposon_type(self) -> Optional[str]:
        if not self.retrotransposon:
            return None
        return self.retrotransposon.feature_type

    @property
    def internal_len(self) -> Optional[int]:
        if self.ltr5 and self.ltr3:
            return max(0, self.ltr3.start - self.ltr5.end - 1)
        return None

    @property
    def ltr5_len(self) -> Optional[int]:
        return self.ltr5.length if self.ltr5 else None

    @property
    def ltr3_len(self) -> Optional[int]:
        return self.ltr3.length if self.ltr3 else None

    @property
    def has_both_tsd(self) -> bool:
        return self.tsd5 is not None and self.tsd3 is not None

    @property
    def tsd_len(self) -> Optional[int]:
        if self.tsd and self.tsd != "NA":
            return len(self.tsd)
        if self.tsd5_span and self.tsd3_span and self.tsd5_span.length == self.tsd3_span.length:
            return self.tsd5_span.length
        return None

    @property
    def qc_status(self) -> str:
        if self.validation_errors:
            return "FAIL"
        if self.validation_warnings:
            return "WARN"
        return "PASS"

    @property
    def is_intact(self) -> bool:
        return (
            not self.validation_errors
            and self.retrotransposon is not None
            and self.ltr5 is not None
            and self.ltr3 is not None
        )

    @property
    def validation_warning_text(self) -> str:
        return "; ".join(self.validation_warnings)

    @property
    def validation_error_text(self) -> str:
        return "; ".join(self.validation_errors)
