"""Data structures used across parsing and rendering."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class LTRChild:
    """Basic coordinate span for an LTR child feature."""

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
    superfamily: Optional[str]
    ltr_identity: Optional[float]
    motif: Optional[str]
    tsd: Optional[str]
    ltr5: Optional[LTRChild]
    ltr3: Optional[LTRChild]
    tsd5: Optional[int]
    tsd3: Optional[int]
    n_children_warn: int = 0

    @property
    def length_bp(self) -> int:
        return self.end - self.start + 1

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
