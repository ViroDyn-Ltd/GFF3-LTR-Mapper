"""Public package exports for gff3_ltr_map."""

from importlib.metadata import version, PackageNotFoundError

__all__ = [
    "parser",
    "model",
    "render_ascii",
    "render_svg",
    "summary",
    "aggregates",
    "average_map",
]

try:
    __version__ = version("gff3-ltr-map")
except PackageNotFoundError:  # pragma: no cover - during local source usage
    __version__ = "0.0.0"
