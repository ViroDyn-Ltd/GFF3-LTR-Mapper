"""Console entry-point for gff3-ltr-map."""

from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

from .aggregates import AggregateRow, compute_aggregates, summarize_cohort, write_aggregate_json, write_aggregate_tsv
from .average_map import average_ascii_map, average_svg_map, build_average_profile
from .batch import (
    BatchSampleRow,
    BatchSuperfamilyRow,
    summarize_sample,
    summarize_sample_superfamilies,
    write_batch_sample_tsv,
    write_batch_superfamily_tsv,
)
from .model import RepeatRegion
from .parser import load_gff, to_repeat_region_object
from .render_ascii import ascii_map
from .render_svg import svg_map
from .scientist_view import (
    ScientistRegionRow,
    build_scientist_region_row,
    build_scientist_element_rows,
    print_scientist_region_table,
    write_scientist_region_tsv,
    write_scientist_element_tsv,
)
from .summary import write_summary, write_validation_report


@dataclass(frozen=True)
class ScopeSpec:
    mode: str  # "genome" or "chrom"
    value: Optional[str]

    @property
    def label(self) -> str:
        return self.value if self.mode == "chrom" and self.value else "genome"

    def matches(self, elem: RepeatRegion) -> bool:
        if self.mode == "chrom" and self.value:
            return elem.scaffold == self.value
        return True


@dataclass(frozen=True)
class RegionFilter:
    scaffold: str
    start: int
    end: int

    @property
    def label(self) -> str:
        return f"{self.scaffold}:{self.start}-{self.end}"

    def matches(self, elem: RepeatRegion) -> bool:
        return elem.scaffold == self.scaffold and elem.end >= self.start and elem.start <= self.end


@dataclass(frozen=True)
class IdentityBin:
    label: str
    min_identity: Optional[float]
    max_identity: Optional[float]

    def matches(self, elem: RepeatRegion) -> bool:
        if elem.ltr_identity is None:
            return False
        if self.min_identity is not None and elem.ltr_identity < self.min_identity:
            return False
        if self.max_identity is not None and elem.ltr_identity > self.max_identity:
            return False
        return True

    @property
    def range_tuple(self) -> Tuple[Optional[float], Optional[float]]:
        return self.min_identity, self.max_identity

    @property
    def slug(self) -> str:
        return _slugify(self.label)


@dataclass
class IdentityReport:
    scope: ScopeSpec
    bin: IdentityBin
    elements: List[RepeatRegion]
    aggregate: AggregateRow
    profile: Optional["AverageMapProfile"]


@dataclass(frozen=True)
class SamplePaths:
    outdir: str
    summary: str
    validation: str
    identity_tsv: str
    identity_json: Optional[str]
    group_tsv: Optional[str]
    group_json: Optional[str]
    scientist_tsv: str
    scientist_elements_tsv: str
    bed: Optional[str]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="EDTA intact GFF3 -> QC-validated LTR summaries, aggregates, and postcards",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input_path", help="Input EDTA intact GFF3, .gff3.gz, or directory of such files")
    parser.add_argument("--outdir", default="out", help="Output directory for map and report artifacts")
    parser.add_argument("--ascii", action="store_true", help="Emit ASCII maps (*.txt) for intact elements")
    parser.add_argument("--svg", action="store_true", help="Emit SVG maps (*.svg) for intact elements")
    parser.add_argument("--summary", default="summary.tsv", help="Single-file mode: path to summary TSV")
    parser.add_argument("--validation-report", help="Single-file mode: path to validation TSV")
    parser.add_argument(
        "--scientist-tsv",
        help=(
            "Simplified scientist-facing region TSV. Single-file mode defaults to <outdir>/scientist_regions.tsv; "
            "directory mode defaults to <outdir>/scientist_batch.tsv."
        ),
    )
    parser.add_argument(
        "--scientist-cli-max-rows",
        type=int,
        default=25,
        help="Maximum scientist-facing rows to print directly to the CLI",
    )
    parser.add_argument(
        "--region",
        help=(
            "Optional coordinate filter. Use scaffold:start-end, or start-end when --chrom names one scaffold. "
            "If no scaffold is supplied in directory mode, the tool infers one primary scaffold per sample. "
            "Elements overlapping the interval are retained."
        ),
    )
    parser.add_argument("--width", type=int, default=800, help="SVG canvas width (px) or ASCII columns")
    parser.add_argument("--height", type=int, default=80, help="SVG canvas height (px)")
    parser.add_argument(
        "--palette",
        choices=["classic", "mono", "protanopia"],
        default="classic",
        help="Color palette for SVG maps",
    )
    parser.add_argument("--bed", help="Single-file mode: also emit BED with intact repeat_region spans")
    parser.add_argument("--index-html", action="store_true", help="Generate HTML index embedding SVGs")
    parser.add_argument("--ruler", action="store_true", help="Add coordinate ruler to ASCII/SVG outputs")
    parser.add_argument("--workers", type=int, default=1, help="Parallel workers for per-element rendering")
    parser.add_argument("--max-elements", type=int, help="Stop after N elements (debug)")
    parser.add_argument("--limit-files", type=int, help="Directory mode: stop after N input files")
    parser.add_argument("--chrom", help="Chromosome scope (default: genome + each chromosome)")
    parser.add_argument("--scope", help=argparse.SUPPRESS)
    parser.add_argument("--scaffold", help=argparse.SUPPRESS)
    parser.add_argument(
        "--identity",
        default="auto",
        help="Identity bins (auto|all|>=X|label=X..Y|bins=0.90..0.94,0.94..0.98,>=0.98)",
    )
    parser.add_argument(
        "--visual",
        choices=["none", "postcard", "postcard+quantiles"],
        default="none",
        help="Aggregate visual output flavor",
    )
    parser.add_argument(
        "--out",
        choices=["text", "text+svg"],
        default="text",
        dest="visual_output",
        help="Aggregate postcard artifact types",
    )
    parser.add_argument("--top-k", type=int, default=3, help="Top-K motifs/TSDs to report")
    parser.add_argument(
        "--min-n",
        type=int,
        default=20,
        help="Warn when cohorts contain fewer elements than this threshold",
    )
    parser.add_argument(
        "--postcard-ascii-width",
        type=int,
        default=100,
        help="ASCII width for aggregate postcards",
    )
    parser.add_argument(
        "--postcard-svg-width",
        type=int,
        default=800,
        help="SVG width for aggregate postcards",
    )
    parser.add_argument(
        "--postcard-svg-height",
        type=int,
        default=120,
        help="SVG height for aggregate postcards",
    )
    parser.add_argument(
        "--aggregate-tsv",
        help="Single-file mode: path for identity-bin aggregate TSV (defaults to <outdir>/identity_bins.tsv)",
    )
    parser.add_argument(
        "--aggregate-json",
        help="Single-file mode: optional path for identity-bin aggregate JSON output",
    )
    parser.add_argument(
        "--group-aggregates",
        default="genome,scaffold,superfamily",
        help="Comma-separated cohort aggregate groups to write (none|genome,scaffold,superfamily)",
    )
    parser.add_argument(
        "--group-aggregate-tsv",
        help="Single-file mode: path for cohort aggregate TSV (defaults to <outdir>/cohort_aggregates.tsv)",
    )
    parser.add_argument(
        "--group-aggregate-json",
        help="Single-file mode: optional path for cohort aggregate JSON output",
    )
    parser.add_argument(
        "--batch-sample-summary",
        help="Directory mode: path for batch sample summary TSV (defaults to <outdir>/batch_samples.tsv)",
    )
    parser.add_argument(
        "--batch-superfamily-summary",
        help="Directory mode: path for batch per-sample superfamily TSV (defaults to <outdir>/batch_superfamilies.tsv)",
    )
    parser.add_argument(
        "--validation",
        choices=["error", "warn", "off"],
        default="error",
        help="How to handle repeat_region blocks that fail intact-LTR QC",
    )
    parser.add_argument(
        "--substitution-rate",
        type=float,
        help="Substitution rate (subs/site/year) for age estimates from ltr_identity",
    )
    return parser.parse_args()


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _ensure_parent(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def _discover_input_files(input_path: str, limit_files: Optional[int]) -> List[str]:
    path = Path(input_path)
    if path.is_dir():
        files = sorted(str(item) for item in [*path.glob("*.gff3"), *path.glob("*.gff3.gz")])
        if limit_files:
            files = files[:limit_files]
        if not files:
            raise SystemExit(f"No .gff3 or .gff3.gz files found under {input_path}")
        return files
    if not path.exists():
        raise SystemExit(f"Input path does not exist: {input_path}")
    return [str(path)]


def _sample_name(path: str) -> str:
    name = Path(path).name
    for suffix in (".gff3.gz", ".gff3"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


def _load_elements(gff_path: str, max_elements: Optional[int]) -> List[RepeatRegion]:
    regions, by_parent = load_gff(gff_path)
    elements: List[RepeatRegion] = []
    for count, (rid, parent) in enumerate(regions.items(), start=1):
        children = by_parent.get(rid, [])
        elements.append(to_repeat_region_object(rid, parent, children))
        if max_elements and count >= max_elements:
            break
    return elements


def _write_ascii(elem: RepeatRegion, outdir: str, width: int, ruler: bool) -> None:
    ascii_width = width if width <= 200 else 100
    text = ascii_map(elem, width=ascii_width, ruler=ruler)
    path = os.path.join(outdir, f"{elem.id}.txt")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def _write_svg(
    elem: RepeatRegion,
    outdir: str,
    width: int,
    height: int,
    ruler: bool,
    palette: str,
) -> None:
    drawing = svg_map(elem, width=width, height=height, ruler=ruler, palette=palette)
    path = os.path.join(outdir, f"{elem.id}.svg")
    drawing.saveas(path)


def _write_bed(elems: Iterable[RepeatRegion], bed_path: str) -> None:
    _ensure_parent(bed_path)
    with open(bed_path, "w", encoding="utf-8") as handle:
        for elem in elems:
            handle.write(
                f"{elem.scaffold}\t{elem.start - 1}\t{elem.end}\t{elem.id}\t0\t{elem.strand}\n"
            )


def _write_index(elems: Iterable[RepeatRegion], outdir: str, width: int, height: int) -> None:
    html_path = os.path.join(outdir, "index.html")
    with open(html_path, "w", encoding="utf-8") as handle:
        handle.write("<html><head><meta charset='utf-8'><title>GFF3 LTR Maps</title></head><body>\n")
        for elem in elems:
            handle.write(f"<h3>{elem.scaffold}:{elem.id}</h3>\n")
            handle.write(
                f"<object data='{elem.id}.svg' type='image/svg+xml' width='{width}' "
                f"height='{height}'></object>\n"
            )
        handle.write("</body></html>\n")


def _build_identity_reports(
    elements: List[RepeatRegion],
    scope: ScopeSpec,
    bins: List[IdentityBin],
    args: argparse.Namespace,
) -> List[IdentityReport]:
    reports: List[IdentityReport] = []
    scoped_elements = [elem for elem in elements if scope.matches(elem)]
    for bin_spec in bins:
        matched = [elem for elem in scoped_elements if bin_spec.matches(elem)]
        label = f"{scope.label}:{bin_spec.label}"
        aggregate = summarize_cohort(
            label,
            matched,
            substitution_rate=args.substitution_rate,
            top_k=max(1, args.top_k),
            min_n=max(1, args.min_n),
        )
        profile = None
        if args.visual != "none" and matched:
            profile = build_average_profile(
                matched,
                label=bin_spec.label,
                group_label=scope.label,
                identity_range=bin_spec.range_tuple,
                min_n=args.min_n,
            )
        reports.append(
            IdentityReport(
                scope=scope,
                bin=bin_spec,
                elements=matched,
                aggregate=aggregate,
                profile=profile,
            )
        )
    return reports


def _write_identity_summary(
    reports: List[IdentityReport],
    *,
    tsv_path: str,
    json_path: Optional[str],
    emit_stdout: bool,
) -> None:
    if not reports:
        return
    rows = [report.aggregate for report in reports]
    _ensure_parent(tsv_path)
    with open(tsv_path, "w", encoding="utf-8") as handle:
        write_aggregate_tsv(rows, handle)
    if emit_stdout:
        _print_identity_table(rows)
    if json_path:
        _ensure_parent(json_path)
        with open(json_path, "w", encoding="utf-8") as handle:
            write_aggregate_json(rows, handle)


def _write_identity_postcards(reports: List[IdentityReport], args: argparse.Namespace, outdir: str) -> None:
    if args.visual == "none":
        return
    show_quantiles = args.visual == "postcard+quantiles"
    svg_enabled = args.visual_output == "text+svg"
    postcard_dir = os.path.join(outdir, "identity_postcards")
    _ensure_dir(postcard_dir)
    written_files: List[str] = []
    for report in reports:
        profile = report.profile
        if not profile or not profile.has_content:
            continue
        base = f"{report.scope.label}_{report.bin.slug}"
        ascii_path = os.path.join(postcard_dir, f"{base}.txt")
        ascii_txt = average_ascii_map(
            profile,
            width=args.postcard_ascii_width,
            ruler=True,
            show_quantiles=show_quantiles,
        )
        with open(ascii_path, "w", encoding="utf-8") as handle:
            handle.write(ascii_txt)
        written_files.append(f"{base}.txt")
        if svg_enabled:
            drawing = average_svg_map(
                profile,
                width=args.postcard_svg_width,
                height=args.postcard_svg_height,
                palette=args.palette,
                show_quantiles=show_quantiles,
            )
            svg_path = os.path.join(postcard_dir, f"{base}.svg")
            drawing.saveas(svg_path)
            written_files.append(f"{base}.svg")
    if written_files:
        exts = sorted({path.split(".")[-1] for path in written_files})
        print(
            f"[gff3-ltr-map] Identity postcards saved under {postcard_dir} "
            f"({', '.join(exts)} files: {len(written_files)})"
        )


def _resolve_scope_arg(args: argparse.Namespace) -> Optional[str]:
    for key in ("chrom", "scope", "scaffold"):
        value = getattr(args, key, None)
        if value:
            return value
    return None


def _scope_scaffold_name(scope_arg: Optional[str]) -> Optional[str]:
    raw = (scope_arg or "").strip()
    if not raw:
        return None
    lowered = raw.lower()
    if lowered in {"default", "auto", "all", "genome"}:
        return None
    if lowered.startswith("chrom:"):
        name = raw.split(":", 1)[1].strip()
        if not name or name.lower() in {"all", "genome"}:
            return None
        return name
    return raw


def _parse_interval(raw: str) -> Tuple[int, int]:
    cleaned = raw.strip().replace(",", "")
    for delimiter in ("..", "-"):
        if delimiter in cleaned:
            left, right = cleaned.split(delimiter, 1)
            try:
                start = int(left)
                end = int(right)
            except ValueError as exc:
                raise ValueError(f"invalid region coordinates '{raw}'") from exc
            if start <= 0 or end <= 0 or end < start:
                raise ValueError(f"invalid region coordinates '{raw}'")
            return start, end
    raise ValueError(f"invalid region coordinates '{raw}'")


def _infer_primary_scaffold(
    elements: Sequence[RepeatRegion],
    *,
    window_start: int,
    window_end: int,
) -> Optional[str]:
    if not elements:
        return None
    scaffold_stats: dict[str, tuple[int, int, int]] = {}
    for elem in elements:
        overlap_n, best_end, count = scaffold_stats.get(elem.scaffold, (0, 0, 0))
        overlaps_window = elem.end >= window_start and elem.start <= window_end
        scaffold_stats[elem.scaffold] = (
            overlap_n + (1 if overlaps_window else 0),
            max(best_end, elem.end),
            count + 1,
        )
    ranked = sorted(
        scaffold_stats.items(),
        key=lambda item: (-item[1][0], -item[1][1], -item[1][2], item[0]),
    )
    return ranked[0][0] if ranked else None


def _resolve_region_filter(
    region_arg: Optional[str],
    scope_arg: Optional[str],
    elements: Sequence[RepeatRegion],
) -> Tuple[Optional[RegionFilter], Optional[str]]:
    raw = (region_arg or "").strip()
    if not raw:
        return None, None
    scaffold: Optional[str]
    interval_text: str
    note: Optional[str] = None
    if ":" in raw:
        scaffold, interval_text = raw.split(":", 1)
        scaffold = scaffold.strip()
        if not scaffold:
            raise ValueError("region scaffold is empty")
    else:
        start, end = _parse_interval(raw)
        scaffold = _scope_scaffold_name(scope_arg)
        if not scaffold:
            scaffold = _infer_primary_scaffold(elements, window_start=start, window_end=end)
            if not scaffold:
                return (
                    RegionFilter(scaffold="NA", start=start, end=end),
                    "no repeat_region features; primary scaffold unavailable",
                )
            note = f"primary scaffold inferred as {scaffold}"
        return RegionFilter(scaffold=scaffold, start=start, end=end), note
    start, end = _parse_interval(interval_text)
    return RegionFilter(scaffold=scaffold, start=start, end=end), note


def _analysis_region_label(scope_arg: Optional[str], region_filter: Optional[RegionFilter]) -> str:
    if region_filter:
        return region_filter.label
    scaffold = _scope_scaffold_name(scope_arg)
    if scaffold:
        return scaffold
    return "genome"


def _window_bp(region_filter: Optional[RegionFilter]) -> Optional[int]:
    if not region_filter:
        return None
    return region_filter.end - region_filter.start + 1


def _covered_bp_within_window(elements: Sequence[RepeatRegion], region_filter: Optional[RegionFilter]) -> Optional[int]:
    if not region_filter:
        return None
    intervals: List[Tuple[int, int]] = []
    for elem in elements:
        if elem.scaffold != region_filter.scaffold:
            continue
        start = max(elem.start, region_filter.start)
        end = min(elem.end, region_filter.end)
        if start <= end:
            intervals.append((start, end))
    if not intervals:
        return 0
    intervals.sort()
    covered = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
            continue
        covered += cur_end - cur_start + 1
        cur_start, cur_end = start, end
    covered += cur_end - cur_start + 1
    return covered


def _build_scopes(scope_arg: Optional[str], elements: List[RepeatRegion]) -> List[ScopeSpec]:
    scaffolds = sorted({elem.scaffold for elem in elements})
    raw = (scope_arg or "").strip()
    if raw == "" or raw.lower() in ("default", "auto"):
        scopes = [ScopeSpec("genome", None)]
        scopes.extend(ScopeSpec("chrom", name) for name in scaffolds)
        return scopes
    if raw.lower() == "genome":
        return [ScopeSpec("genome", None)]
    if raw.lower() == "all":
        return [ScopeSpec("genome", None)] + [ScopeSpec("chrom", name) for name in scaffolds]
    if raw.lower().startswith("chrom:"):
        name = raw.split(":", 1)[1].strip()
        if not name:
            raise ValueError("chrom scope requires a name (e.g., chrom:chr_2)")
        if name.lower() == "all":
            return _build_scopes("all", elements)
        return [ScopeSpec("chrom", name)]
    return [ScopeSpec("chrom", raw)]


def _parse_identity_spec(spec: str) -> List[IdentityBin]:
    spec = (spec or "auto").strip()
    lowered = spec.lower()
    if lowered in ("auto", "default"):
        return [IdentityBin("all", None, None), IdentityBin(">=0.980", 0.98, None)]
    if lowered == "all":
        return [IdentityBin("all", None, None)]
    bins: List[IdentityBin] = []
    token_source = spec
    if lowered.startswith("bins="):
        token_source = spec.split("=", 1)[1]
    tokens = [chunk.strip() for chunk in token_source.split(",") if chunk.strip()]
    for token in tokens:
        if "=" in token and not token.strip().startswith(">="):
            name, range_str = token.split("=", 1)
            bins.append(_bin_from_range(range_str.strip(), name.strip() or None))
        else:
            bins.append(_bin_from_range(token, None))
    if not bins:
        raise ValueError("no valid identity bins parsed")
    return bins


def _bin_from_range(range_str: str, custom_label: Optional[str]) -> IdentityBin:
    if not range_str:
        raise ValueError("empty identity range")
    text = range_str.strip()
    min_id: Optional[float] = None
    max_id: Optional[float] = None
    if text.startswith(">="):
        try:
            min_id = float(text[2:])
        except ValueError as exc:
            raise ValueError(f"invalid identity threshold '{range_str}'") from exc
    else:
        for delim in ("..", "-", "–"):
            if delim in text:
                left, right = text.split(delim, 1)
                try:
                    min_id = float(left)
                    max_id = float(right)
                except ValueError as exc:
                    raise ValueError(f"invalid identity range '{range_str}'") from exc
                if max_id < min_id:
                    raise ValueError(f"identity range '{range_str}' has max < min")
                break
        else:
            raise ValueError(f"unsupported identity descriptor '{range_str}'")
    label = custom_label or _format_identity_label(min_id, max_id)
    return IdentityBin(label=label, min_identity=min_id, max_identity=max_id)


def _parse_group_aggregates(spec: str) -> List[str]:
    raw = (spec or "").strip()
    if not raw or raw.lower() == "none":
        return []
    items = [item.strip().lower() for item in raw.split(",") if item.strip()]
    allowed = {"genome", "scaffold", "superfamily"}
    invalid = [item for item in items if item not in allowed]
    if invalid:
        raise ValueError(f"unsupported group aggregate(s): {', '.join(invalid)}")
    ordered: List[str] = []
    for item in items:
        if item not in ordered:
            ordered.append(item)
    return ordered


def _format_identity_label(min_id: Optional[float], max_id: Optional[float]) -> str:
    if min_id is not None and max_id is None:
        return f">={min_id:.3f}"
    if min_id is not None and max_id is not None:
        return f"{min_id:.3f}-{max_id:.3f}"
    return "all"


def _slugify(text: str) -> str:
    slug = "".join(ch.lower() if ch.isalnum() else "_" for ch in text)
    slug = slug.strip("_")
    return slug or "bin"


def _print_identity_table(rows: List[AggregateRow]) -> None:
    if not rows:
        return
    headers = ["group", "n", "median_len", "median_id", "median_tsd", "top_motifs", "notes"]
    print("[gff3-ltr-map] Identity bins:")
    print("\t".join(headers))
    for row in rows:
        line = [
            row.group,
            str(row.n_elements),
            _safe_float(row.length_bp_median),
            _safe_float(row.ltr_identity_median),
            _safe_float(row.tsd_len_median),
            row.top_motifs or "NA",
            row.notes or "",
        ]
        print("\t".join(line))


def _safe_float(value: Optional[float]) -> str:
    if value is None:
        return "NA"
    if abs(value) >= 1:
        return f"{value:.1f}"
    return f"{value:.3f}"


def _scientist_element_path(outdir: str, primary_tsv: str) -> str:
    detail_path = os.path.join(outdir, "scientist_elements.tsv")
    if os.path.abspath(detail_path) == os.path.abspath(primary_tsv):
        return os.path.join(outdir, "scientist_elements_detail.tsv")
    return detail_path


def _single_sample_paths(args: argparse.Namespace) -> SamplePaths:
    scientist_tsv = args.scientist_tsv or os.path.join(args.outdir, "scientist_regions.tsv")
    return SamplePaths(
        outdir=args.outdir,
        summary=args.summary,
        validation=args.validation_report or os.path.join(args.outdir, "validation.tsv"),
        identity_tsv=args.aggregate_tsv or os.path.join(args.outdir, "identity_bins.tsv"),
        identity_json=args.aggregate_json,
        group_tsv=args.group_aggregate_tsv or os.path.join(args.outdir, "cohort_aggregates.tsv"),
        group_json=args.group_aggregate_json,
        scientist_tsv=scientist_tsv,
        scientist_elements_tsv=_scientist_element_path(args.outdir, scientist_tsv),
        bed=args.bed,
    )


def _batch_sample_paths(root_outdir: str, sample_name: str, bed_enabled: bool) -> SamplePaths:
    sample_outdir = os.path.join(root_outdir, "samples", sample_name)
    scientist_tsv = os.path.join(sample_outdir, "scientist_regions.tsv")
    return SamplePaths(
        outdir=sample_outdir,
        summary=os.path.join(sample_outdir, "summary.tsv"),
        validation=os.path.join(sample_outdir, "validation.tsv"),
        identity_tsv=os.path.join(sample_outdir, "identity_bins.tsv"),
        identity_json=None,
        group_tsv=os.path.join(sample_outdir, "cohort_aggregates.tsv"),
        group_json=None,
        scientist_tsv=scientist_tsv,
        scientist_elements_tsv=_scientist_element_path(sample_outdir, scientist_tsv),
        bed=os.path.join(sample_outdir, "elements.bed") if bed_enabled else None,
    )


def _write_group_aggregates(
    rows: Sequence[AggregateRow],
    *,
    tsv_path: Optional[str],
    json_path: Optional[str],
) -> None:
    if not rows or not tsv_path:
        return
    _ensure_parent(tsv_path)
    with open(tsv_path, "w", encoding="utf-8") as handle:
        write_aggregate_tsv(rows, handle)
    if json_path:
        _ensure_parent(json_path)
        with open(json_path, "w", encoding="utf-8") as handle:
            write_aggregate_json(rows, handle)


def _render_elements(elements: Sequence[RepeatRegion], paths: SamplePaths, args: argparse.Namespace) -> None:
    if not (args.ascii or args.svg):
        return
    workers = max(1, args.workers)
    if workers > 1:
        futures = []
        with ProcessPoolExecutor(max_workers=workers) as executor:
            for elem in elements:
                if args.ascii:
                    futures.append(executor.submit(_write_ascii, elem, paths.outdir, args.width, args.ruler))
                if args.svg:
                    futures.append(
                        executor.submit(
                            _write_svg,
                            elem,
                            paths.outdir,
                            args.width,
                            args.height,
                            args.ruler,
                            args.palette,
                        )
                    )
            for future in as_completed(futures):
                future.result()
    else:
        for elem in elements:
            if args.ascii:
                _write_ascii(elem, paths.outdir, args.width, args.ruler)
            if args.svg:
                _write_svg(elem, paths.outdir, args.width, args.height, args.ruler, args.palette)


def _run_sample(
    input_path: str,
    args: argparse.Namespace,
    *,
    sample_name: str,
    paths: SamplePaths,
    emit_identity_stdout: bool,
) -> Tuple[List[RepeatRegion], List[RepeatRegion], List[IdentityReport], ScientistRegionRow]:
    _ensure_dir(paths.outdir)
    print(f"[gff3-ltr-map] Reading elements from {input_path} …")
    all_elements = _load_elements(input_path, args.max_elements)
    scope_arg = _resolve_scope_arg(args)
    try:
        region_filter, region_note = _resolve_region_filter(args.region, scope_arg, all_elements)
    except ValueError as exc:
        raise SystemExit(f"--region {exc}") from exc
    region_row_note: Optional[str] = None
    if region_note:
        print(f"[gff3-ltr-map] {region_note}")
        if "unavailable" in region_note:
            region_row_note = region_note
    if region_filter:
        all_elements = [elem for elem in all_elements if region_filter.matches(elem)]
        print(
            f"[gff3-ltr-map] Region filter {region_filter.label} retained "
            f"{len(all_elements)} overlapping repeat_region features"
        )

    intact_elements = [elem for elem in all_elements if elem.is_intact]
    failed_elements = [elem for elem in all_elements if elem.qc_status == "FAIL"]
    warned_elements = [elem for elem in all_elements if elem.qc_status == "WARN"]

    print(
        f"[gff3-ltr-map] Loaded {len(all_elements)} repeat_region features "
        f"({len(intact_elements)} intact, {len(warned_elements)} warn, {len(failed_elements)} fail)"
    )

    _ensure_parent(paths.summary)
    with open(paths.summary, "w", encoding="utf-8") as handle:
        write_summary(all_elements, handle)
    _ensure_parent(paths.validation)
    with open(paths.validation, "w", encoding="utf-8") as handle:
        write_validation_report(all_elements, handle)
    print(f"[gff3-ltr-map] Summary TSV written to {paths.summary}")
    print(f"[gff3-ltr-map] Validation TSV written to {paths.validation}")

    if failed_elements:
        message = (
            f"[gff3-ltr-map] {len(failed_elements)} repeat_region blocks failed intact-LTR QC "
            f"(see {paths.validation})"
        )
        if args.validation == "error":
            raise SystemExit(message)
        if args.validation == "warn":
            print(message)

    if paths.bed:
        _write_bed(intact_elements, paths.bed)
        print(f"[gff3-ltr-map] BED written to {paths.bed}")

    _render_elements(intact_elements, paths, args)
    if args.ascii or args.svg:
        print(
            f"[gff3-ltr-map] Per-element postcard rendering completed "
            f"(ascii={args.ascii}, svg={args.svg})"
        )

    if args.index_html:
        if args.svg:
            _write_index(intact_elements, paths.outdir, args.width, args.height)
        else:
            print("--index-html requested but --svg not enabled; skipping index generation")

    try:
        group_types = _parse_group_aggregates(args.group_aggregates)
    except ValueError as exc:
        raise SystemExit(f"--group-aggregates {exc}") from exc
    group_rows = compute_aggregates(
        intact_elements,
        group_types,
        substitution_rate=args.substitution_rate,
        top_k=max(1, args.top_k),
        min_n=max(1, args.min_n),
    )
    _write_group_aggregates(group_rows, tsv_path=paths.group_tsv, json_path=paths.group_json)
    if group_rows and paths.group_tsv:
        print(f"[gff3-ltr-map] Cohort aggregates written to {paths.group_tsv}")

    try:
        identity_bins = _parse_identity_spec(args.identity)
    except ValueError as exc:
        raise SystemExit(f"--identity {exc}") from exc
    effective_scope_arg = region_filter.scaffold if region_filter else scope_arg
    try:
        scopes = _build_scopes(effective_scope_arg, intact_elements)
    except ValueError as exc:
        raise SystemExit(f"--chrom {exc}") from exc

    reports: List[IdentityReport] = []
    for scope in scopes:
        reports.extend(_build_identity_reports(intact_elements, scope, identity_bins, args))
    _write_identity_summary(
        reports,
        tsv_path=paths.identity_tsv,
        json_path=paths.identity_json,
        emit_stdout=emit_identity_stdout,
    )
    _write_identity_postcards(reports, args, paths.outdir)
    window_bp = _window_bp(region_filter)
    covered_bp = _covered_bp_within_window(intact_elements, region_filter)
    scientist_region_row = build_scientist_region_row(
        sample_name,
        _analysis_region_label(scope_arg, region_filter),
        all_elements,
        intact_elements,
        window_bp=window_bp,
        coverage_pct=((covered_bp / window_bp) * 100.0) if window_bp and covered_bp is not None else None,
        extra_note=region_row_note,
    )
    _ensure_parent(paths.scientist_tsv)
    with open(paths.scientist_tsv, "w", encoding="utf-8") as handle:
        write_scientist_region_tsv([scientist_region_row], handle)
    print(f"[gff3-ltr-map] Scientist region TSV written to {paths.scientist_tsv}")

    scientist_rows = build_scientist_element_rows(intact_elements)
    _ensure_parent(paths.scientist_elements_tsv)
    with open(paths.scientist_elements_tsv, "w", encoding="utf-8") as handle:
        write_scientist_element_tsv(scientist_rows, handle)
    print(f"[gff3-ltr-map] Scientist element TSV written to {paths.scientist_elements_tsv}")
    if emit_identity_stdout:
        print_scientist_region_table([scientist_region_row], max_rows=max(1, args.scientist_cli_max_rows))
    print(f"[gff3-ltr-map] Identity bins processed ({len(reports)} cohorts)")
    return all_elements, intact_elements, reports, scientist_region_row


def _run_batch(args: argparse.Namespace, input_files: Sequence[str]) -> None:
    sample_rows: List[BatchSampleRow] = []
    scientist_batch_rows: List[ScientistRegionRow] = []
    superfamily_rows: List[BatchSuperfamilyRow] = []
    for input_path in input_files:
        sample_name = _sample_name(input_path)
        paths = _batch_sample_paths(args.outdir, sample_name, bed_enabled=bool(args.bed))
        all_elements, intact_elements, _reports, scientist_region_row = _run_sample(
            input_path,
            args,
            sample_name=sample_name,
            paths=paths,
            emit_identity_stdout=False,
        )
        sample_rows.append(
            summarize_sample(
                sample_name,
                input_path,
                all_elements,
                intact_elements,
                top_k=max(1, args.top_k),
            )
        )
        scientist_batch_rows.append(scientist_region_row)
        superfamily_rows.extend(
            summarize_sample_superfamilies(
                sample_name,
                intact_elements,
                top_k=max(1, args.top_k),
            )
        )
        print(f"[gff3-ltr-map] Batch sample completed: {sample_name}")

    sample_summary_path = args.batch_sample_summary or os.path.join(args.outdir, "batch_samples.tsv")
    superfamily_summary_path = args.batch_superfamily_summary or os.path.join(
        args.outdir, "batch_superfamilies.tsv"
    )
    scientist_batch_path = args.scientist_tsv or os.path.join(args.outdir, "scientist_batch.tsv")
    scientist_batch_rows_sorted = sorted(
        scientist_batch_rows,
        key=lambda row: (-row.intact_elements, -(row.coverage_pct or 0.0), row.sample),
    )
    _ensure_parent(sample_summary_path)
    with open(sample_summary_path, "w", encoding="utf-8") as handle:
        write_batch_sample_tsv(sample_rows, handle)
    _ensure_parent(superfamily_summary_path)
    with open(superfamily_summary_path, "w", encoding="utf-8") as handle:
        write_batch_superfamily_tsv(superfamily_rows, handle)
    _ensure_parent(scientist_batch_path)
    with open(scientist_batch_path, "w", encoding="utf-8") as handle:
        write_scientist_region_tsv(scientist_batch_rows_sorted, handle)
    print(f"[gff3-ltr-map] Batch sample summary written to {sample_summary_path}")
    print(f"[gff3-ltr-map] Batch superfamily summary written to {superfamily_summary_path}")
    print(f"[gff3-ltr-map] Scientist batch TSV written to {scientist_batch_path}")
    print_scientist_region_table(
        scientist_batch_rows_sorted,
        max_rows=max(25, args.scientist_cli_max_rows),
    )


def main() -> None:
    args = _parse_args()
    _ensure_dir(args.outdir)

    input_files = _discover_input_files(args.input_path, args.limit_files)
    if len(input_files) == 1 and not Path(args.input_path).is_dir():
        paths = _single_sample_paths(args)
        _run_sample(
            input_files[0],
            args,
            sample_name=_sample_name(input_files[0]),
            paths=paths,
            emit_identity_stdout=True,
        )
        print("[gff3-ltr-map] Done.")
        return

    _run_batch(args, input_files)
    print("[gff3-ltr-map] Done.")


if __name__ == "__main__":  # pragma: no cover
    main()
