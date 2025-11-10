"""Console entry-point for gff3-ltr-map."""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, Iterable, List, Optional, Tuple

from .aggregates import AggregateRow, summarize_cohort, write_aggregate_json, write_aggregate_tsv
from .average_map import average_ascii_map, average_svg_map, build_average_profile
from .model import RepeatRegion
from .parser import load_gff, to_repeat_region_object
from .render_ascii import ascii_map
from .render_svg import svg_map
from .summary import write_summary


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
    aggregate: "AggregateRow"
    profile: Optional["AverageMapProfile"]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="EDTA intact GFF3 -> ASCII/SVG per-element maps + TSV summary",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("gff3", help="Input EDTA intact GFF3")
    parser.add_argument("--outdir", default="out", help="Output directory for map artifacts")
    parser.add_argument("--ascii", action="store_true", help="Emit ASCII maps (*.txt)")
    parser.add_argument("--svg", action="store_true", help="Emit SVG maps (*.svg)")
    parser.add_argument("--summary", default="summary.tsv", help="Path to summary TSV")
    parser.add_argument("--width", type=int, default=800, help="SVG canvas width (px) or ASCII columns")
    parser.add_argument("--height", type=int, default=80, help="SVG canvas height (px)")
    parser.add_argument(
        "--palette",
        choices=["classic", "mono", "protanopia"],
        default="classic",
        help="Color palette for SVG maps",
    )
    parser.add_argument("--bed", help="Also emit BED with repeat_region spans")
    parser.add_argument("--index-html", action="store_true", help="Generate HTML index embedding SVGs")
    parser.add_argument("--ruler", action="store_true", help="Add coordinate ruler to ASCII/SVG outputs")
    parser.add_argument("--workers", type=int, default=1, help="Parallel workers for rendering")
    parser.add_argument("--max-elements", type=int, help="Stop after N elements (debug)")
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
        help="Path for identity-bin aggregate TSV (defaults to <outdir>/identity_bins.tsv)",
    )
    parser.add_argument(
        "--aggregate-json",
        help="Optional path for aggregate statistics JSON output",
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


def _iter_elements(args: argparse.Namespace) -> Iterable[RepeatRegion]:
    regions, by_parent = load_gff(args.gff3)
    count = 0
    for rid, parent in regions.items():
        children = by_parent.get(rid, [])
        elem = to_repeat_region_object(rid, parent, children)
        yield elem
        count += 1
        if args.max_elements and count >= args.max_elements:
            break


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


def _write_index(elems: Iterable[RepeatRegion], args: argparse.Namespace) -> None:
    html_path = os.path.join(args.outdir, "index.html")
    with open(html_path, "w", encoding="utf-8") as handle:
        handle.write("<html><head><meta charset='utf-8'><title>GFF3 LTR Maps</title></head><body>\n")
        for elem in elems:
            handle.write(f"<h3>{elem.scaffold}:{elem.id}</h3>\n")
            handle.write(
                f"<object data='{elem.id}.svg' type='image/svg+xml' width='{args.width}' "
                f"height='{args.height}'></object>\n"
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


def _write_identity_summary(reports: List[IdentityReport], args: argparse.Namespace) -> None:
    if not reports:
        return
    rows = [report.aggregate for report in reports]
    tsv_path = args.aggregate_tsv or os.path.join(args.outdir, "identity_bins.tsv")
    _ensure_parent(tsv_path)
    with open(tsv_path, "w", encoding="utf-8") as handle:
        write_aggregate_tsv(rows, handle)
    _print_identity_table(rows)
    if args.aggregate_json:
        _ensure_parent(args.aggregate_json)
        with open(args.aggregate_json, "w", encoding="utf-8") as handle:
            write_aggregate_json(rows, handle)


def _write_identity_postcards(reports: List[IdentityReport], args: argparse.Namespace) -> None:
    if args.visual == "none":
        return
    show_quantiles = args.visual == "postcard+quantiles"
    svg_enabled = args.visual_output == "text+svg"
    outdir = os.path.join(args.outdir, "identity_postcards")
    _ensure_dir(outdir)
    written_files: List[str] = []
    for report in reports:
        profile = report.profile
        if not profile or not profile.has_content:
            continue
        base = f"{report.scope.label}_{report.bin.slug}"
        ascii_path = os.path.join(outdir, f"{base}.txt")
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
            svg_path = os.path.join(outdir, f"{base}.svg")
            drawing.saveas(svg_path)
            written_files.append(f"{base}.svg")
    if written_files:
        exts = sorted({path.split(".")[-1] for path in written_files})
        print(
            f"[gff3-ltr-map] Identity postcards saved under {outdir} "
            f"({', '.join(exts)} files: {len(written_files)})"
        )


def _resolve_scope_arg(args: argparse.Namespace) -> Optional[str]:
    for key in ("chrom", "scope", "scaffold"):
        value = getattr(args, key, None)
        if value:
            return value
    return None


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
    headers = ["group", "n", "median_len", "median_id", "top_motifs", "notes"]
    print("[gff3-ltr-map] Identity bins:")
    print("\t".join(headers))
    for row in rows:
        line = [
            row.group,
            str(row.n_elements),
            _safe_float(row.length_bp_median),
            _safe_float(row.ltr_identity_median),
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


def main() -> None:
    args = _parse_args()

    outdir = args.outdir or "runs"
    args.outdir = outdir
    _ensure_dir(outdir)

    try:
        identity_bins = _parse_identity_spec(args.identity)
    except ValueError as exc:
        raise SystemExit(f"--identity {exc}") from exc

    print(f"[gff3-ltr-map] Reading elements from {args.gff3} …")
    elements: List[RepeatRegion] = list(_iter_elements(args))
    print(f"[gff3-ltr-map] Loaded {len(elements)} repeat_region features")

    scope_arg = _resolve_scope_arg(args)
    try:
        scopes = _build_scopes(scope_arg, elements)
    except ValueError as exc:
        raise SystemExit(f"--chrom {exc}") from exc

    _ensure_parent(args.summary)
    with open(args.summary, "w", encoding="utf-8") as summary_handle:
        write_summary(elements, summary_handle)
    print(f"[gff3-ltr-map] Summary TSV written to {args.summary}")

    if args.bed:
        _write_bed(elements, args.bed)
        print(f"[gff3-ltr-map] BED written to {args.bed}")

    if args.ascii or args.svg:
        workers = max(1, args.workers)
        if workers > 1:
            futures = []
            with ProcessPoolExecutor(max_workers=workers) as executor:
                for elem in elements:
                    if args.ascii:
                        futures.append(
                            executor.submit(_write_ascii, elem, args.outdir, args.width, args.ruler)
                        )
                    if args.svg:
                        futures.append(
                            executor.submit(
                                _write_svg,
                                elem,
                                args.outdir,
                                args.width,
                                args.height,
                                args.ruler,
                                args.palette,
                            )
                        )
                for _ in as_completed(futures):
                    pass
        else:
            for elem in elements:
                if args.ascii:
                    _write_ascii(elem, args.outdir, args.width, args.ruler)
                if args.svg:
                    _write_svg(elem, args.outdir, args.width, args.height, args.ruler, args.palette)
        print(f"[gff3-ltr-map] Per-element postcard rendering completed (ascii={args.ascii}, svg={args.svg})")

    if args.index_html:
        if args.svg:
            _write_index(elements, args)
        else:
            print("--index-html requested but --svg not enabled; skipping index generation")

    reports: List[IdentityReport] = []
    for scope in scopes:
        reports.extend(_build_identity_reports(elements, scope, identity_bins, args))
    _write_identity_summary(reports, args)
    _write_identity_postcards(reports, args)
    print(f"[gff3-ltr-map] Identity bins processed ({len(reports)} cohorts)")
    print("[gff3-ltr-map] Done.")


if __name__ == "__main__":  # pragma: no cover
    main()
