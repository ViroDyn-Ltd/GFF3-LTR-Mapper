# GFF3-LTR-Mapper

A command-line workflow for summarising EDTA "intact" GFF3 annotations into scientist-friendly tables and postcards (ASCII + SVG) per genome, chromosome, or custom identity cohort.

## Obtain & Run
1. Download the GitHub ZIP archive (or clone the repository) and unzip it anywhere.
2. Ensure Python 3.9+ is available.
3. From the project directory, invoke the CLI via `python -m gff3_ltr_map.cli`:
   ```bash
   python -m gff3_ltr_map.cli /path/to/EDTA.intact.gff3 \
     --outdir runs/default_sweep \
     --visual postcard+quantiles \
     --out text+svg \
     --summary runs/default_sweep/summary.tsv \
     --top-k 5 \
     --min-n 25
   ```
   The default behaviour emits genome-wide aggregates **and** one cohort per chromosome. The CLI prints tables directly to the terminal and reports exactly where postcard files land.

### Alternative scopes
- `--chrom chr_5` : only chr_5.
- `--chrom genome` : genome-wide only.
- `--chrom all` : explicit genome + every chromosome (same as default).

### Identity bins
- `--identity auto` (default) : `all` + `>=0.98`.
- `--identity all` : single cohort.
- `--identity '>=0.97'` : threshold.
- `--identity 'bins=0.90..0.94,>=0.94'` : custom bins.

### Visual controls
- `--visual none|postcard|postcard+quantiles`
- `--out text|text+svg` (postcard files still written silently; CLI only logs their directory and extensions).

Full command cookbook lives in `docs/USAGE.md`, and example outputs (tables, ASCII postcards, SVG preview) sit in `docs/EXAMPLES.md`.

## Outputs
- `summary.tsv` : one row per `repeat_region` with raw attributes.
- `identity_bins.tsv` (path configurable via `--aggregate-tsv`): per-cohort stats (counts, medians, IQR, motif/TSD leaders, consensus warnings, low-N notes).
- `identity_postcards/*.txt|*.svg` : aggregate postcards per scope+identity bin, including the metrics table and annotated bar (LTR labels + Q25/Q75 markers).
- Optional `elements.bed` when `--bed` is provided.
- Legacy per-element postcards remain available via `--ascii/--svg`.

## CLI logging
Every run prints:
1. Input path and element count.
2. Summary/BED destinations.
3. Identity-bin table (tab-separated) for quick inspection.
4. Postcard directory + file-type summary.
5. Completion status.

## Scientific safeguards
- **Quartile markers (Q25/Q75)** visualise the interquartile range of element lengths for each cohort.
- **Consensus warnings** flag motifs/TSDs whose top string occurs &lt;40% of the time.
- **Low-N warnings** appear whenever `n &lt; --min-n` (default 20).
- **Age estimates** stay blank unless a lineage-specific `--substitution-rate` is supplied, and should be treated as rough proxies only.

## Documentation
- `docs/USAGE.md` – flags, scenarios, and CLI behaviour.
- `docs/EXAMPLES.md` – sample tables, ASCII postcard, and an actual SVG postcard screenshot.

## License
MIT (see `LICENSE`).
