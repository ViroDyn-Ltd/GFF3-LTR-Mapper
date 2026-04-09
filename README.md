# GFF3-LTR-Mapper

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19489691.svg)](https://doi.org/10.5281/zenodo.19489691)

A command-line workflow for auditing and summarising EDTA intact LTR GFF3 annotations into QC reports, cohort tables, cross-sample comparison TSVs, and optional ASCII postcards. The parser is strict about intact-LTR structure, keeps non-destructive QC traces for every `repeat_region`, and can run either on a single GFF3 or an entire directory of EDTA outputs. Distributed under the [MIT License](LICENSE).

## Scope
- Input is expected to be **EDTA intact-element GFF3** output only.
- Aggregates and scientist-facing tables are restricted to **intact LTR elements**.
- Cross-species `--region` comparisons are **coordinate-window comparisons**, not homology- or synteny-aware locus comparisons.

## Install
1. Clone the repository or download the ZIP archive.
2. Ensure Python 3.9+ is available.
3. Install the package from the project root:
   ```bash
   pip install .
   ```
4. You can then run either `gff3-ltr-map ...` or `python -m gff3_ltr_map.cli ...`.

## Point The CLI At Data

### Single EDTA intact GFF3
```bash
gff3-ltr-map /absolute/path/sample.intact.gff3 \
  --outdir runs/sample_qc \
  --validation error \
  --visual none
```

### Directory of EDTA intact GFF3 files
```bash
gff3-ltr-map /absolute/path/edta_gff3_directory \
  --outdir runs/batch_qc \
  --chrom genome \
  --identity all \
  --visual none
```

### Cross-species comparison of the same coordinate window
```bash
gff3-ltr-map /absolute/path/edta_gff3_directory \
  --outdir runs/cross_species_window \
  --region 1-3000000 \
  --identity all \
  --visual none
```
This writes per-sample reports under `runs/cross_species_window/samples/<sample>/` plus a root `scientist_batch.tsv` sorted by strongest window signal.

### Single scaffold or interval within one file
```bash
gff3-ltr-map /absolute/path/sample.intact.gff3 \
  --outdir runs/chr2_window \
  --chrom chr_2 \
  --region 1-3000000 \
  --identity all \
  --visual none
```

### Alternative scopes
- `--chrom chr_5` : only chr_5.
- `--chrom genome` : genome-wide only.
- `--chrom all` : explicit genome + every chromosome (same as default).
- `--region chr_5:1-3000000` : restrict outputs to elements overlapping that interval.
- `--chrom chr_5 --region 1-3000000` : same interval syntax when the scaffold is already named.
- `--region 1-3000000` in directory mode : compare the same window on each sample's inferred primary scaffold.

### Identity bins
- `--identity auto` (default) : `all` + `>=0.98`.
- `--identity all` : single cohort.
- `--identity '>=0.97'` : threshold.
- `--identity 'bins=0.90..0.94,>=0.94'` : custom bins.

### QC and cohort reports
- `--validation error|warn|off` controls how hard the CLI should react when a `repeat_region` fails intact-LTR QC.
- `validation.tsv` records per-element pass/warn/fail status plus structural signatures.
- `cohort_aggregates.tsv` summarizes intact elements by `genome`, `scaffold`, and `superfamily` by default.
- `--group-aggregates none|genome,scaffold,superfamily` changes those cohort tables.
- `scientist_regions.tsv` gives a short one-row summary for the selected genome / scaffold / region.
- `scientist_elements.tsv` keeps the per-element detail table for auditability, with the EDTA `repeat_region` ID moved to the last column.
- The CLI prints the simplified region summary directly to the terminal.
- Region summaries now include `window_bp`, `intact_per_mbp`, and `coverage_pct` so cross-species window comparisons are directly sortable.

### ASCII postcard controls
- `--visual none|postcard|postcard+quantiles`
- `--ascii` emits one plain-text map per intact element.

Full command cookbook lives in `docs/USAGE.md`, and example outputs sit in `docs/EXAMPLES.md`.

## Outputs
- `summary.tsv` : one row per `repeat_region` with normalized attributes, retrotransposon type, QC status, and intact/non-intact state.
- `validation.tsv` : compact per-element QC report with warning/error codes and child signatures.
- `scientist_regions.tsv` : simplified one-row summary of the current genome / scaffold / region selection.
- `scientist_elements.tsv` : simplified intact-element detail table, kept as a secondary audit view.
- `cohort_aggregates.tsv` : intact-element summaries across genome/scaffold/superfamily cohorts.
- `identity_bins.tsv` (path configurable via `--aggregate-tsv`): per-cohort stats (counts, medians, TSD lengths, motif/TSD leaders, consensus warnings, low-N notes).
- `identity_postcards/*.txt` : aggregate ASCII postcards per scope+identity bin, including the metrics table and annotated bar (LTR labels + Q25/Q75 markers).
- Optional `elements.bed` when `--bed` is provided.
- Legacy per-element ASCII postcards remain available via `--ascii`.
- In directory mode, each sample gets its own nested output folder with `scientist_regions.tsv` and `scientist_elements.tsv`, and the run root also contains `batch_samples.tsv`, `batch_superfamilies.tsv`, and `scientist_batch.tsv`.

## CLI logging
Every run prints:
1. Input path and total/intact/warn/fail counts.
2. Summary and validation report destinations.
3. Cohort report destinations.
4. Identity-bin table (single-file mode) and postcard directory summary.
5. Completion status.

## Scientific safeguards
- **Quartile markers (Q25/Q75)** visualise the interquartile range of element lengths for each cohort.
- **Consensus warnings** flag motifs/TSDs whose top string occurs &lt;40% of the time.
- **Low-N warnings** appear whenever `n &lt; --min-n` (default 20).
- **QC traceability** keeps every `repeat_region` in `summary.tsv` and `validation.tsv`, even if it is excluded from intact-only aggregates.
- **Age estimates** stay blank unless a lineage-specific `--substitution-rate` is supplied, and should be treated as rough proxies only.

## Documentation
- `docs/USAGE.md` – flags, scenarios, and CLI behaviour.
- `docs/EXAMPLES.md` – sample tables and ASCII postcard examples.
- `CITATION.cff` – citation metadata for reuse and release archiving.

## License
MIT (see `LICENSE`).
