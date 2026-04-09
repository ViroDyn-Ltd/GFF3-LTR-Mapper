# gff3-ltr-map Usage Guide

This document captures the key terminal commands and outputs so you can run the tool entirely from the CLI.

## Getting started
1. Download the GitHub ZIP (or clone) and unzip it somewhere convenient.
2. Ensure Python 3.9+ is on your PATH.
3. Run commands via `python -m gff3_ltr_map.cli ...` from the project directory (examples below continue to use the `gff3-ltr-map` shorthand; feel free to create your own alias or shell script if desired).

## Core command (genome + per-chromosome sweep)
```
gff3-ltr-map /path/to/EDTA.intact.gff3 \
  --outdir runs/default_sweep \
  --validation error \
  --visual postcard+quantiles \
  --summary runs/default_sweep/summary.tsv \
  --top-k 5 \
  --min-n 25
```
- **What happens**: genome + every chromosome produces `summary.tsv`, `validation.tsv`, `cohort_aggregates.tsv`, identity-bin tables, and optional postcards. Only intact elements enter the cohort calculations, but every `repeat_region` still appears in the QC tables.
- **Outputs**:
  - `runs/default_sweep/summary.tsv`
  - `runs/default_sweep/validation.tsv`
  - `runs/default_sweep/scientist_regions.tsv`
  - `runs/default_sweep/scientist_elements.tsv`
  - `runs/default_sweep/cohort_aggregates.tsv`
  - `runs/default_sweep/identity_bins.tsv`
  - `runs/default_sweep/identity_postcards/*.txt`

## Directory / batch mode
```
gff3-ltr-map /path/to/edta_outputs \
  --outdir runs/batch_qc \
  --chrom genome \
  --identity all \
  --visual none
```
- **What happens**: every `.gff3` / `.gff3.gz` under the directory is processed into `runs/batch_qc/samples/<sample>/...`
- **Root outputs**:
  - `runs/batch_qc/batch_samples.tsv`
  - `runs/batch_qc/batch_superfamilies.tsv`
  - `runs/batch_qc/scientist_batch.tsv`

## Focus on a single chromosome
```
gff3-ltr-map /path/to/EDTA.intact.gff3 \
  --outdir runs/chr2_visuals \
  --chrom chr_2 \
  --identity 'bins=0.90..0.94,>=0.94' \
  --visual postcard+quantiles \
  --summary runs/chr2_visuals/summary.tsv
```

## Focus on a coordinate range
```
gff3-ltr-map /path/to/EDTA.intact.gff3 \
  --outdir runs/chr2_window \
  --chrom chr_2 \
  --region 1-3000000 \
  --identity all \
  --visual none
```
- **What happens**: the CLI keeps only elements overlapping `chr_2:1-3000000`, then writes region-level summaries and the detailed intact-element table for that window.

## Cross-species comparison of the same window
```
gff3-ltr-map /path/to/edta_outputs \
  --outdir runs/cross_species_window \
  --region 1-3000000 \
  --identity all \
  --visual none
```
- **What happens**: for each sample, the CLI infers a primary scaffold from the EDTA coordinates, keeps only elements overlapping `1-3000000`, and writes a root `scientist_batch.tsv` with one directly comparable row per species.
- **Key columns**:
  - `window_bp`: requested comparison window length
  - `intact_elements`: intact LTR count in that window
  - `intact_per_mbp`: window-normalized count
  - `coverage_pct`: percent of the window covered by intact element spans

## Genome-only
```
gff3-ltr-map /path/to/EDTA.intact.gff3 --chrom genome --visual postcard
```

## Important flags
- `--identity auto` (default): outputs "all" and ">=0.98" cohorts per scope.
- `--identity all` or `>=0.97` or `bins=0.90..0.94,>=0.94`
- `--validation error|warn|off` controls how the CLI reacts to non-intact structural blocks.
- `--visual none|postcard|postcard+quantiles`
- `--group-aggregates genome,scaffold,superfamily` writes cohort summaries by those groups.
- `--top-k N` controls motif/TSD leaderboards.
- `--min-n N` defines the low-N warning threshold.
- `--outdir PATH` lets you place outputs anywhere (defaults to `runs/`).
- `--summary FILE`, `--validation-report FILE`, `--aggregate-tsv FILE`, `--aggregate-json FILE` override single-file output locations.
- `--scientist-tsv FILE` overrides the simplified scientist-facing TSV location.
- `--region scaffold:start-end` or `--chrom chr_2 --region 1-3000000` restricts outputs to overlapping elements in that interval.
- `--region 1-3000000` in directory mode compares that same window on each inferred primary scaffold.
- `--scientist-cli-max-rows N` controls how many simplified rows are printed directly to the terminal.
- `--bed FILE` writes repeat_region spans (0-based BED).
- `--limit-files N` trims directory mode to the first `N` files.
- `--ascii` still emits per-element ASCII postcards if needed.

## CLI logging
Each run prints:
1. Input being read + element count
2. Any active region filter and where the summary / validation / cohort TSVs were written
3. Identity-bin table (single-file mode) with median TSD length included
4. The simplified region summary directly in the terminal
5. Where postcards were saved (`identity_postcards/` with file counts)

## Example ASCII postcard
```
> AVG chr_2:>=0.940  range:0.940–-
  n: 530
  median len: 7762.0 bp
  IQR: 5074–10784 bp
  identity mean:0.989 median:0.991 range:0.940–-
  LTR5 median: 1047.5 bp   LTR3 median: 1047.5 bp
  Internal median: 5667.0 bp
  TSD pair: 100.0%   median TSD len: 5
  strand +:200 (38%)/-:220 (42%)

=============-------------------------------------------------------------------------==============
LTR5:1048 | INT:6714 | LTR3:7762
Q25:5074bp  Q75:10784bp
```

## Example identity-bin table (printed + TSV)
```
group         n   median_len  median_id  top_motifs
chr_2:0.900-0.940  6   5507.5      0.927   TGCA (5, 83%), TACA (1, 17%)  LOW N (n=6); TSD consensus <40%
chr_2:>=0.940      530 7762.0      0.991   TGCA (485, 92%), TACA (18, 3%), ...   TSD consensus <40%
```

## Notes for scientists
- **Quartile markers** `Q25`/`Q75` visualize the interquartile range of element lengths.
- **Consensus warnings** appear when motifs/TSDs lack a dominant mode (>40%).
- **Low N** warnings highlight cohorts with limited support (`--min-n`).
- **QC reports** preserve pass/warn/fail status for every `repeat_region` and make exclusion rules explicit.
- **Scientist summaries** now lead with genome / scaffold / region rows; EDTA `repeat_region` IDs are kept only in the secondary element-detail TSV.
- **Cross-species windows** are easiest to compare from `scientist_batch.tsv`, which now carries window-normalized counts and coverage metrics.
- **Age estimates** only appear if you pass `--substitution-rate` and should be treated as rough proxies.
