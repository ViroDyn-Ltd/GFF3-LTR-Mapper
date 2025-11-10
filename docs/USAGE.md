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
  --visual postcard+quantiles \
  --out text+svg \
  --summary runs/default_sweep/summary.tsv \
  --top-k 5 \
  --min-n 25
```
- **What happens**: genome + every chromosome produces TSV rows and postcards. CLI prints an identity-bin table and where the postcards live.
- **Outputs**:
  - `runs/default_sweep/summary.tsv`
  - `runs/default_sweep/identity_bins.tsv`
  - `runs/default_sweep/identity_postcards/*.txt` and `.svg`

## Focus on a single chromosome
```
gff3-ltr-map /path/to/EDTA.intact.gff3 \
  --outdir runs/chr2_visuals \
  --chrom chr_2 \
  --identity 'bins=0.90..0.94,>=0.94' \
  --visual postcard+quantiles \
  --out text+svg \
  --summary runs/chr2_visuals/summary.tsv
```

## Genome-only
```
gff3-ltr-map /path/to/EDTA.intact.gff3 --chrom genome --visual postcard --out text
```

## Important flags
- `--identity auto` (default): outputs "all" and ">=0.98" cohorts per scope.
- `--identity all` or `>=0.97` or `bins=0.90..0.94,>=0.94`
- `--visual none|postcard|postcard+quantiles`
- `--out text|text+svg`
- `--top-k N` controls motif/TSD leaderboards.
- `--min-n N` defines the low-N warning threshold.
- `--outdir PATH` lets you place outputs anywhere (defaults to `runs/`).
- `--summary FILE`, `--aggregate-tsv FILE`, `--aggregate-json FILE` override file locations.
- `--bed FILE` writes repeat_region spans (0-based BED).
- `--ascii/--svg` (legacy) still emit per-element postcards if needed.

## CLI logging
Each run prints:
1. Input being read + element count
2. Where the summary TSV (and BED) were written
3. Identity-bin table (group, n, median length/identity, motifs, notes)
4. Where postcards were saved (`identity_postcards/` with file counts)

## Example ASCII postcard
```
> AVG chr_2:>=0.940  range:0.940–-
  n: 530
  median len: 7762.0 bp
  IQR: 5074–10784 bp
  identity mean:0.989 median:0.991 range:0.940–-
  LTR5 median: 1047.5 bp   LTR3 median: 1047.5 bp
  Internal median: 5667.0 bp
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
- **Age estimates** only appear if you pass `--substitution-rate` and should be treated as rough proxies.
