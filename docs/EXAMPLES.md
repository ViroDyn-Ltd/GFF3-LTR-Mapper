# Example Outputs

These captures come from representative EDTA `intact` datasets and document both a focused single-file sweep and a real cross-species comparison.

## Reproducing this example
Set an environment variable (or substitute your own absolute path) pointing to the EDTA file:

```bash
export GFF3=/path/to/your/EDTA.intact.gff3
```

Then run:

```bash
python -m gff3_ltr_map.cli \
  "$GFF3" \
  --chrom chr_2 \
  --outdir runs/example_chr2 \
  --summary runs/example_chr2/summary.tsv \
  --visual postcard+quantiles \
  --identity 'bins=0.90..0.94,>=0.94' \
  --top-k 5 \
  --min-n 25
```
Adjust `GFF3` or `--outdir` as needed for your workstation.

## Cross-species TSV example bundled in the repository

The repository includes one real comparison TSV:

- `verified_output/cross_species_primary_window_1_3000000.tsv`

This file was generated from the ViroDyn EDTA directory using:

```bash
gff3-ltr-map /path/to/edta_gff3_directory \
  --outdir runs/cross_species_window \
  --region 1-3000000 \
  --identity all \
  --visual none
```

It is sorted by strongest window signal first and is intended as the canonical example of the scientist-facing batch output.

## Identity-bin table (CLI & TSV excerpt)
```
group	n	median_len	median_id	top_motifs	notes
chr_2:0.900-0.940	1	4424.0	0.931	TGCA (1, 100%)	LOW N (n=1)
chr_2:>=0.940	49	6134.0	0.979	TGCA (41, 84%), TGTA (3, 6%), TACA (2, 4%), TACT (1, 2%), TATA (1, 2%)	TSD consensus <40%
```

## ASCII postcards

### chr_2 ≥ 0.94 cohort
```
> AVG chr_2:>=0.940  range:0.940–-
  n: 49
  median len: 6134.0 bp
  IQR: 5120–10091 bp
  identity mean:0.977 median:0.979 range:0.940–-
  LTR5 median: 516.0 bp   LTR3 median: 510.0 bp
  Internal median: 5108.0 bp
  strand +:20 (41%)/-:18 (37%)

========-----------------------------------------------------------------------------------=========
LTR5:516 | INT:5624 | LTR3:6134
Q25:5120bp  Q75:10091bp
```

### chr_2 0.90–0.94 cohort
```
> AVG chr_2:0.900-0.940  range:0.900–0.940
  n: 1 (LOW N (n=1))
  median len: 4424.0 bp
  IQR: 4424–4424 bp
  identity mean:0.931 median:0.931 range:0.900–0.940
  LTR5 median: 492.0 bp   LTR3 median: 492.0 bp
  Internal median: 3440.0 bp
  strand +:1 (100%)

===========-----------------------------------------------------------------------------============
LTR5:492 | INT:3932 | LTR3:4424
Q25:4424bp  Q75:4424bp
```

The aggregate postcard output is now ASCII-only. SVG postcard generation has been removed from the tool.
