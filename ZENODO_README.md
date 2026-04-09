# GFF3-LTR-Mapper Zenodo README

This record archives version `0.3.0` of **GFF3-LTR-Mapper**.

## What This Software Does

GFF3-LTR-Mapper is a command-line workflow for:

- auditing **EDTA intact LTR GFF3** annotations
- validating intact-element structure
- writing scientist-facing TSV summaries
- comparing genome, scaffold, and coordinate-window TE patterns across samples

The current release is designed for **EDTA intact-element output only** and reports on **intact LTR elements only**.

## Main Outputs

- `summary.tsv`: one row per EDTA `repeat_region`
- `validation.tsv`: per-element QC pass/warn/fail report
- `scientist_regions.tsv`: simplified region summary
- `scientist_elements.tsv`: simplified intact-element detail table
- `batch_samples.tsv`: per-sample batch summary
- `scientist_batch.tsv`: cross-sample comparison table

## Important Scope Limits

- Input should be EDTA intact LTR GFF3 output.
- Cross-species `--region` comparisons are **coordinate-window comparisons**, not homology-aware or synteny-aware locus comparisons.
- Aggregate postcard output in this release is **ASCII-only**.

## Software Package Contents

The software archive contains:

- source distribution (`sdist`)
- wheel distribution (`wheel`)
- license
- citation metadata
- main project README

The example/documentation archive contains:

- usage documentation
- output examples
- one real example TSV

## Repository

GitHub repository:

<https://github.com/ViroDyn-Ltd/GFF3-LTR-Mapper>
