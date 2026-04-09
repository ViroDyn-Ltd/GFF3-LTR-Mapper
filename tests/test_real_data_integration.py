import csv
import subprocess
import sys

import pytest


@pytest.mark.integration
def test_real_data_smoke(tmp_path, real_data_dir):
    if real_data_dir is None:
        pytest.skip("VIRODYN_EDTA_GFF3_DIR is not configured")

    arabidopsis = real_data_dir / "arabidopsis_thaliana.gff3"
    if not arabidopsis.exists():
        pytest.skip("arabidopsis_thaliana.gff3 is not available")

    outdir = tmp_path / "real_smoke"
    cmd = [
        sys.executable,
        "-m",
        "gff3_ltr_map.cli",
        str(arabidopsis),
        "--outdir",
        str(outdir),
        "--summary",
        str(outdir / "summary.tsv"),
        "--validation-report",
        str(outdir / "validation.tsv"),
        "--chrom",
        "genome",
        "--identity",
        "auto",
        "--visual",
        "none",
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr

    with (outdir / "validation.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 435
    assert all(row["qc_status"] == "PASS" for row in rows)


@pytest.mark.integration
def test_real_data_directory_batch_smoke(tmp_path, real_data_dir):
    if real_data_dir is None:
        pytest.skip("VIRODYN_EDTA_GFF3_DIR is not configured")

    outdir = tmp_path / "real_batch"
    cmd = [
        sys.executable,
        "-m",
        "gff3_ltr_map.cli",
        str(real_data_dir),
        "--outdir",
        str(outdir),
        "--chrom",
        "genome",
        "--identity",
        "all",
        "--visual",
        "none",
        "--limit-files",
        "2",
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr

    with (outdir / "batch_samples.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert all(row["pct_intact"] == "100" for row in rows)
