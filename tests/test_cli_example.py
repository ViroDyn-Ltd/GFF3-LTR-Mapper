import csv
import shutil
import subprocess
import sys


def test_cli_generates_qc_and_identity_outputs(tmp_path, sample_gff_path):
    outdir = tmp_path / "example_chr2"
    summary_path = outdir / "summary.tsv"
    validation_path = outdir / "validation.tsv"
    cmd = [
        sys.executable,
        "-m",
        "gff3_ltr_map.cli",
        str(sample_gff_path),
        "--chrom",
        "chr_2",
        "--outdir",
        str(outdir),
        "--summary",
        str(summary_path),
        "--validation-report",
        str(validation_path),
        "--visual",
        "postcard+quantiles",
        "--identity",
        "bins=0.90..0.94,>=0.94",
        "--top-k",
        "5",
        "--min-n",
        "2",
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr

    assert summary_path.exists()
    assert validation_path.exists()
    assert (outdir / "cohort_aggregates.tsv").exists()
    scientist_region_path = outdir / "scientist_regions.tsv"
    scientist_element_path = outdir / "scientist_elements.tsv"
    assert scientist_region_path.exists()
    assert scientist_element_path.exists()
    identity_path = outdir / "identity_bins.tsv"
    assert identity_path.exists()
    assert "[gff3-ltr-map] Simplified region table:" in result.stdout
    assert "dominant_sf" in result.stdout

    with validation_path.open() as handle:
        validation_rows = {row["element_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert validation_rows["repeat_region_high_1"]["qc_status"] == "PASS"
    assert validation_rows["repeat_region_high_1"]["is_intact"] == "true"

    with identity_path.open() as handle:
        rows = {row["group"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["chr_2:>=0.940"]["n_elements"] == "2"
    assert rows["chr_2:>=0.940"]["top_motifs"] == "TACT (1, 50%), TGCA (1, 50%)"
    assert rows["chr_2:0.900-0.940"]["notes"] == "LOW N (n=1)"
    assert rows["chr_2:>=0.940"]["tsd_len_median"] == "5"

    with scientist_region_path.open() as handle:
        region_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert region_rows == [
        {
            "sample": "sample_edta_chr2",
            "region": "chr_2",
            "window_bp": "",
            "intact_elements": "3",
            "intact_per_mbp": "",
            "qc_warn": "0",
            "median_ltr_identity": "0.979",
            "identity_band": "high",
            "median_length_bp": "1001",
            "coverage_pct": "",
            "pct_with_both_tsd": "100",
            "dominant_superfamily": "Copia",
            "dominant_motif": "TGCA",
            "note": "",
        }
    ]

    with scientist_element_path.open() as handle:
        simple_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert simple_rows[0]["identity_band"] in {"high", "moderate", "very_high"}
    assert "qc_note" in simple_rows[0]
    assert simple_rows[0]["edta_id"].startswith("repeat_region_")

    postcards = outdir / "identity_postcards"
    assert (postcards / "chr_2_0_940.txt").exists()
    assert (postcards / "chr_2_0_900_0_940.txt").exists()


def test_cli_directory_batch_mode(tmp_path, sample_gff_path):
    input_dir = tmp_path / "batch_inputs"
    input_dir.mkdir()
    shutil.copy(sample_gff_path, input_dir / "sample_a.gff3")
    shutil.copy(sample_gff_path, input_dir / "sample_b.gff3")

    outdir = tmp_path / "batch_out"
    cmd = [
        sys.executable,
        "-m",
        "gff3_ltr_map.cli",
        str(input_dir),
        "--outdir",
        str(outdir),
        "--chrom",
        "genome",
        "--identity",
        "all",
        "--visual",
        "none",
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr

    batch_samples = outdir / "batch_samples.tsv"
    batch_superfamilies = outdir / "batch_superfamilies.tsv"
    scientist_batch = outdir / "scientist_batch.tsv"
    assert batch_samples.exists()
    assert batch_superfamilies.exists()
    assert scientist_batch.exists()
    assert (outdir / "samples" / "sample_a" / "summary.tsv").exists()
    assert (outdir / "samples" / "sample_b" / "validation.tsv").exists()
    assert (outdir / "samples" / "sample_a" / "scientist_elements.tsv").exists()
    assert (outdir / "samples" / "sample_a" / "scientist_regions.tsv").exists()
    assert "[gff3-ltr-map] Simplified region table:" in result.stdout

    with batch_samples.open() as handle:
        rows = {row["sample"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["sample_a"]["n_repeat_regions"] == "3"
    assert rows["sample_a"]["n_intact"] == "3"
    assert rows["sample_b"]["pct_intact"] == "100"

    with scientist_batch.open() as handle:
        simple_rows = {row["sample"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert simple_rows["sample_a"]["dominant_superfamily"] == "Copia"
    assert simple_rows["sample_a"]["intact_elements"] == "3"
    assert simple_rows["sample_a"]["region"] == "genome"
    assert simple_rows["sample_a"]["window_bp"] == ""


def test_cli_region_filter_limits_outputs(tmp_path, sample_gff_path):
    outdir = tmp_path / "region_slice"
    summary_path = outdir / "summary.tsv"
    validation_path = outdir / "validation.tsv"
    cmd = [
        sys.executable,
        "-m",
        "gff3_ltr_map.cli",
        str(sample_gff_path),
        "--outdir",
        str(outdir),
        "--summary",
        str(summary_path),
        "--validation-report",
        str(validation_path),
        "--chrom",
        "chr_2",
        "--region",
        "1-4500",
        "--identity",
        "all",
        "--visual",
        "none",
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "Region filter chr_2:1-4500 retained 2 overlapping repeat_region features" in result.stdout

    with (outdir / "scientist_regions.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [
        {
            "sample": "sample_edta_chr2",
            "region": "chr_2:1-4500",
            "window_bp": "4500",
            "intact_elements": "2",
            "intact_per_mbp": "444.444",
            "qc_warn": "0",
            "median_ltr_identity": "0.9805",
            "identity_band": "very_high",
            "median_length_bp": "1101",
            "coverage_pct": "48.9333",
            "pct_with_both_tsd": "100",
            "dominant_superfamily": "Copia",
            "dominant_motif": "TGCA",
            "note": "",
        }
    ]

    with summary_path.open() as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(summary_rows) == 2


def test_cli_directory_batch_region_filter_auto_primary_scaffold(tmp_path, sample_gff_path):
    input_dir = tmp_path / "batch_region_inputs"
    input_dir.mkdir()
    shutil.copy(sample_gff_path, input_dir / "sample_a.gff3")
    shutil.copy(sample_gff_path, input_dir / "sample_b.gff3")

    outdir = tmp_path / "batch_region_out"
    cmd = [
        sys.executable,
        "-m",
        "gff3_ltr_map.cli",
        str(input_dir),
        "--outdir",
        str(outdir),
        "--region",
        "1-4500",
        "--identity",
        "all",
        "--visual",
        "none",
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "primary scaffold inferred as chr_2" in result.stdout

    with (outdir / "scientist_batch.tsv").open() as handle:
        rows = {row["sample"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["sample_a"]["region"] == "chr_2:1-4500"
    assert rows["sample_a"]["window_bp"] == "4500"
    assert rows["sample_a"]["intact_elements"] == "2"
    assert rows["sample_a"]["intact_per_mbp"] == "444.444"
    assert rows["sample_a"]["coverage_pct"] == "48.9333"


def test_cli_directory_batch_region_filter_handles_empty_gff(tmp_path, sample_gff_path):
    input_dir = tmp_path / "batch_empty_inputs"
    input_dir.mkdir()
    shutil.copy(sample_gff_path, input_dir / "sample_a.gff3")
    (input_dir / "sample_empty.gff3").write_text("##gff-version 3\n", encoding="utf-8")

    outdir = tmp_path / "batch_empty_out"
    cmd = [
        sys.executable,
        "-m",
        "gff3_ltr_map.cli",
        str(input_dir),
        "--outdir",
        str(outdir),
        "--region",
        "1-4500",
        "--identity",
        "all",
        "--visual",
        "none",
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "no repeat_region features; primary scaffold unavailable" in result.stdout

    with (outdir / "scientist_batch.tsv").open() as handle:
        rows = {row["sample"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["sample_empty"]["region"] == "NA:1-4500"
    assert rows["sample_empty"]["intact_elements"] == "0"
    assert rows["sample_empty"]["note"] == "no repeat_region features; primary scaffold unavailable; no intact elements"
