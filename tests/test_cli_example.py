import csv
import subprocess
import sys


def test_cli_generates_identity_outputs(tmp_path, sample_gff_path):
    outdir = tmp_path / "example_chr2"
    summary_path = outdir / "summary.tsv"
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
        "--visual",
        "postcard+quantiles",
        "--out",
        "text+svg",
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
    identity_path = outdir / "identity_bins.tsv"
    assert identity_path.exists()

    with identity_path.open() as handle:
        rows = {row["group"]: row for row in csv.DictReader(handle, delimiter="\t")}

    assert rows["chr_2:>=0.940"]["n_elements"] == "2"
    assert rows["chr_2:>=0.940"]["top_motifs"] == "TACT (1, 50%), TGCA (1, 50%)"
    assert rows["chr_2:0.900-0.940"]["notes"] == "LOW N (n=1)"

    postcards = outdir / "identity_postcards"
    assert (postcards / "chr_2_0_940.txt").exists()
    assert (postcards / "chr_2_0_940.svg").exists()
    assert (postcards / "chr_2_0_900_0_940.txt").exists()
