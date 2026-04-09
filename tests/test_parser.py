import gzip

import pytest

from gff3_ltr_map import parser


def test_repeat_region_properties(sample_gff_path):
    regions, children = parser.load_gff(str(sample_gff_path))
    assert set(regions) == {
        "repeat_region_high_1",
        "repeat_region_high_2",
        "repeat_region_low_1",
    }

    rr = parser.to_repeat_region_object(
        "repeat_region_high_1",
        regions["repeat_region_high_1"],
        children["repeat_region_high_1"],
    )

    assert rr.scaffold == "chr_2"
    assert rr.classification == "LTR/Copia"
    assert rr.superfamily == "Copia"
    assert rr.ltr_identity == pytest.approx(0.979)
    assert rr.ltr5_len == 201
    assert rr.ltr3_len == 200
    assert rr.internal_len == 600
    assert rr.motif == "TGCA"
    assert rr.tsd == "AATAT"
    assert rr.tsd_len == 5
    assert rr.has_both_tsd
    assert rr.is_intact
    assert rr.qc_status == "PASS"
    assert rr.child_signature == "LTR_retrotransposon=1; long_terminal_repeat=2; target_site_duplication=2"


def test_parser_supports_gzip(sample_gff_path, tmp_path):
    gz_path = tmp_path / "sample_edta_chr2.gff3.gz"
    with sample_gff_path.open("rb") as src, gzip.open(gz_path, "wb") as dst:
        dst.write(src.read())

    regions, children = parser.load_gff(str(gz_path))
    rr = parser.to_repeat_region_object(
        "repeat_region_high_2",
        regions["repeat_region_high_2"],
        children["repeat_region_high_2"],
    )

    assert rr.superfamily == "Gypsy"
    assert rr.tsd == "TACAT"
    assert rr.qc_status == "PASS"


def test_invalid_repeat_region_is_flagged(tmp_path):
    invalid = tmp_path / "invalid.gff3"
    invalid.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tEDTA\trepeat_region\t100\t300\t.\t+\t.\tID=bad_rr;Classification=LTR/Copia;Method=structural",
                "chr1\tEDTA\tLTR_retrotransposon\t100\t300\t.\t+\t.\tID=bad_rr_rt;Parent=bad_rr;Classification=LTR/Copia;ltr_identity=0.97",
                "chr1\tEDTA\tlong_terminal_repeat\t100\t150\t.\t+\t.\tID=bad_rr_ltr5;Parent=bad_rr",
                "chr1\tEDTA\ttarget_site_duplication\t95\t99\t.\t+\t.\tID=bad_rr_tsd5;Parent=bad_rr;tsd=AAAAA",
                "chr1\tEDTA\ttarget_site_duplication\t301\t305\t.\t+\t.\tID=bad_rr_tsd3;Parent=bad_rr;tsd=AAAAA",
                "",
            ]
        ),
        encoding="utf-8",
    )

    regions, children = parser.load_gff(str(invalid))
    rr = parser.to_repeat_region_object("bad_rr", regions["bad_rr"], children["bad_rr"])

    assert not rr.is_intact
    assert rr.qc_status == "FAIL"
    assert "EXPECTED_TWO_LTRS" in rr.validation_error_text
