import pytest

from gff3_ltr_map import parser


def test_repeat_region_properties(sample_gff_path):
    regions, children = parser.load_gff(str(sample_gff_path))
    # Ensure three repeat regions are detected
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
    assert rr.ltr_identity == pytest.approx(0.979)
    assert rr.ltr5_len == 201
    assert rr.ltr3_len == 200
    assert rr.internal_len == 600
    assert rr.motif == "TGCA"
    assert rr.tsd == "AATAT"
    assert rr.has_both_tsd
