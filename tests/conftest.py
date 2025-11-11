from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def sample_gff_path() -> Path:
    """Path to the bundled miniature EDTA GFF3 file."""

    return Path(__file__).parent / "data" / "sample_edta_chr2.gff3"
