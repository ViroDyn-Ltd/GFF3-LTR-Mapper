import os
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def sample_gff_path() -> Path:
    """Path to the bundled miniature EDTA GFF3 file."""

    return Path(__file__).parent / "data" / "sample_edta_chr2.gff3"


@pytest.fixture(scope="session")
def real_data_dir() -> Path | None:
    """Optional ViroDyn real-data directory used for integration smoke tests."""

    configured = os.getenv("VIRODYN_EDTA_GFF3_DIR")
    if not configured:
        return None
    path = Path(configured)
    if not path.exists():
        return None
    return path
