"""The website must be rebuilt whenever results change: site/index.html has to match a fresh build."""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def test_site_is_up_to_date(tmp_path):
    subprocess.run(
        [sys.executable, str(ROOT / "site" / "build_site.py"), "--out", str(tmp_path)],
        check=True,
        capture_output=True,
    )
    fresh = (tmp_path / "index.html").read_text()
    committed = (ROOT / "site" / "index.html").read_text()
    assert fresh == committed, "site/index.html is stale: run `python site/build_site.py`"
