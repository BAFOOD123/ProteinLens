"""Shared pytest fixtures.

``offline_cache`` points ProteinLens at a temporary copy of the bundled fixture
cache and makes any attempt to use the network fail loudly, so offline tests
prove that cached data is enough.
"""

import shutil
import urllib.request
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def offline_cache(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    shutil.copytree(FIXTURES / "cache", cache)
    monkeypatch.setenv("PROTEINLENS_CACHE", str(cache))

    def no_network(*args, **kwargs):
        raise AssertionError("offline test tried to use the network")

    monkeypatch.setattr(urllib.request, "urlopen", no_network)
    return cache


@pytest.fixture
def online_cache(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    monkeypatch.setenv("PROTEINLENS_CACHE", str(cache))
    return cache
