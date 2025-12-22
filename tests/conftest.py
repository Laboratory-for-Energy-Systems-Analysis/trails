# conftest.py
import os
from pathlib import Path
import re
import json
import sys
import platform
import pytest

from datapackage import Package

import trails.trails as trails_module

BASE = Path(os.environ.get("PYTEST_DEBUG_DIR", ".pytest-debug"))


def _slug(s: str) -> str:
    """Normalize a string into a slug identifier."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


@pytest.fixture(scope="session", autouse=True)
def debug_base():
    """Provide a base temporary directory for debug outputs."""
    BASE.mkdir(parents=True, exist_ok=True)
    # Snapshot environment info once per session
    (BASE / "env.json").write_text(
        json.dumps(
            {
                "python": sys.version,
                "executable": sys.executable,
                "platform": platform.platform(),
                "cwd": os.getcwd(),
                "env": dict(os.environ),
            },
            indent=2,
        )
    )
    return BASE


@pytest.fixture
def test_debug_dir(request, debug_base: Path):
    """Provide a per-test debug directory."""
    d = debug_base / _slug(request.node.nodeid)
    d.mkdir(parents=True, exist_ok=True)
    return d


@pytest.fixture(scope="session")
def example_package() -> Package:
    """Load the example datapackage for tests."""
    package_path = (
        Path(__file__).resolve().parents[1] / "dev" / "example" / "datapackage.json"
    )
    return Package(str(package_path))


@pytest.fixture(scope="session")
def example_trails(example_package: Package):
    """Build a Trails instance for tests."""
    return trails_module.Trails(example_package, interpolate_annual=False)
