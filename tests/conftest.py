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
    """Normalize a string into a slug identifier.

    :param s: Input string to normalize.
    :type s: str
    :returns: Slugified identifier.
    :rtype: str
    """
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


@pytest.fixture(scope="session", autouse=True)
def debug_base() -> Path:
    """Provide a base temporary directory for debug outputs.

    :returns: Base debug directory path.
    :rtype: pathlib.Path
    """
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
def test_debug_dir(request: pytest.FixtureRequest, debug_base: Path) -> Path:
    """Provide a per-test debug directory.

    :param request: Pytest fixture request object.
    :type request: pytest.FixtureRequest
    :param debug_base: Base debug directory.
    :type debug_base: pathlib.Path
    :returns: Path to the per-test debug directory.
    :rtype: pathlib.Path
    """
    d = debug_base / _slug(request.node.nodeid)
    d.mkdir(parents=True, exist_ok=True)
    return d


@pytest.fixture(scope="session")
def example_package() -> Package:
    """Load the example data package datapackage for tests.

    :returns: Example datapackage instance.
    :rtype: datapackage.Package
    """
    package_path = (
        Path(__file__).resolve().parents[1] / "dev" / "example data package" / "datapackage.json"
    )
    return Package(str(package_path))


@pytest.fixture(scope="session")
def example_trails(example_package: Package) -> trails_module.Trails:
    """Build a Trails instance for tests.

    :param example_package: Example datapackage fixture.
    :type example_package: datapackage.Package
    :returns: Trails instance initialized for tests.
    :rtype: trails.trails.Trails
    """
    return trails_module.Trails(example_package, interpolate_annual=False)
