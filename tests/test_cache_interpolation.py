"""Tests for annual-interpolation cache identity."""

from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from trails.cache_interpolation import _cache_key

TRACKED_PATH = "inventories/model/pathway/2030/A_matrix.csv"


def _package(base_path: Path, *, title: str = "same package") -> SimpleNamespace:
    matrix_path = base_path / TRACKED_PATH
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    matrix_path.write_text("0;0;1\n", encoding="utf-8")
    resource = SimpleNamespace(
        name="A",
        descriptor={"path": TRACKED_PATH},
        source=str(matrix_path),
    )
    return SimpleNamespace(
        base_path=str(base_path),
        descriptor={"name": "test-package", "title": title},
        resources=[resource],
    )


def _key(package: SimpleNamespace) -> str:
    return _cache_key(
        package,
        value_dtype="float64",
        index_dtype="int32",
        interpolation_start_year_offset=-20,
        interpolation_end_year_offset=20,
    )


def test_ephemeral_zip_extractions_ignore_resource_mtime(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Equivalent ZIP extractions should resolve to the same cache key."""
    monkeypatch.setattr(
        "trails.cache_interpolation.tempfile.gettempdir",
        lambda: str(tmp_path),
    )
    first = _package(tmp_path / "tmp-first-datapackage")
    second = _package(tmp_path / "tmp-second-datapackage")

    first_source = Path(first.resources[0].source)
    second_source = Path(second.resources[0].source)
    os.utime(first_source, ns=(1_700_000_000_000_000_000,) * 2)
    os.utime(second_source, ns=(1_800_000_000_000_000_000,) * 2)

    assert _key(first) == _key(second)


def test_stable_directory_packages_still_use_resource_mtime(tmp_path: Path) -> None:
    """Changing a normal package resource should still invalidate its cache."""
    package = _package(tmp_path / "stable-package")
    source = Path(package.resources[0].source)
    first_key = _key(package)

    os.utime(source, ns=(1_800_000_000_000_000_000,) * 2)

    assert _key(package) != first_key


def test_package_descriptor_is_part_of_cache_identity(tmp_path: Path) -> None:
    """Package-level metadata must not be shadowed by resource metadata."""
    first = _package(tmp_path / "first", title="first package")
    second = _package(tmp_path / "second", title="second package")
    shared_mtime = 1_700_000_000_000_000_000
    os.utime(Path(first.resources[0].source), ns=(shared_mtime,) * 2)
    os.utime(Path(second.resources[0].source), ns=(shared_mtime,) * 2)

    assert _key(first) != _key(second)
