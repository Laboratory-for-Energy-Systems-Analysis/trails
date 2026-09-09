"""Cache notices and explicit cleanup without touching the user's cache."""

import importlib
from pathlib import Path
import warnings

import pytest

import trails
from trails import cache


@pytest.fixture
def cache_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(cache.platformdirs, "user_data_path", lambda **kwargs: tmp_path)
    monkeypatch.setattr(cache, "_cache_checked", False)
    return tmp_path / "cache"


@pytest.mark.parametrize("size", [0, 5 * 1024**3 - 1, 5 * 1024**3])
def test_no_warning_at_or_below_five_gib(cache_dir, monkeypatch, size):
    monkeypatch.setattr(cache, "_cache_size", lambda path: size)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cache._warn_if_cache_large()
    assert not caught
    assert not cache_dir.exists()


def test_import_warns_only_once_and_preserves_files(cache_dir, monkeypatch):
    cache_dir.mkdir()
    payload = cache_dir / "inventory.npz"
    payload.write_bytes(b"unchanged")
    calls = []

    def size(path):
        calls.append(path)
        return int(48.3 * 1024**3)

    monkeypatch.setattr(cache, "_cache_size", size)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        importlib.reload(trails)
        importlib.reload(trails)
        cache._warn_if_cache_large()
    assert len(caught) == 1
    assert caught[0].category is UserWarning
    assert str(caught[0].message) == (
        "TRAILS interpolation cache uses 48.3 GiB. "
        "When no TRAILS calculations are running, use trails.clear_cache() "
        "to reclaim space. Cached data will be rebuilt when needed."
    )
    assert calls == [cache_dir]
    assert payload.read_bytes() == b"unchanged"


def test_warning_just_above_threshold(cache_dir, monkeypatch):
    monkeypatch.setattr(cache, "_cache_size", lambda path: 5 * 1024**3 + 1)
    with pytest.warns(UserWarning, match="trails.clear_cache"):
        cache._warn_if_cache_large()


@pytest.mark.parametrize("exists", [False, True])
def test_missing_or_empty_cache_is_silent(cache_dir, exists):
    if exists:
        cache_dir.mkdir()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cache._warn_if_cache_large()
    assert not caught
    assert cache_dir.exists() is exists


def test_scan_errors_do_not_break_import(cache_dir, monkeypatch):
    def denied(path):
        raise PermissionError("unreadable cache")

    monkeypatch.setattr(cache, "_cache_size", denied)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        importlib.reload(trails)
    assert not caught


def test_size_scans_metadata_without_reading_payloads(cache_dir, monkeypatch):
    nested = cache_dir / "interp_example"
    nested.mkdir(parents=True)
    payload = nested / "A.npz"
    payload.write_bytes(b"matrix")
    info = payload.stat()
    expected = info.st_blocks * 512 if hasattr(info, "st_blocks") else info.st_size

    def no_read(*args, **kwargs):
        raise AssertionError("cache size check must not open payloads")

    monkeypatch.setattr(Path, "open", no_read)
    assert cache._cache_size(cache_dir) == expected


def test_size_does_not_follow_symlinks(cache_dir, tmp_path):
    cache_dir.mkdir()
    outside = tmp_path / "source"
    outside.mkdir()
    (outside / "matrix").write_bytes(b"source data")
    try:
        (cache_dir / "linked_dir").symlink_to(outside, target_is_directory=True)
        (cache_dir / "linked_file").symlink_to(outside / "matrix")
        alias = tmp_path / "alias"
        alias.symlink_to(outside, target_is_directory=True)
    except OSError:
        pytest.skip("symlinks unavailable on this platform")
    assert cache._cache_size(cache_dir) == 0
    assert cache._cache_size(alias) == 0


def test_clear_cache_preserves_source_data(cache_dir, tmp_path):
    cache_dir.mkdir()
    (cache_dir / "A.npz").write_bytes(b"cache")
    source = tmp_path / "source.zip"
    source.write_bytes(b"source")
    assert cache.clear_cache() == cache_dir
    assert list(cache_dir.iterdir()) == []
    assert source.read_bytes() == b"source"


def test_clear_cache_reports_deletion_failure(cache_dir, monkeypatch):
    cache_dir.mkdir()

    def denied(path):
        raise PermissionError(f"Cannot delete {path}")

    monkeypatch.setattr(cache.shutil, "rmtree", denied)
    with pytest.raises(PermissionError, match="Cannot delete"):
        cache.clear_cache()
