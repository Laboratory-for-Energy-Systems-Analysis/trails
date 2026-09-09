from __future__ import annotations

from pathlib import Path
import os
import shutil
import warnings
import platformdirs

_CACHE_WARNING_BYTES = 5 * 1024**3
_cache_checked = False


def _cache_directory() -> Path:
    return platformdirs.user_data_path(appname="trails", appauthor="pylca") / "cache"


def _cache_size(path: Path) -> int:
    """Estimate allocated bytes using metadata only, without following symlinks."""
    if path.is_symlink():
        return 0
    total = 0
    pending = [path]
    while pending:
        with os.scandir(pending.pop()) as entries:
            for entry in entries:
                if entry.is_dir(follow_symlinks=False):
                    pending.append(Path(entry.path))
                elif entry.is_file(follow_symlinks=False):
                    info = entry.stat(follow_symlinks=False)
                    blocks = getattr(info, "st_blocks", None)
                    total += info.st_size if blocks is None else blocks * 512
    return total


def _warn_if_cache_large() -> None:
    """Check once per session; an unavailable cache must not prevent import."""
    global _cache_checked
    if _cache_checked:
        return
    _cache_checked = True
    try:
        size = _cache_size(_cache_directory())
    except OSError:
        return
    if size > _CACHE_WARNING_BYTES:
        warnings.warn(
            f"TRAILS interpolation cache uses {size / 1024**3:.1f} GiB. "
            "When no TRAILS calculations are running, use trails.clear_cache() "
            "to reclaim space. Cached data will be rebuilt when needed.",
            UserWarning,
            stacklevel=2,
        )


def clear_cache() -> Path:
    """Remove rebuildable interpolation caches when no TRAILS jobs are running.

    :returns: The empty cache directory.
    :rtype: Path
    :raises OSError: If deletion or directory creation fails. Partial deletion
        is possible; errors are not silently ignored.
    """
    cache_dir = _cache_directory()
    if cache_dir.exists():
        shutil.rmtree(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir
