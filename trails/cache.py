from __future__ import annotations

from pathlib import Path
import shutil
import platformdirs


def clear_cache() -> Path:
    """Remove all cached interpolation artifacts and return the cache path."""
    base = platformdirs.user_data_path(appname="trails", appauthor="pylca")
    cache_dir = base / "cache"
    if cache_dir.exists():
        shutil.rmtree(cache_dir, ignore_errors=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir
