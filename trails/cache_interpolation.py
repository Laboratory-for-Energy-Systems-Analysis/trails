from __future__ import annotations

import hashlib
import json
import pickle
import tempfile
from pathlib import Path
from typing import Any, Optional

import platformdirs
import sparse


def _uses_ephemeral_datapackage_base(package: Any) -> bool:
    """Return whether Frictionless unpacked the package into its temp directory."""
    base_path = getattr(package, "base_path", None)
    if not base_path:
        return False
    try:
        resolved = Path(base_path).resolve()
        temp_root = Path(tempfile.gettempdir()).resolve()
    except (OSError, TypeError, ValueError):
        return False
    return resolved.name.endswith("-datapackage") and resolved.parent == temp_root


def _cache_key(
    package: Any,
    *,
    value_dtype: str,
    index_dtype: str,
    interpolate_annual: bool = True,
    interpolation_start_year_offset: int = -1,
    interpolation_end_year_offset: int = 1,
) -> str:
    """cache key.

    :param package: Value for `package`.
    :type package: Any
    :param value_dtype: Value for `value_dtype`.
    :type value_dtype: str
    :param index_dtype: Value for `index_dtype`.
    :type index_dtype: str
    :param interpolate_annual: Value for `interpolate_annual`.
    :type interpolate_annual: bool
    :returns: Return value.
    :rtype: str"""
    try:
        package_desc = getattr(package, "descriptor", {})
    except Exception:
        package_desc = {}
    omit_resource_mtime = _uses_ephemeral_datapackage_base(package)
    file_fingerprint: dict[str, dict[str, int]] = {}
    try:
        tracked_names = {
            "A",
            "B",
            "A_indices",
            "B_indices",
            "A_matrix.csv",
            "B_matrix.csv",
            "A_matrix_index.csv",
            "B_matrix_index.csv",
        }
        tracked_filenames = {
            "A_matrix.csv",
            "B_matrix.csv",
            "A_matrix_index.csv",
            "B_matrix_index.csv",
        }
        for res in getattr(package, "resources", []):
            name = getattr(res, "name", "") or ""
            resource_desc = getattr(res, "descriptor", {}) or {}
            path_label = str(resource_desc.get("path") or "")
            filename = Path(path_label).name
            if name not in tracked_names and filename not in tracked_filenames:
                continue
            try:
                path = Path(res.source)
            except Exception:
                continue
            if not path.exists():
                continue
            stat = path.stat()
            fingerprint = {"size": int(stat.st_size)}
            if not omit_resource_mtime:
                fingerprint["mtime_ns"] = int(stat.st_mtime_ns)
            file_fingerprint[path_label or str(path)] = fingerprint
    except Exception:
        file_fingerprint = {}
    payload = {
        "descriptor": package_desc,
        "file_fingerprint": file_fingerprint,
        "value_dtype": value_dtype,
        "index_dtype": index_dtype,
        "interpolate_annual": interpolate_annual,
        "interpolation_start_year_offset": interpolation_start_year_offset,
        "interpolation_end_year_offset": interpolation_end_year_offset,
    }
    blob = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def cache_dir_for_package(
    package: Any,
    *,
    value_dtype: str,
    index_dtype: str,
    interpolation_start_year_offset: int = -1,
    interpolation_end_year_offset: int = 1,
) -> Path:
    """Cache dir for package.

    :param package: Value for `package`.
    :type package: Any
    :param value_dtype: Value for `value_dtype`.
    :type value_dtype: str
    :param index_dtype: Value for `index_dtype`.
    :type index_dtype: str
    :returns: Return value.
    :rtype: Path"""
    cache_base = platformdirs.user_data_path(appname="trails", appauthor="pylca")
    key = _cache_key(
        package,
        value_dtype=value_dtype,
        index_dtype=index_dtype,
        interpolate_annual=True,
        interpolation_start_year_offset=interpolation_start_year_offset,
        interpolation_end_year_offset=interpolation_end_year_offset,
    )
    short_key = key[:12]
    return cache_base / "cache" / f"interp_{short_key}"


def load_cached_interpolation(
    package: Any,
    *,
    value_dtype: str,
    index_dtype: str,
    interpolation_start_year_offset: int = -1,
    interpolation_end_year_offset: int = 1,
) -> tuple[
    Optional[sparse.COO],
    Optional[sparse.COO],
    list[str],
    list[str],
    Optional[dict],
    Optional[dict],
    Optional[dict],
    Path,
]:
    """Load cached interpolation.

    :param package: Value for `package`.
    :type package: Any
    :param value_dtype: Value for `value_dtype`.
    :type value_dtype: str
    :param index_dtype: Value for `index_dtype`.
    :type index_dtype: str
    :returns: Return value.
    :rtype: tuple[Optional[sparse.COO], Optional[sparse.COO], list[str], list[str], Optional[dict], Optional[dict], Optional[dict], Path]
    """
    cache_dir = cache_dir_for_package(
        package,
        value_dtype=value_dtype,
        index_dtype=index_dtype,
        interpolation_start_year_offset=interpolation_start_year_offset,
        interpolation_end_year_offset=interpolation_end_year_offset,
    )
    meta_path = cache_dir / "meta.json"
    a_path = cache_dir / "A.npz"
    b_path = cache_dir / "B.npz"
    td_path = cache_dir / "temporal.pkl"
    idx_path = cache_dir / "indices.pkl"
    if not (meta_path.exists() and a_path.exists() and b_path.exists()):
        return None, None, [], [], None, None, None, cache_dir
    try:
        with open(meta_path, "r", encoding="utf-8") as f:
            meta = json.load(f)
        if (
            meta.get("version") != 5
            or meta.get("value_dtype") != value_dtype
            or meta.get("index_dtype") != index_dtype
            or meta.get("interpolation_start_year_offset", -1)
            != interpolation_start_year_offset
            or meta.get("interpolation_end_year_offset", 1)
            != interpolation_end_year_offset
        ):
            return None, None, [], [], None, None, None, cache_dir
        A = sparse.load_npz(a_path)
        B = sparse.load_npz(b_path)
        labels = list(meta.get("scenario_labels", []))
        template_labels = list(meta.get("template_labels", []))
        temporal_tech = None
        temporal_bio = None
        indices = None
        if td_path.exists():
            with open(td_path, "rb") as f:
                payload = pickle.load(f)
            temporal_tech = payload.get("temporal_technosphere_exchanges")
            temporal_bio = payload.get("temporal_biosphere_exchanges")
        if idx_path.exists():
            with open(idx_path, "rb") as f:
                indices = pickle.load(f)
        return (
            A,
            B,
            labels,
            template_labels,
            temporal_tech,
            temporal_bio,
            indices,
            cache_dir,
        )
    except Exception:
        return None, None, [], [], None, None, None, cache_dir


def save_cached_interpolation(
    package: Any,
    *,
    value_dtype: str,
    index_dtype: str,
    interpolation_start_year_offset: int = -1,
    interpolation_end_year_offset: int = 1,
    A: sparse.COO,
    B: sparse.COO,
    scenario_labels: list[str],
    template_labels: list[str],
    temporal_technosphere_exchanges: dict,
    temporal_biosphere_exchanges: dict,
    activity_indices: dict,
    biosphere_indices: dict,
) -> Path:
    """Save cached interpolation.

    :param package: Value for `package`.
    :type package: Any
    :param value_dtype: Value for `value_dtype`.
    :type value_dtype: str
    :param index_dtype: Value for `index_dtype`.
    :type index_dtype: str
    :param A: Value for `A`.
    :type A: sparse.COO
    :param B: Value for `B`.
    :type B: sparse.COO
    :param scenario_labels: Value for `scenario_labels`.
    :type scenario_labels: list[str]
    :param template_labels: Value for `template_labels`.
    :type template_labels: list[str]
    :param temporal_technosphere_exchanges: Value for `temporal_technosphere_exchanges`.
    :type temporal_technosphere_exchanges: dict
    :param temporal_biosphere_exchanges: Value for `temporal_biosphere_exchanges`.
    :type temporal_biosphere_exchanges: dict
    :param activity_indices: Value for `activity_indices`.
    :type activity_indices: dict
    :param biosphere_indices: Value for `biosphere_indices`.
    :type biosphere_indices: dict
    :returns: Return value.
    :rtype: Path"""
    cache_dir = cache_dir_for_package(
        package,
        value_dtype=value_dtype,
        index_dtype=index_dtype,
        interpolation_start_year_offset=interpolation_start_year_offset,
        interpolation_end_year_offset=interpolation_end_year_offset,
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(cache_dir / "A.npz", A, compressed=False)
    sparse.save_npz(cache_dir / "B.npz", B, compressed=False)
    meta = {
        "version": 5,
        "value_dtype": value_dtype,
        "index_dtype": index_dtype,
        "interpolation_start_year_offset": interpolation_start_year_offset,
        "interpolation_end_year_offset": interpolation_end_year_offset,
        "scenario_labels": scenario_labels,
        "template_labels": template_labels,
        "compressed": False,
    }
    with open(cache_dir / "meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f)
    td_path = cache_dir / "temporal.pkl"
    with open(td_path, "wb") as f:
        pickle.dump(
            {
                "temporal_technosphere_exchanges": temporal_technosphere_exchanges,
                "temporal_biosphere_exchanges": temporal_biosphere_exchanges,
            },
            f,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    idx_path = cache_dir / "indices.pkl"
    with open(idx_path, "wb") as f:
        pickle.dump(
            {
                "activity_indices": activity_indices,
                "biosphere_indices": biosphere_indices,
            },
            f,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    return cache_dir
