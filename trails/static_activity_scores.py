from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse.linalg import spsolve

from .cache_interpolation import cache_dir_for_package
from .characterization import get_cf_matrix
from .lca import _build_direct_technosphere_for_year
from .lcia import _get_lcia_methods_filepath

if TYPE_CHECKING:
    from .trails import Trails


@dataclass
class StaticActivityScores:
    """Cached static activity impact intensities for scenario years."""

    methods: tuple[str, ...]
    years: np.ndarray
    scores: np.ndarray
    year_index: dict[int, int]
    method_index: dict[str, int]
    cache_path: Path | None = None
    loaded_from_cache: bool = False


def _hash_array(hasher: hashlib._Hash, array: np.ndarray) -> None:
    arr = np.ascontiguousarray(array)
    hasher.update(str(arr.dtype).encode("utf-8"))
    hasher.update(json.dumps(arr.shape).encode("utf-8"))
    hasher.update(memoryview(arr).cast("B"))


def _matrix_fingerprint(trails: Trails) -> str:
    """Build a content fingerprint for the current in-memory A/B matrices."""
    if trails.A is None or trails.B is None:
        raise RuntimeError("Cannot fingerprint activity scores without A and B.")

    cache_state = getattr(trails, "_static_activity_score_fingerprint", None)
    matrix_ids = (
        id(trails.A),
        id(trails.B),
        tuple(int(x) for x in trails.A.shape),
        tuple(int(x) for x in trails.B.shape),
        int(trails.A.nnz),
        int(trails.B.nnz),
    )
    if cache_state is not None and cache_state[0] == matrix_ids:
        return str(cache_state[1])

    hasher = hashlib.blake2b(digest_size=20)
    for name, matrix in (("A", trails.A), ("B", trails.B)):
        hasher.update(name.encode("utf-8"))
        hasher.update(json.dumps(tuple(int(x) for x in matrix.shape)).encode("utf-8"))
        hasher.update(str(matrix.coords.dtype).encode("utf-8"))
        hasher.update(str(matrix.data.dtype).encode("utf-8"))
        _hash_array(hasher, matrix.coords)
        _hash_array(hasher, matrix.data)

    digest = hasher.hexdigest()
    setattr(trails, "_static_activity_score_fingerprint", (matrix_ids, digest))
    return digest


def _lcia_fingerprint(ei_version: str) -> dict[str, int | str]:
    path = _get_lcia_methods_filepath(ei_version)
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def _normalise_methods(methods: str | list[str] | tuple[str, ...]) -> tuple[str, ...]:
    if isinstance(methods, str):
        methods = (methods,)
    out = tuple(str(method) for method in methods)
    if not out:
        raise ValueError("At least one adaptive LCIA method is required.")
    return out


def _score_cache_payload(
    trails: Trails,
    *,
    methods: tuple[str, ...],
    ei_version: str,
    years: np.ndarray,
) -> dict:
    return {
        "version": 1,
        "methods": list(methods),
        "ei_version": str(ei_version),
        "years": [int(year) for year in years],
        "value_dtype": str(getattr(trails, "value_dtype", "")),
        "index_dtype": str(getattr(trails, "index_dtype", "")),
        "matrix_fingerprint": _matrix_fingerprint(trails),
        "lcia_fingerprint": _lcia_fingerprint(str(ei_version)),
    }


def _score_cache_key(payload: dict) -> str:
    blob = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:24]


def _score_cache_dir(trails: Trails) -> Path:
    return (
        cache_dir_for_package(
            trails.package,
            value_dtype=str(getattr(trails, "value_dtype", "")),
            index_dtype=str(getattr(trails, "index_dtype", "")),
            interpolation_start_year_offset=int(
                getattr(trails, "interpolation_start_year_offset", -1)
            ),
            interpolation_end_year_offset=int(
                getattr(trails, "interpolation_end_year_offset", 1)
            ),
        )
        / "static_activity_scores"
    )


def _reference_product_arrays(A_csc) -> tuple[np.ndarray, np.ndarray]:
    n_activities = int(A_csc.shape[1])
    product_indices = np.full(n_activities, -1, dtype=np.int64)
    production_values = np.full(n_activities, np.nan, dtype=np.float64)

    indptr = A_csc.indptr
    indices = A_csc.indices
    data = np.asarray(A_csc.data, dtype=np.float64)
    for activity_index in range(n_activities):
        start = int(indptr[activity_index])
        end = int(indptr[activity_index + 1])
        rows = indices[start:end]
        vals = data[start:end]
        if rows.size == 0:
            continue

        product_index: int | None = None
        production_value = 0.0
        if activity_index < int(A_csc.shape[0]):
            matches = np.where(rows == activity_index)[0]
            if matches.size:
                pos = int(matches[0])
                product_index = int(activity_index)
                production_value = float(vals[pos])

        if product_index is None or production_value == 0.0:
            pos = int(np.argmin(np.abs(np.abs(vals) - 1.0)))
            product_index = int(rows[pos])
            production_value = float(vals[pos])

        product_indices[activity_index] = int(product_index)
        production_values[activity_index] = float(production_value)

    return product_indices, production_values


def _activity_scores_from_product_intensities(
    intensities: np.ndarray,
    product_indices: np.ndarray,
    production_values: np.ndarray,
) -> np.ndarray:
    scores = np.full(
        (int(intensities.shape[1]), int(product_indices.size)),
        np.nan,
        dtype=np.float64,
    )
    valid = (product_indices >= 0) & (product_indices < intensities.shape[0])
    signs = np.where(production_values[valid] < 0.0, -1.0, 1.0)
    scores[:, valid] = intensities[product_indices[valid], :].T * signs[None, :]
    return scores


def _compute_static_activity_scores(
    trails: Trails,
    *,
    methods: tuple[str, ...],
    ei_version: str,
    years: np.ndarray,
) -> np.ndarray:
    if trails.A is None or trails.B is None:
        raise RuntimeError("Cannot compute static activity scores without A and B.")

    cf_matrix = get_cf_matrix(
        trails=trails,
        methods=list(methods),
        char_cache={},
        ei_version=str(ei_version),
    )
    n_methods = int(cf_matrix.shape[0])
    n_activities = int(trails.A.shape[1])
    scores = np.empty((n_methods, int(years.size), n_activities), dtype=np.float64)

    direct_matrix_cache = {}
    for year_pos, year in enumerate(years):
        context = trails._get_scenario_context(int(year))
        if context is None:
            raise RuntimeError(f"No scenario context available for year={int(year)}.")
        scenario_year, _label, t = context
        A_csc, _product_dict, _ref_cache = _build_direct_technosphere_for_year(
            trails=trails,
            year=int(scenario_year),
            cache=direct_matrix_cache,
        )
        B_t = trails.B[int(t), :, :]
        direct_scores = np.vstack(
            [
                np.asarray(B_t @ cf_matrix[method_pos, :], dtype=np.float64).reshape(-1)
                for method_pos in range(n_methods)
            ]
        )
        intensities = np.asarray(
            spsolve(A_csc.T.tocsc(), direct_scores.T),
            dtype=np.float64,
        )
        if intensities.ndim == 1:
            intensities = intensities[:, None]
        product_indices, production_values = _reference_product_arrays(A_csc)
        scores[:, year_pos, :] = _activity_scores_from_product_intensities(
            intensities,
            product_indices,
            production_values,
        )

    return scores


def _load_static_activity_scores(path: Path, payload: dict) -> StaticActivityScores:
    with np.load(path, allow_pickle=False) as data:
        years = np.asarray(data["years"], dtype=np.int64)
        methods = tuple(str(method) for method in data["methods"].tolist())
        scores = np.asarray(data["scores"], dtype=np.float64)

    expected_years = np.asarray(payload["years"], dtype=np.int64)
    expected_methods = tuple(str(method) for method in payload["methods"])
    if methods != expected_methods or not np.array_equal(years, expected_years):
        raise ValueError("Static activity score cache metadata mismatch.")

    return StaticActivityScores(
        methods=methods,
        years=years,
        scores=scores,
        year_index={int(year): pos for pos, year in enumerate(years)},
        method_index={method: pos for pos, method in enumerate(methods)},
        cache_path=path,
        loaded_from_cache=True,
    )


def _save_static_activity_scores(
    path: Path,
    *,
    payload: dict,
    scores: StaticActivityScores,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    meta_path = path.with_suffix(".json")
    np.savez_compressed(
        path,
        methods=np.asarray(scores.methods, dtype=str),
        years=np.asarray(scores.years, dtype=np.int64),
        scores=np.asarray(scores.scores, dtype=np.float64),
    )
    with meta_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def _ensure_static_activity_scores(
    trails: Trails,
    *,
    methods: str | list[str] | tuple[str, ...],
    ei_version: str = "3.11",
    years: np.ndarray | list[int] | None = None,
    use_cache: bool = True,
) -> StaticActivityScores:
    method_tuple = _normalise_methods(methods)
    if years is None:
        years_array = np.asarray(
            sorted(int(label) for label in trails.scenario_labels),
            dtype=np.int64,
        )
    else:
        years_array = np.asarray(sorted({int(year) for year in years}), dtype=np.int64)
    if years_array.size == 0:
        raise ValueError("At least one scenario year is required.")

    payload = _score_cache_payload(
        trails,
        methods=method_tuple,
        ei_version=str(ei_version),
        years=years_array,
    )
    key = _score_cache_key(payload)

    memory_cache = getattr(trails, "_static_activity_score_cache", None)
    if memory_cache is None:
        memory_cache = {}
        setattr(trails, "_static_activity_score_cache", memory_cache)
    if key in memory_cache:
        return memory_cache[key]

    cache_path = _score_cache_dir(trails) / f"activity_scores_{key}.npz"
    if use_cache and cache_path.exists():
        try:
            loaded = _load_static_activity_scores(cache_path, payload)
            memory_cache[key] = loaded
            return loaded
        except Exception:
            pass

    scores_array = _compute_static_activity_scores(
        trails,
        methods=method_tuple,
        ei_version=str(ei_version),
        years=years_array,
    )
    result = StaticActivityScores(
        methods=method_tuple,
        years=years_array,
        scores=scores_array,
        year_index={int(year): pos for pos, year in enumerate(years_array)},
        method_index={method: pos for pos, method in enumerate(method_tuple)},
        cache_path=cache_path if use_cache else None,
        loaded_from_cache=False,
    )

    if use_cache:
        _save_static_activity_scores(cache_path, payload=payload, scores=result)

    memory_cache[key] = result
    return result


def _activity_score_potential(
    scores: StaticActivityScores,
    *,
    year: int,
    activity: int,
    amount: float,
) -> tuple[float, dict[str, float]]:
    year_pos = scores.year_index[int(year)]
    act = int(activity)
    amount_abs = abs(float(amount))
    values = np.nan_to_num(
        scores.scores[:, year_pos, act],
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )
    potentials = amount_abs * np.abs(values)
    by_method = {
        method: float(potentials[pos]) for pos, method in enumerate(scores.methods)
    }
    return float(potentials.max(initial=0.0)), by_method
