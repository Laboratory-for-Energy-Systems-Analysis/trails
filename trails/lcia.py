from __future__ import annotations

from collections.abc import Mapping
import json
import os
from pathlib import Path
import time
import warnings
from typing import Any, Literal, TYPE_CHECKING

import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix

from .filesystem_constants import DATA_DIR

import logging

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    import xarray as xr

    from .trails import Trails

LCIA_METHODS_EI310 = DATA_DIR / "lcia_ei310.json"
LCIA_METHODS_EI311 = DATA_DIR / "lcia_ei311.json"
LCIA_METHODS_EI312 = DATA_DIR / "lcia_ei312.json"


_LCIA_METHODS_CACHE = {}


MethodBackend = Literal["auto", "regular", "edges"]


def _method_label(method: Any) -> str:
    if isinstance(method, Mapping):
        name = method.get("name")
        if isinstance(name, (list, tuple)):
            return " - ".join(str(value) for value in name)
        if name:
            return str(name)
    return str(method)


def resolve_lcia_methods(
    trails: Trails,
    methods: list[Any] | None,
    method_backend: MethodBackend | None,
) -> tuple[list[Any], Literal["regular", "edges"]]:
    """Resolve call-level methods and select a characterization backend."""
    call_has_methods = methods is not None
    selected = list(methods) if methods is not None else None
    backend = method_backend
    if selected is None:
        selected = list(getattr(trails, "default_methods", None) or [])
        if backend is None:
            backend = getattr(trails, "default_method_backend", "auto")
        if not selected:
            legacy_edges = getattr(trails, "default_edges_methods", None)
            if legacy_edges:
                selected = list(legacy_edges)
                backend = "edges"
    elif backend is None:
        backend = "auto"

    if not selected:
        raise ValueError(
            "No LCIA methods configured. Pass methods=... to lcia() or Trails()."
        )
    if backend not in {"auto", "regular", "edges"}:
        raise ValueError("method_backend must be one of {'auto', 'regular', 'edges'}")

    if backend == "auto":
        mappings = [isinstance(method, Mapping) for method in selected]
        strings = [isinstance(method, str) for method in selected]
        if all(mappings):
            backend = "edges"
        elif all(strings):
            backend = "regular"
        else:
            source = "lcia()" if call_has_methods else "Trails()"
            raise ValueError(
                f"Cannot infer one LCIA backend for methods supplied to {source}. "
                "Use a homogeneous method list or set method_backend explicitly."
            )

    if backend == "regular" and not all(isinstance(method, str) for method in selected):
        raise TypeError("Regular LCIA methods must be method-name strings.")
    return selected, backend


def lcia(
    trails: Trails,
    methods: list[Any] | None = None,
    *,
    method_backend: MethodBackend | None = None,
    ei_version: str | None = None,
    additional_topologies: dict[str, Any] | None = None,
    strategies: list[str] | None = None,
    reuse_mappings: bool = True,
    reuse_cached_cfs: bool | None = None,
    show_progress: bool = True,
    result_name: str | None = None,
) -> xr.DataArray:
    """Characterize a finalized temporal inventory without rerunning LCI."""
    if getattr(trails, "inventory", None) is None:
        raise RuntimeError(
            "Temporal inventory not initialized; run trails.lci() before "
            "trails.lcia()."
        )
    if reuse_cached_cfs is not None:
        warnings.warn(
            "reuse_cached_cfs is deprecated; use reuse_mappings. Cached EDGES "
            "templates are reevaluated for each inventory year.",
            FutureWarning,
            stacklevel=2,
        )
        reuse_mappings = bool(reuse_cached_cfs)

    selected, backend = resolve_lcia_methods(trails, methods, method_backend)
    version = str(
        ei_version
        if ei_version is not None
        else getattr(
            trails, "default_ei_version", getattr(trails, "ei_version", "3.11")
        )
    )
    started = time.perf_counter()

    if backend == "regular":
        from .characterization import (
            build_characterized_inventory,
            score_inventory_with_regular_methods,
        )
        from .lca import _CHAR_CACHE

        regular_methods = [str(method) for method in selected]
        scores = score_inventory_with_regular_methods(
            trails,
            regular_methods,
            _CHAR_CACHE,
            ei_version=version,
            show_progress=show_progress,
        )
        characterized = build_characterized_inventory(
            trails,
            regular_methods,
            _CHAR_CACHE,
            ei_version=version,
        )
    else:
        from .edges_matrix import score_inventory_with_edges

        trails.characterized_inventory = None
        scores = score_inventory_with_edges(
            trails,
            selected,
            additional_topologies=additional_topologies,
            strategies=strategies,
            reuse_cached_cfs=bool(reuse_mappings),
            show_progress=show_progress,
        )
        characterized = None

    labels = [_method_label(method) for method in selected]
    key = result_name or f"{backend}:" + " | ".join(labels)
    elapsed = float(time.perf_counter() - started)
    result = {
        "backend": backend,
        "methods": list(selected),
        "method_labels": labels,
        "ei_version": version if backend == "regular" else None,
        "scores": scores,
        "characterized_inventory": characterized,
        "seconds": elapsed,
    }
    trails.lcia_results[key] = result
    trails.current_lcia_result = key
    trails.scores = scores
    trails.characterized_inventory = characterized
    trails.lcia_diagnostics = {
        "backend": backend,
        "methods": labels,
        "total_seconds": elapsed,
        "reuse_mappings": bool(reuse_mappings) if backend == "edges" else None,
    }
    return scores


def _get_lcia_methods_filepath(ei_version: str = "3.11") -> Path:
    """Return the LCIA method JSON file for an ecoinvent version."""
    version = str(ei_version)
    if version == "3.10":
        return LCIA_METHODS_EI310
    if version == "3.11":
        return LCIA_METHODS_EI311
    if version == "3.12":
        env_path = os.environ.get("TRAILS_LCIA_EI312_JSON")
        if env_path:
            path = Path(env_path).expanduser()
            if path.exists():
                return path
            raise FileNotFoundError(
                "TRAILS_LCIA_EI312_JSON points to a missing file: " f"{path}"
            )
        if LCIA_METHODS_EI312.exists():
            return LCIA_METHODS_EI312
        raise FileNotFoundError(
            "No ecoinvent 3.12 LCIA method JSON found. Set "
            "TRAILS_LCIA_EI312_JSON to a compatible lcia_ei312.json file."
        )

    candidate = Path(version).expanduser()
    if candidate.exists():
        return candidate

    raise ValueError(
        f"Unsupported ecoinvent LCIA version: {ei_version!r}. "
        "Supported versions are '3.10', '3.11', and '3.12' "
        "(with TRAILS_LCIA_EI312_JSON set)."
    )


def get_lcia_method_names(ei_version: str = "3.11") -> list[str]:
    """Get lcia method names.

    :param ei_version: Value for `ei_version`.
    :type ei_version: str
    :returns: Return value.
    :rtype: list[str]"""

    filepath = _get_lcia_methods_filepath(ei_version)

    with open(filepath, "r") as f:
        data = json.load(f)

    return [" - ".join(x["name"]) for x in data]


def format_lcia_method_exchanges(
    method: dict,
) -> dict[tuple[str, str, str], float]:
    """Format lcia method exchanges.

    :param method: Value for `method`.
    :type method: dict
    :returns: Return value.
    :rtype: dict[tuple[str, str, str], float]"""

    return {
        (
            x["name"],
            x["categories"][0],
            x["categories"][1] if len(x["categories"]) > 1 else "unspecified",
        ): x["amount"]
        for x in method["exchanges"]
    }


def get_lcia_methods(
    methods: list[str] | None = None, ei_version: str = "3.11"
) -> dict[str, dict[tuple[str, str, str], float]]:
    """Get lcia methods.

    :param methods: Value for `methods`.
    :type methods: list[str] | None
    :param ei_version: Value for `ei_version`.
    :type ei_version: str
    :returns: Return value.
    :rtype: dict[str, dict[tuple[str, str, str], float]]"""
    filepath = _get_lcia_methods_filepath(ei_version)
    key = (str(filepath), tuple(methods) if methods else None)
    if key in _LCIA_METHODS_CACHE:
        return _LCIA_METHODS_CACHE[key]

    with open(filepath, "r") as f:
        data = json.load(f)

    if methods:
        data = [x for x in data if " - ".join(x["name"]) in methods]

    out = {" - ".join(x["name"]): format_lcia_method_exchanges(x) for x in data}
    _LCIA_METHODS_CACHE[key] = out
    return out


def fill_characterization_factors_matrices(
    methods: list[str],
    biosphere_matrix_dict: dict[int, int],
    biosphere_dict: dict[tuple[str, str, str], int],
    debug: bool = False,
    ei_version: str = "3.11",
) -> csr_matrix:
    """Fill characterization factors matrices.

    :param methods: Value for `methods`.
    :type methods: list[str]
    :param biosphere_matrix_dict: Value for `biosphere_matrix_dict`.
    :type biosphere_matrix_dict: dict[int, int]
    :param biosphere_dict: Value for `biosphere_dict`.
    :type biosphere_dict: dict[tuple[str, str, str], int]
    :param debug: Value for `debug`.
    :type debug: bool
    :returns: Return value.
    :rtype: csr_matrix"""

    if debug:
        logger.info(
            "LCIA: building CF matrix for methods=%d biosphere_flows=%d",
            len(methods),
            len(biosphere_matrix_dict),
        )

    lcia_data = get_lcia_methods(methods=methods, ei_version=ei_version)

    # Prepare data for efficient creation of the sparse matrix
    data = []
    rows = []
    cols = []
    cfs = []

    matched, unmatched = 0, 0

    for m, method in enumerate(methods):
        method_data = lcia_data[method]

        for flow_name in method_data:
            if flow_name in biosphere_dict:
                idx = biosphere_dict[flow_name]
                if idx in biosphere_matrix_dict:
                    data.append(method_data[flow_name])
                    rows.append(biosphere_matrix_dict[idx])
                    cols.append(m)
                    cfs.append((method, flow_name, idx, method_data[flow_name]))
                    matched += 1
                else:
                    unmatched += 1
            else:
                unmatched += 1

    if debug:
        logger.info("LCIA: CF nonzeros=%d", len(data))
        if len(data) == 0:
            logger.warning(
                "LCIA: CF matrix has zero entries -> all scores will be zero."
            )

        logger.info("LCIA: matched_flows=%d unmatched_flows=%d", matched, unmatched)

    # Efficiently create the sparse matrix
    matrix = sparse.csr_matrix(
        (data, (cols, rows)),
        shape=(len(methods), len(biosphere_matrix_dict)),
        dtype=np.float64,
    )

    if debug:
        # sort l by method and flow
        cfs = sorted(cfs, key=lambda x: (x[0], x[1]))
        for x in cfs:
            method, flow, f, value = x
            logger.info(
                f"LCIA method: {method}, Flow: {flow}, Index: {f}, Value: {value}"
            )

    return matrix
