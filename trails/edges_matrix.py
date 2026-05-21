"""EDGES integration helpers for edge-level characterization factors."""

from __future__ import annotations

from contextlib import nullcontext, redirect_stdout
from io import StringIO
from types import SimpleNamespace
from typing import Any, Mapping

import numpy as np
from scipy import sparse as sp
import sparse
from tqdm.auto import tqdm
import xarray as xr


class _IdentityMappedDict(dict):
    """Brightway-like matrix dictionary with a ``reversed`` lookup."""

    @property
    def reversed(self) -> dict[int, int]:
        return {int(v): int(k) for k, v in self.items()}


class _TrailsLCAAdapter:
    """Minimal LCA object exposing the attributes EDGES needs for matching."""

    def __init__(
        self,
        *,
        inventory: sp.csr_matrix,
        technosphere_matrix: sp.csr_matrix,
    ) -> None:
        self.inventory = inventory
        self.technosphere_matrix = technosphere_matrix
        self.supply_array = np.ones(inventory.shape[1], dtype=float)

        n_flows, n_activities = inventory.shape
        n_products = technosphere_matrix.shape[0]
        activity = _IdentityMappedDict({i: i for i in range(int(n_activities))})
        biosphere = _IdentityMappedDict({i: i for i in range(int(n_flows))})
        product = _IdentityMappedDict({i: i for i in range(int(n_products))})
        self.dicts = SimpleNamespace(
            activity=activity,
            biosphere=biosphere,
            product=product,
        )


def _to_csr_2d(matrix: Any) -> sp.csr_matrix:
    """Convert a 2D sparse-like matrix to SciPy CSR."""
    if sp.issparse(matrix):
        return matrix.tocsr()
    if isinstance(matrix, sparse.COO):
        return matrix.to_scipy_sparse().tocsr()
    return sp.csr_matrix(np.asarray(matrix))


def _metadata_for_year(
    trails: Any,
    mapping_name: str,
    year: int,
) -> dict[int, dict[str, Any]]:
    """Return activity or biosphere metadata for a scenario/template year."""
    mappings = getattr(trails, mapping_name, None) or {}
    if not mappings:
        return {}

    labels = [str(year)]
    map_template = getattr(trails, "_map_year_to_template_year", None)
    if callable(map_template):
        try:
            labels.append(str(map_template(int(year))))
        except Exception:
            pass

    map_scenario = getattr(trails, "_map_year_to_scenario_year", None)
    if callable(map_scenario):
        try:
            labels.append(str(map_scenario(int(year))))
        except Exception:
            pass

    for label in labels:
        if label in mappings:
            return mappings[label]

    return next(iter(mappings.values()))


def _categories_from_biosphere_metadata(meta: Mapping[str, Any]) -> list[str]:
    categories: list[str] = []
    compartment = meta.get("compartment")
    subcompartment = meta.get("subcompartment")
    if compartment not in (None, "", "unspecified"):
        categories.append(str(compartment))
    if subcompartment not in (None, "", "unspecified"):
        categories.append(str(subcompartment))
    return categories


def _load_activity_classifications(
    trails: Any,
    positions: set[int] | None = None,
) -> dict[int, list[tuple[str, str]]]:
    """Load optional datapackage activity classifications by activity index."""
    package = getattr(trails, "package", None)
    if package is None:
        return {}

    descriptor = getattr(package, "descriptor", None) or {}
    custom = descriptor.get("custom") or {}
    resource_name = custom.get("classification_resource") or "classifications"

    resource = None
    try:
        resource = package.get_resource(resource_name)
    except Exception:
        for candidate in getattr(package, "resources", []) or []:
            desc = getattr(candidate, "descriptor", None) or {}
            if desc.get("name") == resource_name:
                resource = candidate
                break

    if resource is None:
        return {}

    try:
        rows = resource.read(keyed=True)
    except Exception:
        return {}

    out: dict[int, list[tuple[str, str]]] = {}
    for row in rows:
        try:
            idx = int(row["index"])
        except (KeyError, TypeError, ValueError):
            continue
        if positions is not None and idx not in positions:
            continue
        system = row.get("classification_system")
        code = row.get("classification_code")
        if system in (None, "") or code in (None, ""):
            continue
        out.setdefault(idx, []).append((str(system), str(code)))
    return out


def _biosphere_flows_for_edges(
    trails: Any,
    *,
    year: int,
    n_flows: int,
    positions: set[int] | None = None,
) -> list[dict[str, Any]]:
    metadata = _metadata_for_year(trails, "biosphere_indices", year)
    flows = []
    if positions is None:
        iterator = range(int(n_flows))
    else:
        iterator = (pos for pos in sorted(positions) if 0 <= pos < int(n_flows))

    for pos in iterator:
        meta = metadata.get(pos, {})
        flows.append(
            {
                "name": meta.get("name"),
                "categories": _categories_from_biosphere_metadata(meta),
                "unit": meta.get("unit"),
                "location": meta.get("location"),
                "classifications": meta.get("classifications"),
                "position": pos,
            }
        )
    return flows


def _technosphere_flows_for_edges(
    trails: Any,
    *,
    year: int,
    n_activities: int,
    positions: set[int] | None = None,
) -> list[dict[str, Any]]:
    metadata = _metadata_for_year(trails, "activity_indices", year)
    classifications = _load_activity_classifications(trails, positions=positions)

    flows = []
    if positions is None:
        iterator = range(int(n_activities))
    else:
        iterator = (pos for pos in sorted(positions) if 0 <= pos < int(n_activities))

    for pos in iterator:
        meta = metadata.get(pos, {})
        flows.append(
            {
                "name": meta.get("name"),
                "reference product": meta.get("reference product"),
                "unit": meta.get("unit"),
                "location": meta.get("location"),
                "classifications": meta.get(
                    "classifications", classifications.get(pos)
                ),
                "position": pos,
            }
        )
    return flows


def _position_lookup(flows: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    return {
        int(flow["position"]): {k: v for k, v in flow.items() if k != "position"}
        for flow in flows
        if "position" in flow
    }


def _hashable_metadata(value: Any) -> Any:
    if isinstance(value, Mapping):
        return tuple(
            sorted(
                (str(key), _hashable_metadata(val))
                for key, val in value.items()
                if val is not None
            )
        )
    if isinstance(value, np.ndarray):
        return tuple(_hashable_metadata(item) for item in value.tolist())
    if isinstance(value, set):
        return tuple(sorted((_hashable_metadata(item) for item in value), key=repr))
    if isinstance(value, (list, tuple)):
        return tuple(_hashable_metadata(item) for item in value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def _edge_signature(
    edge: tuple[int, int],
    biosphere_lookup: Mapping[int, Mapping[str, Any]],
    technosphere_lookup: Mapping[int, Mapping[str, Any]],
) -> tuple[Any, Any]:
    supplier, consumer = edge
    return (
        _hashable_metadata(biosphere_lookup.get(int(supplier), {})),
        _hashable_metadata(technosphere_lookup.get(int(consumer), {})),
    )


def _cache_cf_entries_by_signature(
    cache: dict[str, Any],
    entries: list[Mapping[str, Any]],
    biosphere_lookup: Mapping[int, Mapping[str, Any]],
    technosphere_lookup: Mapping[int, Mapping[str, Any]],
) -> set[tuple[int, int]]:
    entries_by_signature = cache.setdefault("entries_by_signature", {})
    processed: set[tuple[int, int]] = set()
    for entry in entries:
        if entry.get("direction") != "biosphere-technosphere":
            continue
        for supplier, consumer in entry.get("positions", ()):
            edge = (int(supplier), int(consumer))
            signature = _edge_signature(edge, biosphere_lookup, technosphere_lookup)
            entries_by_signature[signature] = [
                {
                    "direction": entry.get("direction"),
                    "value": entry.get("value"),
                }
            ]
            processed.add(edge)
    return processed


def _resolve_edges_parameters(edge_lca: Any, scenario_year: int) -> dict[str, Any]:
    scenario_name = getattr(edge_lca, "scenario", None)
    if scenario_name is None and isinstance(
        getattr(edge_lca, "parameters", None),
        dict,
    ):
        parameters = getattr(edge_lca, "parameters")
        if parameters:
            scenario_name = next(iter(parameters))

    resolver = getattr(edge_lca, "_resolve_parameters_for_scenario", None)
    if callable(resolver):
        return resolver(scenario_year, scenario_name)
    return {}


def _evaluate_cf_value(
    edge_lca: Any,
    raw_value: Any,
    *,
    resolved_params: Mapping[str, Any],
    scenario_year: int,
    cache: dict[Any, float],
) -> float:
    key = _hashable_metadata(raw_value)
    try:
        hash(key)
    except TypeError:
        key = repr(raw_value)
    if key in cache:
        return cache[key]

    evaluator = getattr(edge_lca, "_evaluate_cf_numeric_value", None)
    if callable(evaluator):
        value = evaluator(
            raw_value,
            resolved_params=resolved_params,
            scenario_idx=scenario_year,
        )
    else:
        try:
            value = float(raw_value)
        except Exception:
            value = 0.0

    cache[key] = float(value)
    return cache[key]


def _build_biosphere_matrix_from_edges_mapping(
    edge_lca: Any,
    *,
    scenario_year: int,
    shape: tuple[int, int],
    cached_entries_by_edge: Mapping[tuple[int, int], list[Mapping[str, Any]]],
    mapped_entries: list[Mapping[str, Any]],
) -> sp.csr_matrix:
    resolved_params = _resolve_edges_parameters(edge_lca, scenario_year)
    value_cache: dict[Any, float] = {}
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []

    def add_positions(entry: Mapping[str, Any], positions: Any) -> None:
        if entry.get("direction") != "biosphere-technosphere":
            return
        value = _evaluate_cf_value(
            edge_lca,
            entry.get("value"),
            resolved_params=resolved_params,
            scenario_year=scenario_year,
            cache=value_cache,
        )
        if value == 0.0:
            return
        for supplier, consumer in positions or ():
            supplier = int(supplier)
            consumer = int(consumer)
            if 0 <= supplier < shape[0] and 0 <= consumer < shape[1]:
                rows.append(supplier)
                cols.append(consumer)
                values.append(value)

    for edge, entries in cached_entries_by_edge.items():
        for entry in entries:
            add_positions(entry, (edge,))

    for entry in mapped_entries:
        add_positions(entry, entry.get("positions", ()))

    if not values:
        return sp.csr_matrix(shape, dtype=float)

    return sp.csr_matrix(
        (
            np.asarray(values, dtype=float),
            (
                np.asarray(rows, dtype=np.int64),
                np.asarray(cols, dtype=np.int64),
            ),
        ),
        shape=shape,
    )


def _edge_sets_for_lookup(edge_lca: Any) -> tuple[set[int], set[int], set[int]]:
    supplier_bio: set[int] = set()
    supplier_tech: set[int] = set()
    consumers: set[int] = set()

    for supplier, consumer in getattr(edge_lca, "biosphere_edges", None) or ():
        supplier_bio.add(int(supplier))
        consumers.add(int(consumer))

    for supplier, consumer in getattr(edge_lca, "technosphere_edges", None) or ():
        supplier_tech.add(int(supplier))
        consumers.add(int(consumer))

    return supplier_bio, supplier_tech, consumers


def _ensure_edges_available() -> type:
    try:
        from edges import EdgeLCIA
    except ImportError as exc:  # pragma: no cover - depends on optional package
        raise ImportError(
            "EDGES methods require the optional 'edges' package. Install EDGES "
            "and its solver dependencies before calling trails.lca(..., "
            "edges_methods=...)."
        ) from exc
    return EdgeLCIA


def get_edges_lcia_method_names() -> list[Any]:
    """Return EDGES LCIA method names, importing EDGES only on demand."""
    try:
        from edges import get_available_methods
    except ImportError as exc:  # pragma: no cover - depends on optional package
        raise ImportError(
            "Listing EDGES methods requires the optional 'edges' package."
        ) from exc
    return get_available_methods()


def _method_supplier_matrices(edge_lca: Any) -> set[str]:
    matrices = set()
    for cf in getattr(edge_lca, "raw_cfs_data", None) or ():
        supplier = cf.get("supplier", {}) or {}
        matrix = str(supplier.get("matrix") or "biosphere").strip().lower()
        matrices.add(matrix)
    return matrices or {"biosphere"}


def _build_edges_characterization_matrices_for_year(
    trails: Any,
    methods: list[Any],
    *,
    year: int,
    additional_topologies: dict[str, Any] | None = None,
    strategies: list[str] | None = None,
    biosphere_edges: set[tuple[int, int]] | None = None,
    mapping_caches: list[dict[str, Any]] | None = None,
    progress: Any | None = None,
    suppress_edges_output: bool = True,
) -> list[sp.csr_matrix]:
    """Build one EDGES CF matrix per method for a Trails scenario year.

    The returned matrices are shaped as EDGES biosphere inventories:
    ``(biosphere flow, consuming activity)``.
    """
    if getattr(trails, "A", None) is None or getattr(trails, "B", None) is None:
        raise ValueError(
            "Cannot build EDGES characterization matrices without A and B."
        )

    EdgeLCIA = _ensure_edges_available()

    context = trails._get_scenario_context(int(year))
    if context is None:
        raise ValueError(f"Scenario year {year!r} is not available in Trails.")
    scenario_year, _label, t = context

    b_shape = trails.B[int(t), :, :].shape
    n_activities, n_flows = int(b_shape[0]), int(b_shape[1])
    a_shape = trails.A[int(t), :, :].shape

    if biosphere_edges is None:
        biosphere_inventory = _to_csr_2d(trails.B[int(t), :, :]).T.tocsr()
        edge_pairs = {
            (int(flow), int(activity))
            for flow, activity in zip(*biosphere_inventory.nonzero())
        }
    else:
        edge_pairs = {
            (int(flow), int(activity))
            for flow, activity in biosphere_edges
            if 0 <= int(flow) < n_flows and 0 <= int(activity) < n_activities
        }
        biosphere_inventory = sp.csr_matrix((n_flows, n_activities), dtype=float)

    technosphere_matrix = sp.csr_matrix((int(a_shape[1]), int(a_shape[0])), dtype=float)
    lca_obj = _TrailsLCAAdapter(
        inventory=biosphere_inventory,
        technosphere_matrix=technosphere_matrix,
    )

    supplier_positions = {flow for flow, _activity in edge_pairs}
    consumer_positions = {activity for _flow, activity in edge_pairs}
    biosphere_flows = _biosphere_flows_for_edges(
        trails,
        year=int(scenario_year),
        n_flows=n_flows,
        positions=supplier_positions if biosphere_edges is not None else None,
    )
    technosphere_flows = _technosphere_flows_for_edges(
        trails,
        year=int(scenario_year),
        n_activities=n_activities,
        positions=consumer_positions if biosphere_edges is not None else None,
    )
    biosphere_lookup = _position_lookup(biosphere_flows)
    technosphere_lookup = _position_lookup(technosphere_flows)

    planes: list[sp.csr_matrix] = []
    for method_idx, method in enumerate(methods):
        mapping_cache = (
            mapping_caches[method_idx]
            if mapping_caches is not None and method_idx < len(mapping_caches)
            else None
        )
        cached_entries_by_edge: dict[tuple[int, int], list[Mapping[str, Any]]] = {}
        edges_to_map = set(edge_pairs)
        if mapping_cache is not None:
            entries_by_signature = mapping_cache.setdefault("entries_by_signature", {})
            miss_signatures = mapping_cache.setdefault("miss_signatures", set())
            edges_to_map = set()
            for edge in edge_pairs:
                signature = _edge_signature(
                    edge,
                    biosphere_lookup,
                    technosphere_lookup,
                )
                if signature in entries_by_signature:
                    cached_entries_by_edge[edge] = entries_by_signature[signature]
                elif signature in miss_signatures:
                    continue
                else:
                    edges_to_map.add(edge)

        if progress is not None:
            progress.set_postfix_str(
                f"year={int(scenario_year)}, method={method_idx + 1}/{len(methods)}, "
                f"edges={len(edge_pairs)}, new={len(edges_to_map)}"
            )
        edge_lca = EdgeLCIA(
            demand={},
            method=method,
            lca=lca_obj,
            additional_topologies=additional_topologies,
        )
        supplier_matrices = _method_supplier_matrices(edge_lca)
        if supplier_matrices != {"biosphere"}:
            raise NotImplementedError(
                "Trails currently supports EDGES methods whose supplier matrix "
                "is only 'biosphere'. Technosphere-supplier and mixed EDGES "
                "methods need technosphere edge inventory support."
            )

        edge_lca.biosphere_edges = edges_to_map
        edge_lca.technosphere_edges = set()
        edge_lca.biosphere_flows = biosphere_flows
        edge_lca.technosphere_flows = technosphere_flows
        edge_lca.position_to_technosphere_flows_lookup = technosphere_lookup
        edge_lca.reversed_activity = lca_obj.dicts.activity.reversed
        edge_lca.reversed_biosphere = lca_obj.dicts.biosphere.reversed
        edge_lca.cfs_mapping = []
        edge_lca._seen_positions = set()

        supplier_bio, supplier_tech, consumers = _edge_sets_for_lookup(edge_lca)
        stdout_context = (
            redirect_stdout(StringIO()) if suppress_edges_output else nullcontext()
        )
        with stdout_context:
            if edges_to_map:
                edge_lca._preprocess_lookups(
                    restrict_supplier_positions_bio=supplier_bio,
                    restrict_supplier_positions_tech=supplier_tech,
                    restrict_consumer_positions=consumers,
                )
                edge_lca.apply_strategies(strategies)
                if mapping_cache is not None:
                    new_entries = edge_lca.cfs_mapping
                    processed_edges = _cache_cf_entries_by_signature(
                        mapping_cache,
                        new_entries,
                        biosphere_lookup,
                        technosphere_lookup,
                    )
                    miss_signatures = mapping_cache.setdefault("miss_signatures", set())
                    for edge in edges_to_map - processed_edges:
                        miss_signatures.add(
                            _edge_signature(
                                edge,
                                biosphere_lookup,
                                technosphere_lookup,
                            )
                        )
        if progress is not None:
            progress.update(1)

        matrix = _build_biosphere_matrix_from_edges_mapping(
            edge_lca,
            scenario_year=int(scenario_year),
            shape=biosphere_inventory.shape,
            cached_entries_by_edge=cached_entries_by_edge,
            mapped_entries=edge_lca.cfs_mapping,
        )

        if matrix.shape != biosphere_inventory.shape:
            raise ValueError(
                "EDGES characterization matrix shape "
                f"{matrix.shape} does not match Trails biosphere inventory "
                f"shape {biosphere_inventory.shape}."
            )
        planes.append(matrix)

    return planes


def _lookup_sparse_values(
    matrix: sp.csr_matrix,
    rows: np.ndarray,
    cols: np.ndarray,
) -> np.ndarray:
    """Vectorized CSR lookup for paired row/column coordinates."""
    if rows.size == 0:
        return np.empty(0, dtype=float)
    matrix = matrix.tocsr(copy=False)
    matrix.sort_indices()

    values = np.zeros(rows.shape[0], dtype=float)
    for row in np.unique(rows):
        mask = rows == row
        start, stop = matrix.indptr[int(row)], matrix.indptr[int(row) + 1]
        row_cols = matrix.indices[start:stop]
        row_data = matrix.data[start:stop]
        if row_cols.size == 0:
            continue
        wanted = cols[mask]
        positions = np.searchsorted(row_cols, wanted)
        found = positions < row_cols.size
        if np.any(found):
            found_idx = np.flatnonzero(mask)[found]
            pos = positions[found]
            exact = row_cols[pos] == wanted[found]
            if np.any(exact):
                values[found_idx[exact]] = row_data[pos[exact]]
    return values


def _method_labels(methods: list[Any]) -> np.ndarray:
    labels: list[str] = []
    for method in methods:
        if isinstance(method, Mapping):
            labels.append(str(method.get("name") or repr(method)))
        elif isinstance(method, tuple):
            labels.append(" | ".join(str(part) for part in method))
        else:
            labels.append(str(method))
    return np.asarray(labels, dtype=object)


def score_inventory_with_edges(
    trails: Any,
    methods: list[Any],
    *,
    additional_topologies: dict[str, Any] | None = None,
    strategies: list[str] | None = None,
    reuse_cached_cfs: bool = True,
    show_progress: bool = True,
) -> xr.DataArray:
    """Score ``trails.inventory`` using EDGES edge-level CF matrices.

    ``reuse_cached_cfs`` reuses EDGES matched CF templates across scenario years
    when supplier and consumer metadata signatures are identical. Numeric CF
    values are still evaluated for each scenario year.
    """
    if not methods:
        raise ValueError("edges_methods must contain at least one EDGES method.")
    if getattr(trails, "inventory", None) is None:
        raise ValueError("EDGES scoring requires a finalized Trails inventory.")

    inventory = trails.inventory
    if not {"activity", "flow", "year"}.issubset(inventory.dims):
        raise ValueError(
            "EDGES scoring expects inventory dimensions including "
            "'activity', 'flow', and 'year'."
        )

    inv = inventory.transpose(
        "activity",
        "flow",
        "year",
        *[dim for dim in inventory.dims if dim not in {"activity", "flow", "year"}],
    )
    has_root = "root activity" in inv.dims
    data = inv.data
    if not isinstance(data, sparse.COO):
        data = sparse.COO.from_numpy(np.asarray(data))

    n_methods = len(methods)
    n_activities = int(inv.sizes["activity"])
    n_years = int(inv.sizes["year"])
    n_roots = int(inv.sizes["root activity"]) if has_root else 0
    years = np.asarray(inv.coords["year"].values, dtype=int)

    if data.nnz == 0:
        if has_root and n_methods > 1:
            arr = sparse.zeros(
                (n_methods, n_activities, n_years, n_roots),
                dtype=getattr(trails, "value_dtype", float),
            )
            scores = xr.DataArray(
                arr,
                dims=("method", "activity", "year", "root activity"),
                coords={
                    "method": _method_labels(methods),
                    "activity": np.arange(n_activities, dtype=int),
                    "year": years,
                    "root activity": inv.coords["root activity"].values,
                },
            )
        elif has_root:
            arr = sparse.zeros(
                (n_activities, n_years, n_roots),
                dtype=getattr(trails, "value_dtype", float),
            )
            scores = xr.DataArray(
                arr,
                dims=("activity", "year", "root activity"),
                coords={
                    "activity": np.arange(n_activities, dtype=int),
                    "year": years,
                    "root activity": inv.coords["root activity"].values,
                },
            )
        elif n_methods > 1:
            arr = sparse.zeros(
                (n_methods, n_activities, n_years),
                dtype=getattr(trails, "value_dtype", float),
            )
            scores = xr.DataArray(
                arr,
                dims=("method", "activity", "year"),
                coords={
                    "method": _method_labels(methods),
                    "activity": np.arange(n_activities, dtype=int),
                    "year": years,
                },
            )
        else:
            arr = sparse.zeros(
                (n_activities, n_years), dtype=getattr(trails, "value_dtype", float)
            )
            scores = xr.DataArray(
                arr,
                dims=("activity", "year"),
                coords={
                    "activity": np.arange(n_activities, dtype=int),
                    "year": years,
                },
            )
        trails.scores = scores
        return scores

    coords = data.coords.astype(np.int64, copy=False)
    values = data.data.astype(getattr(trails, "value_dtype", np.float64), copy=False)
    activities = coords[0]
    flows = coords[1]
    year_indices = coords[2]
    roots = coords[3] if has_root else None

    map_scenario = getattr(trails, "_map_year_to_scenario_year", None)
    if not callable(map_scenario):
        map_scenario = lambda year: int(year)

    output_coords: list[np.ndarray] = []
    output_data: list[np.ndarray] = []
    year_groups: dict[int, list[int]] = {}
    edges_by_scenario_year: dict[int, set[tuple[int, int]]] = {}
    mapping_caches: list[dict[str, Any]] | None = (
        [{} for _method in methods] if reuse_cached_cfs else None
    )

    for year_idx in np.unique(year_indices):
        raw_year = int(years[int(year_idx)])
        scenario_year = int(map_scenario(raw_year))
        year_groups.setdefault(scenario_year, []).append(int(year_idx))

        mask = year_indices == year_idx
        pair_ids = np.unique(
            flows[mask].astype(np.int64, copy=False) * n_activities
            + activities[mask].astype(np.int64, copy=False)
        )
        edges = edges_by_scenario_year.setdefault(scenario_year, set())
        edges.update(
            (int(pair_id // n_activities), int(pair_id % n_activities))
            for pair_id in pair_ids
        )

    progress_total = len(year_groups) * n_methods
    progress_context = (
        tqdm(
            total=progress_total,
            desc="EDGES LCIA",
            unit="method-year",
        )
        if show_progress and progress_total
        else nullcontext()
    )

    with progress_context as progress:
        for scenario_year in sorted(year_groups):
            matrices = _build_edges_characterization_matrices_for_year(
                trails,
                methods,
                year=scenario_year,
                additional_topologies=additional_topologies,
                strategies=strategies,
                biosphere_edges=edges_by_scenario_year[scenario_year],
                mapping_caches=mapping_caches,
                progress=progress,
                suppress_edges_output=True,
            )

            for year_idx in year_groups[scenario_year]:
                mask = year_indices == year_idx
                act_group = activities[mask]
                flow_group = flows[mask]
                value_group = values[mask]
                root_group = roots[mask] if roots is not None else None

                for method_idx, matrix in enumerate(matrices):
                    cf_values = _lookup_sparse_values(matrix, flow_group, act_group)
                    scored = value_group * cf_values.astype(
                        value_group.dtype, copy=False
                    )
                    keep = scored != 0.0
                    if not np.any(keep):
                        continue

                    if has_root and n_methods > 1:
                        output_coords.append(
                            np.vstack(
                                [
                                    np.full(
                                        np.count_nonzero(keep),
                                        method_idx,
                                        dtype=np.int64,
                                    ),
                                    act_group[keep],
                                    np.full(
                                        np.count_nonzero(keep),
                                        year_idx,
                                        dtype=np.int64,
                                    ),
                                    root_group[keep],  # type: ignore[index]
                                ]
                            )
                        )
                    elif has_root:
                        output_coords.append(
                            np.vstack(
                                [
                                    act_group[keep],
                                    np.full(
                                        np.count_nonzero(keep),
                                        year_idx,
                                        dtype=np.int64,
                                    ),
                                    root_group[keep],  # type: ignore[index]
                                ]
                            )
                        )
                    elif n_methods > 1:
                        output_coords.append(
                            np.vstack(
                                [
                                    np.full(
                                        np.count_nonzero(keep),
                                        method_idx,
                                        dtype=np.int64,
                                    ),
                                    act_group[keep],
                                    np.full(
                                        np.count_nonzero(keep),
                                        year_idx,
                                        dtype=np.int64,
                                    ),
                                ]
                            )
                        )
                    else:
                        output_coords.append(
                            np.vstack(
                                [
                                    act_group[keep],
                                    np.full(
                                        np.count_nonzero(keep),
                                        year_idx,
                                        dtype=np.int64,
                                    ),
                                ]
                            )
                        )
                    output_data.append(scored[keep])

            del matrix, matrices

    dtype = getattr(trails, "value_dtype", values.dtype)
    if has_root and n_methods > 1:
        shape = (n_methods, n_activities, n_years, n_roots)
        dims = ("method", "activity", "year", "root activity")
        coords_xr = {
            "method": _method_labels(methods),
            "activity": inv.coords["activity"].values,
            "year": years,
            "root activity": inv.coords["root activity"].values,
        }
    elif has_root:
        shape = (n_activities, n_years, n_roots)
        dims = ("activity", "year", "root activity")
        coords_xr = {
            "activity": inv.coords["activity"].values,
            "year": years,
            "root activity": inv.coords["root activity"].values,
        }
    elif n_methods > 1:
        shape = (n_methods, n_activities, n_years)
        dims = ("method", "activity", "year")
        coords_xr = {
            "method": _method_labels(methods),
            "activity": inv.coords["activity"].values,
            "year": years,
        }
    else:
        shape = (n_activities, n_years)
        dims = ("activity", "year")
        coords_xr = {"activity": inv.coords["activity"].values, "year": years}

    if output_coords:
        arr = sparse.COO(
            np.concatenate(output_coords, axis=1),
            np.concatenate(output_data).astype(dtype, copy=False),
            shape=shape,
        )
    else:
        arr = sparse.zeros(shape, dtype=dtype)

    scores = xr.DataArray(arr, dims=dims, coords=coords_xr)
    trails.scores = scores
    return scores
