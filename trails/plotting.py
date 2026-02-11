from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union
import json
import os

import math
import bisect

from collections import defaultdict
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import sparse
import xarray as xr

from .trails import Trails
from .utils import _format_path_label, _format_path_label_with_years


def _build_activity_label_map(trails: Trails) -> dict[int, str]:
    """Build a mapping of activity indices to display labels.

    For plot_temporal_scores we want labels to show ONLY the reference product
    (used in both legend and hover).
    """
    labels: dict[int, str] = {}
    for scen_label, mapping in trails.activity_indices.items():
        for idx, meta in mapping.items():
            if idx in labels:
                continue

            rp = (meta.get("reference product") or "").strip()
            if rp:
                labels[idx] = rp
            else:
                # Fallback if reference product missing
                name = (meta.get("name") or "").strip()
                labels[idx] = name if name else f"Activity {idx}"

    return labels


def _build_flow_label_map(trails: Trails) -> dict[int, str]:
    labels: dict[int, str] = {}

    # Determine which flow keys are actually used in plotted arrays
    flow_coord = None
    if (
        getattr(trails, "characterized_inventory", None) is not None
        and "flow" in trails.characterized_inventory.coords
    ):
        flow_coord = trails.characterized_inventory.coords["flow"].values
    elif (
        getattr(trails, "inventory", None) is not None
        and "flow" in trails.inventory.coords
    ):
        flow_coord = trails.inventory.coords["flow"].values

    coord_value_set = (
        set(int(v) for v in flow_coord) if flow_coord is not None else None
    )

    for scen_label, mapping in trails.biosphere_indices.items():
        for idx, meta in mapping.items():
            k = int(idx)

            # If we know the plotted flow coordinate values, only keep matching keys
            if coord_value_set is not None and k not in coord_value_set:
                continue

            if k in labels:
                continue

            name = meta.get("name") or f"Flow {k}"
            compartment = meta.get("compartment") or ""
            subcompartment = meta.get("subcompartment") or ""
            unit = meta.get("unit") or ""

            label = name
            if compartment or subcompartment:
                parts = [p for p in (compartment, subcompartment) if p]
                label += " | " + "/".join(parts)
            if unit:
                label += f" ({unit})"
            labels[k] = label

    return labels


def _select_years_from_results(
    results_by_year: dict[int, dict[str, Any]],
    year_range: tuple[int, int] | None,
) -> list[int]:
    """Select years from results with an optional range filter.

    :param results_by_year: Mapping of year to results payload.
    :type results_by_year: dict[int, dict]
    :param year_range: Optional ``(start, end)`` bounds.
    :type year_range: tuple[int, int] | None
    :returns: Sorted list of years to plot.
    :rtype: list[int]
    """
    years_all = sorted(results_by_year.keys())
    if not years_all:
        raise ValueError("results_by_year is empty.")

    if year_range is not None:
        y0, y1 = year_range
        years = [y for y in years_all if y0 <= y <= y1]
    else:
        years = years_all

    if not years:
        raise ValueError("No years available after applying year_range.")

    return years


def _collect_root_scores(
    results_by_year: dict[int, dict[str, Any]],
    years: list[int],
    score_key: str,
) -> list[int]:
    """Collect unique root indices with scores across years.

    :param results_by_year: Mapping of year to results payload.
    :type results_by_year: dict[int, dict]
    :param years: Years to scan.
    :type years: list[int]
    :returns: Sorted list of root indices.
    :rtype: list[int]
    """
    all_roots = sorted(
        {root for year in years for root in results_by_year[year].get(score_key, {})}
    )
    if not all_roots:
        raise ValueError(f"No {score_key} found.")
    return all_roots


def _build_score_matrix(
    results_by_year: dict[int, dict[str, Any]],
    years: list[int],
    all_roots: list[int],
    score_key: str,
) -> np.ndarray:
    """Build a score matrix of shape ``(years, roots)``.

    :param results_by_year: Mapping of year to results payload.
    :type results_by_year: dict[int, dict]
    :param years: Years to include.
    :type years: list[int]
    :param all_roots: Root indices to include.
    :type all_roots: list[int]
    :returns: Score matrix.
    :rtype: numpy.ndarray
    """
    Y = np.zeros((len(years), len(all_roots)), dtype=float)
    for yi, year in enumerate(years):
        spr = results_by_year[year].get(score_key, {})
        for ri, root in enumerate(all_roots):
            Y[yi, ri] = spr.get(root, 0.0)
    return Y


def _add_root_traces(
    fig: go.Figure,
    years: list[int],
    Y: np.ndarray,
    all_roots: list[Any],
    idx_to_label: dict[Any, str],
    method_label: str,
    stacked: bool,
    showlegend_roots: Optional[set[int]] = None,
    showhover_roots: Optional[set[int]] = None,
) -> None:
    """Add root score traces to a Plotly figure.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param years: X-axis years.
    :type years: list[int]
    :param Y: Score matrix with shape ``(years, roots)``.
    :type Y: numpy.ndarray
    :param all_roots: Root indices corresponding to columns in ``Y``.
    :type all_roots: list[int]
    :param idx_to_label: Mapping from activity index to label.
    :type idx_to_label: dict[int, str]
    :param method_label: LCIA method label for hover text.
    :type method_label: str
    :param stacked: Whether to stack traces.
    :type stacked: bool
    :param showlegend_roots: Optional set of root ids to show in the legend.
    :type showlegend_roots: set[int] | None
    :param showhover_roots: Optional set of root ids to show in hover.
    :type showhover_roots: set[int] | None
    """
    alpha = 0.4 if not stacked else 1.0

    def label_for_root(idx: Any) -> str:
        """Resolve a root index to its display label.

        :param idx: Activity index to label.
        :type idx: int
        :returns: Display label for the root.
        :rtype: str
        """
        if idx in idx_to_label:
            return idx_to_label[idx]
        if isinstance(idx, str):
            return idx
        return f"Activity {idx}"

    for ri, root in enumerate(all_roots):
        showlegend = True
        if showlegend_roots is not None:
            showlegend = root in showlegend_roots

        showhover = True
        if showhover_roots is not None:
            showhover = root in showhover_roots

        root_label = label_for_root(root)
        root_label_wrapped = _wrap_hover_label(root_label, max_chars=45)

        fig.add_trace(
            go.Scatter(
                x=years,
                y=Y[:, ri],
                name=root_label,  # keep legend name as before
                meta=root_label_wrapped,  # wrapped version for hover
                showlegend=showlegend,
                mode="lines",
                hoverinfo="skip" if not showhover else None,
                hovertemplate=(
                    (
                        "<b>%{meta}</b><br>"
                        "Year: %{x}<br>"
                        f"{method_label}: %{{y:.6g}}<extra></extra>"
                    )
                    if showhover
                    else None
                ),
                **(
                    {"stackgroup": "one"}
                    if stacked
                    else {
                        "fill": "tozeroy",
                        "line": dict(width=2),
                        "opacity": alpha,
                    }
                ),
            )
        )


def plot_temporal_graph(
    trails: Trails,
    *,
    min_edge_amount: float = 0.0,
    notebook: bool = False,
    filename: str = "trails_graph.html",
    height: str = "800px",
    width: str = "100%",
    physics: bool = False,
    layout_by_year_depth: bool = True,
    year_scale: float = 300.0,
    depth_scale: float = 2000.0,
    max_label_chars: int = 28,
    level0_edge_color: str = "#e45756",
    palette: Optional[list[str]] = None,
    show_year_labels: bool = True,
    year_label_offset: float = 1.2,
    show_band_labels: bool = True,
    band_label_chars: int = 80,
    band_label_offset_px: float | None = None,
    auto_depth_scale: bool = True,
    edge_weight: Literal["amount", "score"] = "amount",
    node_scores: dict[tuple, float] | None = None,
) -> str:
    """Render the routing graph with pyvis.

    :param trails: Trails instance with a populated graph.
    :type trails: Trails
    :param min_edge_amount: Minimum absolute edge value to include (amount or score).
    :type min_edge_amount: float
    :param notebook: Whether to render for Jupyter notebooks.
    :type notebook: bool
    :param filename: Output HTML file name.
    :type filename: str
    :param height: HTML canvas height.
    :type height: str
    :param width: HTML canvas width.
    :type width: str
    :param physics: Enable physics simulation.
    :type physics: bool
    :param layout_by_year_depth: Position nodes by year (x) and depth (y).
    :type layout_by_year_depth: bool
    :param year_scale: Scale factor for year spacing on x-axis.
    :type year_scale: float
    :param depth_scale: Scale factor for depth spacing on y-axis.
    :type depth_scale: float
    :param max_label_chars: Max characters per label line (name/ref).
    :type max_label_chars: int
    :param level0_edge_color: Color for edges from depth-0 nodes.
    :type level0_edge_color: str
    :param palette: Colors for first-level branches (depth=1).
    :type palette: list[str] | None
    :param show_year_labels: Add year labels below the graph.
    :type show_year_labels: bool
    :param year_label_offset: Vertical offset multiplier for year labels.
    :type year_label_offset: float
    :param show_band_labels: Add right-side labels for depth bands.
    :type show_band_labels: bool
    :param band_label_chars: Max characters for band labels.
    :type band_label_chars: int
    :param band_label_offset_px: Pixel offset for band labels (vertical). If None, auto-compute.
    :type band_label_offset_px: float | None
    :param auto_depth_scale: Increase depth spacing based on band count.
    :type auto_depth_scale: bool
    :param edge_weight: Whether to size edges by technosphere amount or score.
    :type edge_weight: Literal["amount", "score"]
    :param node_scores: Optional node score mapping used when edge_weight="score".
    :type node_scores: dict[tuple, float] | None
    :returns: Output filename.
    :rtype: str
    """
    G = getattr(trails, "graph", None)
    if G is None:
        raise RuntimeError(
            "Trails graph is missing; run trails.temporal_routing(...) first."
        )

    try:
        import networkx as nx
    except Exception as exc:  # pragma: no cover - dependency guard
        raise RuntimeError(
            "networkx is required for plot_temporal_graph(). "
            "Install it with `pip install networkx`."
        ) from exc

    try:
        from pyvis.network import Network
        import pyvis
    except Exception as exc:  # pragma: no cover - dependency guard
        raise RuntimeError(
            "pyvis is required for plot_temporal_graph(). "
            "Install it with `pip install pyvis`."
        ) from exc

    def _node_score(node: object) -> float:
        if not node_scores:
            return 0.0
        if node in node_scores:
            return float(node_scores.get(node, 0.0))
        data = G.nodes.get(node, {})
        year = int(data.get("year", -1))
        act = int(data.get("act_idx", -1))
        return float(node_scores.get((year, act), 0.0))

    edge_value_map: dict[tuple, float] = {}
    if edge_weight == "score":
        if node_scores is None:
            raise ValueError(
                "edge_weight='score' requires node_scores (use score_temporal_graph_nodes)."
            )
        incoming_abs: dict[object, float] = defaultdict(float)
        for u, v, d in G.edges(data=True):
            amt = float(d.get("amount", 0.0))
            incoming_abs[v] += abs(amt)
        for u, v, d in G.edges(data=True):
            amt = float(d.get("amount", 0.0))
            denom = float(incoming_abs.get(v, 0.0))
            child_score = _node_score(v)
            if denom > 0.0:
                edge_value_map[(u, v)] = child_score * (abs(amt) / denom)
            else:
                edge_value_map[(u, v)] = 0.0

    H = G.copy()

    # Filter edges by amount
    if min_edge_amount > 0.0:
        if edge_weight == "score":
            edges = [
                (u, v)
                for u, v in H.edges()
                if abs(float(edge_value_map.get((u, v), 0.0))) >= float(min_edge_amount)
            ]
        else:
            edges = [
                (u, v)
                for u, v, d in H.edges(data=True)
                if abs(float(d.get("amount", 0.0))) >= float(min_edge_amount)
            ]
        H = H.edge_subgraph(edges).copy()

    # Relabel tuple node ids to strings for pyvis compatibility
    def _truncate(text: str, limit: int) -> str:
        if len(text) <= limit:
            return text
        return text[:limit]

    def _label_node(node: tuple, data: dict) -> str:
        year = data.get("year", "")
        depth = data.get("depth", "")
        name = data.get("name", "") or ""
        ref = data.get("reference_product", "") or ""
        loc = data.get("location", "") or ""
        act = data.get("act_idx", "")
        name = _truncate(str(name), max_label_chars)
        ref = _truncate(str(ref), max_label_chars)
        label = f"{name}"
        if ref:
            label = f"{label}\n{ref}"
        meta_parts = [str(year), f"d{depth}", str(act), loc]
        return f"{label}\n" + " | ".join([p for p in meta_parts if p])

    mapping = {}
    for n, d in H.nodes(data=True):
        mapping[n] = _label_node(n, d)

    net = Network(height=height, width=width, directed=True, notebook=notebook)
    try:
        net.set_template_dir(os.path.join(os.path.dirname(pyvis.__file__), "templates"))
    except Exception:
        pass

    if palette is None:
        palette = [
            "#4c78a8",
            "#f58518",
            "#54a24b",
            "#e45756",
            "#72b7b2",
            "#b279a2",
            "#ff9da6",
            "#9d755d",
            "#bab0ac",
        ]

    # Assign colors per depth-0 edge, then propagate to descendants of that child.
    depth0_nodes = [n for n, d in H.nodes(data=True) if int(d.get("depth", 0)) == 0]
    edge_colors: dict[tuple, str] = {}
    node_branch_color: dict = {}
    color_idx = 0

    def _branch_key(node: object) -> tuple[str, str, str]:
        data = H.nodes[node]
        return (
            str(data.get("name", "")),
            str(data.get("reference_product", "")),
            str(data.get("location", "")),
        )

    key_colors: dict[tuple[str, str], str] = {}
    for n in depth0_nodes:
        for _, child in H.out_edges(n):
            key = _branch_key(child)
            if key in key_colors:
                color = key_colors[key]
            else:
                color = palette[color_idx % len(palette)]
                color_idx += 1
                key_colors[key] = color
            edge_colors[(n, child)] = color
            node_branch_color[child] = color
            queue = [child]
            while queue:
                cur = queue.pop(0)
                for _, nxt in H.out_edges(cur):
                    if nxt not in node_branch_color:
                        node_branch_color[nxt] = color
                        queue.append(nxt)

    def _edge_color(u: object, v: object) -> str:
        src_depth = int(H.nodes[u].get("depth", 0))
        if src_depth == 0:
            return edge_colors.get((u, v), level0_edge_color)
        return node_branch_color.get(u, palette[0])

    def _node_color(n: object) -> Optional[str]:
        depth = int(H.nodes[n].get("depth", 0))
        if depth == 0:
            return None
        return node_branch_color.get(n)

    initial_depth = 1
    full_nodes: list[dict[str, object]] = []
    full_edges: list[dict[str, object]] = []

    if layout_by_year_depth:
        net.toggle_physics(False)
        years = [float(d.get("year", 0.0)) for _, d in H.nodes(data=True)]
        depths = [float(d.get("depth", 0.0)) for _, d in H.nodes(data=True)]
        min_year = min(years) if years else 0.0
        max_year = max(years) if years else 0.0
        min_depth = min(depths) if depths else 0.0
        max_depth_val = max(depths) if depths else 0.0
        year_mid = (min_year + max_year) / 2.0
        depth_mid = 0.0
        depth_buckets: dict[int, list[tuple[object, dict]]] = {}
        for n, d in H.nodes(data=True):
            depth_buckets.setdefault(int(d.get("depth", 0)), []).append((n, d))
        depth_offsets: dict[object, float] = {}
        band_labels: list[dict[str, object]] = []
        max_bands_per_depth = 1
        for depth, items in depth_buckets.items():
            band_map: dict[tuple[str, str, str], list[tuple[object, dict]]] = {}
            for n, d in items:
                key = (
                    str(d.get("name", "")),
                    str(d.get("reference_product", "")),
                    str(d.get("location", "")),
                )
                band_map.setdefault(key, []).append((n, d))
            bands = sorted(band_map.items(), key=lambda it: it[0])
            if not bands:
                continue
            max_bands_per_depth = max(max_bands_per_depth, len(bands))
            band_step = 32.0
            band_start = -band_step * (len(bands) - 1) / 2.0
            for band_idx, (_key, band_items) in enumerate(bands):
                band_offset = band_start + band_idx * band_step
                max_offset = float(depth_scale) * 0.35
                if band_offset > max_offset:
                    band_offset = max_offset
                if band_offset < -max_offset:
                    band_offset = -max_offset
                name = _truncate(_key[0], int(max_label_chars))
                ref = _truncate(_key[1], int(max_label_chars)) if _key[1] else ""
                loc = _key[2]
                label = name
                if ref:
                    label = f"{label} | {ref}"
                if loc:
                    label = f"{label} | {loc}"
                rep_node = band_items[0][0]
                band_labels.append(
                    {
                        "node": mapping[rep_node],
                        "label": label,
                    }
                )
                for n, _ in band_items:
                    depth_offsets[n] = band_offset
        if auto_depth_scale:
            # Ensure depth spacing dominates band offsets.
            max_band_offset = (24.0 * (max_bands_per_depth - 1)) / 2.0
            depth_scale = max(float(depth_scale), 2.0 * max_band_offset + 200.0)
        for n, d in H.nodes(data=True):
            label = mapping[n]
            year = float(d.get("year", 0.0))
            depth = float(d.get("depth", 0.0))
            band_key = (
                str(d.get("name", "")),
                str(d.get("reference_product", "")),
                str(d.get("location", "")),
            )
            x = (year - year_mid) * float(year_scale)
            y = -(depth - depth_mid) * float(depth_scale) + float(
                depth_offsets.get(n, 0.0)
            )
            color = _node_color(n)
            hidden = int(depth) > int(initial_depth)
            if color:
                node_payload = {
                    "id": label,
                    "label": "",
                    "x": x,
                    "y": y,
                    "physics": False,
                    "color": color,
                    "depth": int(depth),
                    "year": year,
                    "band_key": f"{band_key[0]}|{band_key[1]}|{band_key[2]}",
                    "band_offset": float(depth_offsets.get(n, 0.0)),
                    "shape": "dot",
                    "size": 12,
                }
            else:
                node_payload = {
                    "id": label,
                    "label": "",
                    "x": x,
                    "y": y,
                    "physics": False,
                    "depth": int(depth),
                    "year": year,
                    "band_key": f"{band_key[0]}|{band_key[1]}|{band_key[2]}",
                    "band_offset": float(depth_offsets.get(n, 0.0)),
                    "shape": "dot",
                    "size": 12,
                }
            full_nodes.append(node_payload)
            if not hidden:
                n_id = node_payload["id"]
                node_payload_copy = dict(node_payload)
                node_payload_copy.pop("id", None)
                net.add_node(n_id, **node_payload_copy)

        def _edge_title(
            src_label: str,
            dst_label: str,
            amount: float,
            score: float | None = None,
        ) -> str:
            src = src_label.split("\n", 1)[0]
            dst = dst_label.split("\n", 1)[0]
            if score is None:
                return f"{src} → {dst} | amount={amount:.3g}"
            return f"{src} → {dst} | amount={amount:.3g} | score={score:.3g}"

        if show_year_labels:
            unique_years = sorted({int(y) for y in years})
            bottom_y = min(
                (-(float(d.get("depth", 0.0)) - depth_mid) * float(depth_scale))
                + float(depth_offsets.get(n, 0.0))
                for n, d in H.nodes(data=True)
            )
            # Inject year labels as HTML overlay to avoid vis.js rendering quirks.
            x_positions = {
                year: (float(year) - year_mid) * float(year_scale)
                for year in unique_years
            }
        for u, v, d in H.edges(data=True):
            src = mapping[u]
            dst = mapping[v]
            src_depth = int(H.nodes[u].get("depth", 0))
            dst_depth = int(H.nodes[v].get("depth", 0))
            hidden_edge = (src_depth > initial_depth) or (dst_depth > initial_depth)
            edge_amount = float(d.get("amount", 0.0))
            edge_value = (
                float(edge_value_map.get((u, v), 0.0))
                if edge_weight == "score"
                else edge_amount
            )
            edge_payload = {
                "from": src,
                "to": dst,
                "value": float(edge_value),
                "color": _edge_color(u, v),
                "title": _edge_title(
                    src,
                    dst,
                    float(edge_amount),
                    float(edge_value) if edge_weight == "score" else None,
                ),
                "depth_from": src_depth,
                "depth_to": dst_depth,
            }
            full_edges.append(edge_payload)
            if not hidden_edge:
                src = edge_payload["from"]
                dst = edge_payload["to"]
                edge_payload_copy = dict(edge_payload)
                edge_payload_copy.pop("from", None)
                edge_payload_copy.pop("to", None)
                net.add_edge(src, dst, **edge_payload_copy)
    else:
        color_by_label = {mapping[n]: _node_color(n) for n in H.nodes()}
        H = nx.relabel_nodes(H, mapping, copy=True)
        for n, d in H.nodes(data=True):
            color = color_by_label.get(n)
            depth = int(d.get("depth", 0))
            year = float(d.get("year", 0.0))
            band_key = (
                str(d.get("name", "")),
                str(d.get("reference_product", "")),
                str(d.get("location", "")),
            )
            hidden = int(depth) > int(initial_depth)
            if color:
                node_payload = {
                    "id": n,
                    "label": "",
                    "color": color,
                    "depth": depth,
                    "year": year,
                    "band_key": f"{band_key[0]}|{band_key[1]}|{band_key[2]}",
                    "shape": "dot",
                    "size": 12,
                }
            else:
                node_payload = {
                    "id": n,
                    "label": "",
                    "depth": depth,
                    "year": year,
                    "band_key": f"{band_key[0]}|{band_key[1]}|{band_key[2]}",
                    "shape": "dot",
                    "size": 12,
                }
            full_nodes.append(node_payload)
            if not hidden:
                n_id = node_payload["id"]
                node_payload_copy = dict(node_payload)
                node_payload_copy.pop("id", None)
                net.add_node(n_id, **node_payload_copy)

        def _edge_title(
            src_label: str,
            dst_label: str,
            amount: float,
            score: float | None = None,
        ) -> str:
            src = src_label.split("\n", 1)[0]
            dst = dst_label.split("\n", 1)[0]
            if score is None:
                return f"{src} → {dst} | amount={amount:.3g}"
            return f"{src} → {dst} | amount={amount:.3g} | score={score:.3g}"

        for u, v, d in H.edges(data=True):
            src_depth = int(H.nodes[u].get("depth", 0))
            dst_depth = int(H.nodes[v].get("depth", 0))
            hidden_edge = (src_depth > initial_depth) or (dst_depth > initial_depth)
            edge_amount = float(d.get("amount", 0.0))
            edge_value = (
                float(edge_value_map.get((u, v), 0.0))
                if edge_weight == "score"
                else edge_amount
            )
            edge_payload = {
                "from": u,
                "to": v,
                "value": float(edge_value),
                "color": _edge_color(u, v),
                "title": _edge_title(
                    u,
                    v,
                    float(edge_amount),
                    float(edge_value) if edge_weight == "score" else None,
                ),
                "depth_from": src_depth,
                "depth_to": dst_depth,
            }
            full_edges.append(edge_payload)
            if not hidden_edge:
                src = edge_payload["from"]
                dst = edge_payload["to"]
                edge_payload_copy = dict(edge_payload)
                edge_payload_copy.pop("from", None)
                edge_payload_copy.pop("to", None)
                net.add_edge(src, dst, **edge_payload_copy)
        net.toggle_physics(physics)
    if notebook:
        net.show(filename)
    else:
        net.write_html(filename, open_browser=False)

    if show_year_labels and layout_by_year_depth:
        try:
            with open(filename, "r", encoding="utf-8") as f:
                html = f.read()
            overlay_items = [
                {"year": int(year), "x": float(x)} for year, x in x_positions.items()
            ]
            bottom_y = min(
                (-(float(d.get("depth", 0.0)) - depth_mid) * float(depth_scale))
                + float(depth_offsets.get(n, 0.0))
                for n, d in H.nodes(data=True)
            )
            bottom_y = float(bottom_y) - (float(depth_scale) * float(year_label_offset))
            depth_values = sorted(
                {int(d.get("depth", 0)) for _, d in H.nodes(data=True)}
            )
            band_items = band_labels if show_band_labels else []
            max_depth_in_graph = max(depths) if depths else 0.0
            initial_depth = 1
            overlay_script = f"""
<style>
#mynetwork {{
  position: relative;
}}
#trail-year-overlay {{
  position: absolute;
  left: 0;
  right: 0;
  top: 0;
  bottom: 0;
  pointer-events: none;
  font-size: 12px;
  color: #444;
  z-index: 999;
}}
#trail-year-overlay .trail-year,
#trail-year-overlay .trail-year-label {{
  position: absolute;
  transform: translateX(-50%);
  white-space: nowrap;
}}
#trail-year-overlay .trail-year {{
  width: 1px;
  height: 8px;
  background: #888;
}}
#trail-year-overlay .trail-year-label {{
  top: 10px;
}}
#trail-year-overlay .trail-depth-line {{
  position: absolute;
  left: 0;
  right: 0;
  height: 2px;
  background: rgba(0,0,0,0.2);
}}
#trail-year-overlay .trail-band-label {{
  position: absolute;
  right: 6px;
  transform: translateY(-50%);
  white-space: nowrap;
  text-align: right;
  max-width: none;
}}
</style>
<script>
  window.addEventListener('load', function() {{
    var container = document.getElementById('mynetwork');
    if (!container) return;
    var overlay = document.createElement('div');
    overlay.id = 'trail-year-overlay';
    var items = {overlay_items};
    var depths = {depth_values};
    var depthScale = {float(depth_scale)};
    var bands = {band_items};
    var bandOffsetPx = {float(band_label_offset_px) if band_label_offset_px is not None else 0.0};
    function renderAxis() {{
      overlay.innerHTML = '';
      if (typeof network === 'undefined' || !network) return;
      var view = network.getViewPosition();
      var h = overlay.clientHeight || 0;
      var y_px = h - 18;
      var w = overlay.clientWidth || 0;
      var depthScale = window.trailDepthScale || depthScale;
      var depthMid = window.trailDepthMid || {depth_mid};
      items.forEach(function(item) {{
        var dom = network.canvasToDOM({{x: item.x, y: view.y}});
        var tick = document.createElement('span');
        tick.className = 'trail-year';
        tick.style.left = dom.x.toFixed(1) + 'px';
        tick.style.top = y_px.toFixed(1) + 'px';
        overlay.appendChild(tick);
        var label = document.createElement('span');
        label.className = 'trail-year-label';
        label.style.left = dom.x.toFixed(1) + 'px';
        label.style.top = (y_px + 10).toFixed(1) + 'px';
        label.textContent = String(item.year);
        overlay.appendChild(label);
      }});
      depths.forEach(function(d) {{
        var domY = network.canvasToDOM({{x: view.x, y: -((d + 0.5) - depthMid) * depthScale}});
        var line = document.createElement('div');
        line.className = 'trail-depth-line';
        line.style.top = domY.y.toFixed(1) + 'px';
        overlay.appendChild(line);
        var tick = document.createElement('span');
        tick.className = 'trail-year';
        tick.style.left = '14px';
        tick.style.top = domY.y.toFixed(1) + 'px';
        overlay.appendChild(tick);
        var label = document.createElement('span');
        label.className = 'trail-year-label';
        label.style.left = '28px';
        label.style.top = (domY.y - 6).toFixed(1) + 'px';
        label.textContent = 'd' + String(d);
        overlay.appendChild(label);
      }});
      bands.forEach(function(b) {{
        if (!b.node) return;
        var box = network.getBoundingBox(b.node);
        if (!box) return;
        var centerY = (box.top + box.bottom) / 2.0;
        var domY = network.canvasToDOM({{x: 0, y: centerY}});
        var label = document.createElement('span');
        label.className = 'trail-band-label';
        label.style.top = (domY.y + bandOffsetPx).toFixed(1) + 'px';
        label.textContent = b.label;
        overlay.appendChild(label);
      }});
    }}
    window.renderAxis = renderAxis;
    container.appendChild(overlay);
    var tries = 0;
    function bindAxis() {{
      if (typeof network !== 'undefined' && network) {{
        renderAxis();
        network.on('zoom', renderAxis);
        network.on('dragEnd', renderAxis);
        network.on('stabilized', renderAxis);
      }} else if (tries < 10) {{
        tries += 1;
        setTimeout(bindAxis, 200);
      }}
    }}
    bindAxis();
  }});
</script>
"""
            html = html.replace("</body>", overlay_script + "\n</body>")
            with open(filename, "w", encoding="utf-8") as f:
                f.write(html)
        except Exception:
            pass

    try:
        with open(filename, "r", encoding="utf-8") as f:
            html = f.read()
        # Build branch selector data for depth-0 -> depth-1 edges
        branch_items: list[dict[str, object]] = []
        if H.number_of_nodes() > 0:
            grouped: dict[tuple[str, str], set[str]] = {}
            for u, v in H.edges():
                if (
                    int(H.nodes[u].get("depth", 0)) == 0
                    and int(H.nodes[v].get("depth", 0)) == 1
                ):
                    key = _branch_key(v)
                    label = key[0]
                    if key[1]:
                        label = f"{label} | {key[1]}"
                    if key[2]:
                        label = f"{label} | {key[2]}"
                    nodes_set = grouped.setdefault(key, set())
                    nodes_set.add(mapping[v])
                    nodes_set.add(mapping[u])
                    try:
                        for desc in nx.descendants(H, v):
                            nodes_set.add(mapping[desc])
                    except Exception:
                        pass
            for key, nodes_set in grouped.items():
                label = key[0]
                if key[1]:
                    label = f"{label} | {key[1]}"
                if key[2]:
                    label = f"{label} | {key[2]}"
                branch_items.append(
                    {
                        "label": label,
                        "nodes": sorted(nodes_set),
                    }
                )

        full_nodes_json = json.dumps(full_nodes)
        full_edges_json = json.dumps(full_edges)
        base_depth_scale = float(depth_scale)
        auto_depth_flag = bool(auto_depth_scale)
        year_scale_val = float(year_scale)
        selector_script = f"""
<style>
#trail-depth-selector {{
  position: absolute;
  top: 8px;
  left: 8px;
  z-index: 1000;
  background: rgba(255,255,255,0.9);
  border: 1px solid #ccc;
  padding: 6px 8px;
  font-size: 12px;
  border-radius: 4px;
}}
#trail-branch-selector {{
  position: absolute;
  top: 8px;
  right: 8px;
  z-index: 1000;
  background: rgba(255,255,255,0.9);
  border: 1px solid #ccc;
  padding: 6px 8px;
  font-size: 12px;
  border-radius: 4px;
}}
</style>
<script>
  window.addEventListener('load', function() {{
    var container = document.getElementById('mynetwork');
    if (!container) return;
    var fullNodes = {full_nodes_json};
    var fullEdges = {full_edges_json};
    var baseDepthScale = {base_depth_scale};
    var autoDepthScale = {str(auto_depth_flag).lower()};
    var yearScale = {year_scale_val};
    var sel = document.createElement('select');
    sel.id = 'trail-depth-selector';
    var maxDepth = 0;
    fullNodes.forEach(function(n) {{
      if (n.depth !== undefined && n.depth > maxDepth) maxDepth = n.depth;
    }});
    for (var d = 0; d <= maxDepth; d++) {{
      var opt = document.createElement('option');
      opt.value = d;
      opt.textContent = 'max depth: ' + d;
      sel.appendChild(opt);
    }}
    sel.value = '{initial_depth}';
    var branchItems = {branch_items};
    var currentBranch = 'all';
    function applyFilters(d, branchKey) {{
      if (typeof nodes === 'undefined' || typeof edges === 'undefined') return;
      var visible = {{}};
      var branchSet = null;
      if (branchKey && branchKey !== 'all') {{
        branchItems.forEach(function(b) {{
          if (b.label === branchKey) {{
            branchSet = {{}};
            b.nodes.forEach(function(nid) {{ branchSet[nid] = true; }});
          }}
        }});
      }}
      var nodesFiltered = [];
      fullNodes.forEach(function(n) {{
        var show = (n.depth || 0) <= d;
        if (branchSet) {{
          show = show && !!branchSet[n.id];
        }}
        if (show) {{
          visible[n.id] = true;
          nodesFiltered.push(n);
        }}
      }});
      // recompute depth spacing based on visible bands
      var depthBandCounts = {{}};
      nodesFiltered.forEach(function(n) {{
        var k = String(n.depth) + '|' + String(n.band_key || '');
        depthBandCounts[k] = true;
      }});
      var bandsPerDepth = {{}};
      Object.keys(depthBandCounts).forEach(function(k) {{
        var depth = k.split('|', 1)[0];
        bandsPerDepth[depth] = (bandsPerDepth[depth] || 0) + 1;
      }});
      var maxBands = 1;
      Object.keys(bandsPerDepth).forEach(function(k) {{
        if (bandsPerDepth[k] > maxBands) maxBands = bandsPerDepth[k];
      }});
      var depthScale = baseDepthScale;
      if (autoDepthScale) {{
        var maxBandOffset = (24 * (maxBands - 1)) / 2.0;
        var computed = Math.max(200, 2.0 * maxBandOffset + 200);
        depthScale = (computed > baseDepthScale) ? computed : Math.min(baseDepthScale, computed);
      }}
      // recompute band offsets per depth
      var bandOffsets = {{}};
      var perDepth = {{}};
      nodesFiltered.forEach(function(n) {{
        var depth = n.depth || 0;
        var key = String(n.band_key || '');
        if (!perDepth[depth]) perDepth[depth] = {{}};
        perDepth[depth][key] = true;
      }});
      Object.keys(perDepth).forEach(function(depthStr) {{
        var keys = Object.keys(perDepth[depthStr]).sort();
        var count = keys.length || 1;
        var bandStep = 32;
        var bandStart = -bandStep * (count - 1) / 2.0;
        keys.forEach(function(k, i) {{
          var off = bandStart + i * bandStep;
          var maxOff = depthScale * 0.35;
          if (off > maxOff) off = maxOff;
          if (off < -maxOff) off = -maxOff;
          bandOffsets[depthStr + '|' + k] = off;
        }});
      }});
      // recompute centers
      var minYear = null, maxYear = null;
      nodesFiltered.forEach(function(n) {{
        var y = n.year || 0;
        if (minYear === null || y < minYear) minYear = y;
        if (maxYear === null || y > maxYear) maxYear = y;
      }});
      if (minYear === null) {{ minYear = 0; maxYear = 0; }}
      var yearMid = (minYear + maxYear) / 2.0;
      var depthMid = 0.0;
      window.trailDepthScale = depthScale;
      window.trailDepthMid = depthMid;
      // assign positions
      nodesFiltered.forEach(function(n) {{
        var bandKey = String(n.depth) + '|' + String(n.band_key || '');
        var bandOffset = bandOffsets[bandKey] || 0.0;
        n.x = (n.year - yearMid) * yearScale;
        n.y = -((n.depth || 0) - depthMid) * depthScale + bandOffset;
      }});
      var edgesFiltered = [];
      fullEdges.forEach(function(e) {{
        var show = visible[e.from] && visible[e.to];
        if (show) edgesFiltered.push(e);
      }});
      nodes.clear();
      edges.clear();
      nodes.add(nodesFiltered);
      edges.add(edgesFiltered);
      if (typeof network !== 'undefined' && network) {{
        try {{ network.fit({{animation: false}}); }} catch (e) {{}}
      }}
      if (window.renderAxis) window.renderAxis();
    }}
    sel.addEventListener('change', function() {{
      applyFilters(parseInt(sel.value), currentBranch);
    }});
    container.appendChild(sel);

    var branchSel = document.createElement('select');
    branchSel.id = 'trail-branch-selector';
    var optAll = document.createElement('option');
    optAll.value = 'all';
    optAll.textContent = 'branch: all';
    branchSel.appendChild(optAll);
    branchItems.forEach(function(b) {{
      var opt = document.createElement('option');
      opt.value = b.label;
      opt.textContent = 'branch: ' + b.label;
      branchSel.appendChild(opt);
    }});
    branchSel.addEventListener('change', function() {{
      currentBranch = branchSel.value;
      applyFilters(parseInt(sel.value), currentBranch);
    }});
    container.appendChild(branchSel);
    applyFilters(parseInt(sel.value), currentBranch);
  }});
</script>
"""
        html = html.replace("</body>", selector_script + "\n</body>")
        with open(filename, "w", encoding="utf-8") as f:
            f.write(html)
    except Exception:
        pass
    return filename


def _add_cumulative_trace(
    fig: go.Figure,
    years: list[int],
    total_raw: list[float],
    cumulative_axis_label: str,
    yaxis_type: str,
    log_eps: float,
) -> np.ndarray:
    """Add a cumulative total trace to a Plotly figure.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param years: X-axis years.
    :type years: list[int]
    :param total_raw: Raw total scores per year.
    :type total_raw: list[float]
    :param cumulative_axis_label: Axis label for cumulative scores.
    :type cumulative_axis_label: str
    :param yaxis_type: Axis scale type (e.g. ``"linear"`` or ``"log"``).
    :type yaxis_type: str
    :param log_eps: Small epsilon for log scaling.
    :type log_eps: float
    :returns: Cumulative values.
    :rtype: numpy.ndarray
    """
    cum_vals = np.cumsum(total_raw)
    if yaxis_type == "log":
        cum_vals = np.where(cum_vals > 0, cum_vals, log_eps)

    fig.add_trace(
        go.Scatter(
            x=years,
            y=cum_vals,
            name="Cumulative total",
            showlegend=True,
            mode="lines",
            line=dict(width=2, color="black"),
            yaxis="y2",
            hovertemplate=(
                "<b>Cumulative total</b><br>"
                "Year: %{x}<br>"
                f"{cumulative_axis_label}: %{{y:.6g}}<extra></extra>"
            ),
        )
    )

    return cum_vals


def _add_static_score_trace(
    fig: go.Figure,
    years: list[int],
    static_score: float,
    static_score_label: str,
    static_score_dash: str,
    static_score_color: str,
    method_label: str,
) -> None:
    """Add a horizontal static score trace to a Plotly figure.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param years: X-axis years.
    :type years: list[int]
    :param static_score: Static score value.
    :type static_score: float
    :param static_score_label: Label for the static score trace.
    :type static_score_label: str
    :param static_score_dash: Line dash style.
    :type static_score_dash: str
    :param static_score_color: Line color.
    :type static_score_color: str
    :param method_label: LCIA method label for hover text.
    :type method_label: str
    """
    fig.add_trace(
        go.Scatter(
            x=[years[0], years[-1]],
            y=[static_score, static_score],
            mode="lines",
            name=static_score_label,
            yaxis="y2",
            line=dict(
                dash=static_score_dash,
                color=static_score_color,
                width=2,
            ),
            hovertemplate=(
                f"<b>{static_score_label}</b><br>"
                f"{method_label}: {static_score:.6g}<extra></extra>"
            ),
        )
    )


def _compute_layout_dimensions(
    width: int | None,
    legend_entrywidth: int,
    legend_row_height: int,
    n_items: int,
) -> tuple[int, int]:
    """Compute layout sizing values for the plot.

    :param width: Figure width override.
    :type width: int | None
    :param legend_entrywidth: Legend entry width.
    :type legend_entrywidth: int
    :param legend_row_height: Legend row height.
    :type legend_row_height: int
    :param n_items: Number of legend items.
    :type n_items: int
    :returns: Tuple of entry width and top margin.
    :rtype: tuple[int, int]
    """
    fig_w = int(width) if (width is not None) else 800
    entry_w = max(80, int(legend_entrywidth))
    n_cols = max(1, fig_w // entry_w)
    n_rows = max(1, int(math.ceil(n_items / n_cols)))
    top_margin = 55 + n_rows * int(legend_row_height) + 10
    return entry_w, top_margin


def _apply_base_layout(
    fig: go.Figure,
    width: int,
    height: int,
    title: str,
    legend_y: float,
    entry_w: int,
    top_margin: int,
    method_label: str,
    yaxis_type: str,
    show_cumulative_axis: bool,
    static_score: float | None,
    cumulative_axis_label: str,
) -> None:
    """Apply base layout settings to a Plotly figure.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param width: Figure width.
    :type width: int
    :param height: Figure height.
    :type height: int
    :param title: Figure title.
    :type title: str
    :param legend_y: Legend y position.
    :type legend_y: float
    :param entry_w: Legend entry width.
    :type entry_w: int
    :param top_margin: Top margin height.
    :type top_margin: int
    :param method_label: LCIA method label.
    :type method_label: str
    :param yaxis_type: Axis scale type.
    :type yaxis_type: str
    :param show_cumulative_axis: Whether to show cumulative axis.
    :type show_cumulative_axis: bool
    :param static_score: Optional static score.
    :type static_score: float | None
    :param cumulative_axis_label: Label for cumulative axis.
    :type cumulative_axis_label: str
    """
    fig.update_layout(
        width=width,
        height=height,
        template="plotly_white",
        hovermode="x unified",
        hoverlabel=dict(
            align="left",
            font=dict(size=11),
        ),
        title=dict(
            text=title,
            x=0.5,
            xanchor="center",
        ),
        legend=dict(
            orientation="h",
            x=0,
            xanchor="left",
            y=float(legend_y),
            yanchor="bottom",
            font=dict(size=10),
            entrywidth=entry_w,
            entrywidthmode="pixels",
        ),
        margin=dict(
            l=60,
            r=60 if show_cumulative_axis or static_score is not None else 20,
            t=top_margin,
            b=40,
        ),
        yaxis=dict(
            title=method_label,
            type=yaxis_type,
            showgrid=True,
        ),
    )

    use_y2 = bool(show_cumulative_axis) or (static_score is not None)
    if use_y2:
        fig.update_layout(
            yaxis2=dict(
                title=cumulative_axis_label,
                overlaying="y",
                side="right",
                showgrid=False,
                tickformat=".2e",
                exponentformat="e",
                showexponent="all",
            )
        )


def _apply_linear_yaxis_alignment(
    fig: go.Figure,
    Y: np.ndarray,
    cum_vals: np.ndarray | None,
    static_score: float | None,
    y_min: float | None,
    y_max: float | None,
    y2_max: float | None,
    y2_headroom: float,
    stacked: bool,
) -> None:
    """Align linear y-axes for primary and cumulative traces.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param Y: Score matrix for primary axis.
    :type Y: numpy.ndarray
    :param cum_vals: Cumulative values for secondary axis.
    :type cum_vals: numpy.ndarray | None
    :param static_score: Optional static score value.
    :type static_score: float | None
    :param y_min: Optional min for primary y-axis.
    :type y_min: float | None
    :param y_max: Optional max for primary y-axis.
    :type y_max: float | None
    :param y2_max: Optional max for secondary y-axis.
    :type y2_max: float | None
    :param y2_headroom: Headroom multiplier for secondary axis.
    :type y2_headroom: float
    :param stacked: Whether traces are stacked (affects y-range).
    :type stacked: bool
    """
    if stacked:
        pos_sum = np.sum(np.where(Y > 0, Y, 0.0), axis=1)
        neg_sum = np.sum(np.where(Y < 0, Y, 0.0), axis=1)
        y1_min = float(np.nanmin(neg_sum))
        y1_max_data = float(np.nanmax(pos_sum))
    else:
        y1_min = float(np.nanmin(Y))
        y1_max_data = float(np.nanmax(Y))
    y1_min = min(y1_min, 0.0)
    y1_max_data = max(y1_max_data, 0.0)
    if y1_max_data == y1_min:
        y1_max_data = y1_min + 1.0

    if y_min is not None:
        y1_min = float(y_min)

    y1_max = float(y_max) if (y_max is not None) else y1_max_data
    y1_max = max(y1_max, 0.0)
    if y1_max == y1_min:
        y1_max = y1_min + 1.0

    y2_series = []
    if cum_vals is not None and len(cum_vals) > 0:
        y2_series.append(np.asarray(cum_vals, dtype=float))
    if static_score is not None:
        y2_series.append(np.asarray([float(static_score)], dtype=float))

    if not y2_series:
        fig.update_layout(yaxis=dict(range=[y1_min, y1_max]))
        return

    y2_all = np.concatenate(y2_series)
    y2_min_data = float(np.nanmin(y2_all))
    y2_max_data = float(np.nanmax(y2_all))
    y2_min_data = min(y2_min_data, 0.0)
    y2_max_data = max(y2_max_data, 0.0)
    if y2_max_data == y2_min_data:
        y2_max_data = y2_min_data + 1.0

    y2_max_eff = float(y2_max) if (y2_max is not None) else y2_max_data
    y2_max_eff = max(y2_max_eff, 0.0)
    if y2_max_eff == y2_min_data:
        y2_max_eff = y2_min_data + 1.0

    p = (0.0 - y1_min) / (y1_max - y1_min)

    if y2_max is None:
        y2_max_eff = y2_max_eff * (1.0 + max(0.0, float(y2_headroom)))

    if p <= 0.0:
        y2_range = [0.0, y2_max_eff]
    elif p >= 1.0:
        y2_range = [y2_min_data, 0.0]
    else:
        y2_min_eff = -(p / (1.0 - p)) * y2_max_eff
        y2_min_eff = min(y2_min_eff, y2_min_data)
        y2_range = [y2_min_eff, y2_max_eff]

    fig.update_layout(
        yaxis=dict(range=[y1_min, y1_max]),
        yaxis2=dict(range=y2_range),
    )


def _apply_xaxis_settings(
    fig: go.Figure,
    year_tick: int | None,
    year_range: tuple[int, int] | None,
    years: list[int],
    show_year_grid: bool,
) -> None:
    """Apply x-axis settings for year-based plots.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param year_tick: Tick step for the year axis.
    :type year_tick: int | None
    :param year_range: Optional year range bounds.
    :type year_range: tuple[int, int] | None
    :param years: List of years in the plot.
    :type years: list[int]
    :param show_year_grid: Whether to show grid lines on the year axis.
    :type show_year_grid: bool
    """
    fig.update_xaxes(
        dtick=year_tick,
        tickmode="linear",
        showgrid=show_year_grid,
        tick0=(year_range[0] if year_range else years[0]),
        range=list(year_range) if year_range else None,
    )


def _add_reference_year_line(fig: go.Figure, reference_year: int | None) -> None:
    """Add a vertical reference year line to a Plotly figure.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param reference_year: Year to draw the reference line at.
    :type reference_year: int
    """
    if reference_year is not None:
        fig.add_vline(
            x=reference_year,
            line_width=2,
            line_dash="dash",
            annotation_text="Reference year",
            annotation_position="top",
        )


def to_impact_year_results(
    results: Dict[int, Dict[str, Any]] | Dict[str, Any],
) -> Dict[int, Dict[str, Any]]:
    """
    Normalize different result structures into the "impact-year" format expected by plot_temporal_scores:

        impact_year -> {"scores": float, "scores_per_root": {root: float}}

    Supported inputs:
    1) New explicit structure:
         {"results_by_solve_year": {...}, "results_by_impact_year": {...}}
       -> returns results_by_impact_year

    2) Already impact-year structure:
         {2024: {"scores": ..., "scores_per_root": {...}}, ...}
       -> returned as-is

    3) Legacy per-solve-year duplicated timeline structure:
         {solve_year: {"temporal_scores_by_year": {...}, "temporal_scores_per_root_by_year": {...}, ...}, ...}
       -> aggregated to impact-year (sums across solve years if multiple are present)

    :param results: Results mapping in one of the supported formats.
    :type results: dict
    :returns: Impact-year results mapping.
    :rtype: dict[int, dict[str, Any]]
    """

    # --- Case 1: New explicit structure
    if isinstance(results, dict) and "results_by_impact_year" in results:
        iby = results.get("results_by_impact_year", {})
        if not isinstance(iby, dict):
            raise ValueError("'results_by_impact_year' exists but is not a dict.")
        return iby

    # --- Case 2: Already impact-year structure
    # Heuristic: keys are years (ints) and values have "scores_per_root" or "scores"
    if isinstance(results, dict) and results:
        k0 = next(iter(results.keys()))
        v0 = results[k0]
        if (
            isinstance(k0, int)
            and isinstance(v0, dict)
            and ("scores_per_root" in v0 or "scores" in v0)
        ):
            return results  # already normalized

    # --- Case 3: Legacy duplicated timeline under each solve-year
    # Find one exemplar entry that contains the legacy fields
    exemplar = None
    for _, v in (results or {}).items():
        if isinstance(v, dict) and "temporal_scores_per_root_by_year" in v:
            exemplar = v
            break

    if exemplar is None:
        raise ValueError(
            "Could not interpret results for plotting. Expected either:\n"
            "  (a) {'results_by_impact_year': ...}\n"
            "  (b) {impact_year: {'scores', 'scores_per_root'}}\n"
            "  (c) legacy {solve_year: {'temporal_scores_per_root_by_year', ...}}\n"
        )

    # Aggregate across solve years if needed
    out: Dict[int, Dict[str, Any]] = {}
    for solve_year, entry in results.items():
        if not isinstance(entry, dict):
            continue
        tspr = entry.get(
            "temporal_scores_per_root_by_year", {}
        )  # impact_year -> {root: score}
        tst = entry.get(
            "temporal_scores_by_year", {}
        )  # impact_year -> total score (optional)

        for impact_year, per_root in tspr.items():
            impact_year = int(impact_year)
            out.setdefault(impact_year, {"scores": 0.0, "scores_per_root": {}})

            # per-root sums
            for root, val in (per_root or {}).items():
                out[impact_year]["scores_per_root"][int(root)] = out[impact_year][
                    "scores_per_root"
                ].get(int(root), 0.0) + float(val)

        # total scores (if present)
        for impact_year, total in (tst or {}).items():
            impact_year = int(impact_year)
            out.setdefault(impact_year, {"scores": 0.0, "scores_per_root": {}})
            out[impact_year]["scores"] += float(total)

    return out


def _characterized_inventory_to_results(
    characterized_inventory: xr.DataArray,
    by_flow: bool = False,
) -> Dict[int, Dict[str, Any]]:
    """Convert a characterized inventory DataArray into impact-year results."""
    if "method" in characterized_inventory.dims:
        methods = characterized_inventory.coords["method"].values
        if len(methods) != 1:
            raise ValueError(
                "characterized_inventory has multiple methods; select one before plotting."
            )
        characterized_inventory = characterized_inventory.sel(method=methods[0])
    if "flow" not in characterized_inventory.dims:
        raise ValueError("characterized_inventory must include a 'flow' dimension.")
    if "activity" not in characterized_inventory.dims:
        raise ValueError(
            "characterized_inventory must include an 'activity' dimension."
        )
    if "year" not in characterized_inventory.dims:
        raise ValueError("characterized_inventory must include a 'year' dimension.")

    data = characterized_inventory.data
    if not isinstance(data, sparse.COO):
        data = sparse.COO.from_numpy(np.asarray(data))

    if by_flow:
        summed = data.sum(axis=0)  # flow x year
        score_key = "scores_by_flow"
    else:
        summed = data.sum(axis=1)  # activity x year
        score_key = "scores_by_first_level_child"
    years = characterized_inventory.coords["year"].values

    results: Dict[int, Dict[str, Any]] = {
        int(year): {"scores": 0.0, score_key: {}} for year in years
    }

    if isinstance(summed, sparse.COO):
        act_coords = summed.coords[0]
        year_coords = summed.coords[1]
        values = summed.data
        for act_idx, year_idx, val in zip(act_coords, year_coords, values):
            year = int(years[int(year_idx)])
            results[year][score_key][int(act_idx)] = float(val)
            results[year]["scores"] += float(val)
    else:
        for year_idx, year in enumerate(years):
            col = np.asarray(summed[:, year_idx]).ravel()
            if col.size == 0:
                continue
            year_result = results[int(year)]
            for act_idx, val in enumerate(col):
                if val == 0.0:
                    continue
                year_result[score_key][int(act_idx)] = float(val)
                year_result["scores"] += float(val)

    return results


def _characterized_inventory_to_root_results(
    characterized_inventory: xr.DataArray,
) -> Dict[int, Dict[str, Any]]:
    """Convert a characterized inventory with root activity dimension into results."""
    if "method" in characterized_inventory.dims:
        methods = characterized_inventory.coords["method"].values
        if len(methods) != 1:
            raise ValueError(
                "characterized_inventory has multiple methods; select one before plotting."
            )
        characterized_inventory = characterized_inventory.sel(method=methods[0])
    if "flow" not in characterized_inventory.dims:
        raise ValueError("characterized_inventory must include a 'flow' dimension.")
    if "activity" not in characterized_inventory.dims:
        raise ValueError(
            "characterized_inventory must include an 'activity' dimension."
        )
    if "year" not in characterized_inventory.dims:
        raise ValueError("characterized_inventory must include a 'year' dimension.")
    if "root activity" not in characterized_inventory.dims:
        raise ValueError(
            "characterized_inventory must include a 'root activity' dimension."
        )

    data = characterized_inventory.data
    if not isinstance(data, sparse.COO):
        data = sparse.COO.from_numpy(np.asarray(data))

    summed = data.sum(axis=(0, 1))  # year x root activity
    years = characterized_inventory.coords["year"].values
    roots = characterized_inventory.coords["root activity"].values
    score_key = "scores_by_first_level_child"

    results: Dict[int, Dict[str, Any]] = {
        int(year): {"scores": 0.0, score_key: {}} for year in years
    }

    if isinstance(summed, sparse.COO):
        year_coords = summed.coords[0]
        root_coords = summed.coords[1]
        values = summed.data
        for year_idx, root_idx, val in zip(year_coords, root_coords, values):
            year = int(years[int(year_idx)])
            root = int(roots[int(root_idx)])
            results[year][score_key][root] = float(val)
            results[year]["scores"] += float(val)
    else:
        for year_idx, year in enumerate(years):
            col = np.asarray(summed[year_idx, :]).ravel()
            if col.size == 0:
                continue
            year_result = results[int(year)]
            for root_idx, val in enumerate(col):
                if val == 0.0:
                    continue
                root = int(roots[int(root_idx)])
                year_result[score_key][root] = float(val)
                year_result["scores"] += float(val)

    return results


def _wrap_hover_label(text: str, max_chars: int = 45) -> str:
    """
    Wrap text for Plotly hoverlabels by inserting <br> at word boundaries.

    Plotly hoverlabels don't support fixed pixel widths reliably; wrapping the
    content is the most consistent way to keep hover popups narrow.
    """
    if not text:
        return ""

    # Normalize separators so we can wrap nicely around them
    # (optional but helps LCA-style labels)
    t = str(text).replace(" | ", " | ").replace("|", " | ")

    words = t.split()
    lines: list[str] = []
    current: list[str] = []
    current_len = 0

    for w in words:
        w_len = len(w) + (1 if current else 0)
        if current and (current_len + w_len) > max_chars:
            lines.append(" ".join(current))
            current = [w]
            current_len = len(w)
        else:
            current.append(w)
            current_len += w_len

    if current:
        lines.append(" ".join(current))

    # Extra nicety: encourage breaks after separators
    joined = "<br>".join(lines)
    joined = joined.replace(" | ", " |<br>")
    return joined


def _scores_to_results(
    scores: xr.DataArray,
    *,
    by_flow: bool = False,
) -> Dict[int, Dict[str, Any]]:
    """
    Convert a score DataArray into impact-year results.

    Supported dims (examples):
      - ("year",) totals only
      - ("activity","year") or ("year","activity")
      - ("root activity","year") or ("year","root activity")
      - ("activity","year","root activity") or permutations -> will sum over "activity"
      - ("flow","year") or ("year","flow") if by_flow=True

    Returns:
      impact_year -> {"scores": float, score_key: {idx: float}}
    """
    if "year" not in scores.dims:
        raise ValueError("scores must include a 'year' dimension.")

    years = [int(y) for y in scores.coords["year"].values.tolist()]

    # Pick what we attribute to
    if by_flow:
        attrib_dim = "flow" if "flow" in scores.dims else None
        score_key = "scores_by_flow"
    else:
        if "root activity" in scores.dims:
            attrib_dim = "root activity"
        elif "activity" in scores.dims:
            attrib_dim = "activity"
        elif "flow" in scores.dims:
            attrib_dim = "flow"
        else:
            attrib_dim = None
        score_key = "scores_by_first_level_child"

    data = scores

    # Totals only: sum everything except year
    if attrib_dim is None:
        extra_dims = [d for d in data.dims if d != "year"]
        if extra_dims:
            data = data.sum(dim=extra_dims)
        if data.dims != ("year",):
            data = data.transpose("year")
        results: Dict[int, Dict[str, Any]] = {}
        arr = data.data
        if hasattr(arr, "todense"):
            v = np.asarray(arr.todense(), dtype=float).ravel()
        else:
            v = np.asarray(data.values, dtype=float).ravel()
        for yi, year in enumerate(years):
            results[int(year)] = {"scores": float(v[yi]), score_key: {}}
        return results

    if attrib_dim == "root activity" and "activity" in data.dims:
        data = data.sum(dim="activity")

    if attrib_dim == "activity" and "root activity" in data.dims:
        data = data.sum(dim="root activity")

    extra_dims = [d for d in data.dims if d not in (attrib_dim, "year")]
    if extra_dims:
        data = data.sum(dim=extra_dims)

    if data.dims != (attrib_dim, "year"):
        data = data.transpose(attrib_dim, "year")

    arr = data.data
    if hasattr(arr, "todense"):
        vals = np.asarray(arr.todense(), dtype=float)
    else:
        vals = np.asarray(data.values, dtype=float)

    attrib_ids = data.coords[attrib_dim].values

    results: Dict[int, Dict[str, Any]] = {
        int(y): {"scores": 0.0, score_key: {}} for y in years
    }

    for yi, year in enumerate(years):
        col = vals[:, yi]
        per: Dict[int, float] = {}
        year_total = 0.0
        for ai, v in enumerate(col):
            if v == 0.0:
                continue
            idx = int(attrib_ids[ai])
            fv = float(v)
            per[idx] = fv
            year_total += fv
        results[int(year)][score_key] = per
        results[int(year)]["scores"] = year_total

    return results


def _aggregate_flow_results_by_name(
    results_by_year: Dict[int, Dict[str, Any]],
    trails: Trails,
) -> Dict[int, Dict[str, Any]]:
    """Aggregate per-flow results by biosphere flow name."""
    flow_id_to_name: dict[int, str] = {}
    for _label, meta in trails.biosphere_indices.items():
        for fid, md in meta.items():
            if fid in flow_id_to_name:
                continue
            name = md.get("name") if isinstance(md, dict) else None
            if name:
                flow_id_to_name[int(fid)] = str(name)

    out: Dict[int, Dict[str, Any]] = {}
    for year, payload in results_by_year.items():
        scores_by_flow = payload.get("scores_by_flow", {})
        agg: Dict[str, float] = {}
        total = 0.0
        for fid, val in scores_by_flow.items():
            name = flow_id_to_name.get(int(fid), f"Flow {fid}")
            agg[name] = agg.get(name, 0.0) + float(val)
            total += float(val)
        out[int(year)] = {"scores": total, "scores_by_flow": agg}
    return out


def _plot_results_by_year(
    results_by_year: Union[Dict[int, Dict[str, Any]], Dict[str, Any], xr.DataArray],
    trails: Trails,
    title: str,
    method_label: str,
    cumulative: bool,
    stacked: bool,
    legend_top_n: int,
    show_flow_contributions: bool,
    width: Optional[int],
    height: Optional[int],
    year_tick: int,
    year_range: Optional[Tuple[int, int]],
    show_year_grid: bool,
    yaxis_type: Literal["linear", "log"],
    log_eps: float,
    reference_year: Optional[int],
    show_cumulative_axis: bool,
    cumulative_axis_label: str,
    legend_entrywidth: int,
    legend_row_height: int,
    legend_y: float,
    y2_headroom: float,
    show_cumulative_in_legend: bool,
    static_score: Optional[float],
    static_score_label: str,
    static_score_dash: str,
    static_score_color: str,
    y_min: Optional[float],
    y_max: Optional[float],
    y2_max: Optional[float],
    *,
    flow_groupby_name: bool = False,
) -> go.Figure:

    if isinstance(results_by_year, xr.DataArray):
        # Inventory-style arrays have a "flow" dim; score arrays generally do not.
        if "flow" in results_by_year.dims:
            # Existing behavior: interpret as characterized inventory
            if "root activity" in results_by_year.dims:
                if show_flow_contributions:
                    tmp = results_by_year.sum(dim="root activity")
                    results_by_year = _characterized_inventory_to_results(
                        tmp, by_flow=True
                    )
                else:
                    results_by_year = _characterized_inventory_to_root_results(
                        results_by_year
                    )
            else:
                results_by_year = _characterized_inventory_to_results(
                    results_by_year, by_flow=show_flow_contributions
                )
        else:
            # Score array (e.g., trails.scores)
            if show_flow_contributions:
                # Only works if scores actually has a "flow" dim; otherwise will become totals-only.
                results_by_year = _scores_to_results(results_by_year, by_flow=True)
            else:
                results_by_year = _scores_to_results(results_by_year, by_flow=False)
    else:
        results_by_year = to_impact_year_results(results_by_year)

    if year_tick < 1:
        raise ValueError("year_tick must be >= 1")

    if show_flow_contributions:
        score_key = "scores_by_flow"
    else:
        if isinstance(results_by_year, dict):
            # If caller used a non-root data source that falls back to flow attribution
            # in _scores_to_results, detect that and switch the expected key.
            if any(
                isinstance(results_by_year.get(year), dict)
                and "scores_by_flow" in results_by_year.get(year, {})
                for year in results_by_year
            ):
                score_key = "scores_by_flow"
            else:
                score_key = "scores_by_first_level_child"
        else:
            score_key = "scores_by_first_level_child"
    if show_flow_contributions and flow_groupby_name:
        results_by_year = _aggregate_flow_results_by_name(results_by_year, trails)
    years = _select_years_from_results(results_by_year, year_range)
    if not any(score_key in results_by_year[year] for year in years):
        raise ValueError(f"Expected {score_key} in results_by_year for plotting.")

    all_roots = _collect_root_scores(results_by_year, years, score_key)

    Y_raw = _build_score_matrix(results_by_year, years, all_roots, score_key)
    if cumulative:
        Y = np.cumsum(Y_raw, axis=0)
    else:
        Y = Y_raw

    total_raw = Y.sum(axis=1)

    if yaxis_type == "log":
        Y = np.where(Y > 0, Y, log_eps)
        if static_score is not None:
            static_score = max(static_score, log_eps)

    if show_flow_contributions and flow_groupby_name:
        idx_to_label = {}
    else:
        idx_to_label = (
            _build_flow_label_map(trails)
            if show_flow_contributions
            else _build_activity_label_map(trails)
        )

    if legend_top_n < 0:
        raise ValueError("legend_top_n must be >= 0")
    if legend_top_n == 0:
        legend_roots = set()
    else:
        contributions = np.sum(np.abs(Y_raw), axis=0)
        top_count = min(legend_top_n, len(all_roots))
        top_idx = np.argsort(contributions)[::-1][:top_count]
        legend_roots = {all_roots[i] for i in top_idx}

    fig = go.Figure()
    _add_root_traces(
        fig=fig,
        years=years,
        Y=Y,
        all_roots=all_roots,
        idx_to_label=idx_to_label,
        method_label=method_label,
        stacked=stacked,
        showlegend_roots=legend_roots,
        showhover_roots=legend_roots,
    )

    cum_vals = None
    if show_cumulative_axis:
        cum_vals = _add_cumulative_trace(
            fig=fig,
            years=years,
            total_raw=total_raw,
            cumulative_axis_label=cumulative_axis_label,
            yaxis_type=yaxis_type,
            log_eps=log_eps,
        )

    if static_score is not None:
        _add_static_score_trace(
            fig=fig,
            years=years,
            static_score=static_score,
            static_score_label=static_score_label,
            static_score_dash=static_score_dash,
            static_score_color=static_score_color,
            method_label=method_label,
        )

    n_items = len(legend_roots)
    if show_cumulative_axis and show_cumulative_in_legend:
        n_items += 1
    if static_score is not None:
        n_items += 1

    entry_w, top_margin = _compute_layout_dimensions(
        width=width,
        legend_entrywidth=legend_entrywidth,
        legend_row_height=legend_row_height,
        n_items=n_items,
    )

    _apply_base_layout(
        fig=fig,
        width=width,
        height=height,
        title=title,
        legend_y=legend_y,
        entry_w=entry_w,
        top_margin=top_margin,
        method_label=method_label,
        yaxis_type=yaxis_type,
        show_cumulative_axis=show_cumulative_axis,
        static_score=static_score,
        cumulative_axis_label=cumulative_axis_label,
    )

    if (show_cumulative_axis or (static_score is not None)) and yaxis_type == "linear":
        _apply_linear_yaxis_alignment(
            fig=fig,
            Y=Y,
            cum_vals=cum_vals,
            static_score=static_score,
            y_min=y_min,
            y_max=y_max,
            y2_max=y2_max,
            y2_headroom=y2_headroom,
            stacked=stacked,
        )

    if yaxis_type != "linear" and (y_min is not None or y_max is not None):
        fig.update_layout(yaxis=dict(range=[y_min, y_max]))

    if yaxis_type != "linear" and y2_max is not None:
        use_y2 = bool(show_cumulative_axis) or (static_score is not None)
        if use_y2:
            fig.update_layout(yaxis2=dict(range=[None, float(y2_max)]))

    _apply_xaxis_settings(
        fig=fig,
        year_tick=year_tick,
        year_range=year_range,
        years=years,
        show_year_grid=show_year_grid,
    )

    _add_reference_year_line(fig, reference_year)

    return fig


def plot_temporal_scores(
    trails: Trails,
    title: str = "Temporal impacts by responsible activity",
    method_label: str = "Impact score",
    method: Optional[str] = None,
    cumulative: bool = False,
    stacked: bool = True,
    legend_top_n: int = 5,
    show_flow_contributions: bool = False,
    width: Optional[int] = None,
    height: Optional[int] = None,
    year_tick: int = 1,
    year_range: Optional[Tuple[int, int]] = None,
    show_year_grid: bool = True,
    yaxis_type: Literal["linear", "log"] = "linear",
    log_eps: float = 1e-30,
    reference_year: Optional[int] = None,
    show_cumulative_axis: bool = False,
    cumulative_axis_label: str = "Cumulative impact",
    legend_entrywidth: int = 260,
    legend_row_height: int = 18,
    legend_y: float = 1.02,
    y2_headroom: float = 0.05,
    show_cumulative_in_legend: bool = False,
    flow_groupby_name: bool = False,
    static_score: Optional[float] | dict[str, float] = None,
    static_score_label: str = "Static score",
    static_score_dash: str = "dash",
    static_score_color: str = "black",
    y_min: Optional[float] = None,
    y_max: Optional[float] = None,
    y2_max: Optional[float] = None,
) -> go.Figure | list[go.Figure]:
    """Plot temporal impact scores by responsible activity."""
    if trails.characterized_inventory is not None:
        results_by_year: Union[
            Dict[int, Dict[str, Any]], Dict[str, Any], xr.DataArray
        ] = trails.characterized_inventory
    elif getattr(trails, "scores", None) is not None:
        results_by_year = trails.scores
    else:
        raise ValueError("No characterized inventory or scores available for plotting.")

    if isinstance(results_by_year, xr.DataArray) and "method" in results_by_year.dims:
        methods = results_by_year.coords["method"].values.tolist()
        if method is None and len(methods) > 1:
            figures: list[go.Figure] = []
            for idx, m in enumerate(methods):
                selected = results_by_year.isel(method=idx, drop=True)
                score_for_method: Optional[float] = None
                if isinstance(static_score, dict):
                    if str(m) not in static_score:
                        raise ValueError(
                            f"Static score missing for method '{m}'. "
                            f"Available: {sorted(static_score)}"
                        )
                    score_for_method = float(static_score[str(m)])
                else:
                    score_for_method = static_score
                method_title = ""
                fig = _plot_results_by_year(
                    results_by_year=selected,
                    trails=trails,
                    title=method_title,
                    method_label=method_label,
                    cumulative=cumulative,
                    stacked=stacked,
                    legend_top_n=legend_top_n,
                    show_flow_contributions=show_flow_contributions,
                    width=width,
                    height=height,
                    year_tick=year_tick,
                    year_range=year_range,
                    show_year_grid=show_year_grid,
                    yaxis_type=yaxis_type,
                    log_eps=log_eps,
                    reference_year=reference_year,
                    show_cumulative_axis=show_cumulative_axis,
                    cumulative_axis_label=cumulative_axis_label,
                    legend_entrywidth=legend_entrywidth,
                    legend_row_height=legend_row_height,
                    legend_y=legend_y,
                    y2_headroom=y2_headroom,
                    show_cumulative_in_legend=show_cumulative_in_legend,
                    flow_groupby_name=flow_groupby_name,
                    static_score=score_for_method,
                    static_score_label=static_score_label,
                    static_score_dash=static_score_dash,
                    static_score_color=static_score_color,
                    y_min=y_min,
                    y_max=y_max,
                    y2_max=y2_max,
                )
                figures.append(fig)
            return figures

        if method is None:
            method = str(methods[0])
        if method not in methods:
            raise ValueError(
                f"Requested method '{method}' not found. Available: {methods}"
            )
        if isinstance(static_score, dict):
            if method not in static_score:
                raise ValueError(
                    f"Static score missing for method '{method}'. "
                    f"Available: {sorted(static_score)}"
                )
            static_score = float(static_score[method])
        method_idx = methods.index(method)
        results_by_year = results_by_year.isel(method=method_idx, drop=True)
    elif isinstance(static_score, dict):
        if len(static_score) == 1:
            static_score = float(next(iter(static_score.values())))
        else:
            raise ValueError(
                "Multiple static scores provided; pass method=... to select one."
            )

    return _plot_results_by_year(
        results_by_year=results_by_year,
        trails=trails,
        title=title,
        method_label=method_label,
        cumulative=cumulative,
        stacked=stacked,
        legend_top_n=legend_top_n,
        show_flow_contributions=show_flow_contributions,
        width=width,
        height=height,
        year_tick=year_tick,
        year_range=year_range,
        show_year_grid=show_year_grid,
        yaxis_type=yaxis_type,
        log_eps=log_eps,
        reference_year=reference_year,
        show_cumulative_axis=show_cumulative_axis,
        cumulative_axis_label=cumulative_axis_label,
        legend_entrywidth=legend_entrywidth,
        legend_row_height=legend_row_height,
        legend_y=legend_y,
        y2_headroom=y2_headroom,
        show_cumulative_in_legend=show_cumulative_in_legend,
        flow_groupby_name=flow_groupby_name,
        static_score=static_score,
        static_score_label=static_score_label,
        static_score_dash=static_score_dash,
        static_score_color=static_score_color,
        y_min=y_min,
        y_max=y_max,
        y2_max=y2_max,
    )


def plot_top_paths_for_year(
    provenance: Dict[tuple[int, int], Dict[tuple[int, ...], float]],
    trails: Trails,
    year: int,
    top_n: int = 10,
    title: str = "Top demand paths by amount",
    amount_label: str = "Demand (amount units)",
) -> go.Figure:
    """Visualize the top-N paths contributing to demand for a given year.

    :param provenance: Provenance mapping from temporal traversal.
    :type provenance: dict[tuple[int, int], dict[tuple[int, ...], float]]
    :param trails: Trails instance to resolve activity labels.
    :type trails: Trails
    :param year: Scenario year to visualize.
    :type year: int
    :param top_n: Number of paths to show.
    :type top_n: int
    :param title: Plot title.
    :type title: str
    :param amount_label: Label for the amount axis.
    :type amount_label: str
    :returns: Plotly figure.
    :rtype: plotly.graph_objects.Figure
    """
    from collections import defaultdict
    import numpy as np
    import plotly.graph_objects as go

    # Aggregate demand per path for the given year
    path_amounts = defaultdict(float)
    for (y, act_idx), path_map in provenance.items():
        if y != year:
            continue
        for path, amt in path_map.items():
            if not path:
                continue
            path_amounts[path] += float(amt)

    if not path_amounts:
        raise ValueError(f"No path provenance found for year {year}.")

    # Pick top-N by absolute amount
    sorted_paths = sorted(
        path_amounts.items(),
        key=lambda kv: abs(kv[1]),
        reverse=True,
    )[:top_n]

    idx_to_label = _build_activity_label_map(trails)

    labels = [_format_path_label(path, idx_to_label) for path, _ in sorted_paths]
    values = [amt for _, amt in sorted_paths]

    # Horizontal bar chart
    fig = go.Figure(
        data=[
            go.Bar(
                x=values,
                y=labels,
                orientation="h",
                hovertemplate=(
                    "<b>%{y}</b><br>" f"{amount_label}: %{{x:.6g}}<extra></extra>"
                ),
            )
        ]
    )

    fig.update_layout(
        title=title + f" (year {year})",
        xaxis_title=amount_label,
        yaxis_title="Path (first-level → ... → node)",
        template="plotly_white",
        margin=dict(l=200, r=20, t=60, b=40),
    )

    return fig


def _collect_activity_meta(trails_local: Trails) -> Dict[int, Dict[str, Any]]:
    """Collect activity metadata keyed by activity index.

    :param trails_local: Trails instance with activity indices.
    :type trails_local: Trails
    :returns: Mapping of activity index to metadata dict.
    :rtype: dict[int, dict[str, Any]]
    """
    meta_by_idx: Dict[int, Dict[str, Any]] = {}
    for scen_label, mapping in trails_local.activity_indices.items():
        for idx, meta in mapping.items():
            if idx not in meta_by_idx:
                meta_by_idx[idx] = meta
    return meta_by_idx


def _activity_label_from_meta(act_meta: Dict[int, Dict[str, Any]], act_idx: int) -> str:
    """Build an activity label from metadata for a given index.

    :param act_meta: Mapping of activity index to metadata.
    :type act_meta: dict[int, dict[str, Any]]
    :param act_idx: Activity index to label.
    :type act_idx: int
    :returns: Display label for the activity.
    :rtype: str
    """
    meta = act_meta.get(act_idx, {})
    name = meta.get("name", f"Activity {act_idx}")
    rp = meta.get("reference product") or ""
    loc = meta.get("location") or ""
    label = name
    if rp:
        label += f" | {rp}"
    if loc:
        label += f" ({loc})"
    return label


def _build_full_path_amounts(
    provenance: Dict[tuple[int, int], Dict[tuple[int, ...], float]],
    start_year: int,
    start_act_idx: int,
) -> Dict[Tuple[Tuple[int, int], ...], float]:
    """Build full provenance paths including the root node.

    :param provenance: Provenance mapping from traversal.
    :type provenance: dict[tuple[int, int], dict[tuple, float]]
    :param start_year: Root year.
    :type start_year: int
    :param start_act_idx: Root activity index.
    :type start_act_idx: int
    :returns: Mapping of full paths to amounts.
    :rtype: dict[tuple[tuple[int, int], ...], float]
    """
    full_path_amounts: Dict[Tuple[Tuple[int, int], ...], float] = defaultdict(float)
    root_node = (start_year, start_act_idx)

    for (year_final, act_final), path_map in provenance.items():
        for path, amt in path_map.items():
            if not path:
                continue
            full_path = (root_node,) + path
            full_path_amounts[full_path] += float(amt)

    if not full_path_amounts:
        raise ValueError("No non-empty paths in provenance; nothing to plot.")

    return full_path_amounts


def _select_paths(
    full_path_amounts: Dict[Tuple[Tuple[int, int], ...], float],
    top_n_paths: int | None,
) -> list[tuple[tuple[tuple[int, int], ...], float]]:
    """Select top-N paths from a full path mapping.

    :param full_path_amounts: Mapping of paths to amounts.
    :type full_path_amounts: dict
    :param top_n_paths: Maximum number of paths to return.
    :type top_n_paths: int | None
    :returns: Selected path items.
    :rtype: list[tuple[tuple, float]]
    """
    if top_n_paths is None:
        return list(full_path_amounts.items())
    return sorted(
        full_path_amounts.items(),
        key=lambda kv: abs(kv[1]),
        reverse=True,
    )[:top_n_paths]


def _build_depth_map(
    selected_paths: list[tuple[tuple[tuple[int, int], ...], float]],
) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int]]:
    """Build node depth mapping from selected paths.

    :param selected_paths: Selected path items.
    :type selected_paths: list[tuple[tuple, float]]
    :returns: Tuple of node keys and depth map.
    :rtype: tuple[list[tuple[int, int]], dict[tuple[int, int], int]]
    """
    node_keys: set[Tuple[int, int]] = set()
    depth_map: Dict[Tuple[int, int], int] = {}

    for full_path, _ in selected_paths:
        for depth, node in enumerate(full_path):
            node_keys.add(node)
            if node not in depth_map or depth < depth_map[node]:
                depth_map[node] = depth

    node_keys = sorted(node_keys)
    return node_keys, depth_map


def _aggregate_link_impacts(
    selected_paths: list[tuple[tuple[tuple[int, int], ...], float]],
    depth_map: dict[tuple[int, int], int],
    node_intensity: dict[tuple[int, int], float],
) -> tuple[
    dict[tuple[tuple[int, int], tuple[int, int]], float],
    dict[tuple[tuple[int, int], tuple[int, int]], dict[int, float]],
]:
    """Aggregate link impacts by depth/year nodes.

    :param selected_paths: Selected path items.
    :type selected_paths: list[tuple[tuple, float]]
    :param depth_map: Mapping of node to depth.
    :type depth_map: dict[tuple[int, int], int]
    :param node_intensity: Mapping of node to intensity.
    :type node_intensity: dict[tuple[int, int], float]
    :returns: Aggregated link impacts and activity contribution mapping.
    :rtype: tuple[dict, dict]
    """
    link_impact_agg: Dict[Tuple[Tuple[int, int], Tuple[int, int]], float] = defaultdict(
        float
    )
    edge_activity_contrib: Dict[
        Tuple[Tuple[int, int], Tuple[int, int]], Dict[int, float]
    ] = defaultdict(lambda: defaultdict(float))

    for full_path, amt in selected_paths:
        if len(full_path) < 2:
            continue

        for i in range(len(full_path) - 1):
            parent = full_path[i]
            child = full_path[i + 1]

            parent_depth = depth_map[parent]
            child_depth = depth_map[child]

            parent_year, _ = parent
            child_year, child_act = child

            intensity = node_intensity.get(child, 0.0)
            imp = float(amt) * float(intensity)

            src_agg = (parent_depth, parent_year)
            tgt_agg = (child_depth, child_year)

            link_impact_agg[(src_agg, tgt_agg)] += imp
            edge_activity_contrib[(src_agg, tgt_agg)][child_act] += imp

    if not link_impact_agg:
        raise ValueError("No aggregated link impacts; nothing to plot.")

    return link_impact_agg, edge_activity_contrib


def _build_agg_nodes(
    link_impact_agg: dict[tuple[tuple[int, int], tuple[int, int]], float],
) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int]]:
    """Build aggregated node list and index mapping.

    :param link_impact_agg: Aggregated link impacts.
    :type link_impact_agg: dict
    :returns: Tuple of aggregated nodes and index mapping.
    :rtype: tuple[list[tuple[int, int]], dict[tuple[int, int], int]]
    """
    agg_nodes: set[Tuple[int, int]] = set()
    for src_agg, tgt_agg in link_impact_agg.keys():
        agg_nodes.add(src_agg)
        agg_nodes.add(tgt_agg)

    agg_nodes = sorted(agg_nodes)
    node_index_agg: Dict[Tuple[int, int], int] = {
        key: i for i, key in enumerate(agg_nodes)
    }
    return agg_nodes, node_index_agg


def _compute_node_totals(
    link_impact_agg: dict[tuple[tuple[int, int], tuple[int, int]], float],
) -> dict[tuple[int, int], float]:
    """Compute total incoming impact per aggregated node.

    :param link_impact_agg: Aggregated link impacts.
    :type link_impact_agg: dict
    :returns: Mapping of node to total impact.
    :rtype: dict[tuple[int, int], float]
    """
    node_total_impact: Dict[Tuple[int, int], float] = defaultdict(float)
    for (src_agg, tgt_agg), imp in link_impact_agg.items():
        node_total_impact[tgt_agg] += imp
    return node_total_impact


def _compute_sankey_layout(
    agg_nodes: list[tuple[int, int]],
) -> tuple[list[float], list[float], list[int], list[int]]:
    """Compute Sankey node positions from aggregated nodes.

    :param agg_nodes: Aggregated nodes as ``(depth, year)`` tuples.
    :type agg_nodes: list[tuple[int, int]]
    :returns: Tuple of node coordinates and sorted depth/year lists.
    :rtype: tuple[list[float], list[float], list[int], list[int]]
    """
    depths = sorted({d for (d, y) in agg_nodes})
    years = sorted({y for (d, y) in agg_nodes})

    if len(depths) == 1:
        depth_to_x = {depths[0]: 0.5}
    else:
        depth_to_x = {
            d: 0.05 + 0.9 * (i / (len(depths) - 1)) for i, d in enumerate(depths)
        }

    if len(years) == 1:
        year_to_y = {years[0]: 0.5}
    else:
        year_to_y = {
            y: 0.05 + 0.9 * (i / (len(years) - 1)) for i, y in enumerate(years)
        }

    node_x = [depth_to_x[d] for (d, y) in agg_nodes]
    node_y = [year_to_y[y] for (d, y) in agg_nodes]

    return node_x, node_y, depths, years


def _build_node_labels(
    agg_nodes: list[tuple[int, int]],
    node_total_impact: dict[tuple[int, int], float],
    amount_label: str,
) -> list[str]:
    """Build labels for Sankey nodes.

    :param agg_nodes: Aggregated nodes as ``(depth, year)`` tuples.
    :type agg_nodes: list[tuple[int, int]]
    :param node_total_impact: Total impact per node.
    :type node_total_impact: dict[tuple[int, int], float]
    :param amount_label: Label for impact values.
    :type amount_label: str
    :returns: Node label strings.
    :rtype: list[str]
    """
    node_labels: List[str] = []
    for d, y in agg_nodes:
        total_imp = node_total_impact.get((d, y), 0.0)
        node_labels.append(
            f"Depth {d}, Year {y}<br>" f"Incoming {amount_label}: {total_imp:.3g}"
        )
    return node_labels


def _assign_year_colors(years: list[int]) -> dict[int, str]:
    """Assign colors to years using a sequential palette.

    :param years: Years to assign colors for.
    :type years: list[int]
    :returns: Mapping of year to color string.
    :rtype: dict[int, str]
    """
    year_palette = px.colors.sequential.Viridis
    if len(years) > len(year_palette):
        repeats = (len(years) // len(year_palette)) + 1
        full_year_palette = (year_palette * repeats)[: len(years)]
    else:
        full_year_palette = year_palette[: len(years)]

    year_to_color = {y: col for y, col in zip(years, full_year_palette)}
    return year_to_color


def _build_link_arrays(
    link_impact_agg: dict[tuple[tuple[int, int], tuple[int, int]], float],
    edge_activity_contrib: dict[
        tuple[tuple[int, int], tuple[int, int]], dict[int, float]
    ],
    node_index_agg: dict[tuple[int, int], int],
    year_to_color: dict[int, str],
    amount_label: str,
    activity_label: Callable[[int], str],
) -> tuple[list[int], list[int], list[float], list[str], list[str]]:
    """Build Sankey link arrays from aggregated impact data.

    :param link_impact_agg: Aggregated link impacts.
    :type link_impact_agg: dict
    :param edge_activity_contrib: Activity contributions per link.
    :type edge_activity_contrib: dict
    :param node_index_agg: Mapping of aggregated node to index.
    :type node_index_agg: dict
    :param year_to_color: Mapping of year to color.
    :type year_to_color: dict[int, str]
    :param amount_label: Label for impact values.
    :type amount_label: str
    :param activity_label: Function to label activity indices.
    :type activity_label: callable
    :returns: Link source/target/value/color/hovertemplate arrays.
    :rtype: tuple[list[int], list[int], list[float], list[str], list[str]]
    """
    link_sources: List[int] = []
    link_targets: List[int] = []
    link_values: List[float] = []
    link_colors: List[str] = []
    link_hovertemplates: List[str] = []

    for (src_agg, tgt_agg), imp in link_impact_agg.items():
        src_idx = node_index_agg[src_agg]
        tgt_idx = node_index_agg[tgt_agg]

        src_depth, src_year = src_agg
        tgt_depth, tgt_year = tgt_agg

        color = year_to_color.get(tgt_year, "rgba(150,150,150,0.7)")

        activity_map = edge_activity_contrib[(src_agg, tgt_agg)]
        sorted_acts = sorted(
            activity_map.items(), key=lambda kv: abs(kv[1]), reverse=True
        )
        top_acts = sorted_acts[:8]

        lines = [
            f"<b>Depth {src_depth}, Year {src_year}</b> → "
            f"<b>Depth {tgt_depth}, Year {tgt_year}</b>",
            f"Total {amount_label}: {imp:.6g}",
            "<br><b>Top contributing activities at target:</b>",
        ]
        if not top_acts:
            lines.append("(none)")
        else:
            for act_idx, act_imp in top_acts:
                lines.append(f"- {activity_label(act_idx)}: {act_imp:.6g}")

        hovertemplate = "<br>".join(lines) + "<extra></extra>"

        link_sources.append(src_idx)
        link_targets.append(tgt_idx)
        link_values.append(abs(imp))
        link_colors.append(color)
        link_hovertemplates.append(hovertemplate)

    return link_sources, link_targets, link_values, link_colors, link_hovertemplates


def plot_temporal_sankey(
    provenance: Dict[tuple[int, int], Dict[tuple[tuple[int, int], ...], float]],
    trails: Trails,
    start_year: int,
    start_act_idx: int,
    node_intensity: Optional[Dict[Tuple[int, int], float]] = None,
    top_n_paths: int | None = 30,
    title: str = "Temporal Sankey (impact-weighted, aggregated by year)",
    amount_label: str = "Impact score",
    fig_width: int = 1200,
    fig_height: int = 800,
    node_thickness: int = 20,
    node_pad: int = 15,
    font_size: int = 11,
) -> go.Figure:
    """Build an impact-weighted temporal Sankey plot.

    Nodes are aggregated by ``(depth, year)`` and links are impact-weighted
    sums between those nodes, preserving activity-level contributions in
    hover text.

    :param provenance: Provenance mapping from traversal.
    :type provenance: dict[tuple[int, int], dict[tuple[tuple[int, int], ...], float]]
    :param trails: Trails instance for metadata.
    :type trails: Trails
    :param start_year: Root demand year.
    :type start_year: int
    :param start_act_idx: Root activity index.
    :type start_act_idx: int
    :param node_intensity: Impact intensities keyed by ``(year, act_idx)``.
    :type node_intensity: dict[tuple[int, int], float] | None
    :param top_n_paths: Number of paths to include (``None`` for all).
    :type top_n_paths: int | None
    :param title: Figure title.
    :type title: str
    :param amount_label: Label for impact values.
    :type amount_label: str
    :param fig_width: Figure width in pixels.
    :type fig_width: int
    :param fig_height: Figure height in pixels.
    :type fig_height: int
    :param node_thickness: Node thickness in pixels.
    :type node_thickness: int
    :param node_pad: Node padding in pixels.
    :type node_pad: int
    :param font_size: Font size for labels.
    :type font_size: int
    :returns: Plotly figure.
    :rtype: plotly.graph_objects.Figure
    """

    full_path_amounts = _build_full_path_amounts(provenance, start_year, start_act_idx)
    selected_paths = _select_paths(full_path_amounts, top_n_paths)
    node_keys, depth_map = _build_depth_map(selected_paths)
    if node_intensity is None:
        raise ValueError(
            "node_intensity is required; provide impact intensities keyed by (year, act_idx)."
        )

    act_meta = _collect_activity_meta(trails)
    activity_label = lambda act_idx: _activity_label_from_meta(act_meta, act_idx)

    link_impact_agg, edge_activity_contrib = _aggregate_link_impacts(
        selected_paths=selected_paths,
        depth_map=depth_map,
        node_intensity=node_intensity,
    )

    agg_nodes, node_index_agg = _build_agg_nodes(link_impact_agg)
    node_total_impact = _compute_node_totals(link_impact_agg)

    node_x, node_y, depths, years = _compute_sankey_layout(agg_nodes)
    node_labels = _build_node_labels(agg_nodes, node_total_impact, amount_label)
    year_to_color = _assign_year_colors(years)
    node_colors = [year_to_color[y] for (d, y) in agg_nodes]

    link_sources, link_targets, link_values, link_colors, link_hovertemplates = (
        _build_link_arrays(
            link_impact_agg=link_impact_agg,
            edge_activity_contrib=edge_activity_contrib,
            node_index_agg=node_index_agg,
            year_to_color=year_to_color,
            amount_label=amount_label,
            activity_label=activity_label,
        )
    )

    sankey = go.Sankey(
        domain=dict(x=[0.0, 1.0], y=[0.0, 1.0]),
        arrangement="snap",
        node=dict(
            pad=node_pad,
            thickness=node_thickness,
            label=node_labels,
            x=node_x,
            y=node_y,
            color=node_colors,
        ),
        link=dict(
            source=link_sources,
            target=link_targets,
            value=link_values,
            color=link_colors,
            hovertemplate=link_hovertemplates,  # list of per-link templates
        ),
    )

    fig = go.Figure(sankey)
    fig.update_layout(
        title=title,
        width=fig_width,
        height=fig_height,
        margin=dict(l=40, r=40, t=60, b=60),
        font=dict(size=font_size),
    )

    return fig


def build_sankey_arrays_from_tree(
    tree: dict[str, Any],
    *,
    label_fn: Callable[[dict[str, Any]], str] | None = None,
    edge_weight: Literal["amount", "score"] = "score",
    node_scores: dict[tuple[int, int], float] | None = None,
    trails: Trails | None = None,
    score_years: list[int] | None = None,
) -> dict[str, list]:
    """Convert a nested tree from lca.build_temporal_sankey_tree into Sankey arrays.

    :param tree: Nested tree dict with "node" and "children" entries.
    :type tree: dict
    :param label_fn: Optional function to build labels from a node payload.
    :type label_fn: callable | None
    :param edge_weight: Use "amount" (technosphere) or "score" (characterized).
    :type edge_weight: Literal["amount", "score"]
    :param node_scores: Optional mapping (year, act_idx) -> characterized score.
    :type node_scores: dict[tuple[int, int], float] | None
    :param trails: Trails instance used to derive node_scores from trails.scores.
    :type trails: Trails | None
    :param score_years: Optional list of years to align scores to node years.
    :type score_years: list[int] | None
    :returns: Dict with "labels", "sources", "targets", "values", "node_meta".
    :rtype: dict
    """
    if not tree or "node" not in tree:
        raise ValueError("Tree is empty or missing 'node' key.")

    def _default_label(node: dict[str, Any]) -> str:
        name = node.get("name") or f"Activity {node.get('act_idx')}"
        rp = node.get("reference_product") or ""
        loc = node.get("location") or ""
        label = name
        if rp:
            label += f" | {rp}"
        if loc:
            label += f" ({loc})"
        return label

    label_fn = label_fn or _default_label

    node_index: dict[Any, int] = {}
    node_meta: list[dict[str, Any]] = []
    labels: list[str] = []
    sources: list[int] = []
    targets: list[int] = []
    values: list[float] = []
    edge_scores: list[float] = []

    def _node_key(node_payload: dict[str, Any]) -> Any:
        key = node_payload.get("key")
        if key is None:
            return (
                node_payload.get("year"),
                node_payload.get("depth"),
                node_payload.get("act_idx"),
            )
        return key

    if edge_weight == "score":
        if node_scores is None:
            if trails is None:
                raise ValueError(
                    "edge_weight='score' requires node_scores or a Trails instance."
                )
            node_scores = _node_scores_from_trails(trails)
        node_scores = node_scores or {}

        # Determine if node_scores uses graph node keys or (year, act_idx)
        uses_node_keys = False
        if node_scores:
            first_key = next(iter(node_scores.keys()))
            if isinstance(first_key, tuple) and len(first_key) >= 4:
                uses_node_keys = True

        if not uses_node_keys:
            if score_years is None:
                if trails is None:
                    raise ValueError(
                        "score_years is required when trails is not provided."
                    )
                score_years = _score_years_from_trails(trails)
            score_years = [int(y) for y in score_years]
            score_years.sort()

    incoming_abs: dict[Any, float] = defaultdict(float)
    year_map: dict[int, int] = {}

    def _precompute_incoming(subtree: dict[str, Any]) -> None:
        for child in subtree.get("children", []):
            edge_amt = float(child.get("edge_amount") or 0.0)
            child_payload = child.get("node") or {}
            key = _node_key(child_payload)
            incoming_abs[key] += abs(edge_amt)
            _precompute_incoming(child)

    if edge_weight == "score":
        _precompute_incoming(tree)
        if score_years:
            for node in _iter_tree_nodes(tree):
                y = int(node.get("year"))
                if y in year_map:
                    continue
                year_map[y] = _nearest_year(score_years, y)

    def _get_node_id(node_payload: dict[str, Any]) -> int:
        key = _node_key(node_payload)
        if key in node_index:
            return node_index[key]
        idx = len(node_index)
        node_index[key] = idx
        node_meta.append(node_payload)
        labels.append(label_fn(node_payload))
        return idx

    def _walk(subtree: dict[str, Any]) -> None:
        parent_payload = subtree["node"]
        parent_id = _get_node_id(parent_payload)
        for child in subtree.get("children", []):
            edge_amt = float(child.get("edge_amount") or 0.0)
            child_payload = child.get("node") or {}
            child_id = _get_node_id(child_payload)
            sources.append(parent_id)
            targets.append(child_id)
            if edge_weight == "score":
                year = int(child_payload.get("year"))
                act_idx = int(child_payload.get("act_idx"))
                node_key = child_payload.get("key")
                if node_key in node_scores:
                    child_score = float(node_scores.get(node_key, 0.0))
                else:
                    year_lookup = year_map.get(year, year)
                    child_score = float(node_scores.get((year_lookup, act_idx), 0.0))
                denom = float(incoming_abs.get(_node_key(child_payload), 0.0))
                if denom > 0.0:
                    edge_score = child_score * (abs(edge_amt) / denom)
                else:
                    edge_score = 0.0
                edge_scores.append(edge_score)
                values.append(abs(edge_score))
            else:
                values.append(abs(edge_amt))
                edge_scores.append(float(edge_amt))
            _walk(child)

    _walk(tree)

    return {
        "labels": labels,
        "sources": sources,
        "targets": targets,
        "values": values,
        "node_meta": node_meta,
        "edge_scores": edge_scores,
    }


def _iter_tree_nodes(tree: dict[str, Any]) -> list[dict[str, Any]]:
    """Return a flat list of node payloads from a Sankey tree."""
    nodes: list[dict[str, Any]] = []

    def _walk(subtree: dict[str, Any]) -> None:
        nodes.append(subtree["node"])
        for child in subtree.get("children", []):
            _walk(child)

    _walk(tree)
    return nodes


def _nearest_year(years_sorted: list[int], target: int) -> int:
    """Return nearest year from a sorted list."""
    if not years_sorted:
        return int(target)
    i = bisect.bisect_left(years_sorted, target)
    if i <= 0:
        return years_sorted[0]
    if i >= len(years_sorted):
        return years_sorted[-1]
    before = years_sorted[i - 1]
    after = years_sorted[i]
    if abs(target - before) <= abs(after - target):
        return before
    return after


def _score_years_from_trails(trails: Trails) -> list[int]:
    """Return available score years from trails.scores or characterized inventory."""
    scores = getattr(trails, "scores", None)
    if scores is not None and "year" in scores.coords:
        return [int(y) for y in scores.coords["year"].values.tolist()]
    characterized = getattr(trails, "characterized_inventory", None)
    if characterized is not None and "year" in characterized.coords:
        return [int(y) for y in characterized.coords["year"].values.tolist()]
    return []


def _node_scores_from_trails(trails: Trails) -> dict[tuple[int, int], float]:
    """Build (year, act_idx) -> score mapping from trails.scores."""
    scores = getattr(trails, "scores", None)
    if scores is None:
        characterized = getattr(trails, "characterized_inventory", None)
        if characterized is None:
            raise ValueError(
                "No scores available. Run lca(..., compute_score=True) or "
                "build characterized_inventory first."
            )
        return _node_scores_from_characterized_inventory(characterized)
    if "activity" not in scores.dims or "year" not in scores.dims:
        raise ValueError("trails.scores must include 'activity' and 'year' dims.")

    if "root activity" in scores.dims:
        scores = scores.sum(dim="root activity")

    years = scores.coords["year"].values
    activities = scores.coords["activity"].values

    data = scores.data
    out: dict[tuple[int, int], float] = {}

    if isinstance(data, sparse.COO):
        coords = data.coords
        vals = data.data
        for ai, yi, v in zip(coords[0], coords[1], vals):
            if v == 0.0:
                continue
            year = int(years[int(yi)])
            act = int(activities[int(ai)])
            out[(year, act)] = out.get((year, act), 0.0) + float(v)
        return out

    dense = np.asarray(data)
    if dense.ndim != 2:
        raise ValueError("Expected 2D scores array after root aggregation.")
    idxs = np.nonzero(dense)
    for ai, yi in zip(idxs[0], idxs[1]):
        v = float(dense[ai, yi])
        if v == 0.0:
            continue
        year = int(years[int(yi)])
        act = int(activities[int(ai)])
        out[(year, act)] = out.get((year, act), 0.0) + v

    return out


def _node_scores_from_characterized_inventory(
    characterized_inventory: xr.DataArray,
) -> dict[tuple[int, int], float]:
    """Build (year, act_idx) -> score mapping from characterized inventory."""
    if "method" in characterized_inventory.dims:
        methods = characterized_inventory.coords["method"].values
        if len(methods) != 1:
            raise ValueError(
                "characterized_inventory has multiple methods; select one before plotting."
            )
        characterized_inventory = characterized_inventory.sel(method=methods[0])
    if "activity" not in characterized_inventory.dims:
        raise ValueError(
            "characterized_inventory must include an 'activity' dimension."
        )
    if "year" not in characterized_inventory.dims:
        raise ValueError("characterized_inventory must include a 'year' dimension.")
    if "flow" not in characterized_inventory.dims:
        raise ValueError("characterized_inventory must include a 'flow' dimension.")

    data = characterized_inventory.data
    years = characterized_inventory.coords["year"].values
    activities = characterized_inventory.coords["activity"].values

    out: dict[tuple[int, int], float] = {}

    if "root activity" in characterized_inventory.dims:
        if isinstance(data, sparse.COO):
            coords = data.coords
            vals = data.data
            for ai, fi, yi, ri, v in zip(
                coords[0], coords[1], coords[2], coords[3], vals
            ):
                if v == 0.0:
                    continue
                year = int(years[int(yi)])
                act = int(activities[int(ai)])
                out[(year, act)] = out.get((year, act), 0.0) + float(v)
        else:
            dense = np.asarray(data)
            dense = dense.sum(axis=3)
            idxs = np.nonzero(dense)
            for ai, fi, yi in zip(idxs[0], idxs[1], idxs[2]):
                v = float(dense[ai, fi, yi])
                if v == 0.0:
                    continue
                year = int(years[int(yi)])
                act = int(activities[int(ai)])
                out[(year, act)] = out.get((year, act), 0.0) + v
        return out

    if isinstance(data, sparse.COO):
        coords = data.coords
        vals = data.data
        for ai, fi, yi, v in zip(coords[0], coords[1], coords[2], vals):
            if v == 0.0:
                continue
            year = int(years[int(yi)])
            act = int(activities[int(ai)])
            out[(year, act)] = out.get((year, act), 0.0) + float(v)
        return out

    dense = np.asarray(data)
    if dense.ndim != 3:
        raise ValueError("Expected 3D characterized inventory after root aggregation.")
    idxs = np.nonzero(dense)
    for ai, fi, yi in zip(idxs[0], idxs[1], idxs[2]):
        v = float(dense[ai, fi, yi])
        if v == 0.0:
            continue
        year = int(years[int(yi)])
        act = int(activities[int(ai)])
        out[(year, act)] = out.get((year, act), 0.0) + v

    return out


def _select_depths(
    edges_by_depth: dict[int, dict], depths: list[int] | None
) -> list[int]:
    """Select which traversal depths to plot.

    :param edges_by_depth: Mapping of depth to edges.
    :type edges_by_depth: dict[int, dict]
    :param depths: Optional list of depths to include.
    :type depths: list[int] | None
    :returns: Filtered list of depths.
    :rtype: list[int]
    """
    if depths is None:
        depths_list = sorted(edges_by_depth.keys())
    else:
        depths_list = [d for d in depths if d in edges_by_depth]

    if not depths_list:
        raise ValueError("No depths to plot (depth list empty or no edges).")

    return depths_list


def _collect_activities(
    edges_by_depth: dict[int, dict],
    trails: Trails,
    depths_list: list[int],
    include_all_activities: bool,
) -> list[int]:
    """Collect activity indices referenced by selected edges.

    :param edges_by_depth: Mapping of depth to edges.
    :type edges_by_depth: dict[int, dict]
    :param trails: Trails instance for metadata.
    :type trails: Trails
    :param depths_list: Depths to include.
    :type depths_list: list[int]
    :param include_all_activities: Whether to include all activities from metadata.
    :type include_all_activities: bool
    :returns: Sorted list of activity indices.
    :rtype: list[int]
    """
    if include_all_activities:
        all_activities = set()
        for scen_label, mapping in trails.activity_indices.items():
            all_activities.update(mapping.keys())
    else:
        all_activities = set()
        for d in depths_list:
            for (parent, child), amt in edges_by_depth.get(d, {}).items():
                y_cons, a_cons = parent
                y_sup, a_sup = child
                all_activities.add(int(a_cons))
                all_activities.add(int(a_sup))

    acts = sorted(all_activities)
    if not acts:
        raise ValueError("No activities to plot for the selected depths.")

    return acts


def _collect_global_years(
    edges_by_depth: dict[int, dict], depths_list: list[int]
) -> tuple[int, int, list[int]]:
    """Collect all years referenced by selected depths.

    :param edges_by_depth: Mapping of depth to edges.
    :type edges_by_depth: dict[int, dict]
    :param depths_list: Depths to include.
    :type depths_list: list[int]
    :returns: Tuple of ``(year_min, year_max, years_global)``.
    :rtype: tuple[int, int, list[int]]
    """
    years_global_set = set()
    for d in depths_list:
        for (parent, child), amt in edges_by_depth.get(d, {}).items():
            y_cons, a_cons = parent
            y_sup, a_sup = child
            years_global_set.add(int(y_cons))
            years_global_set.add(int(y_sup))

    if not years_global_set:
        raise ValueError("No years found in edges for selected depths.")

    year_min = min(years_global_set)
    year_max = max(years_global_set)
    years_global = list(range(year_min, year_max + 1))

    return year_min, year_max, years_global


def _init_flow_subplots(panel_labels: list[str], ncols: int) -> tuple[go.Figure, int]:
    """Initialize subplot grid for flow panels.

    :param panel_labels: Labels for each panel.
    :type panel_labels: list[str]
    :param ncols: Number of columns.
    :type ncols: int
    :returns: Plotly figure with subplots.
    :rtype: plotly.graph_objects.Figure
    """
    n_panels = len(panel_labels)
    nrows = int(np.ceil(n_panels / ncols))
    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        subplot_titles=panel_labels,
        shared_yaxes=False,
        shared_xaxes=True,
        horizontal_spacing=0.08,
        vertical_spacing=0.12,
    )
    return fig, nrows


def plot_rf(
    trails: Trails,
    *,
    by: Literal["flow", "root activity"] = "root activity",
    title: str = "Radiative forcing by gas/flow",
    method_label: str = "W/m²",
    quantile: float | None = 50.0,
    show_cumulative_quantile_band: bool = True,
    band_quantiles: tuple[float, float] = (2.5, 97.5),
    cumulative: bool = False,
    stacked: bool = True,
    legend_top_n: int = 5,
    width: Optional[int] = 550,
    height: Optional[int] = 450,
    year_tick: int = 5,
    year_range: Optional[Tuple[int, int]] = None,
    show_year_grid: bool = True,
    yaxis_type: Literal["linear", "log"] = "linear",
    log_eps: float = 1e-30,
    reference_year: Optional[int] = None,
    show_cumulative_axis: bool = True,
    cumulative_axis_label: str = "W·m⁻²·yr",
    legend_entrywidth: int = 260,
    legend_row_height: int = 18,
    legend_y: float = 1.0,
    y2_headroom: float = 0.05,
    show_cumulative_in_legend: bool = False,
    flow_groupby_name: bool = False,
    y_min: Optional[float] = None,
    y_max: Optional[float] = None,
    y2_max: Optional[float] = None,
) -> go.Figure:
    """Plot radiative forcing time series by flow or root activity."""
    rf = getattr(trails, "instant_radiative_forcing", None)
    if rf is None:
        raise ValueError("No radiative forcing data stored on Trails.")

    if "year" not in rf.dims:
        raise ValueError("RF data must include a 'year' dimension.")

    rf_all = rf
    if "quantile" in rf.dims:
        quantiles = [float(q) for q in rf.coords["quantile"].values.tolist()]
        if quantile is None:
            quantile = 50.0 if 50.0 in quantiles else quantiles[0]
        if float(quantile) not in quantiles:
            raise ValueError(
                f"Requested quantile {quantile} not in available {quantiles}."
            )
        rf = rf.sel(quantile=float(quantile), drop=True)
    elif show_cumulative_quantile_band:
        show_cumulative_quantile_band = False

    if by == "flow":
        if "flow" not in rf.dims:
            raise ValueError("RF data must include a 'flow' dimension for by='flow'.")
        if "root activity" in rf.dims:
            rf = rf.sum(dim="root activity")
        if "activity" in rf.dims:
            rf = rf.sum(dim="activity")
        results_by_year = _scores_to_results(rf, by_flow=True)
    elif by == "root activity":
        if "root activity" not in rf.dims:
            raise ValueError(
                "RF data must include 'root activity' for by='root activity'."
            )
        data = rf
        if "flow" in data.dims:
            data = data.sum(dim="flow")
        if "activity" in data.dims:
            data = data.sum(dim="activity")
        results_by_year = _scores_to_results(data, by_flow=False)
    else:
        raise ValueError("by must be 'flow' or 'root activity'.")

    fig = _plot_results_by_year(
        results_by_year=results_by_year,
        trails=trails,
        title=title,
        method_label=method_label,
        cumulative=cumulative,
        stacked=stacked,
        legend_top_n=legend_top_n,
        show_flow_contributions=(by == "flow"),
        width=width,
        height=height,
        year_tick=year_tick,
        year_range=year_range,
        show_year_grid=show_year_grid,
        yaxis_type=yaxis_type,
        log_eps=log_eps,
        reference_year=reference_year,
        show_cumulative_axis=show_cumulative_axis,
        cumulative_axis_label=cumulative_axis_label,
        legend_entrywidth=legend_entrywidth,
        legend_row_height=legend_row_height,
        legend_y=legend_y,
        y2_headroom=y2_headroom,
        show_cumulative_in_legend=show_cumulative_in_legend,
        flow_groupby_name=flow_groupby_name,
        static_score=None,
        static_score_label="Static score",
        static_score_dash="dash",
        static_score_color="black",
        y_min=y_min,
        y_max=y_max,
        y2_max=y2_max,
    )

    if show_cumulative_axis and show_cumulative_quantile_band:
        q_low, q_high = band_quantiles
        quantiles = [float(q) for q in rf_all.coords["quantile"].values.tolist()]
        if q_low not in quantiles or q_high not in quantiles:
            raise ValueError(
                f"band_quantiles {band_quantiles} not in available {quantiles}."
            )
        q_low_data = rf_all.sel(quantile=float(q_low), drop=True)
        q_high_data = rf_all.sel(quantile=float(q_high), drop=True)

        def _results_for(data: xr.DataArray) -> dict[int, dict[str, Any]]:
            if by == "flow":
                if "root activity" in data.dims:
                    data = data.sum(dim="root activity")
                if "activity" in data.dims:
                    data = data.sum(dim="activity")
                return _scores_to_results(data, by_flow=True)
            if "root activity" in data.dims:
                data = data.sum(dim="flow")
            if "activity" in data.dims:
                data = data.sum(dim="activity")
            return _scores_to_results(data, by_flow=False)

        results_low = _results_for(q_low_data)
        results_high = _results_for(q_high_data)
        years = _select_years_from_results(results_by_year, year_range)

        def _totals(res: dict[int, dict[str, Any]], years_seq: list[int]) -> np.ndarray:
            out = np.zeros(len(years_seq), dtype=float)
            for i, y in enumerate(years_seq):
                payload = res.get(int(y), {})
                out[i] = float(payload.get("scores", 0.0))
            return out

        total_low = _totals(results_low, years)
        total_high = _totals(results_high, years)
        if show_cumulative_axis:
            total_low = np.cumsum(total_low)
            total_high = np.cumsum(total_high)
        if yaxis_type == "log":
            total_low = np.where(total_low > 0, total_low, log_eps)
            total_high = np.where(total_high > 0, total_high, log_eps)
        if not np.any(total_low) and not np.any(total_high):
            fig.update_layout(title=dict(text=""))
            fig.add_annotation(
                text=title,
                x=0.5,
                y=-0.2,
                xref="paper",
                yref="paper",
                xanchor="center",
                yanchor="top",
                showarrow=False,
            )
            return fig
        fig.add_trace(
            go.Scatter(
                x=years,
                y=total_low,
                mode="lines",
                line=dict(width=0),
                showlegend=False,
                hoverinfo="skip",
                yaxis="y2",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=years,
                y=total_high,
                mode="lines",
                fill="tonexty",
                fillcolor="rgba(60,60,60,0.35)",
                line=dict(width=0),
                name=f"{q_low:g}-{q_high:g}th percentile",
                yaxis="y2",
            )
        )
    fig.update_layout(title=dict(text=""))
    fig.add_annotation(
        text=title,
        x=0.5,
        y=-0.2,
        xref="paper",
        yref="paper",
        xanchor="center",
        yanchor="top",
        showarrow=False,
    )
    fig.update_yaxes(
        tickformat=".2e",
        exponentformat="e",
        showexponent="all",
    )
    return fig


def plot_temp(
    trails: Trails,
    *,
    by: Literal["flow", "root activity"] = "root activity",
    title: str = "Temperature change by gas/flow",
    method_label: str = "°C",
    quantile: float | None = 50.0,
    stacked: bool = True,
    legend_top_n: int = 5,
    width: Optional[int] = 550,
    height: Optional[int] = 450,
    year_tick: int = 5,
    year_range: Optional[Tuple[int, int]] = None,
    show_year_grid: bool = True,
    yaxis_type: Literal["linear", "log"] = "linear",
    log_eps: float = 1e-30,
    reference_year: Optional[int] = None,
    show_total_axis: bool = True,
    show_total_quantile_band: bool = True,
    total_axis_label: str = "Total temperature change",
    legend_entrywidth: int = 260,
    legend_row_height: int = 18,
    legend_y: float = 1.0,
    y2_headroom: float = 0.05,
    show_cumulative_in_legend: bool = False,
    flow_groupby_name: bool = False,
    y_min: Optional[float] = None,
    y_max: Optional[float] = None,
    y2_max: Optional[float] = None,
) -> go.Figure:
    """Plot delta temperature time series by flow or root activity."""
    delta_t = getattr(trails, "delta_temperature", None)
    if delta_t is None:
        raise ValueError("No delta temperature data stored on Trails.")

    if "year" not in delta_t.dims:
        raise ValueError("Delta temperature data must include a 'year' dimension.")

    if hasattr(delta_t.data, "nnz") and int(delta_t.data.nnz) == 0:
        if getattr(trails, "debug", False):
            print("FAIR debug: delta_temperature is all zeros; proceeding with plot.")

    delta_t_all = delta_t
    if "quantile" in delta_t.dims:
        quantiles = [float(q) for q in delta_t.coords["quantile"].values.tolist()]
        if quantile is None:
            quantile = 50.0 if 50.0 in quantiles else quantiles[0]
        if float(quantile) not in quantiles:
            raise ValueError(
                f"Requested quantile {quantile} not in available {quantiles}."
            )
        delta_t = delta_t.sel(quantile=float(quantile), drop=True)

    if by == "flow":
        if "flow" not in delta_t.dims:
            raise ValueError(
                "Delta temperature data must include a 'flow' dimension for by='flow'."
            )
        if "root activity" in delta_t.dims:
            delta_t = delta_t.sum(dim="root activity")
        if "activity" in delta_t.dims:
            delta_t = delta_t.sum(dim="activity")
        results_by_year = _scores_to_results(delta_t, by_flow=True)
    elif by == "root activity":
        if "root activity" not in delta_t.dims:
            raise ValueError(
                "Delta temperature data must include 'root activity' for by='root activity'."
            )
        data = delta_t
        if "flow" in data.dims:
            data = data.sum(dim="flow")
        if "activity" in data.dims:
            data = data.sum(dim="activity")
        results_by_year = _scores_to_results(data, by_flow=False)
    else:
        raise ValueError("by must be 'flow' or 'root activity'.")

    fig = _plot_results_by_year(
        results_by_year=results_by_year,
        trails=trails,
        title=title,
        method_label=method_label,
        cumulative=False,
        stacked=stacked,
        legend_top_n=legend_top_n,
        show_flow_contributions=(by == "flow"),
        width=width,
        height=height,
        year_tick=year_tick,
        year_range=year_range,
        show_year_grid=show_year_grid,
        yaxis_type=yaxis_type,
        log_eps=log_eps,
        reference_year=reference_year,
        show_cumulative_axis=False,
        cumulative_axis_label=total_axis_label,
        legend_entrywidth=legend_entrywidth,
        legend_row_height=legend_row_height,
        legend_y=legend_y,
        y2_headroom=y2_headroom,
        show_cumulative_in_legend=show_cumulative_in_legend,
        flow_groupby_name=flow_groupby_name,
        static_score=None,
        static_score_label="Static score",
        static_score_dash="dash",
        static_score_color="black",
        y_min=y_min,
        y_max=y_max,
        y2_max=y2_max,
    )

    if show_total_axis:

        years = _select_years_from_results(results_by_year, year_range)

        def _totals(res: dict[int, dict[str, Any]], years_seq: list[int]) -> np.ndarray:
            out = np.zeros(len(years_seq), dtype=float)
            for i, y in enumerate(years_seq):
                payload = res.get(int(y), {})
                out[i] = float(payload.get("scores", 0.0))
            return out

        total_vals = _totals(results_by_year, years)
        if yaxis_type == "log":
            total_vals = np.where(total_vals > 0, total_vals, log_eps)
        fig.add_trace(
            go.Scatter(
                x=years,
                y=total_vals,
                name="Total",
                showlegend=True,
                mode="lines",
                line=dict(width=2, color="black"),
                yaxis="y",
                hovertemplate=(
                    "<b>Total</b><br>"
                    "Year: %{x}<br>"
                    f"{total_axis_label}: %{{y:.6g}}<extra></extra>"
                ),
            )
        )

        if "quantile" in delta_t_all.dims and show_total_quantile_band:
            quantiles = [
                float(q) for q in delta_t_all.coords["quantile"].values.tolist()
            ]
            q_low, q_high = 2.5, 97.5
            if q_low in quantiles and q_high in quantiles:
                q_low_data = delta_t_all.sel(quantile=float(q_low), drop=True)
                q_high_data = delta_t_all.sel(quantile=float(q_high), drop=True)

                def _results_for(data: xr.DataArray) -> dict[int, dict[str, Any]]:
                    if by == "flow":
                        if "root activity" in data.dims:
                            data = data.sum(dim="root activity")
                        if "activity" in data.dims:
                            data = data.sum(dim="activity")
                        return _scores_to_results(data, by_flow=True)
                    if "root activity" in data.dims:
                        data = data.sum(dim="flow")
                    if "activity" in data.dims:
                        data = data.sum(dim="activity")
                    return _scores_to_results(data, by_flow=False)

                results_low = _results_for(q_low_data)
                results_high = _results_for(q_high_data)
                total_low = _totals(results_low, years)
                total_high = _totals(results_high, years)
                if yaxis_type == "log":
                    total_low = np.where(total_low > 0, total_low, log_eps)
                    total_high = np.where(total_high > 0, total_high, log_eps)
                fig.add_trace(
                    go.Scatter(
                        x=years,
                        y=total_low,
                        mode="lines",
                        line=dict(width=0),
                        showlegend=False,
                        hoverinfo="skip",
                        yaxis="y",
                    )
                )
                fig.add_trace(
                    go.Scatter(
                        x=years,
                        y=total_high,
                        mode="lines",
                        fill="tonexty",
                        fillcolor="rgba(60,60,60,0.35)",
                        line=dict(width=0),
                        name="2.5-97.5th percentile",
                        yaxis="y",
                    )
                )

    fig.update_layout(title=dict(text=""))
    fig.add_annotation(
        text=title,
        x=0.5,
        y=-0.2,
        xref="paper",
        yref="paper",
        xanchor="center",
        yanchor="top",
        showarrow=False,
    )
    return fig


def plot_delta_temperature(*args, **kwargs) -> go.Figure:
    """Backward-compatible alias for plot_temp."""
    return plot_temp(*args, **kwargs)


def _configure_flow_axes(
    fig: go.Figure,
    n_panels: int,
    ncols: int,
    acts: list[int],
    year_min: int,
    year_max: int,
) -> None:
    """Configure axes for flow subplots.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param n_panels: Number of panels.
    :type n_panels: int
    :param ncols: Number of columns.
    :type ncols: int
    :param acts: Activity indices used for labels.
    :type acts: list[int]
    :param year_min: Minimum year.
    :type year_min: int
    :param year_max: Maximum year.
    :type year_max: int
    """
    n_acts = len(acts)
    all_row_idx = list(range(n_acts))

    for i in range(n_panels):
        r = i // ncols + 1
        c = i % ncols + 1

        fig.update_yaxes(
            row=r,
            col=c,
            tickmode="array",
            tickvals=all_row_idx,
            ticktext=[str(a) for a in acts],
            range=[-0.5, n_acts - 0.5],
            showgrid=True,
            tick0=0,
            dtick=1,
            automargin=True,
            zeroline=False,
        )

        fig.update_xaxes(
            row=r,
            col=c,
            range=[year_min - 0.5, year_max + 0.5],
            tickangle=-90,
            showgrid=True,
            zeroline=False,
        )


def _merge_all_edges(
    edges_by_depth: dict[int, dict], depths_list: list[int]
) -> dict[tuple[tuple[int, int], tuple[int, int]], float]:
    """Merge edge mappings across depths.

    :param edges_by_depth: Mapping of depth to edges.
    :type edges_by_depth: dict[int, dict]
    :param depths_list: Depths to include.
    :type depths_list: list[int]
    :returns: Mapping of edges to total amounts.
    :rtype: dict[tuple, float]
    """
    merged_edges_all: dict[tuple[tuple[int, int], tuple[int, int]], float] = (
        defaultdict(float)
    )
    for d in depths_list:
        for key, amt in edges_by_depth.get(d, {}).items():
            merged_edges_all[key] += float(amt)
    return merged_edges_all


def _add_flow_panel_traces(
    fig: go.Figure,
    edges: dict[tuple[tuple[int, int], tuple[int, int]], float],
    act_to_row: dict[int, int],
    idx_to_label: dict[int, str],
    dot_size: int,
    row: int,
    col: int,
    show_legend: bool,
) -> None:
    """Add consumer and supplier node traces for a panel.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param edges: Edge mapping for the panel.
    :type edges: dict[tuple[tuple[int, int], tuple[int, int]], float]
    :param act_to_row: Mapping of activity index to row position.
    :type act_to_row: dict[int, int]
    :param idx_to_label: Mapping of activity index to label.
    :type idx_to_label: dict[int, str]
    :param dot_size: Marker size for nodes.
    :type dot_size: int
    :param row: Subplot row index.
    :type row: int
    :param col: Subplot column index.
    :type col: int
    :param show_legend: Whether to show legend entries.
    :type show_legend: bool
    """
    consumer_nodes = set()
    supplier_nodes = set()

    for (parent, child), amt in edges.items():
        y_cons, a_cons = parent
        y_sup, a_sup = child

        consumer_nodes.add((int(y_cons), int(a_cons)))
        supplier_nodes.add((int(y_sup), int(a_sup)))

    cons_x, cons_y, cons_text = [], [], []
    sup_x, sup_y, sup_text = [], [], []

    for year, act in sorted(consumer_nodes):
        cons_x.append(year)
        cons_y.append(act_to_row.get(act, -1))
        cons_text.append(idx_to_label.get(act, f"Activity {act}"))

    for year, act in sorted(supplier_nodes):
        sup_x.append(year)
        sup_y.append(act_to_row.get(act, -1))
        sup_text.append(idx_to_label.get(act, f"Activity {act}"))

    fig.add_trace(
        go.Scatter(
            x=cons_x,
            y=cons_y,
            mode="markers",
            marker=dict(size=dot_size, color="red"),
            name="Consumers",
            text=cons_text,
            hovertemplate=(
                "Year: %{x}<br>"
                "Consumer idx: %{y}<br>"
                "Consumer: %{text}<extra></extra>"
            ),
            showlegend=show_legend,
        ),
        row=row,
        col=col,
    )

    fig.add_trace(
        go.Scatter(
            x=sup_x,
            y=sup_y,
            mode="markers",
            marker=dict(size=dot_size, color="green"),
            name="Suppliers",
            text=sup_text,
            hovertemplate=(
                "Year: %{x}<br>"
                "Supplier idx: %{y}<br>"
                "Supplier: %{text}<extra></extra>"
            ),
            showlegend=show_legend,
        ),
        row=row,
        col=col,
    )

    edge_x = []
    edge_y = []
    head_x = []
    head_y = []

    for (parent, child), amt in edges.items():
        y_cons, a_cons = parent
        y_sup, a_sup = child

        x0 = int(y_sup)
        y0 = act_to_row.get(int(a_sup), -1)
        x1 = int(y_cons)
        y1 = act_to_row.get(int(a_cons), -1)

        if y0 < 0 or y1 < 0:
            continue

        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

        head_x.append(x1)
        head_y.append(y1)

    fig.add_trace(
        go.Scattergl(
            x=edge_x,
            y=edge_y,
            mode="lines",
            line=dict(width=1, color="rgba(100,100,100,0.35)"),
            hoverinfo="skip",
            showlegend=False,
        ),
        row=row,
        col=col,
    )

    fig.add_trace(
        go.Scattergl(
            x=head_x,
            y=head_y,
            mode="markers",
            marker=dict(
                symbol="triangle-left",
                size=6,
                color="rgba(100,100,100,0.7)",
            ),
            hoverinfo="skip",
            showlegend=False,
        ),
        row=row,
        col=col,
    )


def _add_activity_legend(
    fig: go.Figure, acts: list[int], idx_to_label: dict[int, str]
) -> None:
    """Add activity legend entries to the figure.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param acts: Activity indices to label.
    :type acts: list[int]
    :param idx_to_label: Mapping of activity index to label.
    :type idx_to_label: dict[int, str]
    """
    mapping_lines = []
    for a in acts:
        label = idx_to_label.get(a, f"Activity {a}")
        mapping_lines.append(f"<b>{a}</b>: {label}")
    legend_text = "<br>".join(mapping_lines)

    fig.add_annotation(
        x=0.0,
        y=-0.18,
        xref="paper",
        yref="paper",
        xanchor="left",
        yanchor="top",
        align="left",
        showarrow=False,
        text=legend_text,
        font=dict(size=10),
    )


def _apply_flow_layout(
    fig: go.Figure, title: str, base_width: int, base_height: int, nrows: int
) -> None:
    """Apply layout settings for flow subplots.

    :param fig: Plotly figure to update.
    :type fig: plotly.graph_objects.Figure
    :param title: Figure title.
    :type title: str
    :param base_width: Base width in pixels.
    :type base_width: int
    :param base_height: Base height in pixels.
    :type base_height: int
    :param nrows: Number of subplot rows.
    :type nrows: int
    """
    fig.update_layout(
        title=title,
        width=base_width * 2,
        height=base_height * nrows + 160,
        template="plotly_white",
        hovermode="closest",
        margin=dict(l=80, r=40, t=70, b=180),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
        ),
    )


def plot_traversal_grid_flow(
    edges_by_depth: dict[int, dict[tuple[tuple[int, int], tuple[int, int]], float]],
    trails: Trails,
    depths: Optional[List[int]] = None,
    include_all_activities: bool = False,
    title: str = "Traversal flow grid (consumers & suppliers)",
    dot_size: int = 10,
    base_width: int = 550,
    base_height: int = 260,
) -> go.Figure:
    """Plot traversal edges as a grid of flow panels.

    Multi-panel plot of traversal flows on an activity×year grid.

    For each depth d in ``depths`` AND for "All depths":
      - rows  = activity indices
      - cols  = years (global min..max across all depths)
      - red   = consumers (nodes with outgoing edges at that depth)
      - green = suppliers (nodes appearing as children)
      - arrows from supplier -> consumer

    Layout:
      - 2 columns, as many rows as needed.
      - Last subplot is "All depths" (cumulative edges).
      - Y-axis shows activity indices; a legend below maps index → label.

    :param edges_by_depth: Mapping of depth to traversal edges.
    :type edges_by_depth: dict[int, dict[tuple[tuple[int, int], tuple[int, int]], float]]
    :param trails: Trails instance used for metadata.
    :type trails: Trails
    :param depths: Optional list of depths to include.
    :type depths: list[int] | None
    :param include_all_activities: Whether to include all activities.
    :type include_all_activities: bool
    :param title: Figure title.
    :type title: str
    :param dot_size: Marker size for scatter points.
    :type dot_size: int
    :param base_width: Base width in pixels.
    :type base_width: int
    :param base_height: Base height in pixels.
    :type base_height: int
    :returns: Plotly figure with flow panels.
    :rtype: plotly.graph_objects.Figure
    """
    depths_list = _select_depths(edges_by_depth, depths)
    panel_labels: List[str] = [f"Depth {d}" for d in depths_list] + ["All depths"]

    acts = _collect_activities(
        edges_by_depth, trails, depths_list, include_all_activities
    )
    idx_to_label = _build_activity_label_map(trails)
    act_to_row = {act: i for i, act in enumerate(acts)}

    year_min, year_max, years_global = _collect_global_years(
        edges_by_depth, depths_list
    )

    ncols = 2
    fig, nrows = _init_flow_subplots(panel_labels, ncols)
    _configure_flow_axes(fig, len(panel_labels), ncols, acts, year_min, year_max)

    merged_edges_all = _merge_all_edges(edges_by_depth, depths_list)

    panels_depths = depths_list + ["all"]
    for i, depth in enumerate(panels_depths):
        row = i // ncols + 1
        col = i % ncols + 1

        edges = merged_edges_all if depth == "all" else edges_by_depth.get(depth, {})
        if not edges:
            continue

        _add_flow_panel_traces(
            fig=fig,
            edges=edges,
            act_to_row=act_to_row,
            idx_to_label=idx_to_label,
            dot_size=dot_size,
            row=row,
            col=col,
            show_legend=(i == 0),
        )

    _add_activity_legend(fig, acts, idx_to_label)
    _apply_flow_layout(fig, title, base_width, base_height, nrows)

    return fig
