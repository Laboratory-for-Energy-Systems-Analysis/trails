from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union
import json
import os

import math
from pathlib import Path
from textwrap import shorten as _shorten_text

from collections import defaultdict
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import sparse
import xarray as xr

from .trails import Trails
from .utils import _format_path_label
from .chunked_inventory import is_chunked_sparse


def _compute_plot_array(value: xr.DataArray) -> np.ndarray | sparse.COO:
    """Compute only an already-reduced plotting array, sequentially."""
    data = value.data
    if is_chunked_sparse(data):
        return data.compute(scheduler="synchronous")
    return data


def _build_activity_label_map(trails: Trails) -> dict[int, str]:
    """build activity label map.

    :param trails: Value for `trails`.
    :type trails: Trails
    :returns: Return value.
    :rtype: dict[int, str]"""
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
    """build flow label map.

    :param trails: Value for `trails`.
    :type trails: Trails
    :returns: Return value.
    :rtype: dict[int, str]"""
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
    """select years from results.

    :param results_by_year: Value for `results_by_year`.
    :type results_by_year: dict[int, dict[str, Any]]
    :param year_range: Value for `year_range`.
    :type year_range: tuple[int, int] | None
    :returns: Return value.
    :rtype: list[int]
    :raises ValueError: If an error occurs."""
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
    """collect root scores.

    :param results_by_year: Value for `results_by_year`.
    :type results_by_year: dict[int, dict[str, Any]]
    :param years: Value for `years`.
    :type years: list[int]
    :param score_key: Value for `score_key`.
    :type score_key: str
    :returns: Return value.
    :rtype: list[int]
    :raises ValueError: If an error occurs."""
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
    """build score matrix.

    :param results_by_year: Value for `results_by_year`.
    :type results_by_year: dict[int, dict[str, Any]]
    :param years: Value for `years`.
    :type years: list[int]
    :param all_roots: Value for `all_roots`.
    :type all_roots: list[int]
    :param score_key: Value for `score_key`.
    :type score_key: str
    :returns: Return value.
    :rtype: np.ndarray"""
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
    """add root traces.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param years: Value for `years`.
    :type years: list[int]
    :param Y: Value for `Y`.
    :type Y: np.ndarray
    :param all_roots: Value for `all_roots`.
    :type all_roots: list[Any]
    :param idx_to_label: Value for `idx_to_label`.
    :type idx_to_label: dict[Any, str]
    :param method_label: Value for `method_label`.
    :type method_label: str
    :param stacked: Value for `stacked`.
    :type stacked: bool
    :param showlegend_roots: Value for `showlegend_roots`.
    :type showlegend_roots: Optional[set[int]]
    :param showhover_roots: Value for `showhover_roots`.
    :type showhover_roots: Optional[set[int]]"""
    alpha = 0.4 if not stacked else 1.0

    def label_for_root(idx: Any) -> str:
        """Label for root.

        :param idx: Value for `idx`.
        :type idx: Any
        :returns: Return value.
        :rtype: str"""
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
        legend_group = f"root::{root}"

        if stacked:
            y_vals = np.asarray(Y[:, ri], dtype=float)
            y_pos = np.where(y_vals > 0.0, y_vals, 0.0)
            y_neg = np.where(y_vals < 0.0, y_vals, 0.0)
            has_pos = bool(np.any(y_pos != 0.0))
            has_neg = bool(np.any(y_neg != 0.0))

            if has_pos:
                fig.add_trace(
                    go.Scatter(
                        x=years,
                        y=y_pos,
                        name=root_label,
                        meta=root_label_wrapped,
                        legendgroup=legend_group,
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
                        stackgroup="positive",
                    )
                )

            if has_neg:
                fig.add_trace(
                    go.Scatter(
                        x=years,
                        y=y_neg,
                        name=root_label,
                        meta=root_label_wrapped,
                        legendgroup=legend_group,
                        showlegend=(showlegend and not has_pos),
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
                        stackgroup="negative",
                    )
                )
            continue

        fig.add_trace(
            go.Scatter(
                x=years,
                y=Y[:, ri],
                name=root_label,  # keep legend name as before
                meta=root_label_wrapped,  # wrapped version for hover
                legendgroup=legend_group,
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
                fill="tozeroy",
                line=dict(width=2),
                opacity=alpha,
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
    band_label_offset_px: float | None = None,
    auto_depth_scale: bool = True,
    edge_weight: Literal["amount", "score"] = "amount",
    node_scores: dict[tuple, float] | None = None,
) -> str:
    """Plot temporal graph.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param min_edge_amount: Value for `min_edge_amount`.
    :type min_edge_amount: float
    :param notebook: Value for `notebook`.
    :type notebook: bool
    :param filename: Value for `filename`.
    :type filename: str
    :param height: Value for `height`.
    :type height: str
    :param width: Value for `width`.
    :type width: str
    :param physics: Value for `physics`.
    :type physics: bool
    :param layout_by_year_depth: Value for `layout_by_year_depth`.
    :type layout_by_year_depth: bool
    :param year_scale: Value for `year_scale`.
    :type year_scale: float
    :param depth_scale: Value for `depth_scale`.
    :type depth_scale: float
    :param max_label_chars: Value for `max_label_chars`.
    :type max_label_chars: int
    :param level0_edge_color: Value for `level0_edge_color`.
    :type level0_edge_color: str
    :param palette: Value for `palette`.
    :type palette: Optional[list[str]]
    :param show_year_labels: Value for `show_year_labels`.
    :type show_year_labels: bool
    :param year_label_offset: Value for `year_label_offset`.
    :type year_label_offset: float
    :param show_band_labels: Value for `show_band_labels`.
    :type show_band_labels: bool
    :param band_label_offset_px: Value for `band_label_offset_px`.
    :type band_label_offset_px: float | None
    :param auto_depth_scale: Value for `auto_depth_scale`.
    :type auto_depth_scale: bool
    :param edge_weight: Value for `edge_weight`.
    :type edge_weight: Literal['amount', 'score']
    :param node_scores: Value for `node_scores`.
    :type node_scores: dict[tuple, float] | None
    :returns: Return value.
    :rtype: str
    :raises RuntimeError: If an error occurs.
    :raises ValueError: If an error occurs."""
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
        """node score.

        :param node: Value for `node`.
        :type node: object
        :returns: Return value.
        :rtype: float"""
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
        """truncate.

        :param text: Value for `text`.
        :type text: str
        :param limit: Value for `limit`.
        :type limit: int
        :returns: Return value.
        :rtype: str"""
        if len(text) <= limit:
            return text
        return text[:limit]

    def _label_node(node: tuple, data: dict) -> str:
        """label node.

        :param node: Value for `node`.
        :type node: tuple
        :param data: Value for `data`.
        :type data: dict
        :returns: Return value.
        :rtype: str"""
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
        """branch key.

        :param node: Value for `node`.
        :type node: object
        :returns: Return value.
        :rtype: tuple[str, str, str]"""
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
        """edge color.

        :param u: Value for `u`.
        :type u: object
        :param v: Value for `v`.
        :type v: object
        :returns: Return value.
        :rtype: str"""
        src_depth = int(H.nodes[u].get("depth", 0))
        if src_depth == 0:
            return edge_colors.get((u, v), level0_edge_color)
        return node_branch_color.get(u, palette[0])

    def _node_color(n: object) -> Optional[str]:
        """node color.

        :param n: Value for `n`.
        :type n: object
        :returns: Return value.
        :rtype: Optional[str]"""
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
            """edge title.

            :param src_label: Value for `src_label`.
            :type src_label: str
            :param dst_label: Value for `dst_label`.
            :type dst_label: str
            :param amount: Value for `amount`.
            :type amount: float
            :param score: Value for `score`.
            :type score: float | None
            :returns: Return value.
            :rtype: str"""
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
            """edge title.

            :param src_label: Value for `src_label`.
            :type src_label: str
            :param dst_label: Value for `dst_label`.
            :type dst_label: str
            :param amount: Value for `amount`.
            :type amount: float
            :param score: Value for `score`.
            :type score: float | None
            :returns: Return value.
            :rtype: str"""
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
    var yearMidFixed = {float(year_mid)};
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
      // Keep year origin fixed so node x-positions remain aligned with the year overlay.
      var depthMid = 0.0;
      window.trailDepthScale = depthScale;
      window.trailDepthMid = depthMid;
      // assign positions
      nodesFiltered.forEach(function(n) {{
        var bandKey = String(n.depth) + '|' + String(n.band_key || '');
        var bandOffset = bandOffsets[bandKey] || 0.0;
        n.x = (n.year - yearMidFixed) * yearScale;
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
    """add cumulative trace.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param years: Value for `years`.
    :type years: list[int]
    :param total_raw: Value for `total_raw`.
    :type total_raw: list[float]
    :param cumulative_axis_label: Value for `cumulative_axis_label`.
    :type cumulative_axis_label: str
    :param yaxis_type: Value for `yaxis_type`.
    :type yaxis_type: str
    :param log_eps: Value for `log_eps`.
    :type log_eps: float
    :returns: Return value.
    :rtype: np.ndarray"""
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
    """add static score trace.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param years: Value for `years`.
    :type years: list[int]
    :param static_score: Value for `static_score`.
    :type static_score: float
    :param static_score_label: Value for `static_score_label`.
    :type static_score_label: str
    :param static_score_dash: Value for `static_score_dash`.
    :type static_score_dash: str
    :param static_score_color: Value for `static_score_color`.
    :type static_score_color: str
    :param method_label: Value for `method_label`.
    :type method_label: str"""
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
    """compute layout dimensions.

    :param width: Value for `width`.
    :type width: int | None
    :param legend_entrywidth: Value for `legend_entrywidth`.
    :type legend_entrywidth: int
    :param legend_row_height: Value for `legend_row_height`.
    :type legend_row_height: int
    :param n_items: Value for `n_items`.
    :type n_items: int
    :returns: Return value.
    :rtype: tuple[int, int]"""
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
    """apply base layout.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param width: Value for `width`.
    :type width: int
    :param height: Value for `height`.
    :type height: int
    :param title: Value for `title`.
    :type title: str
    :param legend_y: Value for `legend_y`.
    :type legend_y: float
    :param entry_w: Value for `entry_w`.
    :type entry_w: int
    :param top_margin: Value for `top_margin`.
    :type top_margin: int
    :param method_label: Value for `method_label`.
    :type method_label: str
    :param yaxis_type: Value for `yaxis_type`.
    :type yaxis_type: str
    :param show_cumulative_axis: Value for `show_cumulative_axis`.
    :type show_cumulative_axis: bool
    :param static_score: Value for `static_score`.
    :type static_score: float | None
    :param cumulative_axis_label: Value for `cumulative_axis_label`.
    :type cumulative_axis_label: str"""
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
            y=0.98,
            yanchor="top",
            yref="container",
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
            tracegroupgap=0,
            groupclick="togglegroup",
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
    """apply linear yaxis alignment.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param Y: Value for `Y`.
    :type Y: np.ndarray
    :param cum_vals: Value for `cum_vals`.
    :type cum_vals: np.ndarray | None
    :param static_score: Value for `static_score`.
    :type static_score: float | None
    :param y_min: Value for `y_min`.
    :type y_min: float | None
    :param y_max: Value for `y_max`.
    :type y_max: float | None
    :param y2_max: Value for `y2_max`.
    :type y2_max: float | None
    :param y2_headroom: Value for `y2_headroom`.
    :type y2_headroom: float
    :param stacked: Value for `stacked`.
    :type stacked: bool"""
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
    """apply xaxis settings.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param year_tick: Value for `year_tick`.
    :type year_tick: int | None
    :param year_range: Value for `year_range`.
    :type year_range: tuple[int, int] | None
    :param years: Value for `years`.
    :type years: list[int]
    :param show_year_grid: Value for `show_year_grid`.
    :type show_year_grid: bool"""
    fig.update_xaxes(
        dtick=year_tick,
        tickmode="linear",
        showgrid=show_year_grid,
        tick0=(year_range[0] if year_range else years[0]),
        range=list(year_range) if year_range else None,
    )


def _add_reference_year_line(fig: go.Figure, reference_year: int | None) -> None:
    """add reference year line.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param reference_year: Value for `reference_year`.
    :type reference_year: int | None"""
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
    """To impact year results.

    :param results: Value for `results`.
    :type results: Dict[int, Dict[str, Any]] | Dict[str, Any]
    :returns: Return value.
    :rtype: Dict[int, Dict[str, Any]]
    :raises ValueError: If an error occurs."""

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
    """characterized inventory to results.

    :param characterized_inventory: Value for `characterized_inventory`.
    :type characterized_inventory: xr.DataArray
    :param by_flow: Value for `by_flow`.
    :type by_flow: bool
    :returns: Return value.
    :rtype: Dict[int, Dict[str, Any]]
    :raises ValueError: If an error occurs."""
    if "method" in characterized_inventory.dims:
        methods = characterized_inventory.coords["method"].values
        if len(methods) != 1:
            raise ValueError(
                "characterized_inventory has multiple methods; select one before plotting."
            )
        characterized_inventory = characterized_inventory.isel(method=0, drop=True)
    if "flow" not in characterized_inventory.dims:
        raise ValueError("characterized_inventory must include a 'flow' dimension.")
    if "activity" not in characterized_inventory.dims:
        raise ValueError(
            "characterized_inventory must include an 'activity' dimension."
        )
    if "year" not in characterized_inventory.dims:
        raise ValueError("characterized_inventory must include a 'year' dimension.")

    if by_flow:
        summed = _compute_plot_array(
            characterized_inventory.sum(dim="activity")
        )
        score_key = "scores_by_flow"
    else:
        summed = _compute_plot_array(characterized_inventory.sum(dim="flow"))
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
    """characterized inventory to root results.

    :param characterized_inventory: Value for `characterized_inventory`.
    :type characterized_inventory: xr.DataArray
    :returns: Return value.
    :rtype: Dict[int, Dict[str, Any]]
    :raises ValueError: If an error occurs."""
    if "method" in characterized_inventory.dims:
        methods = characterized_inventory.coords["method"].values
        if len(methods) != 1:
            raise ValueError(
                "characterized_inventory has multiple methods; select one before plotting."
            )
        characterized_inventory = characterized_inventory.isel(method=0, drop=True)
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

    summed = _compute_plot_array(
        characterized_inventory.sum(dim=["activity", "flow"])
    )
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
    """wrap hover label.

    :param text: Value for `text`.
    :type text: str
    :param max_chars: Value for `max_chars`.
    :type max_chars: int
    :returns: Return value.
    :rtype: str"""
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
    """scores to results.

    :param scores: Value for `scores`.
    :type scores: xr.DataArray
    :param by_flow: Value for `by_flow`.
    :type by_flow: bool
    :returns: Return value.
    :rtype: Dict[int, Dict[str, Any]]
    :raises ValueError: If an error occurs."""
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
    """aggregate flow results by name.

    :param results_by_year: Value for `results_by_year`.
    :type results_by_year: Dict[int, Dict[str, Any]]
    :param trails: Value for `trails`.
    :type trails: Trails
    :returns: Return value.
    :rtype: Dict[int, Dict[str, Any]]"""
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
    """plot results by year.

    :param results_by_year: Value for `results_by_year`.
    :type results_by_year: Union[Dict[int, Dict[str, Any]], Dict[str, Any], xr.DataArray]
    :param trails: Value for `trails`.
    :type trails: Trails
    :param title: Value for `title`.
    :type title: str
    :param method_label: Value for `method_label`.
    :type method_label: str
    :param cumulative: Value for `cumulative`.
    :type cumulative: bool
    :param stacked: Value for `stacked`.
    :type stacked: bool
    :param legend_top_n: Value for `legend_top_n`.
    :type legend_top_n: int
    :param show_flow_contributions: Value for `show_flow_contributions`.
    :type show_flow_contributions: bool
    :param width: Value for `width`.
    :type width: Optional[int]
    :param height: Value for `height`.
    :type height: Optional[int]
    :param year_tick: Value for `year_tick`.
    :type year_tick: int
    :param year_range: Value for `year_range`.
    :type year_range: Optional[Tuple[int, int]]
    :param show_year_grid: Value for `show_year_grid`.
    :type show_year_grid: bool
    :param yaxis_type: Value for `yaxis_type`.
    :type yaxis_type: Literal['linear', 'log']
    :param log_eps: Value for `log_eps`.
    :type log_eps: float
    :param reference_year: Value for `reference_year`.
    :type reference_year: Optional[int]
    :param show_cumulative_axis: Value for `show_cumulative_axis`.
    :type show_cumulative_axis: bool
    :param cumulative_axis_label: Value for `cumulative_axis_label`.
    :type cumulative_axis_label: str
    :param legend_entrywidth: Value for `legend_entrywidth`.
    :type legend_entrywidth: int
    :param legend_row_height: Value for `legend_row_height`.
    :type legend_row_height: int
    :param legend_y: Value for `legend_y`.
    :type legend_y: float
    :param y2_headroom: Value for `y2_headroom`.
    :type y2_headroom: float
    :param show_cumulative_in_legend: Value for `show_cumulative_in_legend`.
    :type show_cumulative_in_legend: bool
    :param static_score: Value for `static_score`.
    :type static_score: Optional[float]
    :param static_score_label: Value for `static_score_label`.
    :type static_score_label: str
    :param static_score_dash: Value for `static_score_dash`.
    :type static_score_dash: str
    :param static_score_color: Value for `static_score_color`.
    :type static_score_color: str
    :param y_min: Value for `y_min`.
    :type y_min: Optional[float]
    :param y_max: Value for `y_max`.
    :type y_max: Optional[float]
    :param y2_max: Value for `y2_max`.
    :type y2_max: Optional[float]
    :param flow_groupby_name: Value for `flow_groupby_name`.
    :type flow_groupby_name: bool
    :returns: Return value.
    :rtype: go.Figure
    :raises ValueError: If an error occurs."""

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
    if year_range is None and Y_raw.size:
        max_abs = float(np.max(np.abs(Y_raw)))
        tol = np.finfo(float).eps * max(1.0, max_abs) * 10.0
        active_rows = np.any(np.abs(Y_raw) > tol, axis=1)
        if np.any(active_rows):
            first = int(np.argmax(active_rows))
            last = int(len(active_rows) - 1 - np.argmax(active_rows[::-1]))
            years = years[first : last + 1]
            Y_raw = Y_raw[first : last + 1, :]

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
    static_score: Optional[float] | dict[str, float] | list[float] = None,
    static_score_label: str = "Static score",
    static_score_dash: str = "dash",
    static_score_color: str = "black",
    y_min: Optional[float] = None,
    y_max: Optional[float] = None,
    y2_max: Optional[float] = None,
) -> go.Figure | list[go.Figure]:
    """Plot temporal scores.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param title: Value for `title`.
    :type title: str
    :param method_label: Value for `method_label`.
    :type method_label: str
    :param method: Value for `method`.
    :type method: Optional[str]
    :param cumulative: Value for `cumulative`.
    :type cumulative: bool
    :param stacked: Value for `stacked`.
    :type stacked: bool
    :param legend_top_n: Value for `legend_top_n`.
    :type legend_top_n: int
    :param show_flow_contributions: Value for `show_flow_contributions`.
    :type show_flow_contributions: bool
    :param width: Value for `width`.
    :type width: Optional[int]
    :param height: Value for `height`.
    :type height: Optional[int]
    :param year_tick: Value for `year_tick`.
    :type year_tick: int
    :param year_range: Value for `year_range`.
    :type year_range: Optional[Tuple[int, int]]
    :param show_year_grid: Value for `show_year_grid`.
    :type show_year_grid: bool
    :param yaxis_type: Value for `yaxis_type`.
    :type yaxis_type: Literal['linear', 'log']
    :param log_eps: Value for `log_eps`.
    :type log_eps: float
    :param reference_year: Value for `reference_year`.
    :type reference_year: Optional[int]
    :param show_cumulative_axis: Value for `show_cumulative_axis`.
    :type show_cumulative_axis: bool
    :param cumulative_axis_label: Value for `cumulative_axis_label`.
    :type cumulative_axis_label: str
    :param legend_entrywidth: Value for `legend_entrywidth`.
    :type legend_entrywidth: int
    :param legend_row_height: Value for `legend_row_height`.
    :type legend_row_height: int
    :param legend_y: Value for `legend_y`.
    :type legend_y: float
    :param y2_headroom: Value for `y2_headroom`.
    :type y2_headroom: float
    :param show_cumulative_in_legend: Value for `show_cumulative_in_legend`.
    :type show_cumulative_in_legend: bool
    :param flow_groupby_name: Value for `flow_groupby_name`.
    :type flow_groupby_name: bool
    :param static_score: Value for `static_score`.
    :type static_score: Optional[float] | dict[str, float] | list[float]
    :param static_score_label: Value for `static_score_label`.
    :type static_score_label: str
    :param static_score_dash: Value for `static_score_dash`.
    :type static_score_dash: str
    :param static_score_color: Value for `static_score_color`.
    :type static_score_color: str
    :param y_min: Value for `y_min`.
    :type y_min: Optional[float]
    :param y_max: Value for `y_max`.
    :type y_max: Optional[float]
    :param y2_max: Value for `y2_max`.
    :type y2_max: Optional[float]
    :returns: Return value.
    :rtype: go.Figure | list[go.Figure]
    :raises ValueError: If an error occurs."""
    if (
        not show_flow_contributions
        and getattr(trails, "scores", None) is not None
    ):
        # Scores are accumulated incrementally during LCA and are already
        # reduced over elementary flows. Prefer this compact representation for
        # ordinary plots; characterized_inventory remains available for flow
        # contribution analysis and downstream coupling.
        results_by_year = trails.scores
    elif trails.characterized_inventory is not None:
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
                if isinstance(static_score, list):
                    if idx >= len(static_score):
                        raise ValueError(
                            "Static score list length does not match methods."
                        )
                    score_for_method = float(static_score[idx])
                elif isinstance(static_score, dict):
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
        if isinstance(static_score, list):
            if method is None:
                raise ValueError(
                    "Static score list provided but no method was selected."
                )
            methods_list = results_by_year.coords["method"].values.tolist()
            if method not in methods_list:
                raise ValueError(
                    f"Method '{method}' not found in characterized inventory."
                )
            idx = methods_list.index(method)
            if idx >= len(static_score):
                raise ValueError("Static score list length does not match methods.")
            static_score = float(static_score[idx])
        elif isinstance(static_score, dict):
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
    elif isinstance(static_score, list):
        if len(static_score) == 1:
            static_score = float(static_score[0])
        else:
            raise ValueError(
                "Multiple static scores provided but plotting data has no "
                "'method' dimension. Run trails.lci() and trails.lcia(...) "
                "to retain per-method characterized inventory, or pass "
                "method=... with a single static score."
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
    """Plot top paths for year.

    :param provenance: Value for `provenance`.
    :type provenance: Dict[tuple[int, int], Dict[tuple[int, ...], float]]
    :param trails: Value for `trails`.
    :type trails: Trails
    :param year: Value for `year`.
    :type year: int
    :param top_n: Value for `top_n`.
    :type top_n: int
    :param title: Value for `title`.
    :type title: str
    :param amount_label: Value for `amount_label`.
    :type amount_label: str
    :returns: Return value.
    :rtype: go.Figure
    :raises ValueError: If an error occurs."""
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


_ADAPTIVE_SANKEY_PALETTE = [
    "#4c78a8",
    "#f58518",
    "#54a24b",
    "#e45756",
    "#72b7b2",
    "#b279a2",
    "#ff9da6",
    "#9d755d",
    "#bab0ac",
    "#2f4b7c",
    "#a05195",
    "#d45087",
]


def _adaptive_sankey_activity_metadata(
    trails_obj: Trails,
    activity_index: int,
) -> dict[str, Any]:
    """Return activity metadata from the loaded Trails index mappings."""
    for mapping in getattr(trails_obj, "activity_indices", {}).values():
        meta = mapping.get(int(activity_index))
        if isinstance(meta, dict):
            return dict(meta)
    return {}


def _adaptive_sankey_full_label(
    trails_obj: Trails,
    activity_index: int,
) -> str:
    meta = _adaptive_sankey_activity_metadata(trails_obj, int(activity_index))
    label = " | ".join(
        part
        for part in (
            str(meta.get("name") or ""),
            str(meta.get("reference product") or ""),
            str(meta.get("location") or ""),
        )
        if part
    )
    return label or f"Activity {int(activity_index)}"


def _adaptive_sankey_reference_product(
    trails_obj: Trails,
    activity_index: int,
) -> str:
    meta = _adaptive_sankey_activity_metadata(trails_obj, int(activity_index))
    reference_product = str(meta.get("reference product") or "").strip()
    if reference_product:
        return reference_product
    name = str(meta.get("name") or "").strip()
    return name or f"Activity {int(activity_index)}"


def _adaptive_sankey_rgba_color(
    branch: str,
    color_map: dict[str, str],
    alpha: float,
) -> str:
    if branch not in color_map:
        idx = len(color_map) % len(_ADAPTIVE_SANKEY_PALETTE)
        color_map[branch] = _ADAPTIVE_SANKEY_PALETTE[idx]
    color = color_map[branch]
    red = int(color[1:3], 16)
    green = int(color[3:5], 16)
    blue = int(color[5:7], 16)
    return f"rgba({red},{green},{blue},{float(alpha):.3f})"


def _adaptive_sankey_edge_color(
    branch: str,
    color_map: dict[str, str],
) -> str:
    return _adaptive_sankey_rgba_color(branch, color_map, 0.42)


def _adaptive_sankey_root_routed_amounts(
    graph: Any,
    *,
    root_indices: set[int],
) -> tuple[
    dict[object, dict[int, float]],
    dict[tuple[object, object, int], float],
]:
    """Allocate routed edge amounts to first-level root branches."""
    node_root_amounts: dict[object, dict[int, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    edge_root_amounts: dict[tuple[object, object, int], float] = defaultdict(float)

    edges = sorted(
        graph.edges(data=True),
        key=lambda edge: (
            int(graph.nodes[edge[0]].get("depth", 0)),
            int(graph.nodes[edge[1]].get("depth", 0)),
            int(graph.nodes[edge[1]].get("year", 0)),
            repr(edge[0]),
            repr(edge[1]),
        ),
    )
    for source, target, data in edges:
        source_depth = int(graph.nodes[source].get("depth", 0))
        target_depth = int(graph.nodes[target].get("depth", 0))
        if target_depth <= source_depth:
            continue
        amount_abs = abs(float(data.get("amount", 0.0)))
        if amount_abs == 0.0:
            continue

        if source_depth == 0:
            root_idx = int(graph.nodes[target].get("act_idx", -1))
            if root_idx not in root_indices:
                continue
            edge_root_amounts[(source, target, root_idx)] += amount_abs
            node_root_amounts[target][root_idx] += amount_abs
            continue

        source_roots = node_root_amounts.get(source, {})
        source_total = sum(abs(float(value)) for value in source_roots.values())
        if source_total <= 0.0:
            continue
        for root_idx, source_root_amount in source_roots.items():
            share = abs(float(source_root_amount)) / source_total
            root_amount = amount_abs * share
            if root_amount == 0.0:
                continue
            edge_root_amounts[(source, target, int(root_idx))] += root_amount
            node_root_amounts[target][int(root_idx)] += root_amount

    return node_root_amounts, edge_root_amounts


def _adaptive_sankey_edge_rows(
    trails_obj: Trails,
    *,
    method: str | None,
) -> list[dict[str, Any]]:
    graph = getattr(trails_obj, "graph", None)
    if graph is None:
        raise RuntimeError(
            "Trails graph is missing; run trails.temporal_routing(...) first."
        )

    root_indices = {
        int(data.get("act_idx", -1))
        for _node, data in graph.nodes(data=True)
        if int(data.get("depth", -1)) == 1
    }
    if not root_indices:
        return []
    node_root_amounts, edge_root_amounts = _adaptive_sankey_root_routed_amounts(
        graph,
        root_indices=root_indices,
    )

    score_denoms: dict[object, float] = {}
    for node, roots in node_root_amounts.items():
        score_denoms[node] = sum(abs(float(value)) for value in roots.values())

    rows: list[dict[str, Any]] = []
    for u, v, data in graph.edges(data=True):
        raw_amount = float(data.get("amount", 0.0))
        amount_abs = abs(raw_amount)
        if amount_abs == 0.0:
            continue
        child = graph.nodes[v]
        source = graph.nodes[u]
        node_score = abs(float(child.get("score_potential") or 0.0))
        if node_score == 0.0:
            continue
        denom = float(score_denoms.get(v, 0.0))
        if denom <= 0.0:
            continue

        src_idx = int(source.get("act_idx", -1))
        dst_idx = int(child.get("act_idx", -1))
        src_meta = _adaptive_sankey_activity_metadata(trails_obj, src_idx)
        dst_meta = _adaptive_sankey_activity_metadata(trails_obj, dst_idx)
        for root_idx in root_indices:
            edge_amount_abs = float(edge_root_amounts.get((u, v, root_idx), 0.0))
            if edge_amount_abs == 0.0:
                continue
            edge_score = node_score * edge_amount_abs / denom
            if edge_score == 0.0:
                continue
            rows.append(
                {
                    "method": "" if method is None else str(method),
                    "source_node": repr(u),
                    "target_node": repr(v),
                    "source_year": int(source.get("year", -1)),
                    "target_year": int(child.get("year", -1)),
                    "source_depth": int(source.get("depth", -1)),
                    "target_depth": int(child.get("depth", -1)),
                    "source_activity_index": src_idx,
                    "target_activity_index": dst_idx,
                    "source_name": str(src_meta.get("name") or ""),
                    "source_reference_product": str(
                        src_meta.get("reference product") or ""
                    ),
                    "source_location": str(src_meta.get("location") or ""),
                    "target_name": str(dst_meta.get("name") or ""),
                    "target_reference_product": str(
                        dst_meta.get("reference product") or ""
                    ),
                    "target_location": str(dst_meta.get("location") or ""),
                    "target_unit": str(dst_meta.get("unit") or ""),
                    "root_activity_index": int(root_idx),
                    "branch": _adaptive_sankey_full_label(
                        trails_obj,
                        int(root_idx),
                    ),
                    "raw_amount": raw_amount,
                    "raw_amount_abs": amount_abs,
                    "root_routed_amount_abs": edge_amount_abs,
                    "child_node_score": node_score,
                    "child_score_allocation_amount_abs": denom,
                    "edge_score": float(edge_score),
                    "edge_score_abs": abs(float(edge_score)),
                }
            )

    rows.sort(key=lambda row: float(row["edge_score_abs"]), reverse=True)
    total_abs = sum(float(row["edge_score_abs"]) for row in rows)
    running = 0.0
    for rank, row in enumerate(rows, start=1):
        running += float(row["edge_score_abs"])
        row["rank_abs_edge_score"] = rank
        row["cumulative_abs_edge_score"] = running
        row["cumulative_abs_edge_score_share"] = (
            running / total_abs if total_abs else 0.0
        )
    return rows


def _adaptive_sankey_branch_rows(
    edge_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    by_branch: dict[str, dict[str, Any]] = {}
    for row in edge_rows:
        branch = str(row.get("branch") or "unassigned")
        item = by_branch.setdefault(
            branch,
            {
                "branch": branch,
                "edge_count": 0,
                "score_sum": 0.0,
                "score_abs_sum": 0.0,
            },
        )
        item["edge_count"] += 1
        item["score_sum"] += float(row.get("edge_score") or 0.0)
        item["score_abs_sum"] += float(row.get("edge_score_abs") or 0.0)

    rows = list(by_branch.values())
    rows.sort(key=lambda row: abs(float(row["score_abs_sum"])), reverse=True)
    total_abs = sum(float(row["score_abs_sum"]) for row in rows)
    for rank, row in enumerate(rows, start=1):
        row["rank_abs_score"] = rank
        row["score_abs_share"] = (
            float(row["score_abs_sum"]) / total_abs if total_abs else 0.0
        )
    return rows


def _adaptive_sankey_node_potential_rows(
    trails_obj: Trails,
    *,
    method: str | None,
) -> list[dict[str, Any]]:
    """Allocate routed node score potentials to first-level branches."""
    graph = getattr(trails_obj, "graph", None)
    if graph is None:
        raise RuntimeError(
            "Trails graph is missing; run trails.temporal_routing(...) first."
        )
    root_indices = {
        int(data.get("act_idx", -1))
        for _node, data in graph.nodes(data=True)
        if int(data.get("depth", -1)) == 1
    }
    node_root_amounts, _edge_root_amounts = _adaptive_sankey_root_routed_amounts(
        graph,
        root_indices=root_indices,
    )

    rows: list[dict[str, Any]] = []
    for node, data in graph.nodes(data=True):
        depth = int(data.get("depth", -1))
        if depth < 1:
            continue
        node_score = abs(float(data.get("score_potential") or 0.0))
        if node_score == 0.0:
            continue
        roots = {
            int(root): abs(float(amount))
            for root, amount in node_root_amounts.get(node, {}).items()
            if float(amount) != 0.0
        }
        denom = sum(roots.values())
        if denom <= 0.0:
            continue
        year = int(data.get("year", -1))
        act_idx = int(data.get("act_idx", -1))
        activity = _adaptive_sankey_full_label(trails_obj, act_idx)
        for root_idx, amount_abs in roots.items():
            fraction = float(amount_abs) / float(denom)
            allocated_score = node_score * fraction
            if allocated_score == 0.0:
                continue
            rows.append(
                {
                    "method": "" if method is None else str(method),
                    "node": repr(node),
                    "year": year,
                    "depth": depth,
                    "activity_index": act_idx,
                    "activity": activity,
                    "root_activity_index": int(root_idx),
                    "branch": _adaptive_sankey_full_label(
                        trails_obj,
                        int(root_idx),
                    ),
                    "node_score": float(allocated_score),
                    "node_score_abs": abs(float(allocated_score)),
                    "node_score_allocation_fraction": fraction,
                    "node_score_allocation_amount_abs": float(denom),
                    "node_score_source": (
                        "routed_node_score_potential_allocated_by_branch_amount"
                    ),
                }
            )

    rows.sort(key=lambda row: float(row["node_score_abs"]), reverse=True)
    total_abs = sum(float(row["node_score_abs"]) for row in rows)
    running = 0.0
    for rank, row in enumerate(rows, start=1):
        running += float(row["node_score_abs"])
        row["rank_abs_node_score"] = rank
        row["cumulative_abs_node_score"] = running
        row["cumulative_abs_node_score_share"] = (
            running / total_abs if total_abs else 0.0
        )
    return rows


def _adaptive_sankey_allowed_branches(
    branch_rows: list[dict[str, Any]],
    *,
    cutoff: float,
) -> set[str]:
    allowed = {
        str(row.get("branch") or "")
        for row in branch_rows
        if float(row.get("score_abs_share") or 0.0) >= float(cutoff)
    }
    if allowed:
        return allowed
    return {str(row.get("branch") or "") for row in branch_rows}


def _adaptive_sankey_year_depth_positions(
    *,
    years: list[int],
    depths: list[int],
    labels: list[str],
    x_min: float,
    x_max: float,
    max_depth: int | None = None,
    min_year: int | None = None,
    max_year: int | None = None,
) -> tuple[list[float], list[float]]:
    """Return fixed Sankey node positions with depth on x and year on y."""
    if not years:
        return [], []

    min_year = min(int(min_year) if min_year is not None else min(years), min(years))
    max_year = max(int(max_year) if max_year is not None else max(years), max(years))
    depth_axis_max = max(int(max_depth or max(depths)), 1)
    year_span = max(max_year - min_year, 1)
    x_span = max(float(x_max) - float(x_min), 0.01)
    x_values = [
        float(x_min) + (float(depth) / float(depth_axis_max)) * x_span
        for depth in depths
    ]
    y_centers = [(float(year) - float(min_year)) / float(year_span) for year in years]

    unique_centers = sorted(set(y_centers))
    if len(unique_centers) > 1:
        min_spacing = min(
            b - a for a, b in zip(unique_centers[:-1], unique_centers[1:])
        )
        half_band = min(0.045, max(0.0, min_spacing * 0.4))
    else:
        half_band = 0.35

    y_values = list(y_centers)
    indices_by_year: dict[int, list[int]] = defaultdict(list)
    for idx, year in enumerate(years):
        indices_by_year[int(year)].append(idx)

    for year, indices in indices_by_year.items():
        if len(indices) == 1:
            continue
        center = (float(year) - float(min_year)) / float(year_span)
        ordered = sorted(indices, key=lambda idx: (depths[idx], labels[idx]))
        if half_band <= 0.0:
            continue
        step = (2.0 * half_band) / float(max(len(ordered) - 1, 1))
        start = center - half_band
        for offset, idx in enumerate(ordered):
            y_values[idx] = start + float(offset) * step

    x_values = [min(float(x_max), max(float(x_min), value)) for value in x_values]
    y_values = [min(0.98, max(0.02, value)) for value in y_values]
    return x_values, y_values


def _adaptive_sankey_year_label_step(min_year: int, max_year: int) -> int:
    span = max(int(max_year) - int(min_year), 1)
    if span > 80:
        return 10
    if span > 35:
        return 5
    if span > 15:
        return 2
    return 1


def _adaptive_sankey_grid_elements(
    *,
    max_depth: int,
    min_year: int,
    max_year: int,
    x_min: float,
    x_max: float,
    horizontal_x_max: float | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    max_depth = max(int(max_depth), 1)
    min_year = int(min_year)
    max_year = int(max_year)
    year_span = max(max_year - min_year, 1)
    year_step = _adaptive_sankey_year_label_step(min_year, max_year)
    horizontal_x_max = float(
        horizontal_x_max if horizontal_x_max is not None else x_max
    )
    shapes: list[dict[str, Any]] = []
    annotations: list[dict[str, Any]] = []

    def depth_x(depth: int) -> float:
        x_span = max(float(x_max) - float(x_min), 0.01)
        value = float(x_min) + (float(depth) / float(max_depth)) * x_span
        return min(float(x_max), max(float(x_min), value))

    def year_y(year: int) -> float:
        sankey_y = (float(year) - float(min_year)) / float(year_span)
        sankey_y = min(0.98, max(0.02, sankey_y))
        return 1.0 - sankey_y

    for depth in range(max_depth + 1):
        x_pos = depth_x(depth)
        is_major = depth == 0 or depth == max_depth or depth % 5 == 0
        shapes.append(
            {
                "type": "line",
                "xref": "paper",
                "yref": "paper",
                "x0": x_pos,
                "x1": x_pos,
                "y0": 0.02,
                "y1": 0.98,
                "layer": "below",
                "line": {
                    "color": (
                        "rgba(107,114,128,0.34)"
                        if is_major
                        else "rgba(148,163,184,0.20)"
                    ),
                    "width": 1.0 if is_major else 0.7,
                },
            }
        )
        annotations.append(
            {
                "xref": "paper",
                "yref": "paper",
                "x": x_pos,
                "y": 1.025,
                "text": f"d{depth}",
                "showarrow": False,
                "xanchor": "center",
                "yanchor": "bottom",
                "font": {"size": 10, "color": "#4b5563"},
            }
        )

    for year in range(min_year, max_year + 1):
        y_pos = year_y(year)
        is_major = year == min_year or year == max_year or year % year_step == 0
        if not is_major and year_span > 40:
            continue
        shapes.append(
            {
                "type": "line",
                "xref": "paper",
                "yref": "paper",
                "x0": x_min,
                "x1": horizontal_x_max,
                "y0": y_pos,
                "y1": y_pos,
                "layer": "below",
                "line": {
                    "color": (
                        "rgba(107,114,128,0.26)"
                        if is_major
                        else "rgba(148,163,184,0.10)"
                    ),
                    "width": 0.9 if is_major else 0.5,
                },
            }
        )
        if is_major:
            annotations.append(
                {
                    "xref": "paper",
                    "yref": "paper",
                    "x": float(x_min) - 0.02,
                    "y": y_pos,
                    "text": str(year),
                    "showarrow": False,
                    "xanchor": "right",
                    "yanchor": "middle",
                    "font": {"size": 10, "color": "#4b5563"},
                }
            )

    annotations.extend(
        [
            {
                "xref": "paper",
                "yref": "paper",
                "x": (float(x_min) + float(x_max)) / 2.0,
                "y": 1.065,
                "text": "depth",
                "showarrow": False,
                "xanchor": "center",
                "yanchor": "bottom",
                "font": {"size": 11, "color": "#374151"},
            },
            {
                "xref": "paper",
                "yref": "paper",
                "x": -0.055,
                "y": 0.5,
                "text": "year",
                "textangle": -90,
                "showarrow": False,
                "xanchor": "center",
                "yanchor": "middle",
                "font": {"size": 11, "color": "#374151"},
            },
        ]
    )
    return shapes, annotations


def _adaptive_sankey_time_density_elements(
    rows: list[dict[str, Any]],
    *,
    branches: list[str],
    color_map: dict[str, str],
    min_year: int,
    max_year: int,
    label: str,
    scale_max: float | None,
    x0: float,
    x1: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    min_year = int(min_year)
    max_year = int(max_year)
    year_span = max(max_year - min_year, 1)
    panel_width = max(float(x1) - float(x0), 0.01)
    branch_year_scores: dict[str, dict[int, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    for row in rows:
        branch = str(row.get("branch") or "")
        if branch not in branches:
            continue
        year = int(row.get("target_year"))
        if year < min_year or year > max_year:
            continue
        branch_year_scores[branch][year] += float(row.get("edge_score_abs") or 0.0)

    if scale_max is None:
        global_max = max(
            (
                float(value)
                for year_scores in branch_year_scores.values()
                for value in year_scores.values()
            ),
            default=0.0,
        )
    else:
        global_max = float(scale_max)

    shapes: list[dict[str, Any]] = [
        {
            "type": "rect",
            "xref": "paper",
            "yref": "paper",
            "x0": x0,
            "x1": x1,
            "y0": 0.02,
            "y1": 0.98,
            "layer": "below",
            "fillcolor": "rgba(255,255,255,0.00)",
            "line": {"color": "rgba(107,114,128,0.32)", "width": 0.8},
        },
        {
            "type": "line",
            "xref": "paper",
            "yref": "paper",
            "x0": x0,
            "x1": x0,
            "y0": 0.02,
            "y1": 0.98,
            "layer": "above",
            "line": {"color": "rgba(55,65,81,0.55)", "width": 0.8},
        },
    ]
    annotations: list[dict[str, Any]] = [
        {
            "xref": "paper",
            "yref": "paper",
            "x": (float(x0) + float(x1)) / 2.0,
            "y": 0.965,
            "text": label,
            "showarrow": False,
            "xanchor": "center",
            "yanchor": "top",
            "font": {"size": 10, "color": "#374151"},
        }
    ]

    def year_y(year: int) -> float:
        sankey_y = (float(year) - float(min_year)) / float(year_span)
        sankey_y = min(0.98, max(0.02, sankey_y))
        return 1.0 - sankey_y

    years_all = list(range(min_year, max_year + 1))
    for branch in branches:
        year_scores = branch_year_scores.get(branch, {})
        if not year_scores or global_max <= 0.0:
            continue
        points: list[tuple[float, float]] = []
        for year in years_all:
            value = float(year_scores.get(year, 0.0))
            x_pos = float(x0) + panel_width * (value / global_max)
            points.append((x_pos, year_y(year)))
        first_y = points[0][1]
        last_y = points[-1][1]
        path_parts = [f"M {float(x0):.5f},{first_y:.5f}"]
        path_parts.extend(f"L {x:.5f},{y:.5f}" for x, y in points)
        path_parts.append(f"L {float(x0):.5f},{last_y:.5f}")
        path_parts.append("Z")
        shapes.append(
            {
                "type": "path",
                "xref": "paper",
                "yref": "paper",
                "path": " ".join(path_parts),
                "layer": "above",
                "fillcolor": _adaptive_sankey_rgba_color(branch, color_map, 0.16),
                "line": {
                    "color": _adaptive_sankey_rgba_color(branch, color_map, 0.78),
                    "width": 0.9,
                },
            }
        )
    return shapes, annotations


def _adaptive_sankey_depth_density_elements(
    rows: list[dict[str, Any]],
    *,
    branches: list[str],
    color_map: dict[str, str],
    max_depth: int,
    label: str,
    scale_max: float | None,
    x_min: float,
    x_max: float,
    y0: float = -0.245,
    y1: float = -0.140,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    max_depth = max(int(max_depth), 1)
    panel_height = max(float(y1) - float(y0), 0.01)
    x_span = max(float(x_max) - float(x_min), 0.01)
    branch_depth_scores: dict[str, dict[int, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    for row in rows:
        branch = str(row.get("branch") or "")
        if branch not in branches:
            continue
        depth = int(row.get("depth", -1))
        if depth < 0 or depth > max_depth:
            continue
        branch_depth_scores[branch][depth] += float(row.get("node_score_abs") or 0.0)

    if scale_max is None:
        global_max = max(
            (
                float(value)
                for depth_scores in branch_depth_scores.values()
                for value in depth_scores.values()
            ),
            default=0.0,
        )
    else:
        global_max = float(scale_max)

    shapes: list[dict[str, Any]] = [
        {
            "type": "rect",
            "xref": "paper",
            "yref": "paper",
            "x0": x_min,
            "x1": x_max,
            "y0": y0,
            "y1": y1,
            "layer": "below",
            "fillcolor": "rgba(255,255,255,0.00)",
            "line": {"color": "rgba(107,114,128,0.32)", "width": 0.8},
        }
    ]
    annotations: list[dict[str, Any]] = [
        {
            "xref": "paper",
            "yref": "paper",
            "x": (float(x_min) + float(x_max)) / 2.0,
            "y": float(y1) + 0.014,
            "text": label,
            "showarrow": False,
            "xanchor": "center",
            "yanchor": "bottom",
            "font": {"size": 10, "color": "#374151"},
        }
    ]

    def depth_x(depth: int) -> float:
        return float(x_min) + (float(depth) / float(max_depth)) * x_span

    for depth in range(max_depth + 1):
        x_pos = depth_x(depth)
        is_major = depth == 0 or depth == max_depth or depth % 5 == 0
        shapes.append(
            {
                "type": "line",
                "xref": "paper",
                "yref": "paper",
                "x0": x_pos,
                "x1": x_pos,
                "y0": y0,
                "y1": y1,
                "layer": "below",
                "line": {
                    "color": (
                        "rgba(107,114,128,0.28)"
                        if is_major
                        else "rgba(148,163,184,0.16)"
                    ),
                    "width": 0.8 if is_major else 0.5,
                },
            }
        )

    depths_all = list(range(max_depth + 1))
    for branch in branches:
        depth_scores = branch_depth_scores.get(branch, {})
        if not depth_scores or global_max <= 0.0:
            continue
        points: list[tuple[float, float]] = []
        for depth in depths_all:
            value = float(depth_scores.get(depth, 0.0))
            y_pos = float(y0) + panel_height * (value / global_max)
            points.append((depth_x(depth), y_pos))
        first_x = points[0][0]
        last_x = points[-1][0]
        path_parts = [f"M {first_x:.5f},{float(y0):.5f}"]
        path_parts.extend(f"L {x:.5f},{y:.5f}" for x, y in points)
        path_parts.append(f"L {last_x:.5f},{float(y0):.5f}")
        path_parts.append("Z")
        shapes.append(
            {
                "type": "path",
                "xref": "paper",
                "yref": "paper",
                "path": " ".join(path_parts),
                "layer": "above",
                "fillcolor": _adaptive_sankey_rgba_color(branch, color_map, 0.16),
                "line": {
                    "color": _adaptive_sankey_rgba_color(branch, color_map, 0.78),
                    "width": 0.9,
                },
            }
        )
    return shapes, annotations


def _adaptive_sankey_legend_elements(
    entries: list[tuple[str, str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if not entries:
        return [], []

    x0 = 0.615
    x1 = 0.935
    y1 = 1.305
    row_height = 0.022
    header_height = 0.044
    pad = 0.010
    y0 = y1 - header_height - row_height * len(entries) - pad
    shapes: list[dict[str, Any]] = [
        {
            "type": "rect",
            "xref": "paper",
            "yref": "paper",
            "x0": x0,
            "x1": x1,
            "y0": y0,
            "y1": y1,
            "layer": "above",
            "fillcolor": "rgba(255,255,255,0.90)",
            "line": {"color": "rgba(107,114,128,0.55)", "width": 0.8},
        }
    ]
    annotations: list[dict[str, Any]] = [
        {
            "xref": "paper",
            "yref": "paper",
            "x": x0 + 0.014,
            "y": y1 - 0.014,
            "text": "Displayed root activities",
            "showarrow": False,
            "xanchor": "left",
            "yanchor": "top",
            "font": {"size": 11, "color": "#111827"},
        }
    ]

    for idx, (branch, color) in enumerate(entries):
        row_y = y1 - header_height - row_height * idx - 0.016
        square_y = row_y - 0.011
        activity_name = str(branch).split(" | ", maxsplit=1)[0]
        shapes.append(
            {
                "type": "rect",
                "xref": "paper",
                "yref": "paper",
                "x0": x0 + 0.014,
                "x1": x0 + 0.028,
                "y0": square_y,
                "y1": square_y + 0.014,
                "layer": "above",
                "fillcolor": color,
                "line": {"color": "rgba(31,41,55,0.55)", "width": 0.5},
            }
        )
        annotations.append(
            {
                "xref": "paper",
                "yref": "paper",
                "x": x0 + 0.036,
                "y": row_y,
                "text": _shorten_text(
                    activity_name,
                    width=46,
                    placeholder="...",
                ),
                "showarrow": False,
                "xanchor": "left",
                "yanchor": "middle",
                "font": {"size": 10, "color": "#1f2937"},
            }
        )

    return shapes, annotations


def plot_adaptive_sankey(
    trails: Trails,
    *,
    method: str | None = None,
    title: str | None = None,
    adaptive_relative_score_cutoff: float | None = None,
    branch_visual_cutoff: float = 0.001,
    max_sankey_links: int = 0,
    display_score_coverage: float = 1.0,
    width: int = 2100,
    height: int = 1100,
    show_time_density: bool = True,
    show_depth_density: bool = True,
    depth_axis_max: int | None = None,
    year_axis_min: int | None = None,
    year_axis_max: int | None = None,
    node_score_rows: list[dict[str, Any]] | None = None,
    time_density_label: str = "score-potential density<br>over time",
    depth_density_label: str = "score-potential density over depth",
    output_path: str | os.PathLike[str] | None = None,
    png_path: str | os.PathLike[str] | None = None,
    png_scale: int = 3,
) -> go.Figure:
    """Plot explicit routed graph edges from adaptive temporal routing.

    The diagram uses only the explicit routed graph stored on ``trails.graph``.
    Nodes are placed horizontally by routing depth and vertically by calendar
    year. Link width is the child node's adaptive routing score potential,
    allocated to first-level root branches by routed absolute amount. Node
    labels are hidden in the plot body and shown in hover text.

    :param trails: Trails object with an existing routed graph.
    :type trails: Trails
    :param method: Optional method label to show in the title and metadata.
    :type method: str | None
    :param title: Optional figure title. If omitted, a routing summary title is
        generated from ``trails._routing_params``.
    :type title: str | None
    :param adaptive_relative_score_cutoff: Optional cutoff label. If omitted,
        the value stored in ``trails._routing_params`` is used when available.
    :type adaptive_relative_score_cutoff: float | None
    :param branch_visual_cutoff: Minimum branch share of total explicit edge
        score potential to display.
    :type branch_visual_cutoff: float
    :param max_sankey_links: Maximum links to draw. Use ``0`` for no hard cap.
    :type max_sankey_links: int
    :param display_score_coverage: Stop selecting links after this fraction of
        candidate explicit edge score potential is covered.
    :type display_score_coverage: float
    :param width: Figure width in pixels.
    :type width: int
    :param height: Figure height in pixels.
    :type height: int
    :param show_time_density: Include the right-side year density panel.
    :type show_time_density: bool
    :param show_depth_density: Include the bottom depth density panel.
    :type show_depth_density: bool
    :param depth_axis_max: Optional maximum depth shown on the fixed depth
        axis. Use this to make multiple Sankey plots share the same depth
        spacing even when they contain different maximum depths.
    :type depth_axis_max: int | None
    :param year_axis_min: Optional minimum calendar year shown on the fixed
        year axis. The value is expanded when needed to include routed nodes.
    :type year_axis_min: int | None
    :param year_axis_max: Optional maximum calendar year shown on the fixed
        year axis. The value is expanded when needed to include routed nodes.
    :type year_axis_max: int | None
    :param node_score_rows: Optional rows with ``branch``, ``depth``, and
        ``node_score_abs`` keys for the bottom panel. If omitted, routed node
        score potentials are used.
    :type node_score_rows: list[dict[str, Any]] | None
    :param time_density_label: Label for the right-side density panel.
    :type time_density_label: str
    :param depth_density_label: Label for the bottom density panel.
    :type depth_density_label: str
    :param output_path: Optional HTML path to write.
    :type output_path: str | os.PathLike[str] | None
    :param png_path: Optional PNG path to write using Kaleido.
    :type png_path: str | os.PathLike[str] | None
    :param png_scale: PNG scale factor.
    :type png_scale: int
    :returns: Plotly figure.
    :rtype: go.Figure
    :raises RuntimeError: If the routed graph is missing.
    :raises ValueError: If no score-potential graph edges can be plotted."""
    graph = getattr(trails, "graph", None)
    if graph is None:
        raise RuntimeError(
            "Trails graph is missing; run trails.temporal_routing(...) first."
        )

    routing_params = getattr(trails, "_routing_params", {}) or {}
    if adaptive_relative_score_cutoff is None:
        stored_cutoff = routing_params.get("adaptive_relative_score_cutoff")
        adaptive_relative_score_cutoff = (
            None if stored_cutoff is None else float(stored_cutoff)
        )

    sankey_x_min = 0.02
    sankey_x_max = 0.81 if show_time_density else 0.96
    density_x_min = 0.835
    density_x_max = 0.933

    all_edge_rows = _adaptive_sankey_edge_rows(trails, method=method)
    if not all_edge_rows:
        raise ValueError("No routed graph score-potential edges to plot.")
    branch_rows = _adaptive_sankey_branch_rows(all_edge_rows)
    if node_score_rows is None:
        node_score_rows = _adaptive_sankey_node_potential_rows(
            trails,
            method=method,
        )

    node_lookup = {repr(node): node for node in graph.nodes}
    allowed_branches = _adaptive_sankey_allowed_branches(
        branch_rows,
        cutoff=float(branch_visual_cutoff),
    )
    candidate_rows = [
        row for row in all_edge_rows if str(row.get("branch") or "") in allowed_branches
    ]
    candidate_rows.sort(key=lambda row: float(row["edge_score_abs"]), reverse=True)
    if not candidate_rows:
        raise ValueError("No routed graph score-potential edges pass filters.")

    total_abs = sum(float(row["edge_score_abs"]) for row in candidate_rows)
    target = min(max(float(display_score_coverage), 0.0), 1.0)
    max_links = None if int(max_sankey_links) <= 0 else max(1, int(max_sankey_links))

    def select_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
        selected: list[dict[str, Any]] = []
        running = 0.0
        row_total = sum(float(row["edge_score_abs"]) for row in rows)
        for row in rows:
            selected.append(row)
            running += float(row["edge_score_abs"])
            if max_links is not None and len(selected) >= max_links:
                break
            if row_total and target < 1.0 and running / row_total >= target:
                break
        return selected

    selected_rows = select_rows(candidate_rows)

    connected_nodes: set[object] = set()
    for row in selected_rows:
        source = node_lookup.get(str(row["source_node"]))
        target_node = node_lookup.get(str(row["target_node"]))
        if source is None or target_node is None or source == target_node:
            continue
        connected_nodes.add(source)
        connected_nodes.add(target_node)
    if not connected_nodes:
        raise ValueError("No non-self routed graph edges to plot.")

    ordered_graph_nodes = sorted(
        connected_nodes,
        key=lambda node: (
            int(graph.nodes[node].get("depth", 0)),
            int(graph.nodes[node].get("year", 0)),
            _adaptive_sankey_reference_product(
                trails,
                int(graph.nodes[node].get("act_idx", -1)),
            ),
            repr(node),
        ),
    )
    node_ids = {node: idx for idx, node in enumerate(ordered_graph_nodes)}
    labels: list[str] = []
    years: list[int] = []
    depths: list[int] = []
    node_customdata: list[list[Any]] = []
    for node in ordered_graph_nodes:
        data = graph.nodes[node]
        act_idx = int(data.get("act_idx", -1))
        meta = _adaptive_sankey_activity_metadata(trails, act_idx)
        reference_product = str(meta.get("reference product") or "")
        name = str(meta.get("name") or "")
        location = str(meta.get("location") or "")
        year = int(data.get("year", -1))
        depth = int(data.get("depth", -1))
        full_label = " | ".join(
            part for part in (name, reference_product, location) if part
        )
        labels.append(
            _shorten_text(
                reference_product or name or "Activity",
                width=42,
                placeholder="...",
            )
        )
        years.append(year)
        depths.append(depth)
        node_customdata.append(
            [
                full_label or repr(node),
                year,
                depth,
                str(meta.get("unit") or ""),
                "explicit routed graph node",
                float(data.get("score_potential") or 0.0),
                float(data.get("amount") or 0.0),
                float(data.get("frontier_amount") or 0.0),
            ]
        )

    sources: list[int] = []
    targets: list[int] = []
    values: list[float] = []
    colors: list[str] = []
    link_customdata: list[list[Any]] = []
    color_map: dict[str, str] = {}
    for row in selected_rows:
        source = node_lookup.get(str(row["source_node"]))
        target_node = node_lookup.get(str(row["target_node"]))
        if source is None or target_node is None or source == target_node:
            continue
        value = float(row["edge_score_abs"])
        if value <= 0.0:
            continue
        branch = str(row["branch"])
        source_data = graph.nodes[source]
        target_data = graph.nodes[target_node]
        sources.append(node_ids[source])
        targets.append(node_ids[target_node])
        values.append(value)
        colors.append(_adaptive_sankey_edge_color(branch, color_map))
        link_customdata.append(
            [
                _adaptive_sankey_full_label(
                    trails,
                    int(source_data.get("act_idx", -1)),
                ),
                int(source_data.get("year", -1)),
                int(source_data.get("depth", -1)),
                _adaptive_sankey_full_label(
                    trails,
                    int(target_data.get("act_idx", -1)),
                ),
                int(target_data.get("year", -1)),
                int(target_data.get("depth", -1)),
                branch,
                float(row.get("raw_amount_abs") or 0.0),
                float(row.get("root_routed_amount_abs") or 0.0),
                value,
                float(row.get("child_node_score") or 0.0),
                "explicit routed graph edge",
            ]
        )

    if not values:
        raise ValueError("No non-self routed graph links to plot.")

    max_depth = max(max(depths), 1)
    depth_axis_max = max(int(depth_axis_max or max_depth), max_depth, 1)
    data_min_year = min(years)
    data_max_year = max(years)
    year_axis_min = min(
        int(year_axis_min) if year_axis_min is not None else data_min_year,
        data_min_year,
    )
    year_axis_max = max(
        int(year_axis_max) if year_axis_max is not None else data_max_year,
        data_max_year,
    )
    x_values, y_values = _adaptive_sankey_year_depth_positions(
        years=years,
        depths=depths,
        labels=labels,
        x_min=sankey_x_min,
        x_max=sankey_x_max,
        max_depth=depth_axis_max,
        min_year=year_axis_min,
        max_year=year_axis_max,
    )
    min_year = int(year_axis_min)
    max_year = int(year_axis_max)
    grid_shapes, grid_annotations = _adaptive_sankey_grid_elements(
        max_depth=depth_axis_max,
        min_year=min_year,
        max_year=max_year,
        x_min=sankey_x_min,
        x_max=sankey_x_max,
        horizontal_x_max=sankey_x_max,
    )
    legend_entries = [
        (branch, _adaptive_sankey_edge_color(branch, color_map))
        for branch in sorted(
            color_map,
            key=lambda item: next(
                (
                    float(row["score_abs_share"])
                    for row in branch_rows
                    if str(row["branch"]) == item
                ),
                0.0,
            ),
            reverse=True,
        )
    ]
    legend_shapes, legend_annotations = _adaptive_sankey_legend_elements(legend_entries)
    displayed_branches = [branch for branch, _color in legend_entries]

    time_density_scale = max(
        (
            sum(
                float(row.get("edge_score_abs") or 0.0)
                for row in selected_rows
                if str(row.get("branch") or "") == branch
                and int(row.get("target_year")) == year
            )
            for branch in displayed_branches
            for year in range(min_year, max_year + 1)
        ),
        default=0.0,
    )
    branch_depth_scores_for_scale: dict[tuple[str, int], float] = defaultdict(float)
    for row in node_score_rows:
        branch = str(row.get("branch") or "")
        if branch not in displayed_branches:
            continue
        branch_depth_scores_for_scale[
            (branch, int(row.get("depth", -1)))
        ] += float(row.get("node_score_abs") or 0.0)
    depth_density_scale = max(branch_depth_scores_for_scale.values(), default=0.0)

    density_shapes: list[dict[str, Any]] = []
    density_annotations: list[dict[str, Any]] = []
    if show_time_density:
        density_shapes, density_annotations = _adaptive_sankey_time_density_elements(
            selected_rows,
            branches=displayed_branches,
            color_map=color_map,
            min_year=min_year,
            max_year=max_year,
            label=time_density_label,
            scale_max=time_density_scale,
            x0=density_x_min,
            x1=density_x_max,
        )

    depth_density_shapes: list[dict[str, Any]] = []
    depth_density_annotations: list[dict[str, Any]] = []
    if show_depth_density:
        depth_density_shapes, depth_density_annotations = (
            _adaptive_sankey_depth_density_elements(
                node_score_rows,
                branches=displayed_branches,
                color_map=color_map,
                max_depth=depth_axis_max,
                label=depth_density_label,
                scale_max=depth_density_scale,
                x_min=sankey_x_min,
                x_max=sankey_x_max,
            )
        )

    display_abs = sum(values)
    display_share = display_abs / total_abs if total_abs else 0.0
    cutoff_label = (
        "disabled"
        if adaptive_relative_score_cutoff is None
        else f"{float(adaptive_relative_score_cutoff):.0e}"
    )
    if title is None:
        method_part = f"<br>{method}" if method else ""
        title = (
            "Adaptive-routing Sankey"
            f"{method_part}<br>"
            "explicit graph edges only; "
            "width=adaptive routing score potential; "
            f"adaptive cutoff={cutoff_label}; "
            f"displayed potential share={display_share:.1%}"
        )

    node_hovertemplate = (
        "%{customdata[0]}<br>"
        "year=%{customdata[1]} | depth=%{customdata[2]}<br>"
        "node score potential=%{customdata[5]:.3e}<br>"
        "amount=%{customdata[6]:.3e}; "
        "frontier amount=%{customdata[7]:.3e}<br>"
        "%{customdata[4]}<extra></extra>"
    )
    link_hovertemplate = (
        "source=%{customdata[0]} "
        "(%{customdata[1]}, d%{customdata[2]})<br>"
        "target=%{customdata[3]} "
        "(%{customdata[4]}, d%{customdata[5]})<br>"
        "root branch=%{customdata[6]}<br>"
        "%{customdata[11]}<br>"
        "edge score potential=%{customdata[9]:.3e}<br>"
        "target node potential=%{customdata[10]:.3e}<br>"
        "edge amount=%{customdata[7]:.3e}; "
        "root-routed amount=%{customdata[8]:.3e}"
        "<extra></extra>"
    )

    all_trace = go.Sankey(
        arrangement="fixed",
        name="All branches",
        visible=True,
        node=dict(
            label=[""] * len(labels),
            x=x_values,
            y=y_values,
            pad=10,
            thickness=12,
            color="#d1d5db",
            customdata=node_customdata,
            hovertemplate=node_hovertemplate,
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=colors,
            customdata=link_customdata,
            hovertemplate=link_hovertemplate,
        ),
    )

    def build_branch_trace(
        branch: str,
        rows: list[dict[str, Any]],
    ) -> go.Sankey | None:
        branch_nodes: set[object] = set()
        valid_rows: list[dict[str, Any]] = []
        for row in rows:
            source = node_lookup.get(str(row["source_node"]))
            target_node = node_lookup.get(str(row["target_node"]))
            if source is None or target_node is None or source == target_node:
                continue
            valid_rows.append(row)
            branch_nodes.add(source)
            branch_nodes.add(target_node)
        if not valid_rows or not branch_nodes:
            return None

        ordered_nodes = sorted(
            branch_nodes,
            key=lambda node: (
                int(graph.nodes[node].get("depth", 0)),
                int(graph.nodes[node].get("year", 0)),
                _adaptive_sankey_reference_product(
                    trails,
                    int(graph.nodes[node].get("act_idx", -1)),
                ),
                repr(node),
            ),
        )
        branch_node_ids = {node: idx for idx, node in enumerate(ordered_nodes)}
        branch_labels: list[str] = []
        branch_years: list[int] = []
        branch_depths: list[int] = []
        branch_node_customdata: list[list[Any]] = []
        for node in ordered_nodes:
            data = graph.nodes[node]
            act_idx = int(data.get("act_idx", -1))
            meta = _adaptive_sankey_activity_metadata(trails, act_idx)
            reference_product = str(meta.get("reference product") or "")
            name = str(meta.get("name") or "")
            location = str(meta.get("location") or "")
            year = int(data.get("year", -1))
            depth = int(data.get("depth", -1))
            full_label = " | ".join(
                part for part in (name, reference_product, location) if part
            )
            branch_labels.append(
                _shorten_text(
                    reference_product or name or "Activity",
                    width=42,
                    placeholder="...",
                )
            )
            branch_years.append(year)
            branch_depths.append(depth)
            branch_node_customdata.append(
                [
                    full_label or repr(node),
                    year,
                    depth,
                    str(meta.get("unit") or ""),
                    "explicit routed graph node",
                    float(data.get("score_potential") or 0.0),
                    float(data.get("amount") or 0.0),
                    float(data.get("frontier_amount") or 0.0),
                ]
            )

        branch_sources: list[int] = []
        branch_targets: list[int] = []
        branch_values: list[float] = []
        branch_colors: list[str] = []
        branch_link_customdata: list[list[Any]] = []
        for row in valid_rows:
            source = node_lookup.get(str(row["source_node"]))
            target_node = node_lookup.get(str(row["target_node"]))
            if source is None or target_node is None or source == target_node:
                continue
            value = float(row["edge_score_abs"])
            if value <= 0.0:
                continue
            source_data = graph.nodes[source]
            target_data = graph.nodes[target_node]
            branch_sources.append(branch_node_ids[source])
            branch_targets.append(branch_node_ids[target_node])
            branch_values.append(value)
            branch_colors.append(_adaptive_sankey_edge_color(branch, color_map))
            branch_link_customdata.append(
                [
                    _adaptive_sankey_full_label(
                        trails,
                        int(source_data.get("act_idx", -1)),
                    ),
                    int(source_data.get("year", -1)),
                    int(source_data.get("depth", -1)),
                    _adaptive_sankey_full_label(
                        trails,
                        int(target_data.get("act_idx", -1)),
                    ),
                    int(target_data.get("year", -1)),
                    int(target_data.get("depth", -1)),
                    branch,
                    float(row.get("raw_amount_abs") or 0.0),
                    float(row.get("root_routed_amount_abs") or 0.0),
                    value,
                    float(row.get("child_node_score") or 0.0),
                    "explicit routed graph edge",
                ]
            )

        if not branch_values:
            return None
        branch_x, branch_y = _adaptive_sankey_year_depth_positions(
            years=branch_years,
            depths=branch_depths,
            labels=branch_labels,
            x_min=sankey_x_min,
            x_max=sankey_x_max,
            max_depth=depth_axis_max,
            min_year=year_axis_min,
            max_year=year_axis_max,
        )
        return go.Sankey(
            arrangement="fixed",
            name=_shorten_text(branch, width=80, placeholder="..."),
            visible=False,
            node=dict(
                label=[""] * len(branch_labels),
                x=branch_x,
                y=branch_y,
                pad=10,
                thickness=12,
                color="#d1d5db",
                customdata=branch_node_customdata,
                hovertemplate=node_hovertemplate,
            ),
            link=dict(
                source=branch_sources,
                target=branch_targets,
                value=branch_values,
                color=branch_colors,
                customdata=branch_link_customdata,
                hovertemplate=link_hovertemplate,
            ),
        )

    traces: list[Any] = [all_trace]
    trace_labels = ["All branches"]
    rows_by_branch: dict[str, list[dict[str, Any]]] = defaultdict(list)
    node_rows_by_branch: dict[str, list[dict[str, Any]]] = defaultdict(list)
    panel_shapes_by_trace: list[list[dict[str, Any]]] = [
        density_shapes + depth_density_shapes
    ]
    for row in candidate_rows:
        rows_by_branch[str(row["branch"])].append(row)
    for row in node_score_rows:
        node_rows_by_branch[str(row["branch"])].append(row)
    for branch, _color in legend_entries:
        rows = select_rows(rows_by_branch.get(branch, []))
        trace = build_branch_trace(branch, rows)
        if trace is None:
            continue
        traces.append(trace)
        trace_labels.append(_shorten_text(branch, width=80, placeholder="..."))
        branch_panel_shapes: list[dict[str, Any]] = []
        if show_time_density:
            branch_density_shapes, _branch_density_annotations = (
                _adaptive_sankey_time_density_elements(
                    rows,
                    branches=[branch],
                    color_map=color_map,
                    min_year=min_year,
                    max_year=max_year,
                    label=time_density_label,
                    scale_max=time_density_scale,
                    x0=density_x_min,
                    x1=density_x_max,
                )
            )
            branch_panel_shapes.extend(branch_density_shapes)
        if show_depth_density:
            branch_depth_shapes, _branch_depth_annotations = (
                _adaptive_sankey_depth_density_elements(
                    node_rows_by_branch.get(branch, []),
                    branches=[branch],
                    color_map=color_map,
                    max_depth=depth_axis_max,
                    label=depth_density_label,
                    scale_max=depth_density_scale,
                    x_min=sankey_x_min,
                    x_max=sankey_x_max,
                )
            )
            branch_panel_shapes.extend(branch_depth_shapes)
        panel_shapes_by_trace.append(branch_panel_shapes)

    fig = go.Figure(data=traces)
    base_shapes = grid_shapes + legend_shapes
    base_annotations = (
        grid_annotations
        + legend_annotations
        + density_annotations
        + depth_density_annotations
    )
    if len(traces) > 1:
        buttons = []
        for idx, label in enumerate(trace_labels):
            visible = [False] * len(traces)
            visible[idx] = True
            buttons.append(
                {
                    "label": label,
                    "method": "update",
                    "args": [
                        {"visible": visible},
                        {"shapes": base_shapes + panel_shapes_by_trace[idx]},
                    ],
                }
            )
        fig.update_layout(
            updatemenus=[
                {
                    "active": 0,
                    "buttons": buttons,
                    "direction": "down",
                    "showactive": True,
                    "x": 0.01,
                    "xanchor": "left",
                    "y": 1.24,
                    "yanchor": "top",
                }
            ]
        )

    stats = {
        "sankey_display_links": int(len(values)),
        "sankey_display_nodes": int(len(labels)),
        "sankey_display_abs_score": float(display_abs),
        "sankey_total_abs_edge_score": float(total_abs),
        "sankey_display_abs_score_share": float(display_share),
        "sankey_display_min_abs_edge_score": float(min(values) if values else 0.0),
        "sankey_original_edge_count": int(len(candidate_rows)),
        "sankey_selected_seed_edges": int(len(selected_rows)),
        "sankey_aggregation": "explicit_routed_graph_edges",
        "sankey_score_basis": "adaptive_routing_score_potential",
        "graph_nodes": int(graph.number_of_nodes()),
        "graph_edges": int(graph.number_of_edges()),
        "graph_max_depth": int(
            max(
                (int(data.get("depth", 0)) for _node, data in graph.nodes(data=True)),
                default=0,
            )
        ),
        "year_axis_min": int(year_axis_min),
        "year_axis_max": int(year_axis_max),
    }
    fig.update_layout(
        title=dict(
            text=title,
            x=0.01,
            xanchor="left",
            y=0.955,
            yanchor="top",
        ),
        width=int(width),
        height=int(height),
        margin=dict(l=78, r=22, t=205, b=150 if show_depth_density else 30),
        shapes=base_shapes + density_shapes + depth_density_shapes,
        annotations=base_annotations,
        font=dict(size=11),
        meta=stats,
    )

    if output_path is not None:
        html_target = Path(output_path)
        html_target.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(str(html_target))
    if png_path is not None:
        png_target = Path(png_path)
        png_target.parent.mkdir(parents=True, exist_ok=True)
        png_fig = go.Figure(data=[fig.data[0]], layout=fig.layout)
        png_fig.layout.updatemenus = ()
        png_fig.write_image(
            str(png_target),
            format="png",
            scale=int(png_scale),
        )
    return fig


def _node_scores_from_trails(trails: Trails) -> dict[tuple[int, int], float]:
    """node scores from trails.

    :param trails: Value for `trails`.
    :type trails: Trails
    :returns: Return value.
    :rtype: dict[tuple[int, int], float]
    :raises ValueError: If an error occurs."""
    scores = getattr(trails, "scores", None)
    if scores is None:
        characterized = getattr(trails, "characterized_inventory", None)
        if characterized is None:
            raise ValueError(
                "No scores available. Run lcia(...) or "
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

    if is_chunked_sparse(data):
        data = data.compute(scheduler="synchronous")

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
    """node scores from characterized inventory.

    :param characterized_inventory: Value for `characterized_inventory`.
    :type characterized_inventory: xr.DataArray
    :returns: Return value.
    :rtype: dict[tuple[int, int], float]
    :raises ValueError: If an error occurs."""
    if "method" in characterized_inventory.dims:
        methods = characterized_inventory.coords["method"].values
        if len(methods) != 1:
            raise ValueError(
                "characterized_inventory has multiple methods; select one before plotting."
            )
        characterized_inventory = characterized_inventory.isel(method=0, drop=True)
    if "activity" not in characterized_inventory.dims:
        raise ValueError(
            "characterized_inventory must include an 'activity' dimension."
        )
    if "year" not in characterized_inventory.dims:
        raise ValueError("characterized_inventory must include a 'year' dimension.")
    if "flow" not in characterized_inventory.dims:
        raise ValueError("characterized_inventory must include a 'flow' dimension.")

    years = characterized_inventory.coords["year"].values
    activities = characterized_inventory.coords["activity"].values

    out: dict[tuple[int, int], float] = {}

    reduce_dims = ["flow"]
    if "root activity" in characterized_inventory.dims:
        reduce_dims.append("root activity")
    reduced = characterized_inventory.sum(dim=reduce_dims).transpose(
        "activity", "year"
    )
    data = _compute_plot_array(reduced)
    if isinstance(data, sparse.COO):
        for ai, yi, v in zip(data.coords[0], data.coords[1], data.data):
            if v == 0.0:
                continue
            year = int(years[int(yi)])
            act = int(activities[int(ai)])
            out[(year, act)] = out.get((year, act), 0.0) + float(v)
        return out

    dense = np.asarray(data)
    if dense.ndim != 2:
        raise ValueError("Expected activity/year characterized inventory reduction.")
    idxs = np.nonzero(dense)
    for ai, yi in zip(idxs[0], idxs[1]):
        v = float(dense[ai, yi])
        if v == 0.0:
            continue
        year = int(years[int(yi)])
        act = int(activities[int(ai)])
        out[(year, act)] = out.get((year, act), 0.0) + v

    return out


def _select_depths(
    edges_by_depth: dict[int, dict], depths: list[int] | None
) -> list[int]:
    """select depths.

    :param edges_by_depth: Value for `edges_by_depth`.
    :type edges_by_depth: dict[int, dict]
    :param depths: Value for `depths`.
    :type depths: list[int] | None
    :returns: Return value.
    :rtype: list[int]
    :raises ValueError: If an error occurs."""
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
    """collect activities.

    :param edges_by_depth: Value for `edges_by_depth`.
    :type edges_by_depth: dict[int, dict]
    :param trails: Value for `trails`.
    :type trails: Trails
    :param depths_list: Value for `depths_list`.
    :type depths_list: list[int]
    :param include_all_activities: Value for `include_all_activities`.
    :type include_all_activities: bool
    :returns: Return value.
    :rtype: list[int]
    :raises ValueError: If an error occurs."""
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
    """collect global years.

    :param edges_by_depth: Value for `edges_by_depth`.
    :type edges_by_depth: dict[int, dict]
    :param depths_list: Value for `depths_list`.
    :type depths_list: list[int]
    :returns: Return value.
    :rtype: tuple[int, int, list[int]]
    :raises ValueError: If an error occurs."""
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
    """init flow subplots.

    :param panel_labels: Value for `panel_labels`.
    :type panel_labels: list[str]
    :param ncols: Value for `ncols`.
    :type ncols: int
    :returns: Return value.
    :rtype: tuple[go.Figure, int]"""
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
    """Plot rf.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param by: Value for `by`.
    :type by: Literal['flow', 'root activity']
    :param title: Value for `title`.
    :type title: str
    :param method_label: Value for `method_label`.
    :type method_label: str
    :param quantile: Value for `quantile`.
    :type quantile: float | None
    :param show_cumulative_quantile_band: Value for `show_cumulative_quantile_band`.
    :type show_cumulative_quantile_band: bool
    :param band_quantiles: Value for `band_quantiles`.
    :type band_quantiles: tuple[float, float]
    :param cumulative: Value for `cumulative`.
    :type cumulative: bool
    :param stacked: Value for `stacked`.
    :type stacked: bool
    :param legend_top_n: Value for `legend_top_n`.
    :type legend_top_n: int
    :param width: Value for `width`.
    :type width: Optional[int]
    :param height: Value for `height`.
    :type height: Optional[int]
    :param year_tick: Value for `year_tick`.
    :type year_tick: int
    :param year_range: Value for `year_range`.
    :type year_range: Optional[Tuple[int, int]]
    :param show_year_grid: Value for `show_year_grid`.
    :type show_year_grid: bool
    :param yaxis_type: Value for `yaxis_type`.
    :type yaxis_type: Literal['linear', 'log']
    :param log_eps: Value for `log_eps`.
    :type log_eps: float
    :param reference_year: Value for `reference_year`.
    :type reference_year: Optional[int]
    :param show_cumulative_axis: Value for `show_cumulative_axis`.
    :type show_cumulative_axis: bool
    :param cumulative_axis_label: Value for `cumulative_axis_label`.
    :type cumulative_axis_label: str
    :param legend_entrywidth: Value for `legend_entrywidth`.
    :type legend_entrywidth: int
    :param legend_row_height: Value for `legend_row_height`.
    :type legend_row_height: int
    :param legend_y: Value for `legend_y`.
    :type legend_y: float
    :param y2_headroom: Value for `y2_headroom`.
    :type y2_headroom: float
    :param show_cumulative_in_legend: Value for `show_cumulative_in_legend`.
    :type show_cumulative_in_legend: bool
    :param flow_groupby_name: Value for `flow_groupby_name`.
    :type flow_groupby_name: bool
    :param y_min: Value for `y_min`.
    :type y_min: Optional[float]
    :param y_max: Value for `y_max`.
    :type y_max: Optional[float]
    :param y2_max: Value for `y2_max`.
    :type y2_max: Optional[float]
    :returns: Return value.
    :rtype: go.Figure
    :raises ValueError: If an error occurs."""
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
            """results for.

            :param data: Value for `data`.
            :type data: xr.DataArray
            :returns: Return value.
            :rtype: dict[int, dict[str, Any]]"""
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
            """totals.

            :param res: Value for `res`.
            :type res: dict[int, dict[str, Any]]
            :param years_seq: Value for `years_seq`.
            :type years_seq: list[int]
            :returns: Return value.
            :rtype: np.ndarray"""
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
    """Plot temp.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param by: Value for `by`.
    :type by: Literal['flow', 'root activity']
    :param title: Value for `title`.
    :type title: str
    :param method_label: Value for `method_label`.
    :type method_label: str
    :param quantile: Value for `quantile`.
    :type quantile: float | None
    :param stacked: Value for `stacked`.
    :type stacked: bool
    :param legend_top_n: Value for `legend_top_n`.
    :type legend_top_n: int
    :param width: Value for `width`.
    :type width: Optional[int]
    :param height: Value for `height`.
    :type height: Optional[int]
    :param year_tick: Value for `year_tick`.
    :type year_tick: int
    :param year_range: Value for `year_range`.
    :type year_range: Optional[Tuple[int, int]]
    :param show_year_grid: Value for `show_year_grid`.
    :type show_year_grid: bool
    :param yaxis_type: Value for `yaxis_type`.
    :type yaxis_type: Literal['linear', 'log']
    :param log_eps: Value for `log_eps`.
    :type log_eps: float
    :param reference_year: Value for `reference_year`.
    :type reference_year: Optional[int]
    :param show_total_axis: Value for `show_total_axis`.
    :type show_total_axis: bool
    :param show_total_quantile_band: Value for `show_total_quantile_band`.
    :type show_total_quantile_band: bool
    :param total_axis_label: Value for `total_axis_label`.
    :type total_axis_label: str
    :param legend_entrywidth: Value for `legend_entrywidth`.
    :type legend_entrywidth: int
    :param legend_row_height: Value for `legend_row_height`.
    :type legend_row_height: int
    :param legend_y: Value for `legend_y`.
    :type legend_y: float
    :param y2_headroom: Value for `y2_headroom`.
    :type y2_headroom: float
    :param show_cumulative_in_legend: Value for `show_cumulative_in_legend`.
    :type show_cumulative_in_legend: bool
    :param flow_groupby_name: Value for `flow_groupby_name`.
    :type flow_groupby_name: bool
    :param y_min: Value for `y_min`.
    :type y_min: Optional[float]
    :param y_max: Value for `y_max`.
    :type y_max: Optional[float]
    :param y2_max: Value for `y2_max`.
    :type y2_max: Optional[float]
    :returns: Return value.
    :rtype: go.Figure
    :raises ValueError: If an error occurs."""
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
            """totals.

            :param res: Value for `res`.
            :type res: dict[int, dict[str, Any]]
            :param years_seq: Value for `years_seq`.
            :type years_seq: list[int]
            :returns: Return value.
            :rtype: np.ndarray"""
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
                    """results for.

                    :param data: Value for `data`.
                    :type data: xr.DataArray
                    :returns: Return value.
                    :rtype: dict[int, dict[str, Any]]"""
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
    fig.update_yaxes(tickformat=".2e")
    return fig


def plot_delta_temperature(*args: Any, **kwargs: Any) -> go.Figure:
    """Plot delta temperature.

    :param args: Variadic positional arguments.
    :type args: Any
    :param kwargs: Variadic keyword arguments.
    :type kwargs: Any
    :returns: Return value.
    :rtype: go.Figure"""
    return plot_temp(*args, **kwargs)


def _configure_flow_axes(
    fig: go.Figure,
    n_panels: int,
    ncols: int,
    acts: list[int],
    year_min: int,
    year_max: int,
) -> None:
    """configure flow axes.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param n_panels: Value for `n_panels`.
    :type n_panels: int
    :param ncols: Value for `ncols`.
    :type ncols: int
    :param acts: Value for `acts`.
    :type acts: list[int]
    :param year_min: Value for `year_min`.
    :type year_min: int
    :param year_max: Value for `year_max`.
    :type year_max: int"""
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
    """merge all edges.

    :param edges_by_depth: Value for `edges_by_depth`.
    :type edges_by_depth: dict[int, dict]
    :param depths_list: Value for `depths_list`.
    :type depths_list: list[int]
    :returns: Return value.
    :rtype: dict[tuple[tuple[int, int], tuple[int, int]], float]"""
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
    """add flow panel traces.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param edges: Value for `edges`.
    :type edges: dict[tuple[tuple[int, int], tuple[int, int]], float]
    :param act_to_row: Value for `act_to_row`.
    :type act_to_row: dict[int, int]
    :param idx_to_label: Value for `idx_to_label`.
    :type idx_to_label: dict[int, str]
    :param dot_size: Value for `dot_size`.
    :type dot_size: int
    :param row: Value for `row`.
    :type row: int
    :param col: Value for `col`.
    :type col: int
    :param show_legend: Value for `show_legend`.
    :type show_legend: bool"""
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
    """add activity legend.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param acts: Value for `acts`.
    :type acts: list[int]
    :param idx_to_label: Value for `idx_to_label`.
    :type idx_to_label: dict[int, str]"""
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
    """apply flow layout.

    :param fig: Value for `fig`.
    :type fig: go.Figure
    :param title: Value for `title`.
    :type title: str
    :param base_width: Value for `base_width`.
    :type base_width: int
    :param base_height: Value for `base_height`.
    :type base_height: int
    :param nrows: Value for `nrows`.
    :type nrows: int"""
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
    """Plot traversal grid flow.

    :param edges_by_depth: Value for `edges_by_depth`.
    :type edges_by_depth: dict[int, dict[tuple[tuple[int, int], tuple[int, int]], float]]
    :param trails: Value for `trails`.
    :type trails: Trails
    :param depths: Value for `depths`.
    :type depths: Optional[List[int]]
    :param include_all_activities: Value for `include_all_activities`.
    :type include_all_activities: bool
    :param title: Value for `title`.
    :type title: str
    :param dot_size: Value for `dot_size`.
    :type dot_size: int
    :param base_width: Value for `base_width`.
    :type base_width: int
    :param base_height: Value for `base_height`.
    :type base_height: int
    :returns: Return value.
    :rtype: go.Figure"""
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
