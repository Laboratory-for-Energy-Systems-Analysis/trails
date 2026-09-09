"""Render a pyvis/vis-network Trails graph HTML file as a static PNG.

The publication notebook writes ``trails_graph.html`` with fixed node positions
for interactive inspection. This script reuses the embedded node/edge datasets
and redraws the graph with Matplotlib in a style close to manuscript Figure 2.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Any

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "trails-matplotlib-cache"),
)

import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch


DEPTH_BAND_COLORS = [
    "#f1f3f5",
    "#eaf2fb",
    "#f5f0e8",
    "#edf7ef",
    "#f5eef8",
    "#eef6f6",
]

FONT_SIZE_DEPTH = 13
FONT_SIZE_LABEL = 12
FONT_SIZE_TICK = 14
FONT_SIZE_AXIS = 14
FONT_SIZE_NOTE = 10

EDGE_WIDTH_MIN = 0.18
EDGE_WIDTH_MAX = 3.0
EDGE_WIDTH_EXPONENT = 0.5
NODE_DIAMETER_MIN = 5.5
NODE_DIAMETER_MAX = 17.0
NODE_DIAMETER_EXPONENT = 0.5


def _extract_json_array_after(html: str, marker: str) -> list[dict[str, Any]]:
    """Return the JSON array immediately following ``marker``."""
    start = html.find(marker)
    if start < 0:
        raise ValueError(f"Could not find marker {marker!r}.")

    json_start = html.find("[", start + len(marker))
    if json_start < 0:
        raise ValueError(f"Could not find JSON array after marker {marker!r}.")

    decoder = json.JSONDecoder()
    data, _end = decoder.raw_decode(html[json_start:])
    if not isinstance(data, list):
        raise ValueError(f"Data after marker {marker!r} is not a JSON list.")
    return data


def _extract_vis_dataset(html: str, variable: str) -> list[dict[str, Any]]:
    """Return a JSON dataset assigned as ``<variable> = new vis.DataSet(...)``."""
    marker = f"{variable} = new vis.DataSet("
    return _extract_json_array_after(html, marker)


def _extract_graph_data(html: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Extract the full graph when available, otherwise the visible graph."""
    try:
        nodes = _extract_json_array_after(html, "var fullNodes = ")
        edges = _extract_json_array_after(html, "var fullEdges = ")
        if nodes and edges:
            return nodes, edges
    except ValueError:
        pass
    return _extract_vis_dataset(html, "nodes"), _extract_vis_dataset(html, "edges")


def _node_year(node: dict[str, Any]) -> float:
    """Get the node calendar year from explicit metadata or the text label."""
    year = node.get("year")
    if year is not None:
        return float(year)
    label = str(node.get("label", ""))
    match = re.search(r"\b(19|20|21)\d{2}\b", label)
    if not match:
        raise ValueError(f"Node has no year metadata and no year in label: {label!r}")
    return float(match.group(0))


def _band_parts(node: dict[str, Any]) -> tuple[str, str, str]:
    """Return activity name, reference product, and location for a node."""
    band_key = str(node.get("band_key") or "")
    parts = band_key.split("|")
    if len(parts) >= 3:
        return parts[0].strip(), parts[1].strip(), parts[2].strip()

    lines = [line.strip() for line in str(node.get("label", "")).splitlines()]
    name = lines[0] if lines else ""
    ref = lines[1] if len(lines) > 1 else ""
    location = ""
    if lines:
        tail = lines[-1].split("|")
        if tail:
            location = tail[-1].strip()
    return name, ref, location


def _display_label(node: dict[str, Any]) -> str:
    """Build a concise row label similar to the manuscript figure."""
    name, ref, location = _band_parts(node)
    lower_name = name.lower()
    lower_ref = ref.lower()
    if "transport, passenger car" in lower_name:
        vehicle = "ICEV" if "ice" in lower_ref or "icev" in lower_name else "car"
        return f"Passenger car transport\n({vehicle}, {location})"
    return f"{name}\n({location})" if location else name


def _row_key(node: dict[str, Any]) -> tuple[int, str]:
    depth = int(node.get("depth", 0))
    band_key = str(node.get("band_key") or node.get("label") or node.get("id"))
    return depth, band_key


def _build_rows(nodes: list[dict[str, Any]]) -> list[tuple[int, str]]:
    """Order rows by routing depth and the vertical order encoded in the HTML."""
    grouped: dict[tuple[int, str], list[dict[str, Any]]] = defaultdict(list)
    for node in nodes:
        grouped[_row_key(node)].append(node)

    def sort_key(item: tuple[tuple[int, str], list[dict[str, Any]]]) -> tuple[int, float]:
        (depth, _band_key), group_nodes = item
        # In pyvis export, larger y values are visually higher.
        return depth, -mean(float(n.get("y", 0.0)) for n in group_nodes)

    return [key for key, _group in sorted(grouped.items(), key=sort_key)]


def _edge_amount(edge: dict[str, Any]) -> float:
    """Return the exchange amount used to scale flow width."""
    title = str(edge.get("title") or "")
    if "score=" in title:
        match = re.search(r"amount=([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?)", title)
        if match:
            return abs(float(match.group(1)))
    return abs(float(edge.get("value", 0.0) or 0.0))


def _edge_width(edge: dict[str, Any], max_amount: float) -> float:
    value = _edge_amount(edge)
    if max_amount <= 0:
        return EDGE_WIDTH_MIN
    ratio = min(1.0, max(0.0, value / max_amount))
    return EDGE_WIDTH_MIN + (EDGE_WIDTH_MAX - EDGE_WIDTH_MIN) * (
        ratio**EDGE_WIDTH_EXPONENT
    )


def _node_area(weight_coefficient: float) -> float:
    """Return Matplotlib scatter area for a node weighting coefficient."""
    coeff = min(1.0, max(0.0, float(weight_coefficient)))
    diameter = NODE_DIAMETER_MIN + (NODE_DIAMETER_MAX - NODE_DIAMETER_MIN) * (
        coeff**NODE_DIAMETER_EXPONENT
    )
    return diameter**2


def _infer_node_weight_coefficients(
    nodes: list[dict[str, Any]],
    edges: list[dict[str, Any]],
) -> dict[str, float]:
    """Infer per-node visual weighting coefficients from incoming edges.

    The HTML stores routed exchange amounts but not the original temporal
    weights. Estimate local split coefficients from comparable sibling edges
    and propagate the parent node coefficient through single-child exchanges.
    This keeps capital-input nodes from being drawn as full-weight nodes merely
    because they are the only child of their parent.
    """
    node_by_id = {str(node["id"]): node for node in nodes}
    band_groups: dict[tuple[str, tuple[int, str]], list[dict[str, Any]]] = (
        defaultdict(list)
    )
    fuel_year_groups: dict[tuple[str, float], list[dict[str, Any]]] = defaultdict(list)
    product_year_groups: dict[tuple[str, float, str], list[dict[str, Any]]] = (
        defaultdict(list)
    )

    def is_fuel_production(node: dict[str, Any]) -> bool:
        name, ref_product, _location = _band_parts(node)
        return name.strip().lower() in {"gasoline production", "biofuel production"} and (
            ref_product.strip().lower() in {"gasoline", "biofuel"}
        )

    for edge in edges:
        source_id = str(edge.get("from"))
        target_id = str(edge.get("to"))
        target = node_by_id.get(target_id)
        if target is None:
            continue
        _name, ref_product, _location = _band_parts(target)
        band_groups[(source_id, _row_key(target))].append(edge)
        if is_fuel_production(target):
            fuel_year_groups[(source_id, _node_year(target))].append(edge)
        product_year_groups[
            (source_id, _node_year(target), ref_product.strip().lower())
        ].append(edge)

    fuel_mix_groups: set[tuple[str, float]] = set()
    for group_key, group_edges in fuel_year_groups.items():
        target_fuels = {
            _band_parts(node_by_id[str(edge.get("to"))])[1].strip().lower()
            for edge in group_edges
            if str(edge.get("to")) in node_by_id
        }
        if {"gasoline", "biofuel"}.issubset(target_fuels):
            fuel_mix_groups.add(group_key)

    electricity_mix_groups: set[tuple[str, float, str]] = set()
    for group_key, group_edges in product_year_groups.items():
        _source_id, _year, ref_product = group_key
        if ref_product != "electricity":
            continue
        target_rows = {
            _row_key(node_by_id[str(edge.get("to"))])
            for edge in group_edges
            if str(edge.get("to")) in node_by_id
        }
        if len(target_rows) > 1:
            electricity_mix_groups.add(group_key)

    coefficients = {str(node["id"]): 0.0 for node in nodes}
    for node in nodes:
        if int(node.get("depth", 0)) == 0:
            coefficients[str(node["id"])] = 1.0

    def edge_group(edge: dict[str, Any]) -> tuple[Any, ...]:
        source_id = str(edge.get("from"))
        target_id = str(edge.get("to"))
        target = node_by_id.get(target_id)
        if target is None:
            return ("missing", source_id, target_id)
        _name, ref_product, _location = _band_parts(target)
        fuel_year_key = (source_id, _node_year(target))
        if fuel_year_key in fuel_mix_groups:
            return ("fuel_mix", *fuel_year_key)
        product_year_key = (
            source_id,
            _node_year(target),
            ref_product.strip().lower(),
        )
        if product_year_key in electricity_mix_groups:
            return ("electricity_mix", *product_year_key)
        return ("temporal", source_id, _row_key(target))

    group_totals: dict[tuple[Any, ...], float] = defaultdict(float)
    for edge in edges:
        group_totals[edge_group(edge)] += _edge_amount(edge)

    sorted_edges = sorted(
        edges,
        key=lambda edge: int(node_by_id.get(str(edge.get("to")), {}).get("depth", 0)),
    )
    for edge in sorted_edges:
        source_id = str(edge.get("from"))
        target_id = str(edge.get("to"))
        target = node_by_id.get(target_id)
        source = node_by_id.get(source_id)
        if target is None or source is None:
            continue
        group_key = edge_group(edge)
        total = group_totals[group_key]
        if total <= 0.0:
            continue
        local_coefficient = _edge_amount(edge) / total
        source_coefficient = coefficients.get(
            source_id,
            1.0 if int(source.get("depth", 0)) == 0 else 0.0,
        )
        if group_key[0] in {"fuel_mix", "electricity_mix"}:
            coefficient = local_coefficient
        else:
            coefficient = source_coefficient * local_coefficient
        coefficients[target_id] = max(coefficients.get(target_id, 0.0), coefficient)

    return coefficients


def _draw_curved_edge(ax: plt.Axes, x0: float, y0: float, x1: float, y1: float, **kwargs: Any) -> None:
    """Draw a smooth cubic Bezier edge."""
    dx = x1 - x0
    c1 = (x0 + 0.45 * dx, y0)
    c2 = (x1 - 0.45 * dx, y1)
    path = MplPath(
        [(x0, y0), c1, c2, (x1, y1)],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    ax.add_patch(PathPatch(path, fill=False, **kwargs))


def render_graph(
    html_path: Path,
    output_path: Path,
    *,
    width: float = 13.5,
    height: float | None = None,
    dpi: int = 300,
    max_depth: int | None = None,
    figure2_style: bool = False,
) -> None:
    html = html_path.read_text(encoding="utf-8", errors="replace")
    nodes, edges = _extract_graph_data(html)

    if not nodes:
        raise ValueError(f"No nodes found in {html_path}.")
    if figure2_style and max_depth is None:
        max_depth = 3
    if max_depth is not None:
        max_depth = int(max_depth)
        kept_ids = {
            str(node["id"])
            for node in nodes
            if int(node.get("depth", 0)) <= max_depth
        }
        nodes = [node for node in nodes if str(node["id"]) in kept_ids]
        edges = [
            edge
            for edge in edges
            if str(edge.get("from")) in kept_ids and str(edge.get("to")) in kept_ids
        ]

    rows = _build_rows(nodes)
    row_position = {row: -idx for idx, row in enumerate(rows)}
    node_by_id = {str(node["id"]): node for node in nodes}
    years = [_node_year(node) for node in nodes]
    year_min = min(years)
    year_max = max(years)
    year_range = max(year_max - year_min, 1.0)
    left_pad = max(20.0, 0.42 * year_range)
    right_pad = max(3.0, 0.08 * year_range)
    x_min = math.floor(year_min - left_pad)
    x_max = math.ceil(year_max + right_pad)

    if height is None:
        if figure2_style:
            height = max(4.6, 0.39 * len(rows) + 0.8)
        else:
            height = max(4.8, 0.48 * len(rows) + 0.9)

    node_weight_coefficients = _infer_node_weight_coefficients(nodes, edges)

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Depth bands.
    rows_by_depth: dict[int, list[float]] = defaultdict(list)
    for depth, band_key in rows:
        rows_by_depth[depth].append(row_position[(depth, band_key)])

    for depth, ys in rows_by_depth.items():
        ymin = min(ys) - 0.43
        ymax = max(ys) + 0.43
        color = DEPTH_BAND_COLORS[depth % len(DEPTH_BAND_COLORS)]
        ax.axhspan(ymin, ymax, color=color, zorder=0)
        ax.text(
            x_min + 0.25,
            (ymin + ymax) / 2,
            f"d{depth}",
            ha="left",
            va="center",
            fontsize=FONT_SIZE_DEPTH,
            fontweight="bold",
            color="#5f6b7a",
            zorder=5,
        )

    # Vertical year grid.
    tick_step = 5 if year_range <= 60 else 10
    tick_start = int(math.ceil(x_min / tick_step) * tick_step)
    ticks = list(range(tick_start, int(x_max) + 1, tick_step))
    for tick in ticks:
        ax.axvline(tick, color="#d8dee6", linewidth=0.8, zorder=1)

    max_edge_amount = max((_edge_amount(edge) for edge in edges), default=1.0)
    for edge in edges:
        source = node_by_id.get(str(edge.get("from")))
        target = node_by_id.get(str(edge.get("to")))
        if source is None or target is None:
            continue
        x0 = _node_year(source)
        y0 = row_position[_row_key(source)]
        x1 = _node_year(target)
        y1 = row_position[_row_key(target)]
        color = str(edge.get("color") or target.get("color") or "#7f8c8d")
        _draw_curved_edge(
            ax,
            x0,
            y0,
            x1,
            y1,
            linewidth=_edge_width(edge, max_edge_amount),
            edgecolor=color,
            alpha=0.24 if figure2_style else 0.26,
            zorder=2,
        )

    # Activity labels, one per row.
    label_x = x_min + 3.8
    row_examples: dict[tuple[int, str], dict[str, Any]] = {}
    for node in nodes:
        row_examples.setdefault(_row_key(node), node)
    for row in rows:
        node = row_examples[row]
        y = row_position[row]
        color = str(node.get("color") or "#2f3b45")
        ax.text(
            label_x,
            y,
            _display_label(node),
            ha="left",
            va="center",
            fontsize=FONT_SIZE_LABEL,
            color=color,
            linespacing=0.96,
            zorder=5,
        )

    # Nodes.
    for node in nodes:
        x = _node_year(node)
        y = row_position[_row_key(node)]
        depth = int(node.get("depth", 0))
        area = _node_area(node_weight_coefficients.get(str(node["id"]), 1.0))
        color = str(node.get("color") or "#4c78a8")
        if depth == 0:
            ax.scatter(
                [x],
                [y],
                s=area * 1.25,
                marker="D",
                c=["#1f2933"],
                edgecolors="white",
                linewidths=0.9,
                zorder=4,
            )
        else:
            ax.scatter(
                [x],
                [y],
                s=area,
                marker="o",
                c=[color],
                edgecolors="white",
                linewidths=0.8,
                zorder=4,
            )

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(min(row_position.values()) - 0.8, 0.8)
    ax.set_yticks([])
    ax.set_xticks(ticks)
    ax.tick_params(axis="x", labelsize=FONT_SIZE_TICK, colors="#5f6b7a")
    ax.set_xlabel("Calendar year", fontsize=FONT_SIZE_AXIS)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#aab3bd")
    ax.spines["bottom"].set_linewidth(0.8)
    ax.text(
        x_max,
        min(row_position.values()) - 0.45,
        "Diameters: distribution weights (power: annual mix shares); lines: exchange amounts",
        ha="right",
        va="top",
        fontsize=FONT_SIZE_NOTE,
        color="#6b7280",
    )

    fig.tight_layout(pad=0.4)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "html",
        nargs="?",
        default="dev/publication/trails_graph.html",
        type=Path,
        help="Input pyvis/vis-network HTML file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=Path("dev/publication/trails_graph.png"),
        type=Path,
        help="Output PNG path.",
    )
    parser.add_argument("--width", default=13.5, type=float, help="Figure width in inches.")
    parser.add_argument("--height", default=None, type=float, help="Figure height in inches.")
    parser.add_argument("--dpi", default=300, type=int, help="Output resolution.")
    parser.add_argument(
        "--max-depth",
        default=None,
        type=int,
        help="Optional maximum routing depth to render.",
    )
    parser.add_argument(
        "--figure2-style",
        action="store_true",
        help="Use a compact manuscript Figure 2 style; defaults to max depth 3.",
    )
    args = parser.parse_args()

    render_graph(
        args.html,
        args.output,
        width=args.width,
        height=args.height,
        dpi=args.dpi,
        max_depth=args.max_depth,
        figure2_style=args.figure2_style,
    )
    print(args.output)


if __name__ == "__main__":
    main()
