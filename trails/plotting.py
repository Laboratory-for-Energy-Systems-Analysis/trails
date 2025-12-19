from typing import Dict, Any, List, Tuple, Optional, Literal

import math

from collections import defaultdict
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

from .trails import Trails
from .lca import compute_node_impact_intensities
from .utils import _format_path_label, _format_path_label_with_years


def _build_activity_label_map(trails: Trails):
    labels = {}
    for scen_label, mapping in trails.activity_indices.items():
        for idx, meta in mapping.items():
            if idx not in labels:
                name = meta.get("name") or f"Activity {idx}"
                rp = meta.get("reference product")
                loc = meta.get("location")

                label = name
                if rp:
                    label += f" | {rp}"
                if loc:
                    label += f" ({loc})"
                labels[idx] = label
    return labels


def _activity_label_map(trails: Trails) -> dict[int, str]:
    labels = {}
    for scen_label, mapping in trails.activity_indices.items():
        for idx, meta in mapping.items():
            if idx not in labels:
                name = meta.get("name") or f"Activity {idx}"
                rp = meta.get("reference product")
                loc = meta.get("location")
                label = name
                if rp:
                    label += f" | {rp}"
                if loc:
                    label += f" ({loc})"
                labels[idx] = label
    return labels


def to_impact_year_results(results: Dict[int, Dict[str, Any]] | Dict[str, Any]) -> Dict[int, Dict[str, Any]]:
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
        if isinstance(k0, int) and isinstance(v0, dict) and ("scores_per_root" in v0 or "scores" in v0):
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
        tspr = entry.get("temporal_scores_per_root_by_year", {})  # impact_year -> {root: score}
        tst  = entry.get("temporal_scores_by_year", {})           # impact_year -> total score (optional)

        for impact_year, per_root in tspr.items():
            impact_year = int(impact_year)
            out.setdefault(impact_year, {"scores": 0.0, "scores_per_root": {}})

            # per-root sums
            for root, val in (per_root or {}).items():
                out[impact_year]["scores_per_root"][int(root)] = (
                    out[impact_year]["scores_per_root"].get(int(root), 0.0) + float(val)
                )

        # total scores (if present)
        for impact_year, total in (tst or {}).items():
            impact_year = int(impact_year)
            out.setdefault(impact_year, {"scores": 0.0, "scores_per_root": {}})
            out[impact_year]["scores"] += float(total)

    return out

def plot_temporal_scores(
    results_by_year: Dict[int, Dict[str, Any]],
    trails,
    title: str = "Temporal impacts by responsible activity",
    method_label: str = "Impact score",
    cumulative: bool = False,
    stacked: bool = True,
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
    static_score: Optional[float] = None,
    static_score_label: str = "Static score",
    static_score_dash: str = "dash",
    static_score_color: str = "black",
    y_max: Optional[float] = None,
    y2_max: Optional[float] = None,
):
    results_by_year = to_impact_year_results(results_by_year)

    if year_tick < 1:
        raise ValueError("year_tick must be >= 1")

    # ------------------------------------------------------------------
    # 1. Years selection
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # 2. Roots and score matrix
    # ------------------------------------------------------------------
    all_roots = sorted({
        root
        for year in years
        for root in results_by_year[year].get("scores_per_root", {})
    })
    if not all_roots:
        raise ValueError("No scores_per_root found.")

    Y = np.zeros((len(years), len(all_roots)), dtype=float)
    for yi, year in enumerate(years):
        spr = results_by_year[year].get("scores_per_root", {})
        for ri, root in enumerate(all_roots):
            Y[yi, ri] = spr.get(root, 0.0)

    if cumulative:
        Y = np.cumsum(Y, axis=0)

    total_raw = Y.sum(axis=1)

    if yaxis_type == "log":
        Y = np.where(Y > 0, Y, log_eps)
        if static_score is not None:
            static_score = max(static_score, log_eps)

    # ------------------------------------------------------------------
    # 3. Labels
    # ------------------------------------------------------------------
    idx_to_label = _build_activity_label_map(trails)

    def label_for_root(idx: int) -> str:
        return idx_to_label.get(idx, f"Activity {idx}")

    # ------------------------------------------------------------------
    # 4. Figure + main traces
    # ------------------------------------------------------------------
    fig = go.Figure()
    alpha = 0.4 if not stacked else 1.0

    for ri, root in enumerate(all_roots):
        fig.add_trace(
            go.Scatter(
                x=years,
                y=Y[:, ri],
                name=label_for_root(root),
                mode="lines",
                hovertemplate=(
                    "<b>%{fullData.name}</b><br>"
                    "Year: %{x}<br>"
                    f"{method_label}: %{{y:.6g}}<extra></extra>"
                ),
                **({"stackgroup": "one"} if stacked else {
                    "fill": "tozeroy",
                    "line": dict(width=2),
                    "opacity": alpha,
                }),
            )
        )

    # ------------------------------------------------------------------
    # 5. Optional secondary cumulative axis
    # ------------------------------------------------------------------
    cum_vals = None
    if show_cumulative_axis:
        cum_vals = np.cumsum(total_raw)
        if yaxis_type == "log":
            cum_vals = np.where(cum_vals > 0, cum_vals, log_eps)

        fig.add_trace(
            go.Scatter(
                x=years,
                y=cum_vals,
                name="Cumulative total",
                showlegend=True,  # <-- force legend handle
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

    # ------------------------------------------------------------------
    # 5b. NEW: static score horizontal line on y2
    # ------------------------------------------------------------------
    if static_score is not None:
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

    # ------------------------------------------------------------------
    # 6. Layout
    # ------------------------------------------------------------------
    fig_w = int(width) if (width is not None) else 800
    entry_w = max(80, int(legend_entrywidth))
    n_items = len(all_roots)
    if show_cumulative_axis and show_cumulative_in_legend:
        n_items += 1
    if static_score is not None:
        n_items += 1

    n_cols = max(1, fig_w // entry_w)
    n_rows = max(1, int(math.ceil(n_items / n_cols)))
    top_margin = 55 + n_rows * int(legend_row_height) + 10

    fig.update_layout(
        width=width,
        height=height,
        template="plotly_white",
        hovermode="x unified",
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

    # Ensure y2 exists if we use it (cumulative and/or static score)
    use_y2 = bool(show_cumulative_axis) or (static_score is not None)
    if use_y2:
        fig.update_layout(
            yaxis2=dict(
                title=cumulative_axis_label,
                overlaying="y",
                side="right",
                showgrid=False,
            )
        )

    # ------------------------------------------------------------------
    # 6b. Optional fixed y maxima + zero alignment between y and y2 (linear)
    # ------------------------------------------------------------------
    if (show_cumulative_axis or (static_score is not None)) and yaxis_type == "linear":
        # -------- y1 data-derived min/max (include 0)
        y1_min = float(np.nanmin(Y))
        y1_max_data = float(np.nanmax(Y))
        y1_min = min(y1_min, 0.0)
        y1_max_data = max(y1_max_data, 0.0)
        if y1_max_data == y1_min:
            y1_max_data = y1_min + 1.0

        # Apply user-specified y_max if provided
        y1_max = float(y_max) if (y_max is not None) else y1_max_data
        y1_max = max(y1_max, 0.0)
        if y1_max == y1_min:
            y1_max = y1_min + 1.0

        # -------- collect y2 data
        y2_series = []
        if show_cumulative_axis and (cum_vals is not None) and len(cum_vals) > 0:
            y2_series.append(np.asarray(cum_vals, dtype=float))
        if static_score is not None:
            y2_series.append(np.asarray([float(static_score)], dtype=float))

        if y2_series:
            y2_all = np.concatenate(y2_series)

            y2_min_data = float(np.nanmin(y2_all))
            y2_max_data = float(np.nanmax(y2_all))
            y2_min_data = min(y2_min_data, 0.0)
            y2_max_data = max(y2_max_data, 0.0)
            if y2_max_data == y2_min_data:
                y2_max_data = y2_min_data + 1.0

            # Apply user-specified y2_max if provided
            y2_max_eff = float(y2_max) if (y2_max is not None) else y2_max_data
            y2_max_eff = max(y2_max_eff, 0.0)
            if y2_max_eff == y2_min_data:
                y2_max_eff = y2_min_data + 1.0

            # Relative position of 0 on y1
            p = (0.0 - y1_min) / (y1_max - y1_min)

            # Headroom guard (apply to y2 when auto-scaling; if user fixes y2_max, we respect it)
            hr = max(0.0, float(y2_headroom))
            if y2_max is None:
                y2_max_eff = y2_max_eff * (1.0 + hr)

            # Compute y2 min so that 0 aligns at fraction p, while respecting fixed y2_max if given
            if p <= 0.0:
                # 0 at/below bottom on y1 -> keep y2 starting at 0
                y2_range = [0.0, y2_max_eff]
            elif p >= 1.0:
                # 0 at/above top on y1
                y2_range = [y2_min_data, 0.0]
            else:
                # Require y2_min such that (0 - y2_min)/(y2_max - y2_min) = p
                # => y2_min = - (p/(1-p)) * y2_max
                y2_min_eff = - (p / (1.0 - p)) * y2_max_eff

                # Ensure we don't clip actual y2 data below y2_min_eff unless user forced y2_max and data demands more span.
                # If data has more negative values than y2_min_eff, extend downward to fit.
                y2_min_eff = min(y2_min_eff, y2_min_data)

                y2_range = [y2_min_eff, y2_max_eff]

            fig.update_layout(
                yaxis=dict(range=[y1_min, y1_max]),
                yaxis2=dict(range=y2_range),
            )

    # If user provided y_max in non-linear mode, apply it directly.
    if yaxis_type != "linear" and y_max is not None:
        fig.update_layout(yaxis=dict(range=[None, float(y_max)]))

    # If user provided y2_max in non-linear mode and y2 exists, apply it directly.
    if yaxis_type != "linear" and y2_max is not None:
        use_y2 = bool(show_cumulative_axis) or (static_score is not None)
        if use_y2:
            fig.update_layout(yaxis2=dict(range=[None, float(y2_max)]))

    # ------------------------------------------------------------------
    # 7. X-axis
    # ------------------------------------------------------------------
    fig.update_xaxes(
        dtick=year_tick,
        tickmode="linear",
        showgrid=show_year_grid,
        tick0=(year_range[0] if year_range else years[0]),
        range=list(year_range) if year_range else None,
    )

    # ------------------------------------------------------------------
    # 8. Reference year
    # ------------------------------------------------------------------
    if reference_year is not None:
        fig.add_vline(
            x=reference_year,
            line_width=2,
            line_dash="dash",
            annotation_text="Reference year",
            annotation_position="top",
        )

    return fig

def plot_top_paths_for_year(
    provenance: Dict[tuple[int, int], Dict[tuple[int, ...], float]],
    trails: Trails,
    year: int,
    top_n: int = 10,
    title: str = "Top demand paths by amount",
    amount_label: str = "Demand (amount units)",
):
    """
    Visualize the top-N *paths* (chains of suppliers) contributing to demand
    in a given year, based on path-wise provenance.

    Parameters
    ----------
    provenance : dict[(year, act_idx) -> dict[path_tuple -> amount]]
        Returned by `lca(..., return_provenance=True)` via temporal_traversal.
    trails : Trails
        To resolve activity indices to labels.
    year : int
        Scenario year to visualize.
    top_n : int
        Number of paths to show.
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
                    "<b>%{y}</b><br>"
                    f"{amount_label}: %{{x:.6g}}<extra></extra>"
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


def plot_temporal_sankey(
    provenance,
    trails: Trails,
    start_year: int,
    start_act_idx: int,
    methods: List[str],
    top_n_paths: int | None = 30,
    title: str = "Temporal Sankey (impact-weighted, aggregated by year)",
    amount_label: str = "Impact score",
    fig_width: int = 1200,
    fig_height: int = 800,
    node_thickness: int = 20,
    node_pad: int = 15,
    font_size: int = 11,
    remove_uncertainty: bool = True,
    debug: bool = False,
):
    """
    Impact-weighted temporal Sankey, simplified:

    - Nodes are aggregated by (depth, year), NOT by activity.
      At each depth level, there is at most one node per year.
    - Links are impact-weighted sums of all activities between these
      (depth, year) nodes, but we preserve activity-level contributions
      in the hover text.

    Layout:
      - x-axis = depth (0 = root, 1 = first level, etc.), left→right.
      - y-axis = year (earlier years towards top, later towards bottom),
        aligned across depths (same year → same vertical position).

    Parameters
    ----------
    provenance : dict[(int, int) -> dict[path_tuple -> amount]]
        Returned by lca(..., return_provenance=True).
        path_tuple is a tuple of (year, act_idx) pairs from first-level onward.
    trails : Trails
        Trails wrapper.
    start_year : int
        Root demand year.
    start_act_idx : int
        Root activity index.
    methods : list[str]
        Single LCIA method to use for impact intensities.
    top_n_paths : int | None
        If None, use all paths. If int, keep only the top-N by |amount|.
    """

    # --- helper: minimal activity meta for labels and hover ---
    def _collect_activity_meta(trails_local: Trails) -> Dict[int, Dict[str, Any]]:
        meta_by_idx: Dict[int, Dict[str, Any]] = {}
        for scen_label, mapping in trails_local.activity_indices.items():
            for idx, meta in mapping.items():
                if idx not in meta_by_idx:
                    meta_by_idx[idx] = meta
        return meta_by_idx

    def _activity_label(act_idx: int) -> str:
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

    # ------------------------------------------------------------------
    # 1. Aggregate amounts per full path including root
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # 2. Top-N pruning (on amount, not impact)
    # ------------------------------------------------------------------
    if top_n_paths is None:
        selected_paths = list(full_path_amounts.items())
    else:
        selected_paths = sorted(
            full_path_amounts.items(),
            key=lambda kv: abs(kv[1]),
            reverse=True,
        )[:top_n_paths]

    # ------------------------------------------------------------------
    # 3. Original node set & depth map (per (year, act))
    # ------------------------------------------------------------------
    node_keys: set[Tuple[int, int]] = set()
    depth_map: Dict[Tuple[int, int], int] = {}

    for full_path, _ in selected_paths:
        for depth, node in enumerate(full_path):
            node_keys.add(node)
            if node not in depth_map or depth < depth_map[node]:
                depth_map[node] = depth

    node_keys = sorted(node_keys)

    # ------------------------------------------------------------------
    # 4. Impact intensities for all (year, act) nodes in selected paths
    # ------------------------------------------------------------------
    nodes_for_intensity = list(node_keys)
    node_intensity = compute_node_impact_intensities(
        trails=trails,
        nodes=nodes_for_intensity,
        methods=methods,
        remove_uncertainty=remove_uncertainty,
        debug=debug,
    )

    # activity metadata for hover labels
    act_meta = _collect_activity_meta(trails)

    # ------------------------------------------------------------------
    # 5. Aggregate link impacts between (depth, year) nodes,
    #    but keep per-activity contributions for hover.
    # ------------------------------------------------------------------
    # aggregated nodes are (depth, year)
    link_impact_agg: Dict[Tuple[Tuple[int, int], Tuple[int, int]], float] = defaultdict(float)
    # per-activity contributions on each aggregated edge:
    #   edge_activity_contrib[(src_agg, tgt_agg)][child_act_idx] = impact
    edge_activity_contrib: Dict[Tuple[Tuple[int, int], Tuple[int, int]], Dict[int, float]] = defaultdict(
        lambda: defaultdict(float)
    )

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

            # impact along this edge from this child activity
            intensity = node_intensity.get(child, 0.0)
            imp = float(amt) * float(intensity)

            src_agg = (parent_depth, parent_year)
            tgt_agg = (child_depth, child_year)

            link_impact_agg[(src_agg, tgt_agg)] += imp
            edge_activity_contrib[(src_agg, tgt_agg)][child_act] += imp

    if not link_impact_agg:
        raise ValueError("No aggregated link impacts; nothing to plot.")

    # Collect aggregated node set
    agg_nodes: set[Tuple[int, int]] = set()
    for src_agg, tgt_agg in link_impact_agg.keys():
        agg_nodes.add(src_agg)
        agg_nodes.add(tgt_agg)

    agg_nodes = sorted(agg_nodes)  # sort by (depth, year)
    node_index_agg: Dict[Tuple[int, int], int] = {
        key: i for i, key in enumerate(agg_nodes)
    }

    # ------------------------------------------------------------------
    # 6. Compute total impact per aggregated node (for labels)
    # ------------------------------------------------------------------
    node_total_impact: Dict[Tuple[int, int], float] = defaultdict(float)
    for (src_agg, tgt_agg), imp in link_impact_agg.items():
        node_total_impact[tgt_agg] += imp

    # ------------------------------------------------------------------
    # 7. Layout: x = depth, y = year
    # ------------------------------------------------------------------
    depths = sorted({d for (d, y) in agg_nodes})
    years = sorted({y for (d, y) in agg_nodes})

    if len(depths) == 1:
        depth_to_x = {depths[0]: 0.5}
    else:
        depth_to_x = {
            d: 0.05 + 0.9 * (i / (len(depths) - 1))
            for i, d in enumerate(depths)
        }

    if len(years) == 1:
        year_to_y = {years[0]: 0.5}
    else:
        year_to_y = {
            y: 0.05 + 0.9 * (i / (len(years) - 1))
            for i, y in enumerate(years)
        }

    node_x = [depth_to_x[d] for (d, y) in agg_nodes]
    node_y = [year_to_y[y] for (d, y) in agg_nodes]

    # ------------------------------------------------------------------
    # 8. Node labels and colors (by year)
    # ------------------------------------------------------------------
    node_labels: List[str] = []
    for (d, y) in agg_nodes:
        total_imp = node_total_impact.get((d, y), 0.0)
        node_labels.append(
            f"Depth {d}, Year {y}<br>"
            f"Incoming {amount_label}: {total_imp:.3g}"
        )

    unique_years = years
    year_palette = px.colors.sequential.Viridis
    if len(unique_years) > len(year_palette):
        repeats = (len(unique_years) // len(year_palette)) + 1
        full_year_palette = (year_palette * repeats)[: len(unique_years)]
    else:
        full_year_palette = year_palette[: len(unique_years)]

    year_to_color = {
        y: col for y, col in zip(unique_years, full_year_palette)
    }
    node_colors = [year_to_color[y] for (d, y) in agg_nodes]

    # ------------------------------------------------------------------
    # 9. Build link arrays and activity-rich hover text
    # ------------------------------------------------------------------
    link_sources: List[int] = []
    link_targets: List[int] = []
    link_values: List[float] = []
    link_colors: List[str] = []
    link_hovertemplates: List[str] = []

    for (src_agg, tgt_agg), imp in link_impact_agg.items():
        src_idx = node_index_agg[src_agg]
        tgt_idx = node_index_agg[tgt_agg]

        (src_depth, src_year) = src_agg
        (tgt_depth, tgt_year) = tgt_agg

        # color link by child (target) year
        color = year_to_color.get(tgt_year, "rgba(150,150,150,0.7)")

        # build hover text with per-activity breakdown (top contributors)
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
                lines.append(f"- {_activity_label(act_idx)}: {act_imp:.6g}")

        hovertemplate = "<br>".join(lines) + "<extra></extra>"

        link_sources.append(src_idx)
        link_targets.append(tgt_idx)
        link_values.append(abs(imp))
        link_colors.append(color)
        link_hovertemplates.append(hovertemplate)

    # ------------------------------------------------------------------
    # 10. Sankey figure
    # ------------------------------------------------------------------
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

def plot_traversal_grid_flow(
    edges_by_depth: dict[int, dict[tuple[tuple[int, int], tuple[int, int]], float]],
    trails: Trails,
    depths: Optional[List[int]] = None,
    include_all_activities: bool = False,
    title: str = "Traversal flow grid (consumers & suppliers)",
    dot_size: int = 10,
    base_width: int = 550,
    base_height: int = 260,
):
    """
    Multi-panel plot of traversal flows on an activity×year grid.

    For each depth d in `depths` AND for "All depths":
      - rows  = activity indices
      - cols  = years (global min..max across all depths)
      - red   = consumers (nodes with outgoing edges at that depth)
      - green = suppliers (nodes appearing as children)
      - arrows from supplier -> consumer

    Layout:
      - 2 columns, as many rows as needed.
      - Last subplot is "All depths" (cumulative edges).
      - Y-axis shows activity indices; a legend below maps index → label.
    """
    # ------------------------------------------------------------------
    # 1. Determine depths to plot
    # ------------------------------------------------------------------
    if depths is None:
        depths_list = sorted(edges_by_depth.keys())
    else:
        depths_list = [d for d in depths if d in edges_by_depth]

    if not depths_list:
        raise ValueError("No depths to plot (depth list empty or no edges).")

    # We'll add one extra panel for "All depths"
    panel_labels: List[str] = [f"Depth {d}" for d in depths_list] + ["All depths"]

    # ------------------------------------------------------------------
    # 2. Activities to show (rows)
    # ------------------------------------------------------------------
    if include_all_activities:
        all_activities = set()
        for scen_label, mapping in trails.activity_indices.items():
            all_activities.update(mapping.keys())
    else:
        all_activities = set()
        for d in depths_list:
            for (parent, child), amt in edges_by_depth.get(d, {}).items():
                (y_cons, a_cons) = parent
                (y_sup, a_sup) = child
                all_activities.add(int(a_cons))
                all_activities.add(int(a_sup))

    acts = sorted(all_activities)
    if not acts:
        raise ValueError("No activities to plot for the selected depths.")

    idx_to_label = _activity_label_map(trails)

    act_to_row = {act: i for i, act in enumerate(acts)}
    y_tickvals = list(range(len(acts)))
    y_ticktext = [str(a) for a in acts]  # show indices only

    # ------------------------------------------------------------------
    # 3. Global time axis: min/max year across ALL selected depths
    # ------------------------------------------------------------------
    years_global_set = set()
    for d in depths_list:
        for (parent, child), amt in edges_by_depth.get(d, {}).items():
            (y_cons, a_cons) = parent
            (y_sup, a_sup) = child
            years_global_set.add(int(y_cons))
            years_global_set.add(int(y_sup))

    if not years_global_set:
        raise ValueError("No years found in edges for selected depths.")

    year_min = min(years_global_set)
    year_max = max(years_global_set)
    years_global = list(range(year_min, year_max + 1))

    # ------------------------------------------------------------------
    # 4. Prepare subplot grid: 2 columns, as many rows as needed
    # ------------------------------------------------------------------
    n_panels = len(panel_labels)  # depths + "All depths"
    ncols = 2
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

    # ------------------------------------------------------------------
    # 4b. Force every subplot to show ALL activity rows (even if empty)
    # ------------------------------------------------------------------
    n_acts = len(acts)

    # y positions are row indices 0..n_acts-1
    all_row_idx = list(range(n_acts))

    for i in range(n_panels):
        r = i // ncols + 1
        c = i % ncols + 1

        fig.update_yaxes(
            row=r,
            col=c,
            tickmode="array",
            tickvals=all_row_idx,
            ticktext=[str(a) for a in acts],  # indices only (as you want)
            range=[-0.5, n_acts - 0.5],  # prevents autoscale cropping
            showgrid=True,  # show horizontal lines
            tick0=0,
            dtick=1,
            automargin=True,
            zeroline=False,
        )

        fig.update_xaxes(
            row=r,
            col=c,
            range=[year_min - 0.5, year_max + 0.5],  # consistent across panels
            tickangle=-90,
            showgrid=True,
            zeroline=False,
        )

    # Pre-build merged edges for "All depths"
    merged_edges_all: dict[tuple[tuple[int, int], tuple[int, int]], float] = defaultdict(float)
    for d in depths_list:
        for key, amt in edges_by_depth.get(d, {}).items():
            merged_edges_all[key] += float(amt)

    # ------------------------------------------------------------------
    # 5. Add traces & arrows per panel
    # ------------------------------------------------------------------
    panels_depths = depths_list + ["all"]  # last panel = cumulative
    for i, depth in enumerate(panels_depths):
        row = i // ncols + 1
        col = i % ncols + 1
        axis_id = i + 1  # x1,y1 for first subplot, etc.

        if depth == "all":
            edges = merged_edges_all
        else:
            edges = edges_by_depth.get(depth, {})

        if not edges:
            continue

        consumer_nodes = set()
        supplier_nodes = set()

        for (parent, child), amt in edges.items():
            (y_cons, a_cons) = parent
            (y_sup, a_sup) = child

            consumer_nodes.add((int(y_cons), int(a_cons)))
            supplier_nodes.add((int(y_sup), int(a_sup)))

        # Build scatter data
        cons_x, cons_y, cons_text = [], [], []
        sup_x, sup_y, sup_text = [], [], []

        for (year, act) in sorted(consumer_nodes):
            cons_x.append(year)
            cons_y.append(act_to_row.get(act, -1))
            cons_text.append(idx_to_label.get(act, f"Activity {act}"))

        for (year, act) in sorted(supplier_nodes):
            sup_x.append(year)
            sup_y.append(act_to_row.get(act, -1))
            sup_text.append(idx_to_label.get(act, f"Activity {act}"))

        # Consumers: red
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
                showlegend=(i == 0),
            ),
            row=row,
            col=col,
        )

        # Suppliers: green
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
                showlegend=(i == 0),
            ),
            row=row,
            col=col,
        )

        # Edges as batched segments (FAST)
        edge_x = []
        edge_y = []
        head_x = []
        head_y = []

        for (parent, child), amt in edges.items():
            (y_cons, a_cons) = parent
            (y_sup, a_sup) = child

            x0 = int(y_sup)
            y0 = act_to_row.get(int(a_sup), -1)
            x1 = int(y_cons)
            y1 = act_to_row.get(int(a_cons), -1)

            if y0 < 0 or y1 < 0:
                continue

            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

            # direction cue at the target (consumer)
            head_x.append(x1)
            head_y.append(y1)

        # Use WebGL for speed
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

        # Arrowheads as markers at consumer end
        fig.add_trace(
            go.Scattergl(
                x=head_x,
                y=head_y,
                mode="markers",
                marker=dict(
                    symbol="triangle-left",  # see note below
                    size=6,
                    color="rgba(100,100,100,0.7)",
                ),
                hoverinfo="skip",
                showlegend=False,
            ),
            row=row,
            col=col,
        )

    # ------------------------------------------------------------------
    # 6. Activity legend below: index -> full label
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # 7. Layout
    # ------------------------------------------------------------------
    fig.update_layout(
        title=title,
        width=base_width * 2,
        height=base_height * nrows + 160,  # extra space for bottom legend
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

    return fig
