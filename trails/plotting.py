from typing import Dict, Any, List, Tuple, Optional
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


def plot_temporal_scores(
    results_by_year: Dict[int, Dict[str, Any]],
    trails: Trails,
    title: str = "Temporal impacts by responsible activity",
    method_label: str = "Impact score",
    cumulative: bool = False,
    stacked: bool = True,
):
    """
    Plot LCIA scores over time, split by responsible first-level activity.

    Features:
    - stacked or unstacked mode
    - filled areas even when unstacked (semi-transparent)
    - optional cumulative view
    """
    import numpy as np
    import plotly.graph_objects as go

    # --- 1. Collect years and root activities ---
    years = sorted(results_by_year.keys())
    all_roots = sorted({
        root
        for year in years
        for root in results_by_year[year].get("scores_per_root", {})
    })

    if not all_roots:
        raise ValueError("No scores_per_root found.")

    # --- 2. Build score matrix (year × root) ---
    Y = np.zeros((len(years), len(all_roots)), dtype=float)
    for yi, year in enumerate(years):
        spr = results_by_year[year].get("scores_per_root", {})
        for ri, root in enumerate(all_roots):
            Y[yi, ri] = spr.get(root, 0.0)

    # --- 3. Cumulative transform if requested ---
    if cumulative:
        Y = np.cumsum(Y, axis=0)

    # --- 4. Lookup labels ---
    idx_to_label = _build_activity_label_map(trails)
    label_for_root = lambda idx: idx_to_label.get(idx, f"Activity {idx}")

    # --- 5. Create figure ---
    fig = go.Figure()

    # Precompute transparency for unstacked mode
    alpha = 0.4 if not stacked else 1.0

    for ri, root in enumerate(all_roots):
        label = label_for_root(root)
        y_vals = Y[:, ri]

        trace_kwargs = dict(
            x=years,
            y=y_vals,
            name=label,
            mode="lines",
            hovertemplate=(
                "<b>%{fullData.name}</b><br>"
                "Year: %{x}<br>"
                f"{method_label}: %{{y:.6g}}<extra></extra>"
            ),
        )

        if stacked:
            # Normal stacked area
            trace_kwargs.update(stackgroup="one")
        else:
            # Filled but NOT stacked: semi-transparent areas
            trace_kwargs.update(
                fill="tozeroy",
                line=dict(width=2),
                opacity=alpha,
            )

        fig.add_trace(go.Scatter(**trace_kwargs))

    # --- 6. Total overlay (only for stacked mode) ---
    if stacked:
        total_scores = Y.sum(axis=1)
        fig.add_trace(
            go.Scatter(
                x=years,
                y=total_scores,
                name="Total",
                mode="lines",
                line=dict(width=2, dash="dot"),
                hovertemplate=(
                    "<b>Total</b><br>"
                    "Year: %{x}<br>"
                    f"{method_label}: %{{y:.6g}}<extra></extra>"
                ),
                showlegend=True,
            )
        )

    # --- 7. Layout ---
    fig.update_layout(
        title=title,
        xaxis_title="Year",
        yaxis_title=method_label,
        template="plotly_white",
        hovermode="x unified",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            font=dict(size=10),
        ),
        margin=dict(l=60, r=20, t=60, b=40),
    )

    fig.update_xaxes(dtick=1, showgrid=True)
    fig.update_yaxes(showgrid=True)

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
