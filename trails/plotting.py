from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

import math

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
    """Build a mapping of biosphere flow indices to display labels.

    :param trails: Trails instance with biosphere metadata.
    :type trails: Trails
    :returns: Mapping from flow index to label.
    :rtype: dict[int, str]
    """
    labels = {}
    for scen_label, mapping in trails.biosphere_indices.items():
        for idx, meta in mapping.items():
            if idx not in labels:
                name = meta.get("name") or f"Flow {idx}"
                compartment = meta.get("compartment") or ""
                subcompartment = meta.get("subcompartment") or ""
                unit = meta.get("unit") or ""

                label = name
                if compartment or subcompartment:
                    parts = [p for p in (compartment, subcompartment) if p]
                    label += " | " + "/".join(parts)
                if unit:
                    label += f" ({unit})"
                labels[idx] = label
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
    all_roots: list[int],
    idx_to_label: dict[int, str],
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

    def label_for_root(idx: int) -> str:
        """Resolve a root index to its display label.

        :param idx: Activity index to label.
        :type idx: int
        :returns: Display label for the root.
        :rtype: str
        """
        return idx_to_label.get(idx, f"Activity {idx}")

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
            )
        )


def _apply_linear_yaxis_alignment(
    fig: go.Figure,
    Y: np.ndarray,
    cum_vals: np.ndarray | None,
    static_score: float | None,
    y_max: float | None,
    y2_max: float | None,
    y2_headroom: float,
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
    :param y_max: Optional max for primary y-axis.
    :type y_max: float | None
    :param y2_max: Optional max for secondary y-axis.
    :type y2_max: float | None
    :param y2_headroom: Headroom multiplier for secondary axis.
    :type y2_headroom: float
    """
    y1_min = float(np.nanmin(Y))
    y1_max_data = float(np.nanmax(Y))
    y1_min = min(y1_min, 0.0)
    y1_max_data = max(y1_max_data, 0.0)
    if y1_max_data == y1_min:
        y1_max_data = y1_min + 1.0

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
        if "flow" in scores.dims:
            attrib_dim = "flow"
        else:
            # If user asked by_flow but there is no flow dim, fall back to totals
            attrib_dim = None
        score_key = "scores_by_flow"
    else:
        attrib_dim = "root activity" if "root activity" in scores.dims else "activity"
        score_key = "scores_by_first_level_child"

    # Start from input
    data = scores

    # Collapse to (attrib_dim, year) if possible
    if attrib_dim is None:
        # totals only: sum everything except year
        extra_dims = [d for d in data.dims if d != "year"]
        if extra_dims:
            data = data.sum(dim=extra_dims)
        # ensure ("year",) ordering
        if data.dims != ("year",):
            data = data.transpose("year")
        # build results
        results: Dict[int, Dict[str, Any]] = {}
        arr = data.data
        if hasattr(arr, "todense"):
            v = np.asarray(arr.todense(), dtype=float).ravel()
        else:
            v = np.asarray(data.values, dtype=float).ravel()
        for yi, year in enumerate(years):
            results[int(year)] = {"scores": float(v[yi]), score_key: {}}
        return results

    # If we have exchange-attributed tensor activity x year x root activity,
    # and we're plotting by root activity, reduce over activity.
    if attrib_dim == "root activity" and "activity" in data.dims:
        data = data.sum(dim="activity")

    # If we are plotting by activity but root activity exists, reduce it away.
    if attrib_dim == "activity" and "root activity" in data.dims:
        data = data.sum(dim="root activity")

    # If any leftover dims exist, collapse them too (e.g. "method", "scenario")
    extra_dims = [d for d in data.dims if d not in (attrib_dim, "year")]
    if extra_dims:
        data = data.sum(dim=extra_dims)

    # Ensure (attrib_dim, year) order
    if data.dims != (attrib_dim, "year"):
        data = data.transpose(attrib_dim, "year")

    # Safely densify (compatible with sparse.AUTO_DENSIFY=False)
    arr = data.data
    if hasattr(arr, "todense"):
        vals = np.asarray(arr.todense(), dtype=float)
    else:
        vals = np.asarray(data.values, dtype=float)

    attrib_ids = data.coords[attrib_dim].values

    results: Dict[int, Dict[str, Any]] = {
        int(y): {"scores": 0.0, score_key: {}} for y in years
    }

    # vals shape: (attrib, year)
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


def plot_temporal_scores(
    results_by_year: Union[
        Dict[int, Dict[str, Any]], Dict[str, Any], xr.DataArray, None
    ],
    trails: Trails,
    title: str = "Temporal impacts by responsible activity",
    method_label: str = "Impact score",
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
    static_score: Optional[float] = None,
    static_score_label: str = "Static score",
    static_score_dash: str = "dash",
    static_score_color: str = "black",
    y_max: Optional[float] = None,
    y2_max: Optional[float] = None,
) -> go.Figure:
    """Plot temporal impact scores by responsible activity.

    :param results_by_year: Impact-year results mapping or characterized inventory.
    :type results_by_year: dict[int, dict[str, Any]] | xarray.DataArray | None
    :param trails: Trails instance used for labels.
    :type trails: Trails
    :param title: Plot title.
    :type title: str
    :param method_label: Label for the impact score axis.
    :type method_label: str
    :param cumulative: Whether to accumulate scores over time.
    :type cumulative: bool
    :param stacked: Whether to stack root traces.
    :type stacked: bool
    :param legend_top_n: Number of top-contributing items to show in legend.
    :type legend_top_n: int
    :param show_flow_contributions: Whether to plot flow contributions instead of activities.
    :type show_flow_contributions: bool
    :param width: Figure width in pixels.
    :type width: int | None
    :param height: Figure height in pixels.
    :type height: int | None
    :param year_tick: Tick step for the year axis.
    :type year_tick: int
    :param year_range: Optional ``(start, end)`` year bounds.
    :type year_range: tuple[int, int] | None
    :param show_year_grid: Whether to show grid lines on the year axis.
    :type show_year_grid: bool
    :param yaxis_type: Axis scale type.
    :type yaxis_type: Literal["linear", "log"]
    :param log_eps: Epsilon for log scaling.
    :type log_eps: float
    :param reference_year: Optional reference year line.
    :type reference_year: int | None
    :param show_cumulative_axis: Whether to show cumulative axis.
    :type show_cumulative_axis: bool
    :param cumulative_axis_label: Label for cumulative axis.
    :type cumulative_axis_label: str
    :param legend_entrywidth: Legend entry width in pixels.
    :type legend_entrywidth: int
    :param legend_row_height: Legend row height in pixels.
    :type legend_row_height: int
    :param legend_y: Legend y position.
    :type legend_y: float
    :param y2_headroom: Headroom multiplier for secondary axis.
    :type y2_headroom: float
    :param show_cumulative_in_legend: Whether cumulative trace is shown in legend.
    :type show_cumulative_in_legend: bool
    :param static_score: Optional static score value.
    :type static_score: float | None
    :param static_score_label: Label for static score trace.
    :type static_score_label: str
    :param static_score_dash: Dash style for static score trace.
    :type static_score_dash: str
    :param static_score_color: Color for static score trace.
    :type static_score_color: str
    :param y_max: Optional max for primary y-axis.
    :type y_max: float | None
    :param y2_max: Optional max for secondary y-axis.
    :type y2_max: float | None
    :returns: Plotly figure.
    :rtype: plotly.graph_objects.Figure
    """
    if results_by_year is None:
        # Prefer characterized inventory if present; otherwise fall back to scores
        if trails.characterized_inventory is not None:
            results_by_year = trails.characterized_inventory
        elif getattr(trails, "scores", None) is not None:
            results_by_year = trails.scores
        else:
            raise ValueError(
                "No characterized inventory or scores available for plotting."
            )

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
            # New behavior: interpret as score array (e.g., trails.scores)
            if show_flow_contributions:
                # Only works if scores actually has a "flow" dim; otherwise will become totals-only.
                results_by_year = _scores_to_results(results_by_year, by_flow=True)
            else:
                results_by_year = _scores_to_results(results_by_year, by_flow=False)
    else:
        results_by_year = to_impact_year_results(results_by_year)

    if year_tick < 1:
        raise ValueError("year_tick must be >= 1")

    score_key = (
        "scores_by_flow" if show_flow_contributions else "scores_by_first_level_child"
    )
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
            y_max=y_max,
            y2_max=y2_max,
            y2_headroom=y2_headroom,
        )

    if yaxis_type != "linear" and y_max is not None:
        fig.update_layout(yaxis=dict(range=[None, float(y_max)]))

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

        (src_depth, src_year) = src_agg
        (tgt_depth, tgt_year) = tgt_agg

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
                (y_cons, a_cons) = parent
                (y_sup, a_sup) = child
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
            (y_cons, a_cons) = parent
            (y_sup, a_sup) = child
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
        (y_cons, a_cons) = parent
        (y_sup, a_sup) = child

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
