from typing import Dict, Any
import numpy as np
import plotly.graph_objects as go

from .a3 import A3  # adjust import as needed


def _build_activity_label_map(a3: A3):
    labels = {}
    for scen_label, mapping in a3.activity_indices.items():
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
    a3: A3,
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
    idx_to_label = _build_activity_label_map(a3)
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

