"""Generate the manuscript Figure 1 graph-matrix flow diagram."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle


OUTPUT = Path(__file__).resolve().with_name("algorithm_flow_diagram")


BOX = "#eaf2fb"
DECISION = "#fff3d9"
FRONTIER = "#eaf6ec"
SOLVE = "#f1edf8"
TEXT = "#17202a"
EDGE = "#374151"
LOOP = "#f8fafc"


def add_box(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    text: str,
    *,
    facecolor: str = BOX,
    fontsize: float = 10.5,
    weight: str = "normal",
) -> None:
    ax.add_patch(
        Rectangle(
            (x, y),
            w,
            h,
            facecolor=facecolor,
            edgecolor=EDGE,
            linewidth=1.4,
            joinstyle="round",
        )
    )
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color=TEXT,
        fontweight=weight,
        linespacing=1.1,
    )


def add_diamond(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    text: str,
    *,
    fontsize: float = 10.0,
) -> None:
    points = [
        (x + w / 2, y + h),
        (x + w, y + h / 2),
        (x + w / 2, y),
        (x, y + h / 2),
    ]
    ax.add_patch(
        Polygon(
            points,
            closed=True,
            facecolor=DECISION,
            edgecolor=EDGE,
            linewidth=1.4,
            joinstyle="round",
        )
    )
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color=TEXT,
        linespacing=1.05,
    )


def add_arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    rad: float = 0.0,
    label: str | None = None,
    label_xy: tuple[float, float] | None = None,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=1.4,
            color=EDGE,
            connectionstyle=f"arc3,rad={rad}",
            shrinkA=2,
            shrinkB=2,
        )
    )
    if label:
        if label_xy is None:
            label_xy = ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2)
        ax.text(
            label_xy[0],
            label_xy[1],
            label,
            ha="center",
            va="center",
            fontsize=8.7,
            color=EDGE,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.7},
        )


def add_elbow_arrow(
    ax,
    points: list[tuple[float, float]],
) -> None:
    """Draw a polyline with an arrowhead on the last segment."""
    for start, end in zip(points[:-2], points[1:-1]):
        ax.plot(
            [start[0], end[0]],
            [start[1], end[1]],
            color=EDGE,
            linewidth=1.4,
        )
    add_arrow(ax, points[-2], points[-1])


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": ["DejaVu Sans"],
            "font.size": 10,
            "svg.fonttype": "none",
        }
    )

    fig, ax = plt.subplots(figsize=(7.5, 9.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 15)
    ax.axis("off")

    add_box(
        ax,
        0.8,
        13.75,
        8.4,
        0.85,
        (
            "Inputs: FU product demand $(a_0, y_0, Q_0)$; $A_y$, $B_y$;\n"
            "temporal metadata; LCIA method(s); relative cutoff $\\tau_{rel}$"
        ),
        fontsize=10.5,
    )
    add_box(
        ax,
        2.05,
        12.6,
        5.9,
        0.75,
        "Initialise graph with root activity-scaling node",
    )

    ax.add_patch(
        Rectangle(
            (0.55, 3.85),
            8.9,
            8.3,
            facecolor=LOOP,
            edgecolor="#9ca3af",
            linewidth=1.1,
        )
    )
    ax.text(
        0.9,
        11.9,
        "routing loop",
        ha="left",
        va="center",
        fontsize=9.5,
        color="#4b5563",
        fontweight="bold",
    )

    add_box(ax, 2.05, 11.0, 5.9, 0.75, "Pop queued state $(y, a, x, d)$")
    add_box(
        ax,
        1.75,
        9.95,
        6.5,
        0.8,
        "Read $A_y$ technosphere row\nand flag routed direct $B_y$ flows",
    )
    add_diamond(ax, 2.5, 8.65, 5.0, 0.95, "Temporal distribution\non exchange?")
    add_box(ax, 0.95, 7.6, 3.2, 0.8, "No: child\n$(j, y)$")
    add_box(ax, 5.85, 7.6, 3.2, 0.8, "Yes: pulses\n$(j, y + offset)$")
    add_box(
        ax,
        1.55,
        6.5,
        6.9,
        0.8,
        "Child scaling $x_{child}$: ported or target-year coefficient,\nnormalised by production amount",
    )
    add_box(
        ax,
        1.55,
        5.35,
        6.9,
        0.8,
        "Convert to reference-product demand $Q_{child}$\n"
        "and estimate $p_{child}=|Q_{child}|\\max_m|s_m(j,y)|$",
    )
    add_diamond(
        ax,
        2.25,
        4.0,
        5.5,
        1.0,
        "Expand child?\n$p_{child} > \\tau_{rel} p_{FU}$",
    )
    add_box(
        ax,
        0.95,
        2.9,
        3.35,
        0.85,
        "Record frontier demand\n(cap, amount, or cutoff)",
        facecolor=FRONTIER,
    )
    add_box(
        ax,
        5.7,
        2.9,
        3.35,
        0.85,
        "Queue child\nfor expansion",
        facecolor=FRONTIER,
    )
    add_box(
        ax,
        1.55,
        1.6,
        6.9,
        0.8,
        "Aggregate frontier demands by year",
        facecolor="#e8f6f4",
    )
    add_box(
        ax,
        0.9,
        0.35,
        8.2,
        0.9,
        (
            "Year-wise solves: $A_y x_y = d_y$ for frontier; add routed $B_y$ flows;\n"
            "characterise/store inventory; optional dynamic impacts"
        ),
        facecolor=SOLVE,
        fontsize=9.4,
    )

    add_arrow(ax, (5.0, 13.75), (5.0, 13.35))
    add_arrow(ax, (5.0, 12.6), (5.0, 11.75))
    add_arrow(ax, (5.0, 11.0), (5.0, 10.75))
    add_arrow(ax, (5.0, 9.95), (5.0, 9.6))
    add_arrow(ax, (4.0, 8.65), (2.85, 8.4), label="no", label_xy=(3.35, 8.36))
    add_arrow(ax, (6.0, 8.65), (7.15, 8.4), label="yes", label_xy=(6.65, 8.36))
    add_arrow(ax, (2.55, 7.6), (4.6, 7.3), rad=-0.08)
    add_arrow(ax, (7.45, 7.6), (5.4, 7.3), rad=0.08)
    add_arrow(ax, (5.0, 6.5), (5.0, 6.15))
    add_arrow(ax, (5.0, 5.35), (5.0, 5.0))
    add_arrow(ax, (3.75, 4.0), (2.7, 3.75), label="no", label_xy=(3.35, 3.98))
    add_arrow(ax, (6.25, 4.0), (7.3, 3.75), label="yes", label_xy=(6.65, 3.98))
    add_elbow_arrow(
        ax,
        [
            (9.05, 3.34),
            (9.65, 3.34),
            (9.65, 11.38),
            (7.96, 11.38),
        ],
    )
    add_arrow(ax, (2.62, 2.9), (4.45, 2.4), rad=-0.1)
    add_arrow(ax, (5.0, 1.6), (5.0, 1.25))

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
