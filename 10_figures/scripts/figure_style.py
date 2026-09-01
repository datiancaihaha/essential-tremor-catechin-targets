from __future__ import annotations

import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap


MM = 1 / 25.4

COLORS = ["#4F587D", "#C68DC0", "#C2E0EE", "#4A6FD3", "#6EA8E6", "#FFE5A3", "#FFB07A", "#FF7F7A", "#B7D8F7", "#FF9BC5"]


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "axes.titlesize": 9,
            "axes.labelsize": 7.2,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "text.color": "black",
            "axes.labelcolor": "black",
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def figure_palette() -> dict[str, str]:
    p = COLORS
    return {
        "microbe": p[0],
        "metabolite": p[1],
        "genetics": p[3],
        "bulk": p[4],
        "purkinje": p[6],
        "cell": p[8],
        "positive": p[7],
        "negative": p[3],
        "support": p[9],
        "light": p[2],
        "neutral": "#B8B8B8",
        "pale": "#F2F2F2",
        "dark": "#333333",
    }


def diverging_cmap(colors: dict[str, str]):
    return LinearSegmentedColormap.from_list("signed", [colors["negative"], "#F7F7F7", colors["positive"]])


def panel_label(ax, label: str, x: float = -0.08, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=18, fontweight="bold", va="bottom", ha="left", clip_on=False)


def clean_axis(ax, grid: bool = False) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis="x", color="#D9D9D9", lw=0.5, zorder=0)
