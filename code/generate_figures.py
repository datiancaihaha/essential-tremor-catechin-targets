from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.ticker import NullFormatter
import numpy as np
import pandas as pd
from PIL import Image, ImageOps

from figure_style import (
    MM,
    clean_axis,
    diverging_cmap,
    panel_label,
    figure_palette,
    setup_style,
)


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data"
ASSETS = ROOT / "assets"
FORMAL = ROOT / "output"
REGISTRY = ROOT / "figure_documentation"
FIG_WIDTH_MM = 183

mpl.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})


def save_formal(fig, stem: Path, dpi: int = 600) -> dict:
    """Export the submission pair directly from the plotting backend."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    pdf = stem.with_suffix(".pdf")
    tiff = stem.with_suffix(".tiff")
    with tempfile.TemporaryDirectory() as temporary_directory:
        editable_svg = Path(temporary_directory) / f"{stem.name}.svg"
        fig.savefig(editable_svg, format="svg", bbox_inches=None, facecolor="white")
        if "<text" not in editable_svg.read_text(encoding="utf-8"):
            raise RuntimeError("The vector export does not contain editable text.")
    fig.savefig(
        pdf,
        format="pdf",
        bbox_inches=None,
        facecolor="white",
        metadata={"Creator": "", "Producer": ""},
    )
    fig.savefig(
        tiff,
        format="tiff",
        dpi=dpi,
        pil_kwargs={"compression": "tiff_lzw"},
        bbox_inches=None,
        facecolor="white",
    )
    with Image.open(tiff) as source:
        rgb = source.convert("RGB")
        rgb.save(tiff, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
    with Image.open(tiff) as image:
        raw_dpi = image.info.get("dpi", (None, None))
        dpi_values = [float(value) if value is not None else None for value in raw_dpi]
        return {
            "pdf": str(pdf),
            "tiff": str(tiff),
            "pixels": [int(value) for value in image.size],
            "mode": image.mode,
            "dpi": dpi_values,
        }


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "yes", "1"})


def neglog10(values) -> np.ndarray:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    return -np.log10(np.clip(arr, 1e-300, None))


def short_taxon(value: str) -> str:
    value = str(value)
    replacements = {
        "LachnospiraceaeUCG001": "Lachnospiraceae UCG-001",
        "RuminococcaceaeUCG011": "Ruminococcaceae UCG-011",
        "Prevotella7": "Prevotella 7",
        "Methanobacteria": "Methanobacteria lineage",
    }
    return replacements.get(value, value)


def load_data() -> dict[str, pd.DataFrame | dict]:
    data: dict[str, pd.DataFrame | dict] = {
        "mr": pd.read_csv(SRC / "mendelian_randomization_nonduplicate_instrument_sets.csv"),
        "mr10": pd.read_csv(SRC / "mendelian_randomization_nominal_associations.csv"),
        "methods": pd.read_csv(SRC / "mendelian_randomization_focal_estimators.csv"),
        "ivs": pd.read_csv(SRC / "mendelian_randomization_instruments.csv"),
        "loo": pd.read_csv(SRC / "mendelian_randomization_leave_one_out.csv"),
        "all_methods": pd.read_csv(SRC / "mendelian_randomization_all_estimators.csv"),
        "sensitivity": pd.read_csv(SRC / "mendelian_randomization_sensitivity_statistics.csv"),
        "swiss": pd.read_csv(SRC / "swisstargetprediction_targets.csv"),
        "sea": pd.read_csv(SRC / "similarity_ensemble_approach_targets.csv"),
        "union": pd.read_csv(SRC / "predicted_target_union.csv"),
        "magma": pd.read_csv(SRC / "magma_predicted_targets.csv"),
        "locus": pd.read_csv(SRC / "carbonic_anhydrase_locus.csv"),
        "host": pd.read_csv(SRC / "predicted_target_annotations.csv"),
        "cell": pd.read_csv(SRC / "cerebellar_cell_type_expression.csv"),
        "host23": pd.read_csv(SRC / "predicted_targets_with_human_evidence.csv"),
        "author36": pd.read_csv(SRC / "purkinje_published_differential_expression.csv"),
        "published_comparison": pd.read_csv(SRC / "purkinje_differential_expression_comparison.csv"),
        "idmap": pd.read_csv(SRC / "gene_identifier_mapping.csv"),
        "mroast": pd.read_csv(SRC / "mroast_gene_set_results.csv"),
        "fry": pd.read_csv(SRC / "fry_gene_set_results.csv"),
        "weighted_roast": pd.read_csv(SRC / "prediction_probability_weighted_roast.csv"),
        "target_expression": pd.read_csv(SRC / "predicted_target_expression_results.csv"),
        "ca3_coloc_grid": pd.read_csv(SRC / "ca3_colocalization_prior_grid.csv"),
        "ca_eqtl": pd.read_csv(SRC / "carbonic_anhydrase_cell_type_cis_eqtl.csv"),
        "umap": pd.read_csv(SRC / "scp3177_umap_coordinates.csv.gz"),
        "pca": pd.read_csv(SRC / "gse197345_rlog_pca_coordinates.csv"),
        "pca_variance": pd.read_csv(SRC / "gse197345_rlog_pca_variance.csv"),
        "purkinje_heatmap": pd.read_csv(SRC / "gse197345_selected_gene_rlog_zscores.csv"),
        "purkinje_heatmap_order": pd.read_csv(SRC / "gse197345_heatmap_sample_order.csv"),
    }
    data["identity"] = json.loads((SRC / "metabolite_identity.json").read_text(encoding="utf-8"))
    return data


def _panel(ax, label: str, title: str, x: float = -0.08, y: float = 1.02, pad: float = 7) -> None:
    panel_label(ax, label, x, y)
    ax.set_title(title, loc="left", pad=pad)


def _cell_abbreviation(value: str) -> str:
    lookup = {
        "Bergmann glial cell": "BG",
        "oligodendrocyte": "OL",
        "microglia": "MG",
        "endothelial cell": "EC",
        "granule cell": "GC",
        "Purkinje cell": "PC",
        "astrocyte": "AS",
        "molecular layer interneuron": "MLI",
        "oligodendrocyte precursor cell": "OPC",
        "pericyte": "PER",
        "Golgi cell": "GOL",
        "unipolar brush cell": "UBC",
        "choroid plexus epithelial cell": "CPE",
    }
    return lookup.get(str(value), str(value)[:4].upper())


def _load_structure(path: Path):
    image = Image.open(path).convert("L")
    return ImageOps.autocontrast(image)


def _structure(ax, path: Path, title: str, subtitle: str, title_size: float = 7.2) -> None:
    ax.imshow(_load_structure(path), cmap="gray", aspect="equal")
    ax.set_axis_off()
    ax.text(0.5, -0.04, title, transform=ax.transAxes, ha="center", va="top",
            fontsize=title_size, linespacing=0.92)
    subtitle_y = -0.14 - 0.08 * title.count("\n")
    ax.text(0.5, subtitle_y, subtitle, transform=ax.transAxes, ha="center", va="top",
            fontsize=6.0, color="#333333", linespacing=0.92)


def _ribbon(ax, x0: float, x1: float, y0: float, y1: float, width0: float, width1: float, color: str, alpha: float = 0.25) -> None:
    verts = [
        (x0, y0 - width0 / 2),
        (x0 + (x1 - x0) * 0.35, y0 - width0 / 2),
        (x0 + (x1 - x0) * 0.65, y1 - width1 / 2),
        (x1, y1 - width1 / 2),
        (x1, y1 + width1 / 2),
        (x0 + (x1 - x0) * 0.65, y1 + width1 / 2),
        (x0 + (x1 - x0) * 0.35, y0 + width0 / 2),
        (x0, y0 + width0 / 2),
        (x0, y0 - width0 / 2),
    ]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts, codes), fc=color, ec="none", alpha=alpha, transform=ax.transAxes))


def figure1(d, c):
    # A single full-width panel follows the study-overview hierarchy used in
    # Skuladottir et al. Fig. 1, with aligned quantitative result tracks in
    # the style of their Fig. 3.
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 175 * MM))
    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.set_xlim(0, 1); canvas.set_ylim(0, 1); canvas.set_axis_off()

    canvas.text(0.055, 0.954, "Study overview", fontsize=10.2, fontweight="bold", va="center", ha="left")

    # Exposure and outcome data.
    canvas.add_patch(Rectangle((0.055, 0.813), 0.505, 0.117, fc=c["microbe"], ec="none", alpha=0.10))
    canvas.add_patch(Rectangle((0.055, 0.922), 0.505, 0.008, fc=c["microbe"], ec="none", alpha=0.95))
    canvas.text(0.070, 0.902, "Exposure GWAS", fontsize=7.6, va="center")
    canvas.text(0.070, 0.874, "MiBioGen microbiome GWAS", fontsize=6.8, va="center")
    stage_x = [0.290, 0.395, 0.500]
    stage_n = ["211", "10", "8"]
    stage_label = ["taxa", "nominal\nassociations", "nonduplicate\ninstrument sets"]
    canvas.plot([stage_x[0], stage_x[-1]], [0.875, 0.875], color=c["microbe"], lw=1.0, alpha=0.55)
    for x, n, label in zip(stage_x, stage_n, stage_label):
        canvas.scatter([x], [0.875], s=54, fc="white", ec=c["microbe"], lw=1.1, zorder=3)
        canvas.text(x, 0.875, n, fontsize=6.6, ha="center", va="center")
        canvas.text(x, 0.857, label, fontsize=5.0, ha="center", va="top", linespacing=0.90)
    canvas.text(0.070, 0.821, "Nominal P < 0.05; identical instrument sets were collapsed", fontsize=5.8, va="center")

    canvas.add_patch(Rectangle((0.585, 0.813), 0.360, 0.117, fc=c["genetics"], ec="none", alpha=0.10))
    canvas.add_patch(Rectangle((0.585, 0.922), 0.360, 0.008, fc=c["genetics"], ec="none", alpha=0.95))
    canvas.text(0.600, 0.902, "Outcome GWAS", fontsize=7.6, va="center")
    canvas.text(0.600, 0.875, "Essential tremor GWAS meta-analysis", fontsize=6.8, va="center")
    canvas.text(0.630, 0.859, "16,480", fontsize=9.0, ha="center", va="center", color=c["genetics"])
    canvas.text(0.630, 0.842, "cases", fontsize=6.1, ha="center", va="top")
    canvas.text(0.790, 0.859, "1,936,173", fontsize=9.0, ha="center", va="center", color=c["genetics"])
    canvas.text(0.790, 0.842, "controls", fontsize=6.1, ha="center", va="top")

    # Mendelian randomization estimates.
    canvas.text(0.055, 0.795, "Mendelian randomization", fontsize=8.1, va="center")
    canvas.plot([0.055, 0.480], [0.782, 0.782], color="#CFCFCF", lw=0.6)
    mr = d["mr"].copy().sort_values("p_value", ascending=True).reset_index(drop=True)
    ax_mr = fig.add_axes([0.195, 0.500, 0.290, 0.265])
    y = np.arange(len(mr))
    focal = mr.taxon_name.isin(["Faecalibacterium", "Flavonifractor"])
    for i in np.where(focal)[0]:
        ax_mr.axhspan(i - 0.43, i + 0.43, color=c["metabolite"], alpha=0.10, zorder=0)
    for i, row in mr.iterrows():
        color = c["metabolite"] if focal.iloc[i] else c["microbe"]
        ax_mr.errorbar(row.odds_ratio, i,
                       xerr=[[row.odds_ratio - row.or_ci_lower], [row.or_ci_upper - row.odds_ratio]],
                       fmt="o", color=color, ecolor="#555555", elinewidth=0.75,
                       capsize=1.8, ms=4.2, zorder=3)
        ax_mr.text(1.55, i, f"P={row.p_value:.3f}", ha="right", va="center", fontsize=6.0)
    ax_mr.axvline(1, color="#666666", lw=0.75, ls="--")
    ax_mr.set_xscale("log"); ax_mr.set_xlim(0.56, 1.58)
    ax_mr.set_xticks([0.6, 0.8, 1.0, 1.2, 1.4], ["0.6", "0.8", "1.0", "1.2", "1.4"])
    ax_mr.xaxis.set_minor_formatter(NullFormatter())
    mr_labels = [short_taxon(value) for value in mr.taxon_name]
    ax_mr.set_yticks(y, mr_labels)
    for tick, label in zip(ax_mr.get_yticklabels(), mr_labels):
        tick.set_fontstyle("normal" if label == "Methanobacteria lineage" else "italic")
    ax_mr.tick_params(axis="y", labelsize=6.1, length=0, pad=3)
    ax_mr.tick_params(axis="x", labelsize=6.1)
    ax_mr.set_xlabel("Odds ratio for essential tremor (95% CI)", fontsize=6.5, labelpad=2)
    clean_axis(ax_mr, grid=True); ax_mr.invert_yaxis()

    # Evidence lanes preserve the distinction between association and pure-culture transformation.
    canvas.text(0.515, 0.795, "Microbe–metabolite evidence", fontsize=8.1, va="center")
    canvas.plot([0.515, 0.810], [0.782, 0.782], color="#CFCFCF", lw=0.6)
    lane_y = [0.738, 0.625]
    canvas.plot([0.485, 0.510], [lane_y[0], lane_y[0]], color=c["metabolite"], lw=1.1)
    canvas.scatter([0.520], [lane_y[0]], s=38, fc="white", ec=c["metabolite"], lw=1.2, zorder=3)
    canvas.text(0.535, lane_y[0] + 0.013, "Faecalibacterium", fontsize=6.7, fontstyle="italic", va="center")
    canvas.text(0.535, lane_y[0] - 0.018, "Human association\n(correlational)",
                fontsize=5.9, va="center", linespacing=1.0)
    canvas.plot([0.655, 0.681], [lane_y[0], 0.690], color=c["metabolite"], lw=1.0, alpha=0.75)

    canvas.plot([0.485, 0.510], [lane_y[1], lane_y[1]], color=c["support"], lw=1.1)
    canvas.scatter([0.520], [lane_y[1]], s=38, fc=c["support"], ec="black", lw=0.45, zorder=3)
    canvas.text(0.535, lane_y[1] + 0.032, "Flavonifractor", fontsize=6.7, fontstyle="italic", va="center")
    canvas.text(0.535, lane_y[1] + 0.004, "F. plautii", fontsize=6.1, fontstyle="italic", va="center")
    canvas.text(0.535, lane_y[1] - 0.022, "Pure-culture metabolism", fontsize=6.1, va="center")
    canvas.text(0.535, lane_y[1] - 0.048, "Catechin ring-fission products", fontsize=6.0, va="center")
    canvas.plot([0.655, 0.704], [lane_y[1], 0.548], color=c["support"], lw=1.0, alpha=0.70)
    for xx in [0.708, 0.724, 0.740]:
        canvas.scatter([xx], [0.542], s=22, fc="white", ec=c["support"], lw=0.9)
    canvas.text(0.724, 0.518, "Related products", fontsize=6.0, ha="center", va="center")

    ax_structure = fig.add_axes([0.690, 0.650, 0.115, 0.105])
    ax_structure.imshow(_load_structure(ASSETS / "PubChem_CID49831816.png"), cmap="gray", aspect="equal")
    ax_structure.set_axis_off()
    canvas.text(0.748, 0.634, "5-(3,4-dihydroxyphenyl)\npentanoic acid", fontsize=6.0,
                ha="center", va="top", linespacing=1.05)
    canvas.text(0.748, 0.595, "CID 49831816", fontsize=6.0, ha="center", va="top", color=c["metabolite"])

    # Target prediction is connected only to the exact metabolite association.
    canvas.text(0.835, 0.795, "Target prediction", fontsize=8.1, va="center")
    canvas.plot([0.835, 0.945], [0.782, 0.782], color="#CFCFCF", lw=0.6)
    canvas.plot([0.805, 0.835], [0.690, 0.690], color=c["metabolite"], lw=4.5, alpha=0.20)
    canvas.text(0.890, 0.724, "105", fontsize=13.0, ha="center", va="center", color=c["support"])
    canvas.text(0.890, 0.699, "predicted targets", fontsize=6.2, ha="center", va="center")
    bar_x, bar_y, bar_w, bar_h = 0.835, 0.655, 0.110, 0.025
    swiss_w = bar_w * 98 / 105
    canvas.add_patch(Rectangle((bar_x, bar_y), swiss_w, bar_h, fc=c["support"], ec="white", lw=0.5, alpha=0.78))
    canvas.add_patch(Rectangle((bar_x + swiss_w, bar_y), bar_w - swiss_w, bar_h, fc=c["genetics"], ec="white", lw=0.5, alpha=0.95))
    canvas.text(bar_x, 0.642, "SwissTargetPrediction 98", fontsize=6.0, ha="left", va="top")
    canvas.text(bar_x + bar_w, 0.620, "SEA 7", fontsize=6.0, ha="right", va="top")
    canvas.text(0.890, 0.592, "Observed overlap = 0", fontsize=6.0, ha="center", va="center")

    canvas.text(0.515, 0.487,
                "F. plautii pure-culture products are chemically distinct from CID 49831816;\n"
                "genus-level and species-level results are displayed separately.",
                fontsize=6.1, va="top", ha="left", color="#333333", linespacing=1.15)

    # Aligned quantitative tracks summarize gene association and cerebellar expression.
    canvas.text(0.055, 0.430, "Gene association and cerebellar expression", fontsize=8.1, va="center")
    canvas.plot([0.055, 0.945], [0.417, 0.417], color="#CFCFCF", lw=0.6)
    ax_host = fig.add_axes([0.235, 0.150, 0.545, 0.225])
    yhost = np.array([3, 2, 1, 0])
    labels = ["MAGMA gene-based association", "GSE134878 whole cerebellum",
              "GSE197345 Purkinje cells", "SCP3177 cell-type expression"]
    track_specs = [
        [(8, c["neutral"], "//"), (86, "#EEEEEE", None), (11, c["genetics"], None)],
        [(41, c["neutral"], "//"), (59, "#EEEEEE", None), (5, c["bulk"], None)],
        [(43, c["neutral"], "//"), (54, "#EEEEEE", None), (8, c["purkinje"], None)],
        [(105, c["cell"], None)],
    ]
    for yy, segments in zip(yhost, track_specs):
        left = 0
        for width, color, hatch in segments:
            ax_host.barh(yy, width, left=left, height=0.56, color=color, edgecolor="white", lw=0.55, hatch=hatch)
            if width >= 5:
                ax_host.text(left + width / 2, yy, str(width), fontsize=6.0, ha="center", va="center")
            left += width
    ax_host.set_xlim(0, 105); ax_host.set_ylim(-0.65, 3.65)
    ax_host.set_xticks([0, 25, 50, 75, 105])
    ax_host.set_yticks(yhost, labels)
    ax_host.tick_params(axis="y", labelsize=6.1, length=0, pad=4)
    ax_host.tick_params(axis="x", labelsize=6.1)
    ax_host.set_xlabel("Predicted targets (n)", fontsize=6.5, labelpad=2)
    clean_axis(ax_host, grid=True)

    summary_x = 0.780
    summary_rows = [
        (0.345, c["genetics"], "97 tested | 11 nominal"),
        (0.292, c["bulk"], "64 testable | 5 nominal"),
        (0.239, c["purkinje"], "62 testable | mroast P = 0.048"),
        (0.186, c["cell"], "105 genes | 13 cell types"),
    ]
    for yy, color, text_value in summary_rows:
        canvas.add_patch(Rectangle((summary_x, yy - 0.008), 0.010, 0.016, fc=color, ec="none", alpha=0.9))
        canvas.text(summary_x + 0.016, yy, text_value, fontsize=6.0, va="center", ha="left")

    canvas.text(0.235, 0.105,
                "Gray hatch, unavailable or excluded by independent expression filtering; pale gray, evaluated with nominal P ≥ 0.05;\n"
                "colored segments, nominal P < 0.05.",
                fontsize=6.0, va="top", ha="left", linespacing=1.15)
    canvas.text(0.055, 0.043,
                "Nominal P < 0.05 defined the screening threshold. Gene association, differential expression and cell-type expression were evaluated separately.",
                fontsize=6.2, va="center", ha="left")
    return fig


def figure2(d, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 250 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.94, 1.08, 0.92], width_ratios=[1.08, 0.92],
                          hspace=0.43, wspace=0.42, left=0.225, right=0.985, top=0.965, bottom=0.060)
    mr = d["mr"].copy()
    mr["display_order"] = mr.taxon_name.map({"Faecalibacterium": 0, "Flavonifractor": 1}).fillna(2)
    mr = mr.sort_values(["display_order", "p_value", "taxon_name"]).reset_index(drop=True)
    forest_grid = gs[0, 0].subgridspec(1, 2, width_ratios=[0.78, 0.22], wspace=0.02)
    ax = fig.add_subplot(forest_grid[0]); _panel(ax, "A", "IVW estimates", -0.39, 1.02)
    y = np.arange(len(mr))
    ax.errorbar(mr.odds_ratio, y, xerr=[mr.odds_ratio - mr.or_ci_lower, mr.or_ci_upper - mr.odds_ratio],
                fmt="o", color=c["microbe"], ecolor="#3D3D3D", elinewidth=0.9, capsize=2.2, ms=5, zorder=3)
    ax.axvline(1, color="#666666", lw=0.8, ls="--")
    labels = [short_taxon(x) for x in mr.taxon_name]
    ax.set_yticks(y, labels=labels)
    for tick, label in zip(ax.get_yticklabels(), labels):
        tick.set_fontstyle("normal" if label == "Methanobacteria lineage" else "italic")
    ax.set_xscale("log"); ax.set_xlim(0.58, 1.50)
    ax.set_xticks([0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4], ["0.6", "0.7", "0.8", "0.9", "1.0", "1.2", "1.4"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Odds ratio for essential tremor (95% CI)")
    clean_axis(ax, grid=True)
    ax.invert_yaxis()
    pax = fig.add_subplot(forest_grid[1], sharey=ax)
    pax.set_xlim(0, 1); pax.set_ylim(ax.get_ylim()); pax.set_axis_off()
    pax.text(0.02, 1.035, "P value", transform=pax.transAxes, fontsize=7.0, fontweight="bold", va="bottom")
    for yi, row in mr.iterrows():
        pax.text(0.02, yi, f"{row.p_value:.3f}", ha="left", va="center", fontsize=6.5)

    sub = gs[0, 1].subgridspec(2, 1, height_ratios=[0.62, 0.38], hspace=0.35)
    ax = fig.add_subplot(sub[0]); _panel(ax, "B", "MR estimators and diagnostics", -0.19, 1.02)
    methods = d["methods"].copy()
    method_order = ["Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode"]
    markers = ["o", "s", "D", "^"]
    cols = [c["microbe"], c["genetics"], c["bulk"], c["purkinje"]]
    offsets = [-0.18, -0.06, 0.06, 0.18]
    focal = ["Faecalibacterium", "Flavonifractor"]
    for j, method in enumerate(method_order):
        q = methods[methods.method.eq(method)].copy()
        q["taxon"] = q.taxon_id.str.extract(r"genus\.([^.]+)")[0]
        q = q.set_index("taxon")
        for i, taxon in enumerate(focal):
            if taxon not in q.index:
                continue
            row = q.loc[taxon]
            ax.errorbar(row.odds_ratio, i + offsets[j],
                        xerr=[[row.odds_ratio - row.or_ci_lower], [row.or_ci_upper - row.odds_ratio]],
                        fmt=markers[j], color=cols[j], ecolor="#555555", ms=4, capsize=2, lw=0.8)
    ax.axvline(1, color="#666666", lw=0.8, ls="--")
    ax.set_xscale("log"); ax.set_xlim(0.12, 1.68)
    ax.set_xticks([0.2, 0.5, 1.0, 1.5], ["0.2", "0.5", "1.0", "1.5"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(range(2), focal, fontstyle="italic"); ax.invert_yaxis(); ax.set_xlabel("Odds ratio (95% CI)")
    clean_axis(ax, grid=True)
    handles = [Line2D([], [], marker=mk, color=co, ls="", label=lab.replace("Inverse variance weighted", "IVW"))
               for mk, co, lab in zip(markers, cols, method_order)]
    axq = fig.add_subplot(sub[1]); axq.set_axis_off()
    axq.legend(handles=handles, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=2,
               columnspacing=0.9, handletextpad=0.35, fontsize=5.8)
    mr_index = d["mr"].set_index("taxon_name")
    diag_rows = []
    for taxon in focal:
        r = mr_index.loc[taxon]
        diag_rows.append((taxon, int(r.nsnp), float(r.ivw_q_p), float(r.egger_intercept_p), float(r.mr_presso_global_p), bool(r.leave_one_out_direction_concordant)))
    headers = ["Taxon", "SNPs", "Q P", "Egger P", "MR-PRESSO\nP", "LOO"]
    x = [0.00, 0.36, 0.48, 0.61, 0.80, 0.96]
    for index, (xx, h) in enumerate(zip(x, headers)):
        axq.text(xx, 0.46, h, transform=axq.transAxes,
                 fontsize=5.8, va="center",
                 ha="left" if index == 0 else "center", linespacing=0.9)
    for i, row in enumerate(diag_rows):
        yy = 0.23 - i * 0.20
        vals = [row[0], str(row[1]), f"{row[2]:.3f}", f"{row[3]:.3f}", f"{row[4]:.3f}", "Retained" if row[5] else "Changed"]
        for index, (xx, val) in enumerate(zip(x, vals)):
            axq.text(xx, yy, val, transform=axq.transAxes, fontsize=5.8, va="center",
                     ha="left" if index == 0 else "center")
        if i == 0:
            axq.plot([0, 1], [0.13, 0.13], transform=axq.transAxes, color="#E0E0E0", lw=0.5)

    focal_specs = [("Faecalibacterium", "genus.Faecalibacterium.id.2057", c["microbe"]),
                   ("Flavonifractor", "genus.Flavonifractor.id.2059", c["metabolite"])]

    axbg = fig.add_subplot(gs[1, 0]); axbg.set_axis_off()
    _panel(axbg, "C", "SNP association estimates", -0.39, 1.02)
    for k, (taxon, taxon_id, color) in enumerate(focal_specs):
        inset = axbg.inset_axes([0.02, 0.06 + (1 - k) * 0.50, 0.96, 0.36])
        q = _focal_iv_data(d, taxon_id)
        inset.errorbar(q.beta_exposure, q.beta_outcome_aligned, xerr=q.se_exposure, yerr=q.se_outcome,
                       fmt="o", ms=3.5, color=color, ecolor="#888888", lw=0.55, alpha=0.9)
        method_rows = d["methods"][d["methods"].taxon_id.eq(taxon_id)].copy()
        xline = np.linspace(min(0, q.beta_exposure.min() * 1.15), q.beta_exposure.max() * 1.15, 100)
        styles = [("Inverse variance weighted", "-", c["microbe"]), ("MR Egger", "--", c["genetics"]),
                  ("Weighted median", ":", c["bulk"]), ("Weighted mode", "-.", c["purkinje"])]
        sens = d["sensitivity"].set_index("taxon_id").loc[taxon_id]
        for method, line_style, method_color in styles:
            result = method_rows[method_rows.method.eq(method)].iloc[0]
            intercept = float(sens.egger_intercept) if method == "MR Egger" else 0.0
            inset.plot(xline, intercept + float(result.beta) * xline, ls=line_style,
                       color=method_color, lw=0.9,
                       label=method.replace("Inverse variance weighted", "IVW"))
        inset.set_title(taxon, loc="left", fontsize=6.4, fontstyle="italic", pad=1)
        if k == 1:
            inset.set_xlabel("SNP effect on microbial abundance", fontsize=5.6)
        inset.set_ylabel("SNP effect on ET", fontsize=5.6)
        inset.tick_params(labelsize=5.2)
        clean_axis(inset)
        if k == 0:
            inset.legend(frameon=False, fontsize=5.2, ncol=2, loc="lower right",
                         columnspacing=0.7, handlelength=1.5)

    axbg = fig.add_subplot(gs[1, 1]); axbg.set_axis_off()
    _panel(axbg, "D", "Leave-one-out estimates", -0.19, 1.02)
    for k, (taxon, taxon_id, color) in enumerate(focal_specs):
        inset = axbg.inset_axes([0.20, 0.06 + (1 - k) * 0.50, 0.77, 0.36])
        loo = d["loo"][d["loo"].taxon_id.eq(taxon_id)].copy().sort_values("beta")
        yv = np.arange(len(loo))
        inset.errorbar(loo.beta, yv, xerr=[loo.beta - loo.ci_lower, loo.ci_upper - loo.beta],
                       fmt="o", color=color, ecolor="#555555", lw=0.65, capsize=1.5, ms=3.5)
        full = float(d["mr"].set_index("taxon_id").loc[taxon_id, "beta"])
        inset.axvline(full, color="#555555", lw=0.8, ls="--")
        inset.set_yticks(yv, loo.omitted_rsid, fontsize=5.2)
        inset.invert_yaxis()
        inset.set_title(taxon, loc="left", fontsize=6.4, fontstyle="italic", pad=1)
        if k == 1:
            inset.set_xlabel("IVW estimate after SNP exclusion", fontsize=5.6)
        inset.tick_params(axis="x", labelsize=5.2)
        clean_axis(inset, grid=True)

    ax = fig.add_subplot(gs[2, :]); _panel(ax, "E", "Sensitivity analyses and microbial evidence", -0.14, 1.02, pad=34)
    rows = [short_taxon(x) for x in mr.taxon_name]
    col_labels = ["IVW\nP < 0.05", "Concordant sensitivity\nestimates", "Leave-one-out\ndirection", "ET case-control\nassociation", "gutMGene\nassociation", "Pure-culture\nmetabolism"]
    mat = np.zeros((8, 6), int); mat[:, 0] = 2
    mat[:, 1] = np.where(as_bool(mr.sensitivity_estimates_direction_concordant), 2, 0)
    mat[:, 2] = np.where(as_bool(mr.leave_one_out_direction_concordant), 2, 0)
    for i, taxon in enumerate(mr.taxon_name):
        if taxon == "Faecalibacterium": mat[i, 3], mat[i, 4] = 2, 1
        if taxon == "Flavonifractor": mat[i, 5] = 2
    ax.set_xlim(-0.5, 5.5); ax.set_ylim(7.5, -0.5)
    for i in range(8):
        for j in range(6):
            ax.add_patch(Rectangle((j - 0.48, i - 0.45), 0.96, 0.90, fc="#F5F5F5", ec="white"))
            if mat[i, j] == 2:
                ax.scatter(j, i, s=70, fc=c["support"], ec="black", lw=0.45, zorder=3)
            elif mat[i, j] == 1:
                ax.scatter(j, i, s=70, fc="white", ec=c["support"], lw=1.4, zorder=3)
    ax.set_xticks(range(6), col_labels); ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, length=0)
    ax.set_yticks(range(8), rows); ax.tick_params(axis="y", length=0)
    for tick, label in zip(ax.get_yticklabels(), rows):
        tick.set_fontstyle("normal" if label == "Methanobacteria lineage" else "italic")
    ax.spines[:].set_visible(False)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=c["support"], mec="black", label="Evidence present"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec=c["support"], label="Correlational association")],
              loc="lower center", bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False)
    return fig


def figure3(d, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 205 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.72, 0.50, 1.18], width_ratios=[1.0, 1.0],
                          hspace=0.18, wspace=0.35, left=0.075, right=0.985, top=0.96, bottom=0.055)

    axbg = fig.add_subplot(gs[0, :]); axbg.set_axis_off(); _panel(axbg, "A", "Microbial associations and compound identity", -0.055, 1.02)
    axbg.text(0.055, 0.69, "Faecalibacterium", transform=axbg.transAxes, ha="center", va="center",
              fontsize=6.6, fontstyle="italic")
    axbg.add_patch(FancyArrowPatch((0.12, 0.69), (0.285, 0.69), transform=axbg.transAxes,
                                   arrowstyle="<->", mutation_scale=8, lw=1.0, ls="--", color=c["microbe"]))
    axbg.text(0.202, 0.80, "Correlational human association", transform=axbg.transAxes,
              ha="center", va="center", fontsize=6.0)
    axbg.text(0.36, 0.69, "5-(3,4-Dihydroxyphenyl)\npentanoic acid\nCID 49831816",
              transform=axbg.transAxes, ha="center", va="center", fontsize=6.4, linespacing=1.05)
    axbg.add_patch(FancyArrowPatch((0.445, 0.69), (0.555, 0.69), transform=axbg.transAxes,
                                   arrowstyle="-|>", mutation_scale=8, lw=1.0, ls=":", color=c["support"]))
    axbg.text(0.50, 0.80, "In silico target prediction", transform=axbg.transAxes,
              ha="center", va="center", fontsize=6.0)
    axbg.text(0.61, 0.69, "105 human targets", transform=axbg.transAxes,
              ha="center", va="center", fontsize=6.6)

    axbg.text(0.055, 0.25, "Catechin ring-fission\nintermediate", transform=axbg.transAxes,
              ha="center", va="center", fontsize=6.0)
    axbg.add_patch(FancyArrowPatch((0.12, 0.25), (0.225, 0.25), transform=axbg.transAxes,
                                   arrowstyle="-|>", mutation_scale=8, lw=1.0, color=c["metabolite"]))
    axbg.text(0.29, 0.25, "Flavonifractor plautii", transform=axbg.transAxes, ha="center", va="center",
              fontsize=6.4, fontstyle="italic")
    axbg.add_patch(FancyArrowPatch((0.365, 0.25), (0.445, 0.25), transform=axbg.transAxes,
                                   arrowstyle="-|>", mutation_scale=8, lw=1.0, color=c["metabolite"]))
    axbg.text(0.565, 0.25,
              "5-(3,4-Dihydroxyphenyl)-γ-valerolactone\nand 4-hydroxy-5-(3,4-dihydroxyphenyl)\nvaleric acid",
              transform=axbg.transAxes, ha="center", va="center", fontsize=5.7, linespacing=1.05)
    axbg.text(0.29, 0.065, "Pure-culture metabolism", transform=axbg.transAxes,
              ha="center", va="center", fontsize=6.0)

    acid_ax = axbg.inset_axes([0.73, 0.22, 0.11, 0.60])
    lactone_ax = axbg.inset_axes([0.88, 0.22, 0.11, 0.60])
    _structure(acid_ax, ASSETS / "PubChem_CID49831816.png", "Open-chain acid", "CID 49831816\nC11H14O4", 5.6)
    _structure(lactone_ax, ASSETS / "PubChem_CID45093080_valerolactone.png", "Cyclic lactone", "CID 45093080\nC11H12O4", 5.6)
    axbg.text(0.86, 0.54, "≠", transform=axbg.transAxes, ha="center", va="center", fontsize=14, fontweight="bold")
    axbg.text(0.86, 0.075, "Chemically distinct", transform=axbg.transAxes, ha="center", va="center", fontsize=6.0)

    ax = fig.add_subplot(gs[1, :]); ax.set_axis_off(); _panel(ax, "B", "Predicted targets by method", -0.055, 1.04)
    ax.add_patch(Rectangle((0.08, 0.24), 0.20, 0.54, transform=ax.transAxes, fc=c["support"], ec="none", alpha=0.72))
    ax.add_patch(Rectangle((0.35, 0.24), 0.20, 0.54, transform=ax.transAxes, fc=c["genetics"], ec="none", alpha=0.82))
    ax.add_patch(Rectangle((0.68, 0.24), 0.24, 0.54, transform=ax.transAxes, fc="#F5F5F5", ec="#777777", lw=0.7))
    ax.text(0.18, 0.55, "98", transform=ax.transAxes, ha="center", va="center", fontsize=14)
    ax.text(0.18, 0.36, "SwissTargetPrediction", transform=ax.transAxes, ha="center", fontsize=6.7)
    ax.text(0.45, 0.55, "7", transform=ax.transAxes, ha="center", va="center", fontsize=14)
    ax.text(0.45, 0.36, "Similarity Ensemble Approach", transform=ax.transAxes, ha="center", fontsize=6.2)
    ax.text(0.80, 0.55, "105", transform=ax.transAxes, ha="center", va="center", fontsize=14)
    ax.text(0.80, 0.36, "Unique predicted targets", transform=ax.transAxes, ha="center", fontsize=6.7)
    ax.add_patch(FancyArrowPatch((0.28, 0.55), (0.68, 0.55), transform=ax.transAxes,
                                 arrowstyle="-|>", mutation_scale=9, lw=1.2, color=c["support"]))
    ax.add_patch(FancyArrowPatch((0.55, 0.39), (0.68, 0.48), transform=ax.transAxes,
                                 arrowstyle="-|>", mutation_scale=9, lw=1.2, color=c["genetics"]))
    ax.text(0.50, 0.08, "Observed overlap = 0; the union contained 105 unique targets.",
            transform=ax.transAxes, ha="center", fontsize=6.3)

    ax = fig.add_subplot(gs[2, 0]); _panel(ax, "C", "Highest-ranked SwissTargetPrediction targets", -0.16, 1.02)
    swiss = d["swiss"].copy()
    swiss["rank"] = pd.to_numeric(swiss["rank"]); swiss["probability"] = pd.to_numeric(swiss["probability"])
    gene_swiss = swiss.sort_values(["rank", "probability"], ascending=[True, False]).drop_duplicates("gene_symbol")
    top15 = gene_swiss.nsmallest(15, "rank").sort_values("rank")
    y = np.arange(len(top15))
    bar_colors = [c["metabolite"] if str(name).startswith("Carbonic anhydrase") else c["support"]
                  for name in top15.target_name]
    ax.barh(y, top15.probability, color=bar_colors, alpha=0.86, edgecolor="none", height=0.68)
    ax.set_yticks(y, top15.gene_symbol, fontsize=5.8)
    for yi, row in enumerate(top15.itertuples()):
        ax.text(row.probability + 0.012, yi, f"{row.probability:.2f}", va="center", ha="left", fontsize=5.4)
    ax.invert_yaxis()
    ax.set_xlabel("Predicted probability")
    ax.set_xlim(0, max(top15.probability) + 0.11); clean_axis(ax, grid=True)

    ax = fig.add_subplot(gs[2, 1]); _panel(ax, "D", "SEA score and molecular similarity", -0.16, 1.02)
    sea = d["sea"].copy()
    ax.scatter(sea.tanimoto_score, sea.z_score, s=42, fc=c["genetics"], ec="black", lw=0.45, zorder=3)
    sea_label_offsets = {
        "FFAR1": (5, 5), "FFAR4": (5, 5), "TDP1": (-5, -11), "HSD17B10": (-5, 7),
        "MMP8": (20, -13), "MMP12": (29, 1), "MMP10": (38, 14),
    }
    for row in sea.itertuples():
        dx, dy = sea_label_offsets[row.gene_symbol]
        ax.annotate(row.gene_symbol, (row.tanimoto_score, row.z_score), xytext=(dx, dy),
                    textcoords="offset points", fontsize=5.8,
                    ha="left" if dx >= 0 else "right", va="bottom" if dy >= 0 else "top",
                    arrowprops=dict(arrowstyle="-", color="#777777", lw=0.4))
    ax.set_xlabel("Tanimoto coefficient"); ax.set_ylabel("SEA z score")
    ax.set_xlim(0.28, 0.62); ax.set_ylim(5.0, 32.0); clean_axis(ax, grid=True)
    return fig


def figure4(d, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 244 * MM))
    gs = fig.add_gridspec(4, 1, height_ratios=[1.00, 1.06, 0.76, 0.92], hspace=0.56,
                          left=0.15, right=0.975, top=0.965, bottom=0.07)

    ax = fig.add_subplot(gs[0]); _panel(ax, "A", "Gene-based association of predicted targets", -0.085, 1.02)
    magma = d["magma"].copy()
    magma = magma[as_bool(magma.included_in_autosomal_magma)].copy()
    for col in ["magma_chr", "magma_start_build37", "magma_pvalue"]:
        magma[col] = pd.to_numeric(magma[col], errors="coerce")
    magma = magma.sort_values(["magma_chr", "magma_start_build37"]).reset_index(drop=True)
    magma["x"] = np.arange(len(magma)); magma["y"] = neglog10(magma.magma_pvalue)
    nominal = as_bool(magma.magma_nominal_p_lt_0_05)
    parity = magma.magma_chr.astype(int) % 2 == 0
    ax.scatter(magma.loc[~nominal & parity, "x"], magma.loc[~nominal & parity, "y"], s=17,
               fc=c["light"], ec="none", alpha=0.85)
    ax.scatter(magma.loc[~nominal & ~parity, "x"], magma.loc[~nominal & ~parity, "y"], s=17,
               fc=c["neutral"], ec="none", alpha=0.72)
    ax.scatter(magma.loc[nominal, "x"], magma.loc[nominal, "y"], s=34, fc=c["genetics"],
               ec="black", lw=0.35, zorder=3)
    threshold = -math.log10(0.05)
    ax.axhline(threshold, ls="--", lw=0.9, color="#666666")
    ax.text(96, threshold + 0.16, "Nominal P = 0.05", ha="right", fontsize=6.5)
    offsets = {
        "CA3": (-16, -15), "CA2": (14, -19), "CA1": (-4, 8), "ADORA3": (-12, 10),
        "HSD17B3": (8, 6), "NR1H3": (-7, 10), "IGF1R": (0, 9), "TTR": (8, 6),
        "AKR1C3": (-9, 8), "PTPN1": (8, 6), "SIRT1": (0, 9),
    }
    for row in magma[nominal].itertuples():
        dx, dy = offsets.get(row.gene_symbol, (0, 7))
        ax.annotate(row.gene_symbol, (row.x, row.y), xytext=(dx, dy), textcoords="offset points",
                    ha="center", va="bottom" if dy >= 0 else "top", fontsize=6.2,
                    arrowprops=dict(arrowstyle="-", color="#777777", lw=0.35) if abs(dx) > 6 else None)
    centers = magma.groupby("magma_chr").x.mean()
    chromosome_labels = [str(int(chromosome)) if int(chromosome) % 2 == 1 or int(chromosome) == 22 else ""
                         for chromosome in centers.index]
    ax.set_xticks(centers.values, chromosome_labels); ax.set_xlabel("Chromosome", labelpad=1)
    ax.set_ylabel("−log10(MAGMA gene P)"); clean_axis(ax)
    ax.text(0.995, 0.97, "97 tested | 11 nominal", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.2)
    ax.text(0.995, 0.86, "Competitive gene-set P = 0.50068", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.2)

    ax = fig.add_subplot(gs[1]); _panel(ax, "B", "Genetic and transcriptomic evidence for genes with nominal MAGMA P < 0.05", -0.085, 1.02, pad=32)
    host = d["host"].copy()
    q = host[as_bool(host.magma_nominal_p_lt_0_05)].copy().sort_values("magma_pvalue").reset_index(drop=True)
    columns = ["MAGMA\n−log10 P", "Swiss", "SEA", "Published ET\ngenetics",
               "Whole\ncerebellum", "Purkinje\ncells", "Highest-expression\ncell type"]
    ax.set_xlim(-0.5, len(columns) - 0.5); ax.set_ylim(len(q) - 0.5, -0.5)
    cell_colors = {}
    palette_cycle = [c["cell"], c["bulk"], c["microbe"], c["metabolite"], c["purkinje"], c["support"], c["light"]]
    for i, row in q.iterrows():
        for j in range(len(columns)):
            ax.add_patch(Rectangle((j - 0.48, i - 0.43), 0.96, 0.86, fc="#F6F6F6", ec="white"))
        strength = -math.log10(float(row.magma_pvalue))
        ax.scatter(0, i, s=24 + strength * 8, fc=c["genetics"], ec="black", lw=0.35)
        binary = [(1, row.in_swiss, c["support"]), (2, row.in_sea, c["genetics"]),
                  (3, row.published_et_genetic_evidence, c["metabolite"]),
                  (4, row.bulk_nominal_p_lt_0_05, c["bulk"]), (5, row.purkinje_nominal_p_lt_0_05, c["purkinje"])]
        for j, value, color in binary:
            if bool(value): ax.scatter(j, i, s=45, fc=color, ec="black", lw=0.35)
        ct = str(row.highest_mean_expression_cell_type)
        if ct not in cell_colors:
            cell_colors[ct] = palette_cycle[len(cell_colors) % len(palette_cycle)]
        ax.scatter(6, i, s=74, fc=cell_colors[ct], ec="black", lw=0.35)
        ax.text(6, i, _cell_abbreviation(ct), ha="center", va="center", fontsize=6.0)
    ax.set_xticks(range(len(columns)), columns)
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, length=0)
    ax.set_yticks(range(len(q)), q.gene_symbol); ax.tick_params(axis="y", length=0); ax.spines[:].set_visible(False)
    legend_text = "  ".join(f"{_cell_abbreviation(k)}={k}" for k in cell_colors)
    ax.text(0.0, -0.12, legend_text,
            transform=ax.transAxes, fontsize=5.8, va="top")

    ax = fig.add_subplot(gs[2]); _panel(ax, "C", "Carbonic anhydrase locus", -0.095, 1.02)
    locus = d["locus"].copy()
    numeric_cols = ["gene_start_grch37", "gene_stop_grch37", "magma_window_start_grch37", "magma_window_stop_grch37",
                    "top_snp_position_grch37_1kg_eur", "magma_gene_pvalue"]
    for col in numeric_cols: locus[col] = pd.to_numeric(locus[col], errors="coerce")
    locus = locus.sort_values("gene_start_grch37").reset_index(drop=True)
    xmin = locus.magma_window_start_grch37.min(); xmax = locus.magma_window_stop_grch37.max()
    ypos = {"CA1": 2, "CA3": 1, "CA2": 0}
    for row in locus.itertuples():
        y0 = ypos.get(row.gene_symbol, 0)
        ax.hlines(y0, row.magma_window_start_grch37, row.magma_window_stop_grch37, color=c["light"], lw=9, zorder=1)
        ax.hlines(y0, row.gene_start_grch37, row.gene_stop_grch37, color=c["genetics"], lw=3.5, zorder=2)
        ax.scatter(row.top_snp_position_grch37_1kg_eur, y0, s=45, marker="v", fc=c["support"], ec="black", lw=0.45, zorder=3)
        status = "reported putative causal gene" if row.gene_symbol == "CA3" else "same-locus gene"
        if row.gene_symbol == "CA2":
            tx, ha = row.gene_stop_grch37, "right"
        elif row.gene_symbol == "CA3":
            tx, ha = row.gene_start_grch37, "right"
        else:
            tx, ha = row.gene_start_grch37, "left"
        ax.text(tx, y0 + 0.20, f"{row.gene_symbol}  P = {row.magma_gene_pvalue:.2g}\n{status}", fontsize=5.6, va="bottom", ha=ha)
    ax.plot([xmin, xmax], [-0.52, -0.52], color="#444444", lw=0.7)
    ax.plot([xmin, xmin], [-0.52, -0.40], color="#444444", lw=0.7)
    ax.plot([xmax, xmax], [-0.52, -0.40], color="#444444", lw=0.7)
    ax.text((xmin + xmax) / 2, -0.78,
            "CA3 and CA2 share rs955007 in overlapping MAGMA windows; conditional independence was not tested.",
            ha="center", fontsize=5.8)
    ax.annotate("rs955007", (float(locus.loc[locus.gene_symbol.eq("CA3"), "top_snp_position_grch37_1kg_eur"].iloc[0]), 1),
                xytext=(8, 10), textcoords="offset points", fontsize=5.8, va="bottom",
                arrowprops=dict(arrowstyle="-", lw=0.45, color="#555555"))
    ax.set_xlim(xmin - 3000, xmax + 3000); ax.set_ylim(-1.00, 2.60); ax.set_yticks([])
    ax.ticklabel_format(style="plain", axis="x"); ax.set_xlabel("GRCh37 position on chromosome 8 (bp)")
    clean_axis(ax)
    ax.legend(handles=[Line2D([], [], color=c["genetics"], lw=3.5, label="Gene body"),
                       Line2D([], [], color=c["light"], lw=8, label="MAGMA window"),
                       Line2D([], [], marker="v", ls="", mfc=c["support"], mec="black", label="Top SNP in window")],
              frameon=False, ncol=3, loc="upper right")

    regulatory_grid = gs[3].subgridspec(1, 2, width_ratios=[0.96, 1.04], wspace=0.43)
    ax = fig.add_subplot(regulatory_grid[0]); _panel(ax, "D", "CA3 cis-eQTL evidence", -0.22, 1.04)
    eqtl = d["ca_eqtl"].copy()
    eqtl = eqtl[eqtl.gene.eq("CA3")].copy()
    eqtl["minimum_nominal_p"] = pd.to_numeric(eqtl.minimum_nominal_p, errors="coerce")
    eqtl["bonferroni_p_global"] = (eqtl["minimum_nominal_p"] * 18888).clip(upper=1.0)
    eqtl = eqtl.sort_values("minimum_nominal_p").reset_index(drop=True)
    eqtl["x"] = neglog10(eqtl.minimum_nominal_p)
    passed = eqtl.bonferroni_p_global < 0.05
    ax.scatter(eqtl.loc[~passed, "x"], np.arange(len(eqtl))[~passed], s=34, fc=c["light"], ec="black", lw=0.45)
    ax.scatter(eqtl.loc[passed, "x"], np.arange(len(eqtl))[passed], s=54, fc=c["genetics"], ec="black", lw=0.8)
    threshold = -math.log10(0.05 / 18888.0)
    ax.axvline(threshold, color="#666666", ls="--", lw=0.9)
    labels = ["Oligodendrocytes" if value == "Oligodendrocytes" else value.replace("_", " ") for value in eqtl.cell_type]
    ax.set_yticks(np.arange(len(eqtl)), labels)
    ax.set_xlabel("−log10(cis-eQTL P)")
    ax.set_ylabel("Cell type")
    ax.set_xlim(0, max(8.0, float(eqtl.x.max()) + 0.7))
    ax.set_ylim(len(eqtl) - 0.5, -0.5)
    clean_axis(ax, grid=True)
    ax.text(0.98, 0.05, "Filled: Bonferroni P < 0.05 across 18,888 tests", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.8)

    ax = fig.add_subplot(regulatory_grid[1]); _panel(ax, "E", "Prior sensitivity of CA3 colocalization", -0.19, 1.04)
    coloc = d["ca3_coloc_grid"].copy()
    coloc["p12"] = pd.to_numeric(coloc.p12, errors="coerce")
    coloc["PPH4"] = pd.to_numeric(coloc.PPH4, errors="coerce")
    colors = {"Granule": c["genetics"], "Oligodendrocytes": c["purkinje"]}
    for cell_type, q in coloc.groupby("cell_type", sort=False):
        q = q.sort_values("p12")
        ax.plot(q.p12, q.PPH4, marker="o", ms=4.4, lw=1.5, color=colors[cell_type], label=cell_type)
        default = q[np.isclose(q.p12, 1e-5)]
        if not default.empty:
            row = default.iloc[0]
            ax.annotate(f"{row.PPH4:.3f}", (row.p12, row.PPH4), xytext=(4, 5),
                        textcoords="offset points", fontsize=5.8, color=colors[cell_type])
    ax.axvline(1e-5, color="#666666", ls="--", lw=0.9)
    ax.set_xscale("log")
    ax.set_xticks([1e-6, 5e-6, 1e-5, 5e-5], ["10^-6", "5×10^-6", "10^-5", "5×10^-5"])
    ax.set_ylim(0, 0.86)
    ax.set_xlabel("Prior probability p12")
    ax.set_ylabel("PPH4")
    clean_axis(ax)
    ax.legend(frameon=False, loc="upper left")
    ax.text(0.98, 0.05, "Dashed line: default p12", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.8)
    return fig


def _volcano(ax, host: pd.DataFrame, c, xcol: str, pcol: str, sigcol: str, title: str,
             labels: list[str], label_offsets: dict[str, tuple[float, float]] | None = None):
    x = pd.to_numeric(host[xcol], errors="coerce")
    p = pd.to_numeric(host[pcol], errors="coerce")
    valid = x.notna() & p.notna()
    sig = as_bool(host[sigcol]) & valid
    y = -np.log10(p.clip(lower=1e-300))
    ax.scatter(x[valid & ~sig], y[valid & ~sig], s=14, fc=c["neutral"], ec="none", alpha=0.66)
    ax.scatter(x[sig], y[sig], s=31, c=np.where(x[sig] >= 0, c["positive"], c["negative"]), ec="black", lw=0.35)
    ax.axhline(-math.log10(0.05), ls="--", lw=0.8, color="#666666")
    ax.axvline(0, lw=0.7, color="#999999")
    ymax = float(y[valid].max()) if valid.any() else 1.0
    ax.set_ylim(0, max(1.0, ymax * 1.20))
    positions = {gene: idx for idx, gene in enumerate(labels)}
    for row in host[sig].itertuples():
        gene = row.gene_symbol
        if gene not in positions: continue
        xx = float(getattr(row, xcol)); yy = -math.log10(float(getattr(row, pcol)))
        idx = positions[gene]
        default_dx = 4 if xx >= 0 else -4
        dx, dy = (label_offsets or {}).get(gene, (default_dx, 5 + (idx % 3) * 5))
        ha = "left" if dx >= 0 else "right"
        ax.annotate(gene, (xx, yy), xytext=(dx, dy), textcoords="offset points", fontsize=6.2,
                    ha=ha, va="bottom" if dy >= 0 else "top",
                    arrowprops=dict(arrowstyle="-", color="#777777", lw=0.35))
    ax.set_xlabel("log2 fold change (ET versus control)"); ax.set_ylabel("−log10(nominal P)")
    ax.set_title(title, loc="left"); clean_axis(ax)
    ax.text(0.02, 0.03, f"Predicted targets with estimates: n = {int(valid.sum())}", transform=ax.transAxes, fontsize=6.2)


def figure5(d, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 260 * MM))
    gs = fig.add_gridspec(4, 12, height_ratios=[1.28, 0.92, 1.08, 1.12], hspace=0.58, wspace=0.90,
                          left=0.16, right=0.975, top=0.972, bottom=0.082)
    host = d["host"].copy()
    target_expression = d["target_expression"].copy()
    target_expression["pvalue"] = pd.to_numeric(target_expression.pvalue, errors="coerce")
    target_expression["effect_et_vs_control"] = pd.to_numeric(target_expression.effect_et_vs_control, errors="coerce")
    target_expression["bh_adjusted_p"] = pd.to_numeric(target_expression.bh_adjusted_p, errors="coerce")
    bulk = target_expression[target_expression.dataset.str.startswith("Whole") & as_bool(target_expression.included_in_expression_analysis)].copy()
    purk = target_expression[target_expression.dataset.str.startswith("Purkinje") & as_bool(target_expression.included_in_expression_analysis)].copy()
    purk_genes = purk.loc[purk.pvalue < 0.05, "gene_symbol"].tolist()

    cell_order = ["Purkinje cell", "Bergmann glial cell", "granule cell", "Golgi cell", "molecular layer interneuron",
                  "interneuron", "unipolar brush cell", "astrocyte", "microglia", "oligodendrocyte",
                  "oligodendrocyte precursor cell", "endothelial cell", "pericyte"]
    cell_abbrev = {cell_type: _cell_abbreviation(cell_type) for cell_type in cell_order}
    cell_abbrev["interneuron"] = "IN"
    cell_palette = [c["purkinje"], c["cell"], c["microbe"], c["metabolite"], c["genetics"], c["bulk"],
                    c["support"], c["positive"], c["light"], "#7A6DA8", "#74A58C", "#D69A64", "#8E8E8E"]
    cell_color = dict(zip(cell_order, cell_palette))

    ax = fig.add_subplot(gs[0, 0:6]); _panel(ax, "A", "Cerebellar cell-type atlas", -0.27, 1.02)
    umap = d["umap"].copy()
    for cell_type in cell_order:
        q = umap[umap.cell_type.eq(cell_type)]
        ax.scatter(q.UMAP1, q.UMAP2, s=0.35, c=cell_color[cell_type], alpha=0.48,
                   linewidths=0, rasterized=True)
        center = q[["UMAP1", "UMAP2"]].median()
        ax.text(center.UMAP1, center.UMAP2, cell_abbrev[cell_type], ha="center", va="center",
                fontsize=5.7, bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.76})
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[:].set_visible(False)
    ax.text(0.02, 0.02, "100,000-nucleus portal visualization subsample", transform=ax.transAxes,
            fontsize=5.4, ha="left", va="bottom")

    selected = ["CA3", "ADCY10", "CA7", "TET3", "IGF1R", "SIRT1", "KDM6B", "PTGES", "PTGS1", "FYN", "PDE4B", "NR3C1"]
    ax = fig.add_subplot(gs[0, 6:12]); _panel(ax, "B", "Cell-type expression", -0.24, 1.02)
    cell = d["cell"].copy()
    cq = cell[cell.gene_symbol.isin(selected)].copy()
    vmax = max(0.01, cq.mean_expression.quantile(0.98))
    expr_cmap = mpl.colors.LinearSegmentedColormap.from_list("expr", ["#F7F7F7", c["cell"], c["genetics"]])
    highest = host.set_index("gene_symbol").highest_mean_expression_cell_type.to_dict()
    for row in cq.itertuples():
        if row.annotation_value not in cell_order:
            continue
        x = cell_order.index(row.annotation_value); y = selected.index(row.gene_symbol)
        ax.scatter(x, y, s=max(5, row.fraction_expressing * 225), c=[row.mean_expression], cmap=expr_cmap,
                   vmin=0, vmax=vmax, ec="black" if row.annotation_value == highest.get(row.gene_symbol) else "none", lw=0.45)
    ax.set_xlim(-0.6, 12.6); ax.set_ylim(len(selected) - 0.5, -0.5)
    ax.set_xticks(range(13), [cell_abbrev[value] for value in cell_order], rotation=55, ha="right", fontsize=5.8)
    ax.set_yticks(range(len(selected)), selected, fontsize=5.8)
    ax.tick_params(length=0); ax.spines[:].set_visible(False)
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, vmax), cmap=expr_cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.032, pad=0.025)
    cb.set_label("Mean expression", fontsize=5.8); cb.ax.tick_params(labelsize=5.2)
    size_handles = [ax.scatter([], [], s=value * 225, fc=c["neutral"], ec="black", lw=0.3, label=f"{value:.1f}")
                    for value in [0.1, 0.3, 0.5]]
    ax.legend(handles=size_handles, title="Fraction expressing", title_fontsize=5.8, fontsize=5.4,
              frameon=False, ncol=3, loc="upper right", bbox_to_anchor=(1.0, 1.15),
              columnspacing=0.7, handletextpad=0.25)

    ax = fig.add_subplot(gs[1, 0:5]); _panel(ax, "C", "Purkinje-cell PCA", -0.34, 1.03)
    pca = d["pca"].copy()
    variance = d["pca_variance"].set_index("component").variance_percent
    pca_styles = {"Control": (c["neutral"], "o"), "ET": (c["purkinje"], "^")}
    for condition, (color, marker) in pca_styles.items():
        q = pca[pca.condition.eq(condition)]
        ax.scatter(q.PC1, q.PC2, s=30, marker=marker, fc=color, ec="black", lw=0.45,
                   alpha=0.90, label=f"{condition} (n={len(q)})")
    ax.axhline(0, color="#E4E4E4", lw=0.55, zorder=0)
    ax.axvline(0, color="#E4E4E4", lw=0.55, zorder=0)
    ax.set_xlabel(f"PC1 ({variance.loc['PC1']:.1f}%)")
    ax.set_ylabel(f"PC2 ({variance.loc['PC2']:.1f}%)")
    ax.margins(x=0.12, y=0.08)
    clean_axis(ax)
    ax.legend(frameon=False, loc="best", fontsize=5.5)

    ax = fig.add_subplot(gs[1, 5:12]); _panel(ax, "D", "Purkinje-cell target expression", -0.18, 1.03)
    heat = d["purkinje_heatmap"].set_index("gene_symbol")
    sample_order = d["purkinje_heatmap_order"].sort_values("display_order")
    heat = heat.loc[:, sample_order.sample_id]
    heat_cmap = diverging_cmap(c)
    bound = max(2.0, float(np.nanquantile(np.abs(heat.to_numpy(dtype=float)), 0.98)))
    image = ax.imshow(heat.to_numpy(dtype=float), aspect="auto", interpolation="nearest",
                      cmap=heat_cmap, vmin=-bound, vmax=bound)
    ax.set_yticks(range(len(heat)), heat.index, fontsize=5.8)
    ax.set_xticks([])
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    strip = ax.inset_axes([0.0, 1.01, 1.0, 0.06])
    condition_codes = sample_order.condition.map({"Control": 0, "ET": 1}).to_numpy()[None, :]
    condition_cmap = mpl.colors.ListedColormap([c["neutral"], c["purkinje"]])
    strip.imshow(condition_codes, aspect="auto", interpolation="nearest", cmap=condition_cmap, vmin=0, vmax=1)
    strip.set_axis_off()
    cb = fig.colorbar(image, ax=ax, fraction=0.032, pad=0.025)
    cb.set_label("Row z score", fontsize=5.8); cb.ax.tick_params(labelsize=5.2)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color=c["neutral"], label="Control"),
                       Line2D([], [], marker="s", ls="", color=c["purkinje"], label="ET")],
              frameon=False, ncol=2, fontsize=5.4, loc="upper right", bbox_to_anchor=(1.0, 1.18))

    ax = fig.add_subplot(gs[2, 0:6]); panel_label(ax, "E", -0.27, 1.03)
    _volcano(ax, purk, c, "effect_et_vs_control", "pvalue", "nominal_p_lt_0_05",
             "Purkinje-cell target expression", purk_genes,
             {"CA7": (-3, 14), "ADCY10": (-6, 12), "TET3": (-18, 18), "KDM6B": (20, 2),
              "SIRT1": (-20, 8), "IGF1R": (10, -15), "PTGES": (-28, 10), "PTGS1": (12, 12)})

    ax = fig.add_subplot(gs[2, 7:12]); _panel(ax, "F", "Rotation gene-set tests", -0.31, 1.03)
    mroast = d["mroast"].copy(); fry = d["fry"].copy(); weighted = d["weighted_roast"].copy()
    rows = []
    for dataset, short_name, color in [("Whole cerebellum (GSE134878)", "Whole cerebellum", c["bulk"]),
                                       ("Purkinje cells (GSE197345)", "Purkinje cells", c["purkinje"])]:
        mr = mroast[(mroast.dataset == dataset) & (mroast.gene_set == "Combined predicted targets")].iloc[0]
        fr = fry[(fry.dataset == dataset) & (fry.gene_set == "Combined predicted targets")].iloc[0]
        wr = weighted[(weighted.dataset == dataset) & (weighted.gene_set == "SwissTargetPrediction targets")].iloc[0]
        rows.extend([(short_name, "mroast", float(mr["PValue.Mixed"]), int(mr.NGenes), color, "o"),
                     (short_name, "fry", float(fr["PValue.Mixed"]), int(fr.NGenes), color, "D"),
                     (short_name, "probability-weighted mroast", float(wr.PValue), int(wr.NGenes), color, "^")])
    summary = pd.DataFrame(rows, columns=["dataset", "method", "p", "n", "color", "marker"])
    y = np.arange(len(summary))
    for i, row in summary.iterrows():
        ax.scatter(-math.log10(row.p), i, s=46, marker=row.marker, fc=row.color, ec="black", lw=0.55, zorder=3)
        ax.text(-math.log10(row.p) + 0.035, i, f"P={row.p:.3f}", va="center", fontsize=5.3)
    ax.axvline(-math.log10(0.05), color="#666666", ls="--", lw=0.8)
    display_labels = ["Whole | mroast\n(n=64)", "Whole | fry\n(n=64)", "Whole | weighted mroast\n(n=62)",
                      "Purkinje | mroast\n(n=62)", "Purkinje | fry\n(n=62)", "Purkinje | weighted mroast\n(n=60)"]
    ax.set_yticks(y, display_labels, fontsize=5.0)
    ax.set_ylim(len(summary) - 0.5, -0.5)
    ax.set_xlim(0, max(1.75, float((-np.log10(summary.p)).max()) + 0.34))
    ax.set_xlabel("−log10(P)")
    clean_axis(ax, grid=True)
    ax.text(0.99, 0.97, "Nominal P=0.05", transform=ax.transAxes, ha="right", va="top", fontsize=5.3)

    ax = fig.add_subplot(gs[3, :]); _panel(ax, "G", "Gene association and cerebellar transcription", -0.105, 1.03, pad=28)
    q = host.set_index("gene_symbol").loc[selected].reset_index()
    bulk_map = bulk.set_index("gene_symbol")
    purk_map = purk.set_index("gene_symbol")
    columns = ["MAGMA", "Whole cerebellum", "Purkinje cells", "CA3 cell-type\ncis-eQTL", "Highest-expression\ncell type"]
    ax.set_xlim(-0.5, len(columns) - 0.5); ax.set_ylim(len(q) - 0.5, -0.5)
    cmap = diverging_cmap(c)
    for i, row in q.iterrows():
        for j in range(len(columns)):
            ax.add_patch(Rectangle((j - 0.48, i - 0.43), 0.96, 0.86, fc="#F6F6F6", ec="white"))
        mp = pd.to_numeric(pd.Series([row.magma_pvalue]), errors="coerce").iloc[0]
        if np.isfinite(mp):
            ax.scatter(0, i, s=18 + min(90, -math.log10(mp) * 20), fc=c["genetics"], ec="black", lw=0.4)
        for j, data_map in [(1, bulk_map), (2, purk_map)]:
            if row.gene_symbol in data_map.index:
                expr = data_map.loc[row.gene_symbol]
                fc = float(expr.effect_et_vs_control); pv = float(expr.pvalue)
                color = cmap(np.clip((fc + 1.5) / 3.0, 0, 1))
                ax.scatter(j, i, s=18 + min(90, -math.log10(max(pv, 1e-300)) * 20), fc=color, ec="black", lw=0.3)
        if row.gene_symbol == "CA3":
            ax.scatter(3, i, s=58, marker="D", fc=c["support"], ec="black", lw=0.8)
        ct = str(row.highest_mean_expression_cell_type)
        ax.scatter(4, i, s=76, fc=cell_color.get(ct, c["neutral"]), ec="black", lw=0.4)
        ax.text(4, i, cell_abbrev.get(ct, _cell_abbreviation(ct)), ha="center", va="center", fontsize=5.8)
    ax.set_xticks(range(len(columns)), columns)
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, length=0)
    ax.set_yticks(range(len(q)), q.gene_symbol); ax.tick_params(axis="y", length=0); ax.spines[:].set_visible(False)
    ax.text(0.0, -0.12, "Point size: −log10(P); expression color: log2 fold change; MAGMA P values are unsigned.",
            transform=ax.transAxes, fontsize=5.8, va="top")
    fig.text(0.16, 0.034,
             "AS, astrocyte; BG, Bergmann glia; EC, endothelial cell; GC, granule cell; GOL, Golgi cell; IN, interneuron; MG, microglia; MLI, molecular layer interneuron;",
             fontsize=5.3, va="bottom")
    fig.text(0.16, 0.018,
             "OL, oligodendrocyte; OPC, oligodendrocyte precursor cell; PC, Purkinje cell; PER, pericyte; UBC, unipolar brush cell.",
             fontsize=5.3, va="bottom")
    fig.text(0.16, 0.004,
        "SCP3177 UMAP coordinates and cell labels were obtained from the Single Cell Portal; cell-type expression is descriptive and is not an ET-versus-control comparison.",
             fontsize=5.3, va="bottom")
    return fig


def _focal_iv_data(d, taxon_id: str) -> pd.DataFrame:
    q = d["ivs"].copy()
    q = q[(q.taxon_id == taxon_id) & as_bool(q.included_in_mr)].copy()
    numeric = ["beta_exposure", "se_exposure", "beta_outcome_aligned", "se_outcome"]
    for col in numeric: q[col] = pd.to_numeric(q[col], errors="coerce")
    complete_rows = q[numeric].notna().all(axis=1)
    q = q.loc[complete_rows].copy()
    q["wald_beta"] = q.beta_outcome_aligned / q.beta_exposure
    q["wald_se"] = q.se_outcome / q.beta_exposure.abs()
    q["wald_low"] = q.wald_beta - 1.96 * q.wald_se
    q["wald_high"] = q.wald_beta + 1.96 * q.wald_se
    q["precision"] = 1.0 / q.wald_se
    return q.sort_values("wald_beta").reset_index(drop=True)


def supplementary_figure_s1(d, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 195 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.10, 0.92], hspace=0.34, wspace=0.38,
                          left=0.14, right=0.985, top=0.955, bottom=0.105)
    focal = [("Faecalibacterium", "genus.Faecalibacterium.id.2057", c["microbe"]),
             ("Flavonifractor", "genus.Flavonifractor.id.2059", c["metabolite"])]

    ax = fig.add_subplot(gs[0, 0]); _panel(ax, "A", "Single-variant estimates (9 and 4 instruments)", -0.27, 1.02)
    rows = []
    for taxon, taxon_id, color in focal:
        q = _focal_iv_data(d, taxon_id)
        for row in q.itertuples(): rows.append((taxon, row.variant_id, row.wald_beta, row.wald_low, row.wald_high, color))
    y = np.arange(len(rows))
    for yi, (taxon, rsid, beta, lo, hi, color) in enumerate(rows):
        ax.errorbar(beta, yi, xerr=[[beta - lo], [hi - beta]], fmt="o", ms=3.8, color=color, ecolor="#555555", lw=0.7, capsize=1.8)
    ax.axvline(0, color="#777777", ls="--", lw=0.8)
    ax.set_yticks(y, [f"{'Fae.' if r[0] == 'Faecalibacterium' else 'Flav.'} | {r[1]}" for r in rows], fontsize=5.6); ax.invert_yaxis()
    ax.set_xlabel("Wald-ratio estimate (log odds)"); clean_axis(ax, grid=True)

    axbg = fig.add_subplot(gs[0, 1]); axbg.set_axis_off(); _panel(axbg, "B", "Funnel plots", -0.18, 1.02)
    for k, (taxon, taxon_id, color) in enumerate(focal):
        inset = axbg.inset_axes([0.04, 0.06 + (1 - k) * 0.51, 0.92, 0.32])
        q = _focal_iv_data(d, taxon_id)
        ivw = float(d["mr"].set_index("taxon_id").loc[taxon_id, "beta"])
        inset.scatter(q.wald_beta, q.precision, s=20, fc=color, ec="black", lw=0.25)
        inset.axvline(ivw, color="#555555", lw=0.8, ls="--")
        inset.set_title(taxon, loc="left", fontsize=6.4, fontstyle="italic", pad=1)
        if k == 1: inset.set_xlabel("Wald-ratio estimate", fontsize=5.6)
        inset.set_ylabel("Precision (1/SE)", fontsize=5.6)
        inset.tick_params(labelsize=5.2); clean_axis(inset)

    ax = fig.add_subplot(gs[1, 0]); _panel(ax, "C", "Sensitivity estimators: six additional sets", -0.27, 1.02)
    mr = d["mr"].copy()
    other = mr[~mr.taxon_name.isin(["Faecalibacterium", "Flavonifractor"])].copy()
    method_order = ["Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode"]
    method_colors = [c["microbe"], c["genetics"], c["bulk"], c["purkinje"]]
    method_markers = ["o", "s", "D", "^"]; offsets = [-0.18, -0.06, 0.06, 0.18]
    for i, row in other.reset_index(drop=True).iterrows():
        q = d["all_methods"][d["all_methods"].taxon_id.eq(row.taxon_id)].set_index("method")
        for j, method in enumerate(method_order):
            rr = q.loc[method]
            ax.errorbar(rr.odds_ratio, i + offsets[j], xerr=[[rr.odds_ratio - rr.or_ci_lower], [rr.or_ci_upper - rr.odds_ratio]],
                        fmt=method_markers[j], ms=3.4, color=method_colors[j], ecolor="#666666", lw=0.6, capsize=1.4)
    ax.axvline(1, color="#777777", ls="--", lw=0.8)
    short_other = {"LachnospiraceaeUCG001": "Lachnosp. UCG-001", "RuminococcaceaeUCG011": "Ruminococc. UCG-011",
                   "Methanobrevibacter": "Methanobrev.", "Butyrivibrio": "Butyrivibrio", "Methanobacteria": "Methanobacteria", "Prevotella7": "Prevotella 7"}
    other_labels = [short_other.get(x, x) for x in other.taxon_name]
    ax.set_yticks(range(len(other)), other_labels, fontsize=5.6); ax.invert_yaxis()
    for tick, label in zip(ax.get_yticklabels(), other_labels):
        tick.set_fontstyle("normal" if label == "Methanobacteria" else "italic")
    ax.set_xscale("log"); ax.xaxis.set_minor_formatter(NullFormatter()); ax.set_xlabel("Odds ratio (95% CI)")
    clean_axis(ax, grid=True)
    fig.legend(handles=[Line2D([], [], marker=m, ls="", color=co, label=lab.replace("Inverse variance weighted", "IVW"))
                        for m, co, lab in zip(method_markers, method_colors, method_order)],
               frameon=False, fontsize=5.6, ncol=2, loc="lower left", bbox_to_anchor=(0.145, 0.012))

    ax = fig.add_subplot(gs[1, 1]); _panel(ax, "D", "Methanobacteria lineage (10 shared SNPs)", -0.18, 1.02)
    lineage = d["mr10"][d["mr10"]["rank"].isin(["class", "order", "family"])].copy()
    lineage["rank"] = pd.Categorical(lineage["rank"], categories=["class", "order", "family"], ordered=True)
    lineage = lineage.sort_values("rank").reset_index(drop=True)
    yv = np.arange(len(lineage))
    ax.errorbar(lineage.odds_ratio, yv, xerr=[lineage.odds_ratio - lineage.or_ci_lower, lineage.or_ci_upper - lineage.odds_ratio],
                fmt="o", ms=5, color=c["microbe"], ecolor="#555555", capsize=2, lw=0.8)
    ax.axvline(1, color="#777777", ls="--", lw=0.8)
    ax.set_yticks(yv, [str(rank).capitalize() for rank in lineage["rank"]], fontsize=5.8)
    ax.invert_yaxis(); ax.set_xlabel("Odds ratio for essential tremor (95% CI)"); clean_axis(ax, grid=True)
    ax.text(0.02, 0.84, "Shared 10-SNP instrument set", transform=ax.transAxes, fontsize=5.8)
    ax.text(0.98, 0.70, "Identical IVW estimates", transform=ax.transAxes, fontsize=5.8, ha="right")
    return fig


def supplementary_figure_s2(d, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 280 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.82, 1.00, 1.72], hspace=0.46, wspace=0.36,
                          left=0.075, right=0.985, top=0.965, bottom=0.040)
    identity = d["identity"]

    ax = fig.add_subplot(gs[0, 0]); ax.set_axis_off(); _panel(ax, "A", "Compound identity and nomenclature", -0.15, 1.02)
    fields = [
        ("Preferred name", identity["compound"]),
        ("Common synonym", identity["preferred_synonym"]),
        ("PubChem CID", str(identity["pubchem_cid"])),
        ("HMDB", identity["hmdb"]),
        ("ChEBI", identity["chebi"]),
        ("Molecular formula", identity["molecular_formula"]),
        ("Molecular weight", f"{identity['molecular_weight']} g mol−1"),
        ("Canonical SMILES", identity["canonical_smiles"]),
    ]
    for i, (label, value) in enumerate(fields):
        y = 0.88 - i * 0.105
        ax.add_patch(Rectangle((0.0, y - 0.043), 1.0, 0.086, transform=ax.transAxes,
                               fc="#F7F7F7" if i % 2 == 0 else "white", ec="none"))
        ax.text(0.02, y, label, transform=ax.transAxes, fontsize=5.8, va="center")
        ax.text(0.36, y, value, transform=ax.transAxes, fontsize=5.5 if label == "Canonical SMILES" else 5.8,
                va="center", wrap=True)
    ax.text(0.02, 0.01, "Additional searchable synonym: 5-(3,4-dihydroxyphenyl)valeric acid",
            transform=ax.transAxes, fontsize=5.5)

    axbg = fig.add_subplot(gs[0, 1]); axbg.set_axis_off(); _panel(axbg, "B", "Open-chain acid and γ-valerolactone", -0.16, 1.02)
    left = axbg.inset_axes([0.03, 0.26, 0.44, 0.66]); right = axbg.inset_axes([0.53, 0.26, 0.44, 0.66])
    _structure(left, ASSETS / "PubChem_CID49831816.png", "Open-chain acid", "CID 49831816 | C11H14O4")
    _structure(right, ASSETS / "PubChem_CID45093080_valerolactone.png", "γ-Valerolactone", "CID 45093080 | C11H12O4")
    axbg.text(0.50, 0.54, "≠", transform=axbg.transAxes, ha="center", va="center", fontsize=14, fontweight="bold")
    axbg.text(0.50, 0.075, "Chemically distinct", transform=axbg.transAxes, ha="center", fontsize=6.0)
    axbg.text(0.50, 0.025, "Different molecular formulae, ring states and database identifiers",
              transform=axbg.transAxes, ha="center", fontsize=5.6)

    ax = fig.add_subplot(gs[1, 0]); _panel(ax, "C", "SwissTargetPrediction: all 98 unique genes", -0.15, 1.02)
    swiss = d["swiss"].copy(); swiss["rank"] = pd.to_numeric(swiss["rank"]); swiss["probability"] = pd.to_numeric(swiss["probability"])
    swiss = swiss.sort_values(["rank", "probability"], ascending=[True, False]).drop_duplicates("gene_symbol")
    ax.plot(swiss["rank"], swiss["probability"], color=c["support"], lw=1.0)
    ax.scatter(swiss["rank"], swiss["probability"], s=15, fc=c["support"], ec="black", lw=0.25)
    top8 = swiss.nsmallest(8, "rank").gene_symbol.tolist()
    top_lines = [" | ".join(top8[i:i + 4]) for i in range(0, len(top8), 4)]
    ax.text(0.77, 0.96, "Top-ranked genes\n" + "\n".join(top_lines), transform=ax.transAxes,
            ha="right", va="top", fontsize=5.5,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#D0D0D0", lw=0.4, alpha=0.92))
    for gene in ["CA3", "IGF1R", "SIRT1"]:
        q = swiss[swiss.gene_symbol.eq(gene)]
        if len(q):
            row = q.iloc[0]
            ax.scatter(row["rank"], row["probability"], s=34, fc="white", ec=c["metabolite"], lw=1.2, zorder=4)
            ax.annotate(gene, (row["rank"], row["probability"]), xytext=(3, 6), textcoords="offset points", fontsize=5.7)
    ax.set_xlim(0, 102); ax.set_xlabel("Rank"); ax.set_ylabel("Probability"); clean_axis(ax)

    ax = fig.add_subplot(gs[1, 1]); _panel(ax, "D", "Similarity Ensemble Approach: all seven genes", -0.16, 1.02)
    sea = d["sea"].copy().sort_values("z_score", ascending=True)
    y = np.arange(len(sea)); ax.barh(y, sea.z_score, color=c["genetics"], alpha=0.80)
    ax.set_yticks(y, sea.gene_symbol, fontsize=5.8)
    ax.set_xlabel("SEA z-score")
    for yi, row in enumerate(sea.itertuples()):
        ax.text(row.z_score - 0.12, yi, f"T={row.tanimoto_score:.2f}",
                va="center", ha="right", fontsize=5.5)
    ax.set_xlim(0, max(sea.z_score) + 1.0); clean_axis(ax, grid=True)

    axbg = fig.add_subplot(gs[2, :]); axbg.set_axis_off()
    panel_label(axbg, "E", -0.075, 1.00)
    axbg.text(0.0, 1.00, "Target-source membership and MAGMA coverage",
              transform=axbg.transAxes, fontsize=9, fontweight="bold", va="bottom")
    mapping = d["union"][["gene_symbol", "in_swiss", "in_sea"]].copy()
    mapping = mapping.merge(d["magma"][["gene_symbol", "included_in_autosomal_magma"]], on="gene_symbol", how="left").sort_values("gene_symbol").reset_index(drop=True)
    columns = ["Swiss", "SEA", "MAGMA"]
    blocks = [mapping.iloc[i:i + 35].copy() for i in range(0, len(mapping), 35)]
    for bidx, block in enumerate(blocks):
        inset = axbg.inset_axes([0.01 + bidx * 0.335, 0.09, 0.32, 0.73])
        inset.set_xlim(-2.4, len(columns) - 0.5); inset.set_ylim(len(block) - 0.5, -0.5)
        for i, row in block.reset_index(drop=True).iterrows():
            inset.text(-0.60, i, row.gene_symbol, ha="right", va="center", fontsize=5.2)
            vals = [bool(row.in_swiss), bool(row.in_sea), bool(row.included_in_autosomal_magma)]
            for j, val in enumerate(vals):
                inset.add_patch(Rectangle((j - 0.48, i - 0.43), 0.96, 0.86, fc="#F6F6F6", ec="white"))
                if val:
                    inset.scatter(j, i, s=15, fc=[c["support"], c["genetics"], c["genetics"]][j],
                                  ec="black", lw=0.22)
                elif j == 2 and not bool(row.included_in_autosomal_magma):
                    inset.text(j, i, "X", ha="center", va="center", fontsize=5.5, color="#777777")
        inset.set_yticks([])
        inset.set_xticks(range(len(columns)), columns, rotation=90, ha="left", fontsize=5.5)
        inset.tick_params(length=0, top=True, labeltop=True, bottom=False, labelbottom=False)
        inset.spines[:].set_visible(False)
    counts = {
        "SwissTargetPrediction": int(as_bool(mapping.in_swiss).sum()),
        "SEA": int(as_bool(mapping.in_sea).sum()),
        "MAGMA": int(as_bool(mapping.included_in_autosomal_magma).sum()),
    }
    axbg.text(0.01, 0.015,
              "Filled circles indicate target-source membership or MAGMA coverage; X, not tested in autosomal MAGMA.",
              transform=axbg.transAxes, ha="left", fontsize=5.7)
    axbg.text(0.99, 0.075, " | ".join(f"{k}: {v}/105" for k, v in counts.items()), transform=axbg.transAxes,
              ha="right", fontsize=5.7)
    return fig


def supplementary_figure_s3(d, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 280 * MM))
    gs = fig.add_gridspec(5, 1, height_ratios=[0.56, 1.10, 0.72, 0.72, 1.50], hspace=0.42,
                          left=0.095, right=0.985, top=0.970, bottom=0.042)
    magma = d["magma"].copy()
    tested = magma[as_bool(magma.included_in_autosomal_magma)].copy()
    tested["magma_pvalue"] = pd.to_numeric(tested.magma_pvalue, errors="coerce")

    axbg = fig.add_subplot(gs[0]); axbg.set_axis_off(); _panel(axbg, "A", "Q–Q plot of MAGMA gene P values", -0.085, 1.02)
    ax = axbg.inset_axes([0.05, 0.03, 0.36, 0.85])
    valid_magma = tested.magma_pvalue.notna()
    p = np.sort(tested.loc[valid_magma, "magma_pvalue"].to_numpy())
    n = len(p); expected = -np.log10((np.arange(1, n + 1) - 0.5) / n); observed = -np.log10(p)
    ax.scatter(expected, observed, s=18, fc=c["genetics"], ec="black", lw=0.25)
    limit = max(expected.max(), observed.max()) + 0.35
    ax.plot([0, limit], [0, limit], color="#777777", ls="--", lw=0.8)
    ax.set_xlim(0, limit); ax.set_ylim(0, limit); ax.set_xlabel("Expected −log10 P"); ax.set_ylabel("Observed −log10 P")
    clean_axis(ax)
    axbg.text(0.50, 0.69, "97 autosomal predicted targets", transform=axbg.transAxes, fontsize=7.0)
    axbg.text(0.50, 0.47, "11 genes with nominal P < 0.05", transform=axbg.transAxes, fontsize=6.5)
    axbg.text(0.50, 0.27, "3 genes with BH-adjusted P < 0.05", transform=axbg.transAxes, fontsize=6.5)
    axbg.text(0.50, 0.08, "Competitive gene-set P = 0.50068", transform=axbg.transAxes, fontsize=6.5)

    axbg = fig.add_subplot(gs[1]); axbg.set_axis_off(); _panel(axbg, "B", "MAGMA gene P values for 97 predicted targets", -0.085, 1.01)
    ranked = tested.sort_values(["magma_chr", "magma_start_build37"]).reset_index(drop=True)
    blocks = [ranked.iloc[i:i + 33].copy() for i in range(0, len(ranked), 33)]
    for bidx, block in enumerate(blocks):
        inset = axbg.inset_axes([0.02 + bidx * 0.34, 0.03, 0.28, 0.89])
        vals = -np.log10(block.magma_pvalue.to_numpy())
        y = np.arange(len(block))
        colors = np.where(block.magma_pvalue < 0.05, c["genetics"], c["neutral"])
        inset.hlines(y, 0, vals, color=colors, lw=1.1, alpha=0.8)
        inset.scatter(vals, y, s=np.where(block.magma_pvalue < 0.05, 22, 12), fc=colors,
                      ec=np.where(as_bool(block.magma_bh_adjusted_p_lt_0_05), "black", colors),
                      lw=np.where(as_bool(block.magma_bh_adjusted_p_lt_0_05), 0.9, 0.2))
        inset.axvline(-math.log10(0.05), color="#777777", ls="--", lw=0.65)
        inset.set_yticks(y, block.gene_symbol, fontsize=5.5)
        inset.invert_yaxis(); inset.set_xlabel("−log10 P", fontsize=5.8); inset.tick_params(axis="x", labelsize=5.5)
        clean_axis(inset, grid=True)

    target_expression = d["target_expression"].copy()
    target_expression["pvalue"] = pd.to_numeric(target_expression.pvalue, errors="coerce")
    target_expression["effect_et_vs_control"] = pd.to_numeric(target_expression.effect_et_vs_control, errors="coerce")
    bulk = target_expression[target_expression.dataset.str.startswith("Whole") & as_bool(target_expression.included_in_expression_analysis)].copy()
    bulk_genes = bulk.loc[bulk.pvalue < 0.05, "gene_symbol"].tolist()
    ax = fig.add_subplot(gs[2]); panel_label(ax, "C", -0.085, 1.02)
    _volcano(ax, bulk, c, "effect_et_vs_control", "pvalue", "nominal_p_lt_0_05",
             "Whole-cerebellum target expression", bulk_genes,
             {"NR3C1": (-4, 12), "DUSP3": (4, -2), "FYN": (4, 9), "PDE4B": (4, 13), "HTR2B": (4, -5)})

    ax = fig.add_subplot(gs[3]); _panel(ax, "D", "Published and reanalyzed Purkinje-cell estimates", -0.085, 1.02)
    rec = d["published_comparison"].copy()
    x = pd.to_numeric(rec.published_log2FoldChange, errors="coerce")
    y = pd.to_numeric(rec.reanalyzed_log2FoldChange, errors="coerce")
    ok = x.notna() & y.notna()
    recovered = as_bool(rec.reanalyzed_bh_adjusted_p_lt_0_05) & ok
    ax.scatter(x[ok & ~recovered], y[ok & ~recovered], s=34, fc="white", ec=c["purkinje"], lw=1.2)
    ax.scatter(x[recovered], y[recovered], s=25, fc=c["purkinje"], ec="black", lw=0.3)
    lo = min(x[ok].min(), y[ok].min()) - 0.08; hi = max(x[ok].max(), y[ok].max()) + 0.08
    ax.plot([lo, hi], [lo, hi], color="#777777", ls="--", lw=0.8)
    label_offsets = {"SPTBN5": (-5, 4), "GNB3": (-5, -10), "LAMP5": (4, 5)}
    for row in rec[ok & rec.gene_symbol.isin(label_offsets)].itertuples():
        dx, dy = label_offsets[row.gene_symbol]
        ax.annotate(row.gene_symbol, (row.published_log2FoldChange, row.reanalyzed_log2FoldChange),
                    xytext=(dx, dy), textcoords="offset points", fontsize=5.8,
                    ha="left" if dx >= 0 else "right",
                    va="bottom" if dy >= 0 else "top")
    not_fdr_recovered = rec.loc[ok & ~recovered, "gene_symbol"].tolist()
    corr = np.corrcoef(x[ok], y[ok])[0, 1]
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_xlabel("Published log2 fold change"); ax.set_ylabel("Reanalyzed log2 fold change")
    clean_axis(ax)
    ax.text(0.03, 0.95, f"Estimates available: {int(ok.sum())}/36 | BH-adjusted P < 0.05: {int(recovered.sum())}/36 | Pearson r = {corr:.4f}",
            transform=ax.transAxes, va="top", fontsize=6.0)
    ax.text(0.03, 0.84, "BH-adjusted P ≥ 0.05 in reanalysis: " + (", ".join(not_fdr_recovered) if not_fdr_recovered else "none"),
            transform=ax.transAxes, va="top", fontsize=6.0)

    axbg = fig.add_subplot(gs[4]); axbg.set_axis_off(); _panel(axbg, "E", "Cell-type expression of 105 predicted genes", -0.085, 1.01)
    cell = d["cell"].copy()
    genes = sorted(cell.gene_symbol.unique())
    cell_order = ["Purkinje cell", "Bergmann glial cell", "granule cell", "Golgi cell", "molecular layer interneuron",
                  "interneuron", "unipolar brush cell", "astrocyte", "microglia", "oligodendrocyte",
                  "oligodendrocyte precursor cell", "endothelial cell", "pericyte"]
    labels = ["PC", "BG", "GC", "GOL", "MLI", "IN", "UBC", "AS", "MG", "OL", "OPC", "EC", "PER"]
    vmax = max(0.01, cell.mean_expression.quantile(0.98)); cmap = mpl.colors.LinearSegmentedColormap.from_list("full_expr", ["#F7F7F7", c["cell"], c["genetics"]])
    highest = d["host"].set_index("gene_symbol").highest_mean_expression_cell_type.to_dict()
    gene_blocks = [genes[i:i + 35] for i in range(0, len(genes), 35)]
    for bidx, block in enumerate(gene_blocks):
        inset = axbg.inset_axes([0.005 + bidx * 0.33, 0.04, 0.295, 0.89])
        q = cell[cell.gene_symbol.isin(block)]
        for row in q.itertuples():
            if row.annotation_value not in cell_order: continue
            xx = cell_order.index(row.annotation_value); yy = block.index(row.gene_symbol)
            inset.scatter(xx, yy, s=max(3, row.fraction_expressing * 135), c=[row.mean_expression], cmap=cmap,
                          vmin=0, vmax=vmax, ec="black" if row.annotation_value == highest.get(row.gene_symbol) else "none", lw=0.35)
        inset.set_xlim(-0.55, 12.55); inset.set_ylim(len(block) - 0.5, -0.5)
        inset.set_xticks(range(13), labels, rotation=67, ha="left", fontsize=5.5)
        inset.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, length=0)
        inset.set_yticks(range(len(block)), block, fontsize=5.5); inset.tick_params(axis="y", length=0); inset.spines[:].set_visible(False)
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, vmax), cmap=cmap); sm.set_array([])
    cax = axbg.inset_axes([0.968, 0.19, 0.012, 0.52])
    cb = fig.colorbar(sm, cax=cax); cb.set_label("Mean expression", fontsize=5.8); cb.ax.tick_params(labelsize=5.5)
    axbg.text(0.99, -0.005, "Point size: fraction expressing | Black outline: highest mean expression cell type",
              transform=axbg.transAxes, ha="right", fontsize=5.7)
    return fig


def validate_data_dimensions(d) -> dict[str, int]:
    swiss_unique = int(d["swiss"].gene_symbol.nunique())
    sea_unique = int(d["sea"].gene_symbol.nunique())
    tested = as_bool(d["magma"].included_in_autosomal_magma)
    nominal = as_bool(d["magma"].magma_nominal_p_lt_0_05)
    fdr = as_bool(d["magma"].magma_bh_adjusted_p_lt_0_05)
    recovery = as_bool(d["published_comparison"].reanalyzed_bh_adjusted_p_lt_0_05)
    counts = {
        "nonduplicate_mr_sets": len(d["mr"]),
        "nominal_mr_rows": len(d["mr10"]),
        "swiss_unique_genes": swiss_unique,
        "sea_unique_genes": sea_unique,
        "predicted_target_union": int(d["union"].gene_symbol.nunique()),
        "swiss_sea_overlap": len(set(d["swiss"].gene_symbol) & set(d["sea"].gene_symbol)),
        "magma_tested_autosomal": int(tested.sum()),
        "magma_not_tested": int((~tested).sum()),
        "magma_nominal": int((tested & nominal).sum()),
        "magma_bh_fdr": int((tested & fdr).sum()),
        "whole_cerebellum_nominal": int(as_bool(d["host"].bulk_nominal_p_lt_0_05).sum()),
        "purkinje_nominal": int(as_bool(d["host"].purkinje_nominal_p_lt_0_05).sum()),
        "celltype_genes": int(d["cell"].gene_symbol.nunique()),
        "cell_types": int(d["cell"].annotation_value.nunique()),
        "published_purkinje_fdr_genes": len(d["author36"]),
        "published_genes_local_bh_fdr_recovered": int(recovery.sum()),
        "host_evidence_genes": len(d["host23"]),
    }
    expected = {
        "nonduplicate_mr_sets": 8, "nominal_mr_rows": 10, "swiss_unique_genes": 98, "sea_unique_genes": 7,
        "predicted_target_union": 105, "swiss_sea_overlap": 0, "magma_tested_autosomal": 97, "magma_not_tested": 8,
        "magma_nominal": 11, "magma_bh_fdr": 3, "whole_cerebellum_nominal": 5, "purkinje_nominal": 8,
        "celltype_genes": 105, "cell_types": 13, "published_purkinje_fdr_genes": 36,
        "published_genes_local_bh_fdr_recovered": 35, "host_evidence_genes": 23,
    }
    mismatches = {key: {"observed": counts[key], "expected": value} for key, value in expected.items() if counts[key] != value}
    if mismatches:
        raise RuntimeError(f"Data dimension mismatch: {mismatches}")
    return counts


REFERENCE_ROWS = [
    ("Figure 1", "Full figure", "Skuladottir et al., Communications Biology (2024)", "10.1038/s42003-024-06207-4", "Figure 1",
     "Full-width study overview with quantitative study scales and separate analysis layers",
     "MiBioGen exposure GWAS, essential-tremor outcome GWAS, eight nonduplicate IVW instrument sets, metabolite association, target prediction, MAGMA and cerebellar expression",
     "Skuladottir et al. (2024); Ning et al., Nature Communications (2023), Figure 1c",
     "The graphical hierarchy is borrowed; all numerical values and associations are from the displayed analyses."),
    ("Figure 2", "A", "Dekkers et al., Nature Genetics (2026)", "10.1038/s41588-026-02512-2", "Extended Data Figure 8",
     "Forest plot with aligned effect estimates, 95% confidence intervals and a separate P-value column",
     "Eight IVW odds ratios after collapsing exact duplicate instrument sets",
     "Microbial MR summary statistics",
     "All eight associations met the prespecified nominal P < 0.05 inclusion criterion."),
    ("Figure 2", "B", "Jurgens et al., Nature Genetics (2024)", "10.1038/s41588-024-01975-5", "Figure 4b",
     "Multiple MR estimators displayed on a common odds-ratio scale",
     "IVW, MR-Egger, weighted-median and weighted-mode estimates for Faecalibacterium and Flavonifractor",
     "Microbial MR estimator and sensitivity results",
     "No pooled estimate across MR methods is shown."),
    ("Figure 2", "C", "Wen et al., Nature Communications (2024)", "10.1038/s41467-024-46796-6", "Figure 5D",
     "Paired SNP-exposure and SNP-outcome association plots with estimator slopes",
     "Instrument effects and four MR estimator slopes for Faecalibacterium and Flavonifractor",
     "Microbial MR harmonized instruments and estimator results",
     "The panel visualizes estimator geometry and does not add a separate causal claim."),
    ("Figure 2", "D", "Wen et al., Nature Communications (2024)", "10.1038/s41467-024-46796-6", "Figure 5B",
     "Leave-one-out estimates aligned against the full inverse-variance-weighted estimate",
     "Sequential SNP-exclusion estimates for Faecalibacterium and Flavonifractor",
     "Microbial MR leave-one-out results",
     "The dashed line is the full IVW estimate; no single-variant causal claim is made."),
    ("Figure 2", "E", "Skuladottir et al., Communications Biology (2024)", "10.1038/s42003-024-06207-4", "Figure 3",
     "Taxon-by-analysis matrix",
     "Eight microbial signals with sensitivity-analysis, case-control, gutMGene and pure-culture annotations",
     "Microbial MR and external microbiome studies",
     "Filled and open symbols denote the stated analysis or association; no composite score is calculated."),
    ("Figure 3", "A", "Redan et al., Journal of Agricultural and Food Chemistry (2020)", "10.1021/acs.jafc.0c05890", "Figure 1",
     "Side-by-side chemical structures with names, formulae and identifiers",
     "The open-chain acid CID 49831816 and phenyl-γ-valerolactone CID 45093080",
     "PubChem and HMDB identifiers",
     "The structures are chemically distinct and no direct conversion is asserted."),
    ("Figure 3", "A", "Dong et al., Journal of the American Chemical Society (2025)", "10.1021/jacs.4c09892", "Figure 3F",
     "Organism-aware catechin-metabolism branch",
     "A separate Flavonifractor plautii pure-culture branch with catechin ring-fission products",
     "Kutschera et al., Journal of Applied Microbiology (2011), 10.1111/j.1365-2672.2011.05025.x",
     "The pure-culture branch is not connected to CID 49831816 or the 105 predicted targets."),
    ("Figure 3", "B", "Ning et al., Nature Communications (2023)", "10.1038/s41467-023-42788-0", "Figure 1b-c",
     "Two source-specific inputs connected independently to one combined set",
     "Ninety-eight SwissTargetPrediction genes and seven SEA genes forming a 105-gene union",
     "SwissTargetPrediction and Similarity Ensemble Approach results",
     "Observed overlap was zero."),
    ("Figure 3", "C", "Daina et al., Nucleic Acids Research (2019)", "10.1093/nar/gkz382", "Figure 2",
     "Probability-ranked display of the 15 highest-ranked predicted targets",
     "The 15 highest-ranked SwissTargetPrediction genes with predicted probabilities",
     "SwissTargetPrediction results",
     "Carbonic anhydrases are highlighted as a standard protein family; no probability threshold is imposed."),
    ("Figure 3", "D", "Irwin et al., Journal of Chemical Information and Modeling (2018)", "10.1021/acs.jcim.7b00316", "Figure 3D",
     "Joint presentation of SEA prediction strength and Tanimoto similarity for individual targets",
     "SEA z-scores plotted against reported Tanimoto coefficients for seven genes, with direct labels",
     "Similarity Ensemble Approach results",
     "Both native metrics are retained; no cross-method score or unreported threshold is introduced."),
    ("Figure 4", "A", "Skuladottir et al., Communications Biology (2024)", "10.1038/s42003-024-06207-4", "Figure 2",
     "Chromosome-ordered gene-association display with a nominal threshold and labeled genes",
     "MAGMA P values for 97 autosomal predicted targets",
     "Target-restricted MAGMA analysis",
     "This is a gene-based association plot, not a variant-level Manhattan plot."),
    ("Figure 4", "B", "Skuladottir et al., Communications Biology (2024)", "10.1038/s42003-024-06207-4", "Figure 3",
     "Gene-by-analysis matrix",
     "Eleven nominal MAGMA genes with prediction, published genetics, differential-expression and cell-type annotations",
     "MAGMA, target prediction, published ET genetics and cerebellar expression data",
     "The matrix displays distinct evidence sources without a composite score."),
    ("Figure 4", "C", "Skuladottir et al., Communications Biology (2024)", "10.1038/s42003-024-06207-4", "Figure 3",
     "Locus-level gene bodies, association windows and leading variants",
     "CA1, CA2 and CA3 at the chromosome 8 locus",
     "MAGMA windows and the published ET GWAS",
     "CA3 and CA2 share rs955007 in overlapping windows; conditional independence was not tested."),
    ("Figure 4", "D", "Castonguay et al., Nature Genetics (2026)", "10.1038/s41588-026-02544-8", "Figure 6",
     "Aligned cell-population statistical summary on a common scale",
     "CA3 cis-eQTL P values across eight cerebellar cell types",
     "SCP3177 cell-type cis-eQTL summary statistics",
     "Filled points denote Bonferroni P < 0.05 across 18,888 variant-cell-type tests; association does not establish mediation."),
    ("Figure 4", "E", "Wallace, PLOS Genetics (2020)", "10.1371/journal.pgen.1008720", "Figure 5",
     "Posterior-probability sensitivity curves across p12 with the default prior marked",
     "CA3 PPH4 for Granule cells and Oligodendrocytes across four p12 values",
     "CA3 ET-GWAS and cell-type cis-eQTL colocalization results",
     "Default-prior PPH4 remained below 0.50 and the higher values were prior-sensitive."),
    ("Figure 5", "A", "Castonguay et al., Nature Genetics (2026)", "10.1038/s41588-026-02544-8", "Figure 1A",
     "Cell-type UMAP with compact labels placed within separated atlas clusters",
     "Official SCP3177 UMAP coordinates and cell-type labels for the 100,000-nucleus portal visualization subsample",
     "Single Cell Portal study SCP3177",
     "The embedding was obtained from the Single Cell Portal and is descriptive; it was not recomputed and is not an ET-versus-control separation analysis."),
    ("Figure 5", "B", "Castonguay et al., Nature Genetics (2026)", "10.1038/s41588-026-02544-8", "Figure 1D",
     "Cell-type dot plot using color for mean expression and size for the fraction expressing",
     "Twelve selected genes across 13 SCP3177 cerebellar cell types",
     "SCP3177 summary expression data",
     "Expression is descriptive on the original summary scale; no ET-versus-control test was performed."),
    ("Figure 5", "C", "Gutierrez-Arcelus et al., Nature Communications (2019)", "10.1038/s41467-019-08604-4", "Figure 2A",
     "Sample-level principal-component scores with groups encoded by color and shape",
     "DESeq2 rlog PCA of 40 GSE197345 Purkinje-cell samples using the 500 most variable genes",
     "GSE197345 raw counts; Martuscello et al., Cerebellum (2023), Figure 2 and Methods",
     "The PCA displays sample-level transcriptome structure and is not presented as a significant case-control separation."),
    ("Figure 5", "D", "Gutierrez-Arcelus et al., Nature Communications (2019)", "10.1038/s41467-019-08604-4", "Figure 2C",
     "Row-scaled expression heatmap with a sample-group annotation strip",
     "DESeq2 rlog z scores for the eight nominal Purkinje-cell target transcripts across 16 controls and 24 ET samples",
     "GSE197345 raw counts; Martuscello et al., Cerebellum (2023), Figure 1",
     "Genes were selected by the prespecified target-restricted nominal P < 0.05 criterion; samples are ordered by condition and PC1."),
    ("Figure 5", "E", "Castonguay et al., Nature Genetics (2026)", "10.1038/s41588-026-02544-8", "Figure 6C and 6E",
     "Cell-population differential-expression display with labeled candidate transcripts",
     "Independently expression-filtered GSE197345 Purkinje-cell target estimates",
     "GSE197345",
     "Sixty-two targets were testable and eight met nominal P < 0.05."),
    ("Figure 5", "F", "Castonguay et al., Nature Genetics (2026)", "10.1038/s41588-026-02544-8", "Figure 6D and 6F",
     "Compact aligned transcriptomic summary with a common statistical reference line",
     "mroast, fry and probability-weighted mroast results in whole cerebellum and Purkinje cells",
     "GSE134878 and GSE197345 target-set rotation tests",
     "The probability-weighted result is a sensitivity analysis and no common direction of change is asserted."),
    ("Figure 5", "G", "Ning et al., Nature Communications (2023)", "10.1038/s41467-023-42788-0", "Figure 2E",
     "Aligned gene-by-analysis dot matrix using size for significance and color for signed expression change",
     "Selected genes across MAGMA, whole cerebellum, Purkinje cells, CA3 cis-eQTL and highest-expression cell type",
     "MAGMA, cerebellar expression and SCP3177 results",
     "MAGMA gene P values are unsigned; expression color represents log2 fold change."),
    ("Supplementary Figure S1", "A", "Dekkers et al., Nature Genetics (2026)", "10.1038/s41588-026-02512-2", "Extended Data Figure 8",
     "SNP-specific forest plot", "Wald-ratio estimates for 9 and 4 instruments", "Microbial GWAS instruments", "Variant-specific sensitivity display."),
    ("Supplementary Figure S1", "B", "Wiberg et al., Nature Communications (2019)", "10.1038/s41467-019-08993-6", "Figure 3c-d",
     "Funnel plots of single-variant estimates", "Wald-ratio estimates and precision for two genera", "Microbial MR harmonized instruments", "Visual assessment of estimate symmetry."),
    ("Supplementary Figure S1", "C", "Jurgens et al., Nature Genetics (2024)", "10.1038/s41588-024-01975-5", "Figure 4b",
     "Multiple MR estimators on a common odds-ratio scale", "Four estimators for six additional nonduplicate microbial instrument sets", "Microbial MR estimator results", "Nominal associations only."),
    ("Supplementary Figure S1", "D", "Dekkers et al., Nature Genetics (2026)", "10.1038/s41588-026-02512-2", "Figure 3",
     "Aligned estimates across taxonomic levels", "Class, order and family rows sharing one 10-SNP instrument set", "MiBioGen taxonomic associations", "The three rows represent one nonduplicate instrument set after exact duplicates were collapsed."),
    ("Supplementary Figure S2", "A", "Redan et al., Journal of Agricultural and Food Chemistry (2020)", "10.1021/acs.jafc.0c05890", "Figure 1",
     "Compound identity fields adjacent to chemical nomenclature", "CID, HMDB, ChEBI, formula, molecular mass and SMILES", "PubChem, HMDB and ChEBI", "Identifiers refer to the exact open-chain acid."),
    ("Supplementary Figure S2", "B", "Redan et al., Journal of Agricultural and Food Chemistry (2020)", "10.1021/acs.jafc.0c05890", "Figure 1",
     "Side-by-side chemical structure comparison", "Open-chain acid and phenyl-γ-valerolactone", "PubChem structures", "The compounds are not interchangeable."),
    ("Supplementary Figure S2", "C", "Daina et al., Nucleic Acids Research (2019)", "10.1093/nar/gkz382", "Figure 2",
     "Complete probability-ranked target series", "All 98 unique SwissTargetPrediction genes", "SwissTargetPrediction results", "Native probability values are retained."),
    ("Supplementary Figure S2", "D", "Irwin et al., Journal of Chemical Information and Modeling (2018)", "10.1021/acs.jcim.7b00316", "Figure 3D",
     "Complete target-level SEA display with similarity retained as a second metric", "All seven SEA genes", "Similarity Ensemble Approach results", "Native z-scores and Tanimoto similarities are retained."),
    ("Supplementary Figure S2", "E", "Qi et al., Nucleic Acids Research (2025)", "10.1093/nar/gkae1002", "Figure 1A-B",
     "Gene-by-source membership matrix", "The 105-gene union across SwissTargetPrediction and SEA with autosomal MAGMA coverage", "Target-source membership and MAGMA coverage", "Symbols denote prediction-source membership or analysis coverage, not biological effect."),
    ("Supplementary Figure S3", "A", "Dekkers et al., Nature Genetics (2026)", "10.1038/s41588-026-02512-2", "Extended Data Figure 5",
     "Observed-versus-expected Q-Q plot", "MAGMA P values for 97 autosomal predicted targets", "Target-restricted MAGMA analysis", "The selected target set is not a genome-wide calibration sample."),
    ("Supplementary Figure S3", "B", "Karlsson Linnér et al., Nature Neuroscience (2021)", "10.1038/s41593-021-00908-3", "Extended Data Figure 4",
     "Complete gene-association results arranged in aligned columns", "MAGMA P values for all 97 autosomal predicted targets", "Target-restricted MAGMA analysis", "Nominal and BH-adjusted P-value thresholds are displayed separately."),
    ("Supplementary Figure S3", "C", "Castonguay et al., Nature Genetics (2026)", "10.1038/s41588-026-02544-8", "Figure 6C and 6E",
     "Cerebellar target-expression display with labeled candidate transcripts", "Independently expression-filtered GSE134878 whole-cerebellum target estimates", "GSE134878", "Sixty-four predicted targets were testable and five met nominal P < 0.05."),
    ("Supplementary Figure S3", "D", "Castonguay et al., Nature Genetics (2026)", "10.1038/s41588-026-02544-8", "Figure 6",
     "Cross-analysis comparison of differential-expression estimates", "Published and reanalyzed log2 fold changes for 36 Purkinje-cell genes", "GSE197345 and the published 36-gene table", "Thirty-five of 36 genes had BH-adjusted P < 0.05 in the reanalysis."),
    ("Supplementary Figure S3", "E", "Gabitto et al., Nature Neuroscience (2024)", "10.1038/s41593-024-01774-5", "Figure 7h",
     "Complete cell-type dot matrix using color and size for expression summaries", "All 105 predicted genes across 13 SCP3177 cell types", "SCP3177 summary expression data", "The panel shows descriptive cell-type expression only."),
]


def write_supplementary_tables(d) -> list[str]:
    table_dir = FORMAL / "supplementary_tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    chromosome_x = d["magma"][~as_bool(d["magma"].included_in_autosomal_magma)][
        ["gene_symbol", "magma_chr"]
    ].copy().sort_values("gene_symbol")
    chromosome_x.columns = ["gene_symbol", "chromosome"]
    chromosome_x["chromosome"] = "X"
    chromosome_x["reason_not_analyzed"] = "Chromosome X genes were not included in the autosomal MAGMA analysis."
    chromosome_x_path = table_dir / "Supplementary_Table_chromosome_X_targets_not_tested_by_autosomal_MAGMA.csv"
    chromosome_x.to_csv(chromosome_x_path, index=False, encoding="utf-8-sig")

    annotation_columns = [
        "gene_symbol", "in_swiss", "in_sea", "included_in_autosomal_magma", "magma_pvalue",
        "magma_bh_adjusted_p", "magma_nominal_p_lt_0_05",
        "magma_bh_adjusted_p_lt_0_05", "bulk_log2FoldChange_et_vs_control",
        "bulk_pvalue", "bulk_nominal_p_lt_0_05", "purkinje_log2FoldChange_et_vs_control",
        "purkinje_pvalue", "purkinje_nominal_p_lt_0_05", "in_published_differential_expression_table",
        "published_et_genetic_evidence", "highest_mean_expression_cell_type",
        "highest_mean_expression", "fraction_expressing_in_that_cell_type",
    ]
    annotations = d["host23"][annotation_columns].copy().sort_values("gene_symbol")
    annotations_path = table_dir / "Supplementary_Table_gene_association_and_expression_annotations.csv"
    annotations.to_csv(annotations_path, index=False, encoding="utf-8-sig")

    target_set_columns = [
        "dataset", "gene_set", "method", "set_statistic", "rotations", "NGenes",
        "PropDown", "PropUp", "Direction", "PValue", "FDR", "PValue.Mixed", "FDR.Mixed",
    ]
    target_set = pd.concat([d["mroast"], d["fry"]], ignore_index=True, sort=False)
    target_set = target_set.reindex(columns=target_set_columns)
    target_set_path = table_dir / "Supplementary_Table_target_set_rotation_tests.csv"
    target_set.to_csv(target_set_path, index=False, encoding="utf-8-sig")

    weighted_columns = [
        "dataset", "gene_set", "method", "set_statistic", "rotations", "weight_source",
        "weight_normalization", "NGenes", "PropDown", "PropUp", "Direction", "PValue",
    ]
    weighted_path = table_dir / "Supplementary_Table_probability_weighted_mroast.csv"
    d["weighted_roast"][weighted_columns].to_csv(weighted_path, index=False, encoding="utf-8-sig")

    expression_columns = [
        "dataset", "gene_symbol", "included_in_expression_analysis", "baseMean",
        "effect_et_vs_control", "statistic_et_vs_control", "pvalue", "nominal_p_lt_0_05",
        "bh_adjusted_p", "bh_adjusted_p_lt_0_05",
    ]
    expression_path = table_dir / "Supplementary_Table_predicted_target_expression.csv"
    d["target_expression"][expression_columns].to_csv(expression_path, index=False, encoding="utf-8-sig")

    eqtl_path = table_dir / "Supplementary_Table_CA3_cell_type_cis_eQTL.csv"
    eqtl = d["ca_eqtl"].loc[d["ca_eqtl"].gene.eq("CA3")].copy()
    eqtl["bonferroni_p_across_18888_tests"] = (pd.to_numeric(eqtl["minimum_nominal_p"], errors="coerce") * 18888).clip(upper=1.0)
    eqtl_columns = [
        "cell_type", "gene", "cis_variant_count", "minimum_nominal_p",
        "bonferroni_p_across_18888_tests", "lead_variant", "lead_slope", "lead_slope_se",
    ]
    eqtl[eqtl_columns].to_csv(
        eqtl_path, index=False, encoding="utf-8-sig"
    )

    coloc_path = table_dir / "Supplementary_Table_CA3_colocalization_prior_sensitivity.csv"
    coloc_columns = [
        "cell_type", "gene", "p1", "p2", "p12", "nsnps", "PPH0", "PPH1", "PPH2",
        "PPH3", "PPH4", "top_shared_variant", "top_shared_variant_PPH4", "method_reference",
    ]
    d["ca3_coloc_grid"][coloc_columns].to_csv(coloc_path, index=False, encoding="utf-8-sig")
    return [str(path) for path in [chromosome_x_path, annotations_path, target_set_path, weighted_path,
                                   expression_path, eqtl_path, coloc_path]]


def write_reference_package() -> None:
    REGISTRY.mkdir(parents=True, exist_ok=True)
    columns = ["figure", "panel", "reference_paper", "doi", "reference_figure", "presentation_structure",
               "study_data", "data_source", "interpretation"]
    pd.DataFrame(REFERENCE_ROWS, columns=columns).to_csv(
        REGISTRY / "FIGURE_PANEL_REFERENCE_SOURCES.csv", index=False, encoding="utf-8-sig"
    )
    specifications = """# Figure specifications

- Contents: five main figures and three supplementary figures.
- File formats: paired PDF and 600-dpi RGB TIFF files.
- Typography: black Arial Bold throughout; uppercase panel labels are 18 pt and use a uniform position.
- Figure 1 is a single full-width study overview and therefore has no panel letter.
- The screening threshold is nominal P < 0.05. BH-adjusted values are reported as supplementary multiplicity annotations and are not used for inclusion or exclusion.
- Predicted targets are the union of SwissTargetPrediction and SEA genes: 98 + 7 genes, observed overlap 0, union 105.
- Microbial MR, human microbe-metabolite association, pure-culture metabolism, target prediction, MAGMA, differential expression and cell-type expression are reported as separate analyses.
- The SCP3177 UMAP uses the 100,000-nucleus visualization subsample downloaded from the Single Cell Portal study containing 1,004,112 nuclei; coordinates were not recomputed. SCP3177 expression panels are descriptive and are not ET-versus-control tests.
- GSE197345 principal-component scores were calculated from DESeq2 rlog values for the 500 most variable genes. The selected-gene heatmap shows row z scores across all 40 samples.
- SNP association scatter plots and leave-one-out estimates for the two focal genera are displayed in Figure 2; single-variant forests and funnel plots remain supplementary.
- CA1, CA2 and CA3 are shown at one carbonic anhydrase locus; only CA3 was reported as a putative causal gene in the cited essential-tremor GWAS.
- Mixed-direction self-contained target-set testing allows upregulated and downregulated genes to contribute to the same gene-set statistic; it does not imply a uniform direction of effect.
- CA3 colocalization is displayed across prespecified p12 values because the posterior probability was sensitive to the shared-variant prior.
- Reference figures contribute information organization and layout only; their numerical data and distinctive wording are not reproduced.
"""
    (REGISTRY / "FIGURE_SPECIFICATIONS.md").write_text(specifications, encoding="utf-8")
    legends = """# Figure legends

## Figure 1 | Study overview
The exposure GWAS comprised 211 microbial taxa and the outcome GWAS comprised 16,480 essential-tremor cases and 1,936,173 controls. Ten IVW associations met nominal P < 0.05 and represented eight nonduplicate instrument sets after exact duplicates were collapsed. The human *Faecalibacterium*–5-(3,4-dihydroxyphenyl)pentanoic acid association is correlational. Pure-culture metabolism by *Flavonifractor plautii* concerns related catechin ring-fission products and does not establish production of PubChem CID 49831816. SwissTargetPrediction and SEA contributed 98 and 7 nonoverlapping genes, respectively, to a union of 105 predicted targets. MAGMA tested 97 autosomal targets. Independent expression filtering retained 64 whole-cerebellum and 62 Purkinje-cell targets; the Purkinje-cell union showed mixed-direction differential expression by mroast (P = 0.048).

## Figure 2 | Mendelian randomization and microbial evidence
**A,** IVW odds ratios and 95% confidence intervals for eight nonduplicate microbial instrument sets, with P values shown in a separate column. **B,** IVW, MR-Egger, weighted-median and weighted-mode estimates for *Faecalibacterium* and *Flavonifractor*, together with Cochran Q, MR-Egger intercept, MR-PRESSO global-test and leave-one-out summaries. **C,** SNP-exposure and SNP-outcome associations with IVW, MR-Egger, weighted-median and weighted-mode slopes for the two focal genera. **D,** Leave-one-out IVW estimates after sequential exclusion of each instrument; dashed lines denote the full IVW estimates. **E,** Sensitivity analyses and microbial evidence. Filled circles indicate that the stated criterion or evidence is present; open circles indicate a correlational gutMGene association.

## Figure 3 | Catechin-related metabolite and predicted targets
**A,** The correlational human association between *Faecalibacterium* and 5-(3,4-dihydroxyphenyl)pentanoic acid (PubChem CID 49831816) is shown separately from pure-culture catechin metabolism by *F. plautii*. The latter produces 5-(3,4-dihydroxyphenyl)-γ-valerolactone and 4-hydroxy-5-(3,4-dihydroxyphenyl)valeric acid from a catechin ring-fission intermediate. The open-chain acid CID 49831816 and cyclic lactone CID 45093080 have different formulae and are chemically distinct. **B,** SwissTargetPrediction and SEA contributed nonoverlapping genes to the 105-gene union; observed overlap was zero. **C,** The 15 highest-ranked SwissTargetPrediction genes are ordered by predicted probability; carbonic anhydrases are highlighted. **D,** SEA z-scores are plotted against the reported Tanimoto coefficients for seven genes; each point denotes one predicted target.

## Figure 4 | Gene association and cell-type regulatory evidence at the carbonic anhydrase locus
**A,** Chromosome-ordered MAGMA gene P values for 97 autosomal predicted targets; the dashed line denotes nominal P = 0.05. **B,** Prediction, published genetics, differential-expression and cell-type annotations for 11 genes with nominal MAGMA P < 0.05. **C,** CA1, CA2 and CA3 gene bodies, MAGMA windows and leading variants at the chromosome 8 locus. CA3 and CA2 share rs955007 as the leading variant in overlapping MAGMA windows; conditional independence was not tested. **D,** CA3 cis-eQTL P values across eight cerebellar cell types. Filled points identify Granule cells and Oligodendrocytes, which met Bonferroni P < 0.05 across 18,888 variant-cell-type tests. **E,** Posterior probability of a shared causal variant (PPH4) across four p12 values. The dashed line marks the default p12 of 10−5; default PPH4 was 0.419 in Granule cells and 0.251 in Oligodendrocytes, and higher posterior probabilities at p12 = 5×10−5 were prior-sensitive.

## Figure 5 | Predicted-target transcription in essential-tremor cerebellum
**A,** Cell-type UMAP reconstructed from SCP3177 coordinates and ontology labels downloaded from the Single Cell Portal for the visualization subsample of 100,000 nuclei from the 1,004,112-nucleus atlas. The embedding is descriptive and was not recomputed. **B,** SCP3177 cell-type expression of 12 selected genes. Point size represents the fraction expressing, color represents mean expression and black outlines mark the cell type with the highest mean expression. **C,** Principal-component analysis of DESeq2 rlog-transformed GSE197345 Purkinje-cell counts using the 500 most variable genes; points represent 16 control and 24 essential-tremor samples. **D,** Row-scaled rlog expression of the eight target transcripts meeting nominal P < 0.05 in the target-restricted Purkinje-cell analysis. Samples are grouped by condition and ordered by PC1 within condition. **E,** Purkinje-cell target expression among 62 independently expression-filtered genes; eight met nominal P < 0.05. **F,** Rotation gene-set tests after independent expression filtering. Points show mroast, fry and mroast weighted by SwissTargetPrediction probability for whole cerebellum and Purkinje cells; the dashed line denotes nominal P = 0.05. The Purkinje-cell union test used the mean-square statistic and 19,999 rotations (P = 0.048); the probability-weighted sensitivity analysis yielded P = 0.059. **G,** Gene-based association, whole-cerebellum and Purkinje-cell expression, CA3 cis-eQTL evidence and highest-expression cell type for selected targets. Point size represents −log10(P), expression color represents log2 fold change and MAGMA P values are unsigned. AS, astrocyte; BG, Bergmann glia; EC, endothelial cell; GC, granule cell; GOL, Golgi cell; IN, interneuron; MG, microglia; MLI, molecular layer interneuron; OL, oligodendrocyte; OPC, oligodendrocyte precursor cell; PC, Purkinje cell; PER, pericyte; UBC, unipolar brush cell. SCP3177 expression is shown on the original summary scale and is not an essential-tremor-versus-control test.

## Supplementary Figure S1 | Mendelian randomization diagnostics
**A,** SNP-specific Wald-ratio estimates for 9 *Faecalibacterium* and 4 *Flavonifractor* instruments. **B,** Funnel plots of Wald-ratio estimates against precision. **C,** Four MR estimators for six additional nonduplicate microbial instrument sets. **D,** Class, order and family Methanobacteria-lineage rows sharing the same 10-SNP instrument set and identical IVW estimates.

## Supplementary Figure S2 | Compound identity and target prediction
**A,** Names, identifiers, molecular formula, molecular mass and canonical SMILES for 5-(3,4-dihydroxyphenyl)pentanoic acid. **B,** Structural comparison of the open-chain acid and phenyl-γ-valerolactone. **C,** Complete SwissTargetPrediction ranking. **D,** Complete SEA ranking. **E,** Prediction-source membership and autosomal MAGMA coverage for all 105 predicted targets. Filled circles indicate membership or coverage; X indicates a gene not tested in MAGMA.

## Supplementary Figure S3 | Gene association and cerebellar expression
**A,** Q-Q plot of MAGMA gene P values for 97 autosomal predicted targets. **B,** MAGMA P values for all 97 autosomal predicted targets. **C,** Whole-cerebellum target expression among 64 independently expression-filtered genes. **D,** Published and reanalyzed Purkinje-cell log2 fold changes for 36 genes. **E,** Cell-type expression of all 105 predicted genes across 13 SCP3177 cerebellar cell types. Chromosome-X targets not tested by the autosomal MAGMA analysis and gene-level association and expression annotations are provided in supplementary tables.

"""
    (REGISTRY / "FIGURE_LEGENDS.md").write_text(legends, encoding="utf-8")
    palette = pd.DataFrame([{"visual_category": key, "hex": value} for key, value in figure_palette().items()])
    palette.to_csv(REGISTRY / "FIGURE_COLOR_PALETTE.csv", index=False, encoding="utf-8-sig")


FIGURES = {
        "Figure1": (figure1, FORMAL / "main_figures" / "Figure1_study_overview"),
        "Figure2": (figure2, FORMAL / "main_figures" / "Figure2_mendelian_randomization"),
        "Figure3": (figure3, FORMAL / "main_figures" / "Figure3_metabolite_target_prediction"),
        "Figure4": (figure4, FORMAL / "main_figures" / "Figure4_gene_association_regulatory_evidence"),
        "Figure5": (figure5, FORMAL / "main_figures" / "Figure5_cerebellar_target_transcription"),
        "Supplementary_Figure_S1": (supplementary_figure_s1, FORMAL / "supplementary_figures" / "Supplementary_Figure_S1_MR_diagnostics"),
        "Supplementary_Figure_S2": (supplementary_figure_s2, FORMAL / "supplementary_figures" / "Supplementary_Figure_S2_compound_identity_target_prediction"),
        "Supplementary_Figure_S3": (supplementary_figure_s3, FORMAL / "supplementary_figures" / "Supplementary_Figure_S3_gene_association_cerebellar_expression"),
}


def run() -> None:
    setup_style()
    data = load_data()
    validate_data_dimensions(data)
    write_reference_package()
    write_supplementary_tables(data)
    colors = figure_palette()
    for name, (figure_function, stem) in FIGURES.items():
        print(f"Building {name}...", flush=True)
        fig = figure_function(data, colors)
        save_formal(fig, stem, dpi=600)
        plt.close(fig)


if __name__ == "__main__":
    run()
