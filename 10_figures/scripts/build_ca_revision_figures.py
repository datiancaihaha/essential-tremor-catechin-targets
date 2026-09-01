from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDepictor


ROOT = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1\outputs\Movement_Disorders_CA_Revision_20260828_v1")
CODE = ROOT / "07_reproducible_code"
SCP_PRIMARY = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1\outputs\SCP3177_strengthening_analysis_20260827")
SCP_SENSITIVITY = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1\outputs\SCP3177_completion_analysis_20260828")
MM = 1 / 25.4
FIG_WIDTH_MM = 183


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


sys.path.insert(0, str(CODE))
pub = _load_module("publication_figures", CODE / "generate_publication_figures.py")
base = pub.base
style = _load_module("revision_figure_style", CODE / "figure_style.py")


PALETTES = {
    "palette_09": {
        "microbe": "#008D94",
        "metabolite": "#D60066",
        "genetics": "#7A2D90",
        "bulk": "#4CAF50",
        "purkinje": "#1A237E",
        "cell": "#E53935",
        "positive": "#E53935",
        "negative": "#1A237E",
        "support": "#FFC107",
        "light": "#D8EEEE",
        "neutral": "#BDBDBD",
        "pale": "#F3F3F3",
        "dark": "#2E2E2E",
    },
    "palette_10": {
        "microbe": "#5BBFC3",
        "metabolite": "#E58B89",
        "genetics": "#8478AA",
        "bulk": "#82B99B",
        "purkinje": "#648DB8",
        "cell": "#DFA259",
        "positive": "#E58B89",
        "negative": "#648DB8",
        "support": "#E6C96E",
        "light": "#DDEDEF",
        "neutral": "#B8B8B8",
        "pale": "#F4F4F2",
        "dark": "#343434",
    },
}


FORM_ORDER = pub.FORM_ORDER
FORM_LABELS = pub.FORM_LABELS
FORM_SHORT = ["Acid", "γVL", "(5R)-γVL", "4-OH acid", "3′-S", "4′-S", "3′-G", "4′-G"]


def setup_style() -> None:
    style.setup_style()
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial"],
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel(ax, label: str, x: float = -0.08, y: float = 1.02) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=18, fontweight="bold",
            va="bottom", ha="left", clip_on=False)


def clean(ax, grid: bool = False) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis="x", color="#D8D8D8", lw=0.55, zorder=0)


def muted(color: str, amount: float = 0.38) -> tuple[float, float, float]:
    rgb = np.array(to_rgb(color))
    neutral = np.array(to_rgb("#AFAFAF"))
    return tuple((1 - amount) * rgb + amount * neutral)


def load_data() -> dict:
    data = pub.load_data()
    data["raw_voom"] = pd.read_csv(
        SCP_PRIMARY / "09_raw_count_voom" / "all_target_gene_cell_type_results.csv"
    )
    data["raw_ca"] = pd.read_csv(
        SCP_PRIMARY / "10_raw_count_integrated_results" / "carbonic_anhydrase_family_results.csv"
    )
    data["model_summary"] = pd.read_csv(
        SCP_PRIMARY / "09_raw_count_voom" / "cell_type_model_summary.csv"
    )
    data["composition"] = pd.read_csv(
        SCP_SENSITIVITY / "01_composition_propeller_14_cell_types" / "propeller_blocked_all_14_cell_types.csv"
    )
    data["composition_evidence"] = pd.read_csv(
        SCP_SENSITIVITY / "06_integrated_results" / "composition_evidence_matrix.csv"
    )
    data["ca_thresholds"] = pd.read_csv(
        SCP_SENSITIVITY / "06_integrated_results" / "carbonic_anhydrase_nucleus_threshold_stability.csv"
    )
    data["target_thresholds"] = pd.read_csv(
        SCP_SENSITIVITY / "06_integrated_results" / "target_gene_nucleus_threshold_stability.csv"
    )
    data["camera_thresholds"] = pd.read_csv(
        SCP_SENSITIVITY / "06_integrated_results" / "camera_nucleus_threshold_stability.csv"
    )
    data["loo_composition"] = pd.read_csv(
        SCP_SENSITIVITY / "01_composition_propeller_14_cell_types" / "leave_one_et_donor_out_summary.csv"
    )
    data["chemical_identity"] = pd.read_csv(ROOT / "09_source_data" / "chemical_exposure_classification.csv")
    return data


def _card(ax, xy, wh, color, heading, lines, heading_size=6.4, body_size=5.8):
    x, y = xy
    w, h = wh
    ax.add_patch(FancyBboxPatch((x + 0.003, y - 0.004), w, h,
                                boxstyle="round,pad=0.004,rounding_size=0.012",
                                fc="#000000", ec="none", alpha=0.07, zorder=0))
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.004,rounding_size=0.012",
                                fc="white", ec=color, lw=1.0, zorder=1))
    ax.add_patch(Rectangle((x, y + h - 0.035), w, 0.035, fc=color, ec="none", alpha=0.18, zorder=1.2))
    ax.text(x + 0.010, y + h - 0.0175, heading, fontsize=heading_size, va="center", ha="left")
    ax.text(x + 0.010, y + h - 0.052, lines, fontsize=body_size, va="top", ha="left", linespacing=1.15)


def figure1(data, c):
    """Single integrated study-design canvas; no panel letters."""
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 125 * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    modules = [
        (0.018, 0.145, c["microbe"]),
        (0.185, 0.170, c["genetics"]),
        (0.378, 0.150, c["metabolite"]),
        (0.550, 0.170, c["support"]),
        (0.742, 0.240, c["purkinje"]),
    ]
    headings = [
        "Genome-wide association\nresources",
        "Mendelian\nrandomization",
        "Chemically resolved\nmetabolites",
        "Predicted human\ntargets",
        "Essential-tremor and\ncerebellar evidence",
    ]
    for (x, w, color), title in zip(modules, headings):
        ax.add_patch(FancyBboxPatch((x, 0.075), w, 0.85,
                                    boxstyle="round,pad=0.004,rounding_size=0.018",
                                    fc="white", ec="#B8B8B8", lw=0.8))
        ax.add_patch(FancyBboxPatch((x, 0.845), w, 0.080,
                                    boxstyle="round,pad=0.004,rounding_size=0.018",
                                    fc=color, ec="none", alpha=0.17))
        ax.text(x + w / 2, 0.884, title, ha="center", va="center", fontsize=5.6, linespacing=1.0)
    for left, right in zip(modules[:-1], modules[1:]):
        ax.add_patch(FancyArrowPatch((left[0] + left[1] + 0.004, 0.50), (right[0] - 0.004, 0.50),
                                    arrowstyle="-|>", mutation_scale=9, lw=1.1,
                                    color="#737373", shrinkA=0, shrinkB=0))

    _card(ax, (0.030, 0.650), (0.121, 0.152), c["microbe"], "Microbiome GWAS", "18,340 participants\n211 taxonomic traits")
    _card(ax, (0.030, 0.430), (0.121, 0.152), c["metabolite"], "Species GWAS", "16,017 participants\n17 metagenomic traits")
    _card(ax, (0.030, 0.210), (0.121, 0.152), c["genetics"], "Essential-tremor GWAS", "16,480 cases\n1,936,173 controls")

    mr = data["mr"].copy().sort_values("odds_ratio").reset_index(drop=True)
    mrax = fig.add_axes([0.235, 0.205, 0.095, 0.565])
    y = np.arange(len(mr))
    mrax.errorbar(mr.odds_ratio, y,
                  xerr=[mr.odds_ratio - mr.or_ci_lower, mr.or_ci_upper - mr.odds_ratio],
                  fmt="o", ms=3.6, color=c["microbe"], ecolor="#555555", lw=0.7, capsize=1.5)
    mrax.axvline(1, color="#777777", lw=0.7, ls="--")
    mrax.set_xscale("log")
    mrax.set_xlim(0.60, 1.58)
    taxon_labels = [base.short_taxon(x).replace("Lachnospiraceae", "Lachnosp.").replace("Ruminococcaceae", "Ruminoc.") for x in mr.taxon_name]
    mrax.set_yticks(y, taxon_labels, fontsize=3.9)
    mrax.set_xticks([0.7, 1.0, 1.4], ["0.7", "1.0", "1.4"])
    mrax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    mrax.xaxis.offsetText.set_visible(False)
    mrax.tick_params(length=0, pad=1)
    mrax.set_xlabel("Odds ratio (95% CI)", fontsize=5.0)
    clean(mrax)
    ax.text(0.270, 0.795, "Eight independent microbial exposures", ha="center", fontsize=5.0)

    forms = data["form_summary"].set_index("record_id").loc[FORM_ORDER]
    for i, (label, (_, row)) in enumerate(zip(FORM_SHORT, forms.iterrows())):
        col, rr = i % 2, i // 2
        x = 0.392 + col * 0.064
        y0 = 0.675 - rr * 0.145
        ax.add_patch(FancyBboxPatch((x, y0), 0.052, 0.106,
                                    boxstyle="round,pad=0.003,rounding_size=0.008",
                                    fc=c["metabolite"], ec="none", alpha=0.10))
        ax.text(x + 0.026, y0 + 0.070, label, ha="center", va="center", fontsize=4.7)
        ax.text(x + 0.026, y0 + 0.028, str(int(row.valid_unique_gene_count)),
                ha="center", va="center", fontsize=6.2, color=c["metabolite"])
    ax.text(0.453, 0.135, "targets per exact form", ha="center", fontsize=5.0)

    target_sets = [(198, "Union"), (150, "Conjugate union"), (45, "Conjugate core"), (25, "All-form core")]
    maxw = 0.138
    for i, (n, label) in enumerate(target_sets):
        y0 = 0.735 - i * 0.145
        width = maxw * math.sqrt(n / 198)
        x = 0.635 - width / 2
        ax.add_patch(FancyBboxPatch((x, y0), width, 0.075,
                                    boxstyle="round,pad=0.003,rounding_size=0.012",
                                    fc=c["support"], ec="none", alpha=0.18 + 0.12 * i))
        ax.text(0.635, y0 + 0.050, str(n), ha="center", va="center", fontsize=7.0)
        ax.text(0.635, y0 + 0.020, label, ha="center", va="center", fontsize=4.7)
    ax.text(0.635, 0.145, "CA1  CA2  CA9  CA12: 8/8 forms\nCA3: 5/8 forms",
            ha="center", va="center", fontsize=5.0, linespacing=1.25)

    evidence = [
        ("CA3", "Gene-based association\nGranule-cell cis-eQTL\nOligodendrocyte cis-eQTL", c["genetics"]),
        ("CA12", "Higher expression in\ngranule cells, Bergmann glia\nand astrocytes", c["bulk"]),
        ("CA7", "Purkinje-cell\ndifferential expression", c["purkinje"]),
    ]
    for i, (gene, desc, color) in enumerate(evidence):
        y0 = 0.665 - i * 0.205
        ax.add_patch(FancyBboxPatch((0.760, y0), 0.202, 0.142,
                                    boxstyle="round,pad=0.005,rounding_size=0.012",
                                    fc=color, ec="none", alpha=0.10))
        ax.text(0.777, y0 + 0.095, gene, fontsize=8.2, color=color, va="center")
        ax.text(0.818, y0 + 0.070, desc, fontsize=4.55, va="center", linespacing=1.12)
    ax.text(0.861, 0.125, "109 donors | 1,004,112 nuclei", ha="center", fontsize=5.4)
    return fig


def _draw_structure(ax, smiles: str, label: str, color: str) -> None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {label}")
    rdDepictor.Compute2DCoords(mol)
    conformer = mol.GetConformer()
    coords = np.array([[conformer.GetAtomPosition(i).x, conformer.GetAtomPosition(i).y]
                       for i in range(mol.GetNumAtoms())], dtype=float)
    span = np.ptp(coords, axis=0)
    span[span == 0] = 1.0
    coords = (coords - coords.min(axis=0)) / span
    coords[:, 0] = 0.08 + 0.84 * coords[:, 0]
    coords[:, 1] = 0.12 + 0.78 * coords[:, 1]

    def parallel_line(p1, p2, offset):
        direction = p2 - p1
        length = np.hypot(direction[0], direction[1])
        normal = np.array([-direction[1], direction[0]]) / max(length, 1e-9)
        return p1 + normal * offset, p2 + normal * offset

    for bond in mol.GetBonds():
        p1 = coords[bond.GetBeginAtomIdx()]
        p2 = coords[bond.GetEndAtomIdx()]
        bond_type = bond.GetBondType()
        direction = bond.GetBondDir()
        if direction in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH):
            vec = p2 - p1
            length = np.hypot(vec[0], vec[1])
            normal = np.array([-vec[1], vec[0]]) / max(length, 1e-9)
            if direction == Chem.BondDir.BEGINWEDGE:
                triangle = np.vstack([p1, p2 + normal * 0.018, p2 - normal * 0.018])
                ax.add_patch(Polygon(triangle, closed=True, fc="black", ec="black", lw=0.35))
            else:
                for fraction in np.linspace(0.18, 0.92, 5):
                    centre = p1 + fraction * vec
                    half = 0.002 + 0.012 * fraction
                    ax.plot([centre[0] - normal[0] * half, centre[0] + normal[0] * half],
                            [centre[1] - normal[1] * half, centre[1] + normal[1] * half],
                            color="black", lw=0.55)
            continue
        if bond_type == Chem.BondType.DOUBLE:
            offsets = (-0.009, 0.009)
        elif bond_type == Chem.BondType.TRIPLE:
            offsets = (-0.012, 0.0, 0.012)
        elif bond_type == Chem.BondType.AROMATIC:
            offsets = (-0.007, 0.007)
        else:
            offsets = (0.0,)
        for line_index, offset in enumerate(offsets):
            q1, q2 = parallel_line(p1, p2, offset)
            linestyle = "--" if bond_type == Chem.BondType.AROMATIC and line_index == 1 else "-"
            ax.plot([q1[0], q2[0]], [q1[1], q2[1]], color="black", lw=0.85,
                    linestyle=linestyle, solid_capstyle="round")

    for atom_index, atom in enumerate(mol.GetAtoms()):
        symbol = atom.GetSymbol()
        charge = atom.GetFormalCharge()
        isotope = atom.GetIsotope()
        explicit_h = atom.GetTotalNumHs() if symbol != "C" else 0
        show = symbol != "C" or isotope or charge
        if not show:
            continue
        text = (str(isotope) if isotope else "") + symbol
        if explicit_h:
            text += "H" if explicit_h == 1 else f"H{explicit_h}"
        if charge:
            text += "+" if charge == 1 else "-" if charge == -1 else f"{abs(charge)}{'+' if charge > 0 else '-'}"
        x, y = coords[atom_index]
        ax.text(x, y, text, ha="center", va="center", fontsize=5.3, color="black",
                bbox=dict(boxstyle="round,pad=0.05", fc="white", ec="none"), zorder=5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    ax.text(0.5, -0.03, label, transform=ax.transAxes, ha="center", va="top", fontsize=5.0, color=color)


def _upset(ax, matrix: pd.DataFrame, c, top_n: int = 10) -> None:
    patterns = []
    for _, row in matrix.iterrows():
        patterns.append(tuple(str(row[f"{record}_predicted"]).lower() == "true" for record in FORM_ORDER))
    counts = pd.Series(patterns).value_counts().head(top_n)
    ax.set_axis_off()
    top = ax.inset_axes([0.18, 0.52, 0.80, 0.43])
    bot = ax.inset_axes([0.18, 0.03, 0.80, 0.43], sharex=top)
    x = np.arange(len(counts))
    bars = top.bar(x, counts.values, color=c["support"], ec="black", lw=0.3, width=0.68)
    for bar, value in zip(bars, counts.values):
        top.text(bar.get_x() + bar.get_width() / 2, value + 0.5, str(int(value)), ha="center", fontsize=4.7)
    top.set_ylabel("Targets", fontsize=5.4)
    top.tick_params(axis="x", bottom=False, labelbottom=False)
    clean(top)
    for xi, pattern in enumerate(counts.index):
        active = [i for i, flag in enumerate(pattern) if flag]
        if active:
            bot.plot([xi, xi], [min(active), max(active)], color="#555555", lw=0.75)
        for yi, flag in enumerate(pattern):
            bot.scatter(xi, yi, s=20 if flag else 9, fc=c["metabolite"] if flag else "white",
                        ec="black" if flag else "#BBBBBB", lw=0.3)
    bot.set_yticks(range(8), FORM_SHORT, fontsize=4.4)
    bot.set_ylim(7.6, -0.6)
    bot.set_xticks(x, [str(i + 1) for i in x], fontsize=4.4)
    bot.tick_params(length=0)
    bot.spines[:].set_visible(False)


def figure3(data, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 238 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.02, 0.92, 0.92], width_ratios=[1.02, 0.98],
                          left=0.085, right=0.975, top=0.975, bottom=0.060, hspace=0.30, wspace=0.30)
    identity = data["chemical_identity"].set_index("record_id").loc[FORM_ORDER]
    summary = data["form_summary"].set_index("record_id").loc[FORM_ORDER]
    matrix = data["form_matrix"].copy()
    form_color = muted(c["metabolite"], 0.46)
    support_color = muted(c["support"], 0.48)
    network_color = muted(c["microbe"], 0.44)

    ax = fig.add_subplot(gs[0, 0]); panel(ax, "A", -0.18, 1.04)
    ax.set_axis_off()
    positions = [(0.03, 0.54), (0.27, 0.54), (0.51, 0.54), (0.75, 0.54),
                 (0.03, 0.05), (0.27, 0.05), (0.51, 0.05), (0.75, 0.05)]
    for (record, row), (x, y) in zip(identity.iterrows(), positions):
        sub = ax.inset_axes([x, y, 0.21, 0.40])
        _draw_structure(sub, row.isomeric_smiles, FORM_LABELS[record], form_color)

    ax = fig.add_subplot(gs[0, 1]); panel(ax, "B", -0.18, 1.04)
    y = np.arange(8)
    counts = summary.valid_unique_gene_count.astype(int).to_numpy()
    ax.hlines(y, 88, counts, color=form_color, lw=2.4, alpha=0.60)
    ax.scatter(counts, y, s=55, fc=form_color, ec="black", lw=0.4)
    for yi, value in enumerate(counts):
        ax.text(value + 0.5, yi, str(value), va="center", fontsize=5.1)
    ax.set_yticks(y, FORM_SHORT, fontsize=5.5)
    ax.set_xlim(88, 104)
    ax.set_xlabel("Predicted targets")
    ax.invert_yaxis(); clean(ax, grid=True)

    ax = fig.add_subplot(gs[1, 0]); panel(ax, "C", -0.18, 1.04)
    local = dict(c); local["support"] = support_color; local["metabolite"] = form_color
    _upset(ax, matrix, local, 10)

    ax = fig.add_subplot(gs[1, 1]); panel(ax, "D", -0.18, 1.04)
    ca = ["CA1", "CA2", "CA3", "CA7", "CA9", "CA12"]
    probs = np.full((len(ca), 8), np.nan)
    for i, gene in enumerate(ca):
        row = matrix[matrix.gene_symbol.eq(gene)].iloc[0]
        for j, record in enumerate(FORM_ORDER):
            probs[i, j] = pd.to_numeric(row[f"{record}_probability"], errors="coerce")
    gene_colors = [muted(x, 0.38) for x in
                   [c["microbe"], c["metabolite"], c["genetics"], c["bulk"], c["purkinje"], c["cell"]]]
    x = np.arange(8)
    for i, gene in enumerate(ca):
        valid = np.isfinite(probs[i])
        ax.plot(x[valid], probs[i, valid], marker=["o", "s", "D", "^", "v", "P"][i],
                ms=4.1, lw=1.25, color=gene_colors[i], label=gene)
    ax.set_xticks(x, FORM_SHORT, rotation=38, ha="right", fontsize=5.1)
    ax.set_ylim(0, 0.86)
    ax.set_ylabel("Predicted probability")
    ax.legend(frameon=False, ncol=3, fontsize=5.0, loc="upper right", handlelength=1.6)
    clean(ax)

    ax = fig.add_subplot(gs[2, 0]); panel(ax, "E", -0.18, 1.04)
    sets = pub._form_sets(matrix)
    jac = np.zeros((8, 8))
    for i, a in enumerate(FORM_ORDER):
        for j, b in enumerate(FORM_ORDER):
            jac[i, j] = len(sets[a] & sets[b]) / len(sets[a] | sets[b])
    theta = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    xy = np.c_[np.cos(theta), np.sin(theta)]
    for i in range(8):
        for j in range(i + 1, 8):
            if jac[i, j] < 0.30:
                continue
            ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]],
                    color=network_color, alpha=0.12 + 0.60 * (jac[i, j] - 0.30) / 0.70,
                    lw=0.35 + 2.4 * jac[i, j], zorder=1)
    node_colors = [form_color, support_color, support_color, form_color,
                   muted(c["genetics"], 0.42), muted(c["genetics"], 0.42),
                   muted(c["purkinje"], 0.42), muted(c["purkinje"], 0.42)]
    ax.scatter(xy[:, 0], xy[:, 1], s=135, fc=node_colors, ec="white", lw=1.0, zorder=3)
    for i, label in enumerate(FORM_SHORT):
        ax.text(1.20 * xy[i, 0], 1.20 * xy[i, 1], label, ha="center", va="center", fontsize=5.2)
    ax.set_xlim(-1.38, 1.38); ax.set_ylim(-1.34, 1.34); ax.set_aspect("equal"); ax.set_axis_off()
    ax.text(0.02, 0.02, "Edge width: pairwise target-set Jaccard index", transform=ax.transAxes, fontsize=4.8)

    ax = fig.add_subplot(gs[2, 1]); panel(ax, "F", -0.18, 1.04)
    hierarchy = [(198, "All-form union"), (150, "Circulating-conjugate union"),
                 (45, "Circulating-conjugate core"), (25, "All-form core")]
    max_value = hierarchy[0][0]
    funnel_colors = [muted(c["microbe"], 0.52), muted(c["metabolite"], 0.52),
                     muted(c["genetics"], 0.52), muted(c["support"], 0.52)]
    ax.set_xlim(-1.0, 1.0); ax.set_ylim(-0.15, 4.1); ax.set_axis_off()
    for i, ((value, label), color) in enumerate(zip(hierarchy, funnel_colors)):
        upper = 0.92 * math.sqrt(value / max_value)
        next_value = hierarchy[i + 1][0] if i + 1 < len(hierarchy) else value * 0.72
        lower = 0.92 * math.sqrt(next_value / max_value)
        top = 3.85 - i * 0.95; bottom = top - 0.70
        ax.add_patch(Polygon([(-upper, top), (upper, top), (lower, bottom), (-lower, bottom)],
                             closed=True, fc=color, ec="white", lw=0.8))
        ax.text(0, (top + bottom) / 2 + 0.08, str(value), ha="center", va="center", fontsize=8.0)
        ax.text(0, (top + bottom) / 2 - 0.16, label, ha="center", va="center", fontsize=5.1)
    return fig


def _donor_logcpm(cell_type: str, gene: str) -> pd.DataFrame:
    root = SCP_SENSITIVITY / "05_donor_aggregated_count_sensitivity" / "donor_celltype_raw_counts"
    slug = cell_type.lower().replace(" ", "_")
    metadata = pd.read_csv(root / "sample_metadata" / f"{slug}.csv")
    genes = pd.read_csv(root / "genes.csv")
    gene_row = genes.index[genes.gene_symbol.eq(gene)]
    if len(gene_row) != 1:
        raise ValueError(f"Expected one {gene} row, found {len(gene_row)}")
    n_genes = len(genes)
    n_donors = len(metadata)
    arr = np.memmap(root / "count_matrices" / f"{slug}.gene_by_donor.int32.bin",
                    dtype="<i4", mode="r", shape=(n_genes, n_donors), order="F")
    counts = np.asarray(arr[int(gene_row[0]), :], dtype=float)
    metadata = metadata.copy()
    metadata["logCPM"] = np.log2((counts + 0.5) / (metadata.library_size.astype(float) + 1.0) * 1e6)
    return metadata


def figure5(data, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 244 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.02, 1.00, 0.88], width_ratios=[1.02, 0.98],
                          left=0.095, right=0.975, top=0.975, bottom=0.060, hspace=0.36, wspace=0.32)

    ax = fig.add_subplot(gs[0, 0]); panel(ax, "A", -0.18, 1.03)
    um = data["umap"].copy()
    categories = list(pd.Categorical(um.cell_type).categories)
    base_colors = [c["microbe"], c["metabolite"], c["genetics"], c["bulk"], c["purkinje"], c["cell"], c["support"]]
    cell_colors = [muted(base_colors[i % len(base_colors)], 0.40 + 0.08 * (i // len(base_colors))) for i in range(len(categories))]
    for i, label in enumerate(categories):
        q = um[um.cell_type.eq(label)]
        ax.scatter(q.UMAP1, q.UMAP2, s=0.42, fc=cell_colors[i], ec="none", alpha=0.78, rasterized=True)
        ax.text(q.UMAP1.median(), q.UMAP2.median(), base._cell_abbreviation(label), fontsize=4.6,
                ha="center", va="center", color="#222222")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); clean(ax)

    ax = fig.add_subplot(gs[0, 1]); panel(ax, "B", -0.18, 1.03)
    raw = data["raw_voom"].copy()
    ca_genes = ["CA1", "CA2", "CA3", "CA4", "CA7", "CA9", "CA12", "CA13"]
    cell_order = ["Granule", "Bergmann", "Astrocytes", "Purkinje", "Oligodendrocytes", "OPC", "MLI_1", "MLI_2", "Golgi", "UBC", "Microglia", "Endocytes", "PLI", "Pericytes"]
    heat = raw[raw.gene_symbol.isin(ca_genes)].pivot_table(index="gene_symbol", columns="cell_type", values="logFC", aggfunc="first").reindex(index=ca_genes, columns=cell_order)
    cmap = LinearSegmentedColormap.from_list("signed", [c["negative"], "#F7F7F7", c["positive"]])
    vmax = np.nanmax(np.abs(heat.to_numpy()))
    im = ax.imshow(np.ma.masked_invalid(heat.to_numpy()), aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(len(ca_genes)), ca_genes)
    ax.set_xticks(range(len(cell_order)), [x.replace("Oligodendrocytes", "OL").replace("Astrocytes", "AS").replace("Bergmann", "BG").replace("Purkinje", "PC") for x in cell_order], rotation=48, ha="right", fontsize=4.8)
    sig = raw.set_index(["gene_symbol", "cell_type"])["P.Value"]
    for i, gene in enumerate(ca_genes):
        for j, cell in enumerate(cell_order):
            if (gene, cell) in sig.index and float(sig.loc[(gene, cell)]) < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=7.0)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025); cb.set_label("log2 fold change", fontsize=5.7)

    ax = fig.add_subplot(gs[1, 0]); panel(ax, "C", -0.18, 1.03)
    ca12 = data["raw_ca"].query("gene_symbol == 'CA12'").sort_values("logFC")
    y = np.arange(len(ca12))
    ax.errorbar(ca12.logFC, y,
                xerr=[ca12.logFC - ca12["CI.L"], ca12["CI.R"] - ca12.logFC],
                fmt="o", ms=6, color=c["bulk"], ecolor="#4B4B4B", capsize=2.5, lw=0.9)
    ax.axvline(0, color="#777777", lw=0.8, ls="--")
    ax.set_yticks(y, ca12.cell_type)
    for yi, (_, row) in enumerate(ca12.iterrows()):
        ax.text(1.42, yi, f"P={row['P.Value']:.3g}", ha="right", va="center", fontsize=5.2)
    ax.set_xlim(0.0, 1.45)
    ax.set_xlabel("log2 fold change (95% CI)")
    clean(ax, grid=True)

    ax = fig.add_subplot(gs[1, 1]); panel(ax, "D", -0.18, 1.03)
    rng = np.random.default_rng(20260828)
    positions = [0, 1, 3, 4]
    violin_values = []
    violin_colors = []
    for idx, (cell, label) in enumerate([("granule", "Granule"), ("bergmann", "Bergmann glia")]):
        q = _donor_logcpm(cell, "CA12")
        for group_i, disease in enumerate(["normal", "essential tremor"]):
            sub = q[q.disease.eq(disease)]
            xpos = idx * 3 + group_i
            violin_values.append(sub.logCPM.to_numpy())
            jitter = rng.uniform(-0.16, 0.16, len(sub))
            color = c["neutral"] if disease == "normal" else c["bulk"]
            violin_colors.append(color)
            ax.scatter(np.full(len(sub), xpos) + jitter, sub.logCPM, s=12, fc=color, ec="white", lw=0.25, alpha=0.72)
            mean = sub.logCPM.mean(); se = sub.logCPM.std(ddof=1) / math.sqrt(len(sub))
            ax.errorbar(xpos, mean, yerr=se, fmt="D", ms=5, color="black", mfc=color, capsize=3, lw=1.0, zorder=4)
    violins = ax.violinplot(violin_values, positions=positions, widths=0.72, showextrema=False)
    for body, color in zip(violins["bodies"], violin_colors):
        body.set_facecolor(color); body.set_edgecolor(color); body.set_alpha(0.18)
    ax.set_xticks([0, 1, 3, 4], ["Control", "ET", "Control", "ET"], rotation=25, ha="right", fontsize=5.5)
    ax.text(0.5, 1.02, "Granule cells", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.8)
    ax.text(3.5, 1.02, "Bergmann glia", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.8)
    ax.set_ylabel("CA12 expression (log2 counts per million)")
    clean(ax)

    lower = gs[2, :].subgridspec(1, 3, wspace=0.50)
    ax = fig.add_subplot(lower[0, 0]); panel(ax, "E", -0.26, 1.03)
    purk = data["host"].copy()
    purk = purk[purk.gene_symbol.isin(["CA7", "ADCY10", "TET3"])].copy()
    purk = purk.sort_values("purkinje_log2FoldChange_et_vs_control")
    y = np.arange(len(purk))
    ax.scatter(purk.purkinje_log2FoldChange_et_vs_control, y, s=55,
               fc=[c["purkinje"] if float(p) < 0.05 else c["neutral"] for p in purk.purkinje_pvalue], ec="black", lw=0.4)
    ax.axvline(0, color="#777777", lw=0.8, ls="--")
    ax.set_yticks(y, purk.gene_symbol)
    xvals = pd.to_numeric(purk["purkinje_log2FoldChange_et_vs_control"], errors="coerce")
    x_limit = max(1.5, float(np.ceil(np.nanmax(np.abs(xvals)) * 2.0) / 2)) if np.isfinite(np.nanmax(np.abs(xvals))) else 1.5
    for yi, row in enumerate(purk.itertuples()):
        ax.text(x_limit * 0.96, yi, f"P={row.purkinje_pvalue:.3g}", ha="right", va="center", fontsize=5.2)
    ax.set_xlim(-x_limit, x_limit); ax.set_xlabel("Purkinje-cell log2 fold change")
    clean(ax, grid=True)

    ax = fig.add_subplot(lower[0, 1]); panel(ax, "F", -0.26, 1.03)
    ca = data["ca_thresholds"].query("gene_symbol == 'CA12' and cell_type in ['Granule','Bergmann','Astrocytes','Purkinje']").copy()
    thresholds = [10, 20, 30]
    trace_colors = [muted(x, 0.34) for x in [c["bulk"], c["microbe"], c["genetics"], c["purkinje"]]]
    for color, (_, row) in zip(trace_colors, ca.iterrows()):
        effects = [row[f"logFC_{n}_nuclei"] for n in thresholds]
        ax.plot(thresholds, effects, marker="o", ms=4.0, lw=1.2, color=color)
        ax.text(30.8, effects[-1], row.cell_type, fontsize=4.6, ha="left", va="center", color=color)
    ax.axhline(0, color="#777777", lw=0.75, ls="--")
    ax.set_xticks(thresholds)
    ax.set_xlim(8, 39)
    ax.set_xlabel("Minimum nuclei per pseudobulk sample")
    ax.set_ylabel("CA12 log2 fold change")
    clean(ax)

    ax = fig.add_subplot(lower[0, 2]); panel(ax, "G", -0.26, 1.03)
    comp = data["composition"].copy()
    comp = comp[comp.nominal_p_lt_0_05.astype(str).str.lower().eq("true")].sort_values("logFC")
    y = np.arange(len(comp))
    colors = [c["negative"] if x < 0 else c["positive"] for x in comp.logFC]
    ax.errorbar(comp.logFC, y, xerr=[comp.logFC - comp["CI.L"], comp["CI.R"] - comp.logFC],
                fmt="none", ecolor="#4B4B4B", capsize=2.2, lw=0.9)
    ax.scatter(comp.logFC, y, s=50, fc=colors, ec="black", lw=0.4)
    ax.axvline(0, color="#777777", lw=0.8, ls="--")
    ax.set_yticks(y, comp.cell_type)
    ax.set_xlabel("Logit-scale effect\n(95% CI)")
    clean(ax, grid=True)
    return fig


def supplementary_s2(data, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 238 * MM))
    gs = fig.add_gridspec(3, 2, left=0.10, right=0.975, top=0.975, bottom=0.060,
                          hspace=0.34, wspace=0.32)
    summary = data["form_summary"].set_index("record_id").loc[FORM_ORDER]
    matrix = data["form_matrix"].copy()
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "A", -0.18, 1.03)
    distributions=[]
    for record in FORM_ORDER:
        values=pd.to_numeric(matrix[f"{record}_probability"],errors="coerce").dropna().to_numpy()
        distributions.append(values)
    vp=ax.violinplot(distributions,positions=np.arange(8),widths=.76,showextrema=False,showmedians=True)
    for body in vp["bodies"]:
        body.set_facecolor(muted(c["metabolite"],.46));body.set_edgecolor(muted(c["metabolite"],.30));body.set_alpha(.65)
    vp["cmedians"].set_color("black");vp["cmedians"].set_linewidth(.8)
    ax.set_xticks(range(8),FORM_SHORT,rotation=38,ha="right",fontsize=5);ax.set_ylabel("Predicted probability");clean(ax)
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "B", -0.18, 1.03)
    sets = pub._form_sets(matrix); jac = np.zeros((8,8))
    for i,a in enumerate(FORM_ORDER):
        for j,b in enumerate(FORM_ORDER): jac[i,j]=len(sets[a]&sets[b])/len(sets[a]|sets[b])
    im=ax.imshow(jac,cmap=LinearSegmentedColormap.from_list("j",["#F7F7F7",c["microbe"]]),vmin=.35,vmax=1)
    ax.set_xticks(range(8),FORM_SHORT,rotation=40,ha="right",fontsize=5);ax.set_yticks(range(8),FORM_SHORT,fontsize=5)
    fig.colorbar(im,ax=ax,fraction=.045,pad=.025,label="Jaccard index")
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "C", -0.18, 1.03)
    swiss=data["swiss"].copy().sort_values("probability",ascending=False).reset_index(drop=True)
    ax.scatter(np.arange(1,len(swiss)+1),swiss.probability,s=16,fc=c["genetics"],ec="none",alpha=.75)
    label_offsets={"CA9":(-18,14),"CA2":(-4,-12),"CA1":(12,14),"CA12":(18,-12),"CA3":(4,5),"CA7":(4,-9)}
    for gene in ["CA9","CA2","CA1","CA12","CA3","CA7"]:
        q=swiss[swiss.gene_symbol.eq(gene)]
        if not q.empty:
            ax.annotate(gene,(int(q.index[0])+1,float(q.probability.iloc[0])),xytext=label_offsets[gene],textcoords="offset points",fontsize=4.8,ha="center",va="center",arrowprops=dict(arrowstyle="-",lw=.35,color="#666666"))
    ax.set_xlabel("SwissTargetPrediction rank");ax.set_ylabel("Predicted probability");clean(ax)
    ax=fig.add_subplot(gs[1,1]);panel(ax,"D",-0.18,1.03);_upset(ax,matrix,c,16)
    ax=fig.add_subplot(gs[2,0]);panel(ax,"E",-0.18,1.03)
    core=[]
    for _,row in matrix.iterrows():
        if all(str(row[f"{r}_predicted"]).lower()=="true" for r in FORM_ORDER): core.append(row.gene_symbol)
    core=core[:25]; vals=np.full((len(core),8),np.nan)
    for i,g in enumerate(core):
        row=matrix[matrix.gene_symbol.eq(g)].iloc[0]
        for j,r in enumerate(FORM_ORDER): vals[i,j]=pd.to_numeric(row[f"{r}_probability"],errors="coerce")
    im=ax.imshow(np.ma.masked_invalid(vals),aspect="auto",cmap=LinearSegmentedColormap.from_list("p",["#F7F7F7",c["metabolite"]]))
    ax.set_yticks(range(len(core)),core,fontsize=4.5);ax.set_xticks(range(8),FORM_SHORT,rotation=40,ha="right",fontsize=5)
    fig.colorbar(im,ax=ax,fraction=.035,pad=.02,label="Predicted probability")
    ax=fig.add_subplot(gs[2,1]);panel(ax,"F",-0.18,1.03)
    exp=data["experimental_forms"].set_index("record_id").loc[FORM_ORDER]
    cols=["chembl_exact_molecule_count","chembl_activity_count","pubchem_bioassay_count","bindingdb_affinity_count","chembl_target_resolved_activity_count"]
    labels=["ChEMBL\nrecord","ChEMBL\nactivity","PubChem\nBioAssay","BindingDB\naffinity","Target-resolved\nactivity"]
    vals=exp[cols].apply(pd.to_numeric,errors="coerce").fillna(0).to_numpy(); shown=(vals>0).astype(float)
    ax.imshow(shown,aspect="auto",cmap=LinearSegmentedColormap.from_list("db",["#F2F2F2",c["genetics"]]),vmin=0,vmax=1)
    ax.set_xticks(range(5),labels,fontsize=5);ax.set_yticks(range(8),FORM_SHORT,fontsize=5)
    for i in range(8):
        for j in range(5): ax.text(j,i,str(int(vals[i,j])),ha="center",va="center",fontsize=4.7,color="white" if shown[i,j] else "black")
    return fig


def supplementary_s4(data, c):
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, 260 * MM))
    gs = fig.add_gridspec(4, 2, left=0.105, right=0.975, top=0.975, bottom=0.055,
                          hspace=0.42, wspace=0.34)
    ax=fig.add_subplot(gs[0,0]);panel(ax,"A",-0.20,1.03)
    cell=data["cell"].copy(); genes=["CA1","CA2","CA3","CA4","CA7","CA9","CA12","CA13"]
    q=cell[cell.gene_symbol.isin(genes)].copy(); cell_types=list(dict.fromkeys(q.annotation_value)); xmap={v:i for i,v in enumerate(cell_types)};ymap={v:i for i,v in enumerate(genes)}
    sc=ax.scatter([xmap[v] for v in q.annotation_value],[ymap[v] for v in q.gene_symbol],s=8+q.percent_expressing*1.4,c=q.mean_expression,cmap=LinearSegmentedColormap.from_list("expr",["#F2F2F2",c["purkinje"]]),ec="black",lw=.2)
    ax.set_xticks(range(len(cell_types)),[base._cell_abbreviation(v) for v in cell_types],rotation=45,ha="right",fontsize=4.8);ax.set_yticks(range(len(genes)),genes);fig.colorbar(sc,ax=ax,fraction=.04,pad=.02,label="Mean expression")
    ax=fig.add_subplot(gs[0,1]);panel(ax,"B",-0.20,1.03)
    models=data["model_summary"].copy()
    sc=ax.scatter(models.donors,models.pseudobulk_samples,
                  s=22+75*(models.genes_after_filterByExpr-models.genes_after_filterByExpr.min())/(models.genes_after_filterByExpr.max()-models.genes_after_filterByExpr.min()),
                  c=models.final_consensus_correlation,cmap=LinearSegmentedColormap.from_list("cor",["#ECECEC",muted(c["microbe"],.30)]),ec="black",lw=.35)
    ordered=models.sort_values("pseudobulk_samples").reset_index(drop=True);label_y=[]
    for value in ordered.pseudobulk_samples:
        label_y.append(max(float(value),label_y[-1]+1.75 if label_y else float(value)))
    if label_y[-1]>141.0:
        shift=label_y[-1]-141.0;label_y=[value-shift for value in label_y]
    for (_,row),target_y in zip(ordered.iterrows(),label_y):
        ax.plot([row.donors,109.35],[row.pseudobulk_samples,target_y],color="#8A8A8A",lw=.35)
        ax.text(109.45,target_y,base._cell_abbreviation(row.cell_type),fontsize=3.8,ha="left",va="center")
    ax.set_xlim(101.7,110.5);ax.set_ylim(103.5,142.0)
    ax.set_xlabel("Donors");ax.set_ylabel("Pseudobulk samples");fig.colorbar(sc,ax=ax,fraction=.04,pad=.02,label="Consensus correlation");clean(ax,grid=True)
    ax=fig.add_subplot(gs[1,0]);panel(ax,"C",-0.20,1.03)
    th=data["target_thresholds"].copy(); counts=[]
    for ct,g in th.groupby("cell_type"):
        counts.append((ct,sum(g.nominal_p_lt_0_05_10_nuclei.astype(str).str.lower().eq("true")),sum(g.nominal_p_lt_0_05_20_nuclei.astype(str).str.lower().eq("true")),sum(g.nominal_p_lt_0_05_30_nuclei.astype(str).str.lower().eq("true"))))
    ct=pd.DataFrame(counts,columns=["cell_type","10","20","30"]).sort_values("10")
    y=np.arange(len(ct));ax.hlines(y,ct[["10","20","30"]].min(axis=1),ct[["10","20","30"]].max(axis=1),color="#A8A8A8",lw=1.1)
    for offset,(col,color) in zip([-.14,0,.14],zip(["10","20","30"],[muted(c["microbe"],.32),muted(c["metabolite"],.32),muted(c["support"],.42)])):
        ax.scatter(ct[col],y+offset,s=31,fc=color,ec="black",lw=.3,label=f"{col} nuclei")
    ax.set_yticks(y,ct.cell_type,fontsize=5);ax.set_xlabel("Nominal gene–cell-type results");ax.legend(frameon=False,ncol=3,fontsize=5);clean(ax,grid=True)
    ax=fig.add_subplot(gs[1,1]);panel(ax,"D",-0.20,1.03)
    ca=data["ca_thresholds"].copy();ca_genes=list(dict.fromkeys(ca.gene_symbol));ca_cells=list(dict.fromkeys(ca.cell_type));gx={v:i for i,v in enumerate(ca_cells)};gy={v:i for i,v in enumerate(ca_genes)}
    ca["mean_logFC"]=ca[["logFC_10_nuclei","logFC_20_nuclei","logFC_30_nuclei"]].mean(axis=1)
    vmax=max(abs(ca.mean_logFC.min()),abs(ca.mean_logFC.max()));cm=LinearSegmentedColormap.from_list("signed_ca",[c["negative"],"#F5F5F5",c["positive"]])
    sc=ax.scatter([gx[v] for v in ca.cell_type],[gy[v] for v in ca.gene_symbol],s=12+22*ca.nominal_p_lt_0_05_threshold_count,c=ca.mean_logFC,cmap=cm,vmin=-vmax,vmax=vmax,ec="black",lw=.2)
    ax.set_xticks(range(len(ca_cells)),[base._cell_abbreviation(v) for v in ca_cells],rotation=45,ha="right",fontsize=4.6);ax.set_yticks(range(len(ca_genes)),ca_genes,fontsize=4.8);fig.colorbar(sc,ax=ax,fraction=.04,pad=.02,label="Mean log2 fold change")
    ax=fig.add_subplot(gs[2,0]);panel(ax,"E",-0.20,1.03)
    loo=data["loo_composition"].copy(); nominal=data["composition"].query("nominal_p_lt_0_05 == True").cell_type.tolist();loo=loo[loo.cell_type.isin(nominal)]
    primary=data["composition"][["cell_type","logFC"]].rename(columns={"logFC":"primary_logFC"});loo=loo.merge(primary,on="cell_type",how="left").sort_values("primary_logFC")
    y=np.arange(len(loo));ax.hlines(y,loo.minimum_logFC,loo.maximum_logFC,color=c["negative"],lw=2);ax.scatter(loo.primary_logFC,y,s=45,fc=c["positive"],ec="black",lw=.35)
    ax.axvline(0,color="#777",ls="--",lw=.8);ax.set_yticks(y,loo.cell_type,fontsize=5);ax.set_xlabel("Leave-one-donor-out logit effect range");clean(ax,grid=True)
    ax=fig.add_subplot(gs[2,1]);panel(ax,"F",-0.20,1.03)
    cam=data["camera_thresholds"].copy(); cam=cam[cam.correlation_method.eq("fixed_0.01")]
    sets=list(dict.fromkeys(cam.gene_set)); cells=list(dict.fromkeys(cam.cell_type)); mat=np.full((len(sets),len(cells)),np.nan)
    for i,s in enumerate(sets):
        for j,cty in enumerate(cells):
            q=cam[(cam.gene_set.eq(s))&(cam.cell_type.eq(cty))]
            if not q.empty: mat[i,j]=-math.log10(max(float(q.PValue_10_nuclei.iloc[0]),1e-300))
    im=ax.imshow(mat,aspect="auto",cmap=LinearSegmentedColormap.from_list("cam",["#F7F7F7",c["genetics"]]))
    ax.set_yticks(range(len(sets)),[s.replace("_"," ") for s in sets],fontsize=4.6);ax.set_xticks(range(len(cells)),cells,rotation=48,ha="right",fontsize=4.4);fig.colorbar(im,ax=ax,fraction=.035,pad=.02,label="−log10(P)")
    ax=fig.add_subplot(gs[3,:]);panel(ax,"G",-0.09,1.03)
    ev=data["composition_evidence"].copy(); ev=ev[ev.cell_type.isin(nominal)].copy(); methods=[("primary_nominal_p_lt_0_05","Propeller"),("donor_total_nominal_p_lt_0_05","Donor aggregate"),("matched_nominal_p_lt_0_05","Matched donors"),("sccoda_hmc_da_all_credible_reference_models","scCODA credible models")]
    ax.set_xlim(-.5,len(methods)-.5);ax.set_ylim(len(ev)-.5,-.5)
    for i,row in ev.reset_index(drop=True).iterrows():
        for j,(col,label) in enumerate(methods):
            value=row[col]
            if "credible" in col: size=28+18*float(value); filled=float(value)>0
            else: size=52; filled=str(value).lower()=="true"
            ax.scatter(j,i,s=size,fc=c["cell"] if filled else "white",ec="black",lw=.4)
    ax.set_xticks(range(len(methods)),[x[1] for x in methods]);ax.set_yticks(range(len(ev)),ev.cell_type);ax.tick_params(length=0);ax.spines[:].set_visible(False)
    return fig


def _make_output_dirs(palette_name: str):
    root = ROOT / ("02_figures_palette_09" if palette_name == "palette_09" else "03_figures_palette_10")
    main = root / "Main_Figures"
    supp = root / "Supplementary_Figures"
    main.mkdir(parents=True, exist_ok=True); supp.mkdir(parents=True, exist_ok=True)
    return main, supp


def run() -> None:
    setup_style()
    data = load_data()
    rebuild_stems = {
        "Supplementary_Figure_S4_single_nucleus_sensitivity",
    }
    for palette_name, colors in PALETTES.items():
        main_dir, supp_dir = _make_output_dirs(palette_name)
        figures = [
            (figure1, data, colors, main_dir / "Figure1_integrated_study_design"),
            (pub.figure2, data, colors, main_dir / "Figure2_mendelian_randomization"),
            (figure3, data, colors, main_dir / "Figure3_chemical_forms_and_target_prediction"),
            (pub.figure4, data, colors, main_dir / "Figure4_gene_association_and_regulatory_evidence"),
            (figure5, data, colors, main_dir / "Figure5_cerebellar_case_control_evidence"),
            (base.supplementary_figure_s1, data, colors, supp_dir / "Supplementary_Figure_S1_MR_diagnostics"),
            (supplementary_s2, data, colors, supp_dir / "Supplementary_Figure_S2_target_prediction_sensitivity"),
            (base.supplementary_figure_s3, data, colors, supp_dir / "Supplementary_Figure_S3_gene_association_complete_results"),
            (supplementary_s4, data, colors, supp_dir / "Supplementary_Figure_S4_single_nucleus_sensitivity"),
        ]
        for func, d, c, stem in figures:
            if stem.name not in rebuild_stems and stem.with_suffix(".pdf").exists() and stem.with_suffix(".tiff").exists():
                print(f"Keeping existing {palette_name}: {stem.name}", flush=True)
                continue
            print(f"Building {palette_name}: {stem.name}", flush=True)
            fig = func(d, c)
            pub._remove_panel_titles(fig)
            base.save_formal(fig, stem, dpi=600)
            plt.close(fig)


if __name__ == "__main__":
    run()
