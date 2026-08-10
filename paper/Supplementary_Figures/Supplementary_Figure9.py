#%%
# ---------------------------------------------------------------------------
# UpSet plot per chip: how the CpG-eliminating variants of that chip are shared
# among the six cohorts (AFR, AMR, EAS, EUR, SAS, Korea10K).
#
# Input : Data/Chip_Overlap/<POP>/<CHIP>.CpG_disappeared.AC1.ChipOverlap.tsv
#         Already filtered to AC >= 1, so a variant appears in a cohort's file
#         only if it segregates there.
# Output: Figures/Figure.UpSet.Population_Overlap.png / .pdf  (one panel per chip)
#         Figures/Figure.UpSet.Intersections.tsv   (every intersection size)
#
# Each column of the matrix is one exact combination of cohorts; the bar above it
# counts the variants carried by that combination and by no other cohort. Columns
# are grouped by degree - private (1 cohort) first, then 2, 3, 4, 5 and finally
# the variants common to all 6 - and sorted by size inside each degree.
#
# Caveat carried from the source data: the five 1KGP superpopulations share one
# call set, so absence there means AC = 0. Korea10K is an independent call set, so
# a variant missing from the 1KGP files may simply never have been called there.
# Degree-1 counts for Korea10K therefore mix real privacy with call-set
# differences and should be read as an upper bound.
# ---------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from korea10k.config import PROJECT_DIR

dir_chip_overlap = f"{PROJECT_DIR}/Analysis/Revision/Draw_Figure/Data/Chip_Overlap"
dir_figure_out = f"{PROJECT_DIR}/Analysis/Revision/Draw_Figure/Figures"
os.makedirs(dir_figure_out, exist_ok=True)

path_out_table = os.path.join(dir_figure_out, "Figure.UpSet.Intersections.tsv")

# Directory name -> label. The order fixes the rows of the matrix.
list_pop = [
    ("AFR", "AFR"),
    ("AMR", "AMR"),
    ("EAS", "EAS"),
    ("EUR", "EUR"),
    ("SAS", "SAS"),
    ("Korea10K", "Korea10K"),
]
list_chips = ["HM27", "HM450", "EPICv2", "MSA"]
dict_chip_to_title = {
    "HM27": "Human Methylation 27K (HM27)",
    "HM450": "Human Methylation 450K (HM450)",
    "EPICv2": "Human Methylation EPIC v2.0 (EPICv2)",
    "MSA": "Human Methylation Screening Array (MSA)",
}

# Blue ordinal ramp, light (private) to dark (shared by every cohort).
list_color_degree = ["#86b6ef", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

color_surface = "#ffffff"
color_ink = "#0b0b0b"
color_grid = "#e1e0d9"
color_dot_off = "#dcdbd4"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10.5,
    "axes.labelweight": "bold",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.edgecolor": color_ink,
    "axes.labelcolor": color_ink,
    "text.color": color_ink,
    "xtick.color": color_ink,
    "ytick.color": color_ink,
    "figure.facecolor": color_surface,
    "axes.facecolor": color_surface,
    "savefig.facecolor": color_surface,
    "pdf.fonttype": 42,
})


#%%
def load_variant_keys(dirname, chip):
    """Unique CHROM:POS:REF:ALT keys of the variants this cohort carries."""
    path = os.path.join(dir_chip_overlap, dirname, f"{chip}.CpG_disappeared.AC1.ChipOverlap.tsv")
    table = pd.read_csv(path, sep='\t', usecols=["#CHROM", "POS", "REF", "ALT"])
    series_key = (table["#CHROM"].astype(str) + ":" + table["POS"].astype(str)
                  + ":" + table["REF"].astype(str) + ":" + table["ALT"].astype(str))
    return series_key.unique()


def build_intersections(chip):
    """Size of every non-empty cohort combination, plus each cohort's total."""
    series_mask = pd.Series(dtype="int64")
    dict_pop_to_total = dict()
    for index_pop, (dirname, pop) in enumerate(list_pop):
        array_key = load_variant_keys(dirname, chip)
        dict_pop_to_total[pop] = len(array_key)
        series_mask = series_mask.add(
            pd.Series(1 << index_pop, index=array_key), fill_value=0
        )
    series_mask = series_mask.astype(int)

    series_size = series_mask.value_counts()
    table = pd.DataFrame({"Mask": series_size.index.astype(int), "N_Variant": series_size.to_numpy()})
    table["Degree"] = [bin(mask).count("1") for mask in table["Mask"]]
    table["Combination"] = [
        "+".join(pop for index_pop, (_d, pop) in enumerate(list_pop) if mask >> index_pop & 1)
        for mask in table["Mask"]
    ]
    table = table.sort_values(["Degree", "N_Variant"], ascending=[True, False]).reset_index(drop=True)
    return table, dict_pop_to_total


def format_count(value, _pos=None):
    if value >= 1e6:
        return f"{value / 1e6:g}M"
    if value >= 1e3:
        return f"{value / 1e3:g}k"
    return f"{value:g}"


#%%
def draw_upset_panel(fig, spec, chip, table_inter, dict_pop_to_total):
    """One UpSet panel drawn into an outer gridspec slot."""
    n_combination = len(table_inter)
    n_pop = len(list_pop)

    grid = spec.subgridspec(
        2, 2, width_ratios=[1.15, 4.3], height_ratios=[2.5, 1.7],
        hspace=0.06, wspace=0.20,
    )
    ax_bar = fig.add_subplot(grid[0, 1])
    ax_dot = fig.add_subplot(grid[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(grid[1, 0], sharey=ax_dot)
    ax_note = fig.add_subplot(grid[0, 0])

    array_x = np.arange(n_combination)
    array_color = [list_color_degree[degree - 1] for degree in table_inter["Degree"]]

    # --- intersection sizes ---
    ax_bar.bar(array_x, table_inter["N_Variant"], width=0.72, color=array_color, zorder=3)
    y_max = float(table_inter["N_Variant"].max())
    ax_bar.set_ylim(0, y_max * 1.10)
    ax_bar.set_ylabel("Variants", fontsize=9)
    ax_bar.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(format_count))
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.grid(axis='y', color=color_grid, linewidth=0.6, zorder=0)
    ax_bar.set_axisbelow(True)
    ax_bar.tick_params(axis='x', length=0, labelbottom=False)
    ax_bar.tick_params(axis='y', labelsize=8)

    # --- membership matrix ---
    for index_pop in range(n_pop):
        if index_pop % 2 == 0:
            ax_dot.axhspan(index_pop - 0.5, index_pop + 0.5, color="#f4f3ef", zorder=0)
    ax_dot.scatter(np.repeat(array_x, n_pop), np.tile(np.arange(n_pop), n_combination),
                   s=11, color=color_dot_off, zorder=2)
    for x_value, mask, degree in zip(array_x, table_inter["Mask"], table_inter["Degree"]):
        list_member = [index_pop for index_pop in range(n_pop) if mask >> index_pop & 1]
        ax_dot.scatter([x_value] * len(list_member), list_member, s=11,
                       color=list_color_degree[degree - 1], zorder=3)
        if len(list_member) > 1:
            ax_dot.plot([x_value, x_value], [min(list_member), max(list_member)],
                        color=list_color_degree[degree - 1], linewidth=0.9, zorder=2)

    ax_dot.set_yticks(np.arange(n_pop))
    ax_dot.set_yticklabels([pop for _d, pop in list_pop], fontsize=8)
    ax_dot.set_ylim(n_pop - 0.5, -0.5)
    ax_dot.set_xlim(-0.7, n_combination - 0.3)
    for name_spine in ("top", "right", "bottom", "left"):
        ax_dot.spines[name_spine].set_visible(False)
    ax_dot.tick_params(axis='both', length=0, labelbottom=False)
    ax_dot.tick_params(axis='y', pad=4)

    # --- per-cohort totals ---
    array_total = np.array([dict_pop_to_total[pop] for _d, pop in list_pop], dtype=float)
    ax_set.barh(np.arange(n_pop), array_total, height=0.55, color="#7d7b72", zorder=3)
    for index_pop, total in enumerate(array_total):
        ax_set.text(total * 1.05, index_pop, format_count(total), ha="right", va="center",
                    fontsize=7.5, color=color_ink)
    ax_set.set_xlim(array_total.max() * 1.75, 0)
    ax_set.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(format_count))
    ax_set.spines["top"].set_visible(False)
    ax_set.spines["left"].set_visible(False)
    ax_set.grid(axis='x', color=color_grid, linewidth=0.6, zorder=0)
    ax_set.set_axisbelow(True)
    ax_set.tick_params(axis='y', length=0, labelleft=False)
    ax_set.tick_params(axis='x', labelsize=8)

    # --- totals by degree, in the free corner ---
    ax_note.axis("off")
    table_degree = table_inter.groupby("Degree", as_index=False)["N_Variant"].sum()
    n_union = int(table_inter["N_Variant"].sum())
    ax_note.text(0.0, 1.0, f"{n_union:,} variants\nin the union",
                 transform=ax_note.transAxes, ha="left", va="top",
                 fontsize=8.5, fontweight="bold", linespacing=1.4, color=color_ink)
    for index_line, (degree, n_variant) in enumerate(zip(table_degree["Degree"],
                                                         table_degree["N_Variant"])):
        y_line = 0.66 - 0.132 * index_line
        ax_note.add_patch(plt.Rectangle(
            (0.0, y_line - 0.030), 0.10, 0.060, transform=ax_note.transAxes,
            facecolor=list_color_degree[degree - 1], edgecolor="none", clip_on=False))
        ax_note.text(0.16, y_line, f"{degree}", transform=ax_note.transAxes,
                     ha="left", va="center", fontsize=8, color=color_ink)
        ax_note.text(0.93, y_line, f"{n_variant:,} ({100.0 * n_variant / n_union:.1f}%)",
                     transform=ax_note.transAxes, ha="right", va="center",
                     fontsize=8, color=color_ink)

    return ax_bar, ax_note


#%%
# --- Compose one figure with a panel per chip -------------------------------
list_table = list()
dict_chip_to_result = dict()
for chip in list_chips:
    table_inter, dict_pop_to_total = build_intersections(chip)
    dict_chip_to_result[chip] = (table_inter, dict_pop_to_total)
    table_out = table_inter.copy()
    table_out.insert(0, "Chip", chip)
    list_table.append(table_out[["Chip", "Degree", "Combination", "N_Variant"]])
    table_degree = table_inter.groupby("Degree")["N_Variant"].sum()
    print(f"{chip:7s}: union {int(table_inter['N_Variant'].sum()):>8,}  "
          f"by degree {dict(table_degree.astype(int))}")

n_row = 2
n_col = 2
fig = plt.figure(figsize=(17.5, 10.0))
grid_outer = fig.add_gridspec(
    n_row, n_col,
    left=0.040, right=0.985, top=0.925, bottom=0.055, hspace=0.30, wspace=0.10,
)

for index_chip, chip in enumerate(list_chips):
    row, col = divmod(index_chip, n_col)
    table_inter, dict_pop_to_total = dict_chip_to_result[chip]
    ax_bar, ax_note = draw_upset_panel(
        fig, grid_outer[row, col], chip, table_inter, dict_pop_to_total,
    )
    box_note = ax_note.get_position()
    box_bar = ax_bar.get_position()
    fig.text(box_note.x0 - 0.030, box_bar.y1 + 0.014, "ABCD"[index_chip],
             fontsize=14, fontweight="bold", ha="left", va="bottom", color=color_ink)
    fig.text((box_note.x0 + box_bar.x1) / 2, box_bar.y1 + 0.014, dict_chip_to_title[chip],
             fontsize=12, fontweight="bold", ha="center", va="bottom", color=color_ink)

path_out_png = os.path.join(dir_figure_out, "Figure.UpSet.Population_Overlap.png")
path_out_pdf = os.path.join(dir_figure_out, "Figure.UpSet.Population_Overlap.pdf")
fig.savefig(path_out_png, dpi=300)
fig.savefig(path_out_pdf)
print(f"\nWrote {path_out_png}")
print(f"Wrote {path_out_pdf}")

table_all = pd.concat(list_table, ignore_index=True)
table_all.to_csv(path_out_table, sep='\t', index=False)
print(f"Wrote {path_out_table}")

# %%
