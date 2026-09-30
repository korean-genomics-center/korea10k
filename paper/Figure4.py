#%%
# ---------------------------------------------------------------------------
# Figure: CpG sites lost to sequence variation on methylation arrays,
#         across 1KGP superpopulations and Korea10K.
#         COMMON-FIRST variant of draw_figure_chip_overlap.py.
#
#   A) Per-CpG mean beta against the genotype of the disrupting variant
#   B) Human Methylation 27K       | left : # disrupted CpG sites by AF class
#   C) Human Methylation 450K      | right: cumulative # disrupted CpG sites
#   D) Human Methylation EPIC v2.0 |        vs. allele-frequency cut-off
#   E) Human Methylation MSA       |
#   F) Aging-clock CpG markers lost, per superpopulation and per clock
#
# Difference from draw_figure_chip_overlap.py - both panels of B-E now read from
# common to rare instead of rare to common:
#   - the stacked bars stack the COMMON classes at the BOTTOM, so the singleton
#     class sits on top;
#   - the cumulative curves count sites whose most frequent disrupting variant
#     has AF >= the cut-off, and the x-axis runs from AF 1 down to AF 3e-5, so
#     the curve accumulates common sites first and still ends at the bar total.
# Everything else (panels A and F, palettes, layout) is unchanged.
#
# Run prepare_chip_overlap_figure_data.py first (writes Figure_Source/*.tsv).
#
# Cohort hues (curves of B-E) are categorical slots 1-6 of the reference
# data-viz palette, in fixed order. The four AF-class hues are specified by the
# study and used verbatim; panel F colours its two segments with the classes
# either side of the AF 0.01 cut, and names the cohort under each bar instead.
# ---------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from korea10k.config import PROJECT_DIR

dir_root = f"{PROJECT_DIR}/Analysis/Revision/Draw_Figure_ver20260929"
dir_chip_overlap = os.path.join(dir_root, "Data", "Chip_Overlap")
dir_figure_source = os.path.join(dir_chip_overlap, "Figure_Source")
dir_figure_out = os.path.join(dir_root, "Figures")
os.makedirs(dir_figure_out, exist_ok=True)

path_class = os.path.join(dir_figure_source, "AF_Class_Count.tsv")
path_cumul = os.path.join(dir_figure_source, "AF_Cumulative.tsv")
path_summary_chip = os.path.join(dir_chip_overlap, "Summary.CpG_disappeared.AC1.ChipOverlap.tsv")
path_summary_clock = os.path.join(dir_chip_overlap, "Summary.CpG_disappeared.AC1.AgeClockOverlap.tsv")
path_genotype_beta = os.path.join(dir_root, "Data", "Methylation_Change", "per_cpg_beta_by_genotype.tsv")

path_out_png = os.path.join(dir_figure_out, "Figure4.png")
path_out_pdf = os.path.join(dir_figure_out, "Figure4.pdf")

#%%
# --- Palette / cohort constants --------------------------------------------
# Cohorts are always drawn in this order, with these hues, in every panel.
list_pops = ["AFR", "AMR", "EAS", "EUR", "SAS", "KOREAN"]
dict_pop_to_label = {
    "AFR": "AFR", "AMR": "AMR", "EAS": "EAS",
    "EUR": "EUR", "SAS": "SAS", "KOREAN": "Korea10K",
}
# Cohort hues are deliberately deeper and less saturated than the AF-class hues
# of the bars, and sit in different hue families, so the two encodings never read
# as the same key. Chosen by search under hard gates - among themselves (they
# share one axes) OKLab dE >= 15 normal and >= 8 under simulated deuter-, prot-
# and tritanopia, minimum hue separation 25 deg - while maximising distance to
# the four AF-class colours. Worst cohort-vs-AF distance is now 13.7 (SAS vs
# Very common); the previous slots put Korea10K green 0.6 from the Rare green
# and AMR orange 6.4 from the Singleton orange.
dict_pop_to_color = {
    "AFR": "#448bb9",     # mid blue
    "AMR": "#74301e",     # dark rust
    "EAS": "#145157",     # dark teal
    "EUR": "#926422",     # ochre
    "SAS": "#c45b80",     # rose
    "KOREAN": "#753e92",  # purple - the focal cohort, also drawn thicker
}
# Diploid genotype counts (OBS_CT) observed in the source tables -> cohort size.
# The AF resolution floor is 1/OBS_CT and differs ~20-fold between the 1KGP
# superpopulations and Korea10K, which is why the rarest class is defined by
# allele count instead of a frequency cut. Stated on the figure.
dict_pop_to_n_sample = {
    "AFR": 893, "AMR": 490, "EAS": 585,
    "EUR": 633, "SAS": 601, "KOREAN": 10014,
}

# The rarest class is keyed on allele count, so it means the same thing in every
# cohort; a fixed AF cut there would be unreachable for the smallest cohort.
list_af_class = [
    "Singleton (AC=1)",
    "Rare (AC>1, AF<=0.01)",
    "Common (0.01<AF<=0.05)",
    "Very common (AF>0.05)",
]
dict_af_class_to_label = {
    "Singleton (AC=1)": "Singleton (AC = 1)",
    "Rare (AC>1, AF<=0.01)": "Rare (AC > 1, AF $\\leq$ 0.01)",
    "Common (0.01<AF<=0.05)": "Common (0.01 < AF $\\leq$ 0.05)",
    "Very common (AF>0.05)": "Very common (AF > 0.05)",
}
dict_af_class_to_color = {
    "Singleton (AC=1)": "#ff7b00",
    "Rare (AC>1, AF<=0.01)": "#0B8100",
    "Common (0.01<AF<=0.05)": "#0026ff",
    # Magenta rather than the violet-blue originally specified: against
    # Common #0026ff the violet measured OKLab dE 6.7 (normal vision), far under
    # the 15 floor, so the two top segments of every stacked bar read as one
    # block. #c400a4 lifts that pair to 34.0 normal / 18.1 worst-case CVD, and
    # every other pair of the four classes clears the gates too.
    "Very common (AF>0.05)": "#c400a4",
}

# Bottom-to-top stacking order of the bars in panels B-E: the most common class
# sits on the baseline and the singleton class caps the bar, so the bars read
# common -> rare upward, in step with the cumulative curves beside them.
list_af_class_stack = list(reversed(list_af_class))

# Panel E keeps the single AF 0.01 split carried by the summary table, coloured
# with the class either side of that cut in panels A-D: green below, blue above.
color_af_rare = dict_af_class_to_color["Rare (AC>1, AF<=0.01)"]
color_af_common = dict_af_class_to_color["Common (0.01<AF<=0.05)"]

list_chips = ["HM27", "HM450", "EPICv2", "MSA"]
dict_chip_to_title = {
    "HM27": "Human Methylation 27K (HM27)",
    "HM450": "Human Methylation 450K (HM450)",
    "EPICv2": "Human Methylation EPIC v2.0 (EPICv2)",
    "MSA": "Human Methylation Screening Array (MSA)",
}

list_clocks = ["Horvath", "Horvath SB", "Hannum", "PhenoAge", "GrimAge", "DunedinPACE"]

# Three distinct hues that also read as a progression. Within the panel they
# clear the categorical gates (OKLab dE 25.0 normal, 15.0 worst-case CVD); the
# nearest colour in any other key of the figure is 9.4 away, which is acceptable
# because panel A stands alone and its x-axis already names each genotype.
list_color_genotype = ["#00a6a6", "#e0b000", "#7b1f5c"]
# A CpG site enters panel A only when every genotype class has at least this many
# samples, so the three distributions describe the same set of sites.
n_min_genotype = 10

# Chart chrome.
color_surface = "#ffffff"
color_ink = "#0b0b0b"
# Every piece of text in the figure is black; the greys below are rules only.
color_grid = "#e1e0d9"
color_axis = "#0b0b0b"       # axis spines and tick marks, black
color_guide = "#c3c2b7"      # AF cut-off guide lines, kept recessive

# Colour keys. Reviewer: "the color scale in Figure 4 is quite small and hard to
# read" - the keys were 8 pt with 1.1-em swatches; they are now set well above
# the tick labels (~9 pt once the 12-in figure is printed at 180 mm), with
# swatches and line samples large enough to read the hue.
size_key_text = 15.0
size_key_title = 16.0
size_key_text_panel = 14.0   # key inside panel F, where it shares space with bars

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 10.5,
    "axes.labelweight": "bold",
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 8.5,
    "axes.edgecolor": color_axis,
    "axes.labelcolor": color_ink,
    "text.color": color_ink,
    "xtick.color": color_axis,          # tick marks share the axis line colour
    "ytick.color": color_axis,
    "xtick.labelcolor": color_ink,
    "ytick.labelcolor": color_ink,
    "figure.facecolor": color_surface,
    "axes.facecolor": color_surface,
    "savefig.facecolor": color_surface,
    "pdf.fonttype": 42,
})


#%%
# --- Load ------------------------------------------------------------------
table_class = pd.read_csv(path_class, sep='\t')
table_cumul = pd.read_csv(path_cumul, sep='\t')
table_summary_chip = pd.read_csv(path_summary_chip, sep='\t')
table_summary_clock = pd.read_csv(path_summary_clock, sep='\t')

# Total CpG sites carried by each array -> denominator of the "% of array lost"
# annotation. Constant per chip, so read it off any single superpopulation row.
dict_chip_to_n_site = (
    table_summary_chip.drop_duplicates("Chip").set_index("Chip")["N_Chip_CpG_Site"].to_dict()
)

# Panel A-D bar heights: N_CpG_Site indexed by (chip, pop, class).
pivot_class = table_class.pivot_table(
    index=["Chip", "Superpopulation"], columns="AF_Class", values="N_CpG_Site"
).fillna(0.0)

# Panel E: aging-clock markers lost, split at AF 0.01.
table_clock = table_summary_clock.copy()
table_clock["Pct_Common"] = 100.0 * table_clock["N_Disappeared_Common"] / table_clock["N_Marker_Total"]
table_clock["Pct_Rare"] = 100.0 * table_clock["N_Disappeared_Rare"] / table_clock["N_Marker_Total"]
dict_clock_to_n_marker = (
    table_clock.drop_duplicates("AgeClock").set_index("AgeClock")["N_Marker_Total"].to_dict()
)

# Panel A: per-CpG mean beta by genotype of the disrupting variant.
table_genotype = pd.read_csv(
    path_genotype_beta, sep='\t', na_values=["NA"],
    usecols=["n_gt0", "mean_gt0", "n_gt1", "mean_gt1", "n_gt2", "mean_gt2"],
)
flag_complete = (
    (table_genotype["n_gt0"] >= n_min_genotype)
    & (table_genotype["n_gt1"] >= n_min_genotype)
    & (table_genotype["n_gt2"] >= n_min_genotype)
)
table_genotype = table_genotype[flag_complete]
array_beta_by_gt = {gt: table_genotype[f"mean_gt{gt}"].to_numpy() for gt in (0, 1, 2)}
print(f"Panel A: {len(table_genotype):,} CpG sites with >= {n_min_genotype} samples in all genotypes")
for gt in (0, 1, 2):
    array_beta = array_beta_by_gt[gt]
    print(f"  genotype {gt}: mean {array_beta.mean():.4f}  median {np.median(array_beta):.4f}")


#%%
# --- Small helpers ---------------------------------------------------------
def style_axis(ax, grid_axis='y'):
    """Recessive chrome: hairline grid behind the marks, no top/right spines."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis=grid_axis, color=color_grid, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, width=0.7)


def format_count(value, _pos=None):
    """Axis ticks as 12k / 1.2M rather than raw digits."""
    if value >= 1e6:
        return f"{value / 1e6:g}M"
    if value >= 1e3:
        return f"{value / 1e3:g}k"
    return f"{value:g}"


def place_end_labels(ax, list_label_y_color, x_axes=1.17, min_gap_frac=0.075):
    """Direct-label line ends, nudged apart so no two labels collide, with a
    leader arrow from each label back to the end of the curve it names.

    The arrow is what disambiguates a nudged label: several cohorts end within a
    few hundred sites of each other, so the label a curve owns is often not the
    one closest to it. x is taken in axes fractions (1.0 = the right spine, where
    every curve ends) and y in data units, so the geometry holds at any panel
    size; the arrow is drawn in the curve's own hue and the text stays black.
    The labels share one x, so every arrow runs inside the gutter between the
    spine and that column and no arrow can cross another label's text.
    """
    y_low, y_high = ax.get_ylim()
    min_gap = (y_high - y_low) * min_gap_frac
    list_sorted = sorted(list_label_y_color, key=lambda item: item[1])
    list_placed = list()
    y_prev = -np.inf
    for label, y_value, color in list_sorted:
        y_draw = max(y_value, y_prev + min_gap)
        list_placed.append((label, y_value, y_draw, color))
        y_prev = y_draw

    transform_end = blended_transform_factory(ax.transAxes, ax.transData)
    for label, y_value, y_draw, color in list_placed:
        ax.annotate(
            label,
            xy=(1.0, y_value), xycoords=transform_end,
            xytext=(x_axes, y_draw), textcoords=transform_end,
            ha="left", va="center", fontsize=9.5, color=color_ink,
            annotation_clip=False, zorder=6,
            # The white box is the label's own footprint: matplotlib uses it as
            # patchA, so the tail stops at the edge of the text instead of
            # running under the glyphs, and anything that still reaches it is
            # painted over. Leaving patchA unset (or None) is what made the
            # arrow cross the text. The box sits in the gutter outside the
            # axes, so it hides no data or grid.
            bbox=dict(boxstyle="square,pad=0.18", facecolor=color_surface,
                      edgecolor="none"),
            arrowprops=dict(arrowstyle="-|>", color=color, linewidth=1.0,
                            mutation_scale=8, shrinkA=0.0, shrinkB=0.5,
                            connectionstyle="arc3,rad=0"),
        )
    return list_placed


#%%
# --- Panel painters --------------------------------------------------------
def draw_af_class_bars(ax, chip):
    """Left half of panels B-E: disrupted CpG sites per cohort, split by AF class.
    Stacked most common at the bottom, singletons on top."""
    array_x = np.arange(len(list_pops))
    array_bottom = np.zeros(len(list_pops))
    n_site_chip = dict_chip_to_n_site[chip]

    for name_class in list_af_class_stack:
        array_height = np.array([pivot_class.loc[(chip, pop), name_class] for pop in list_pops])
        ax.bar(array_x, array_height, bottom=array_bottom, width=0.68,
               color=dict_af_class_to_color[name_class],
               edgecolor=color_surface, linewidth=1.2, zorder=3)
        array_bottom = array_bottom + array_height

    # Share of the whole array lost, as text - not a second y-axis.
    for x_value, total in zip(array_x, array_bottom):
        ax.text(x_value, total + array_bottom.max() * 0.025,
                f"{100.0 * total / n_site_chip:.1f}%",
                ha="center", va="bottom", fontsize=8.5, color=color_ink)

    ax.set_xticks(array_x)
    ax.set_xticklabels([dict_pop_to_label[pop] for pop in list_pops], rotation=30, ha="right")
    for tick_label, pop in zip(ax.get_xticklabels(), list_pops):
        tick_label.set_color(color_ink)
        tick_label.set_fontweight("bold" if pop == "KOREAN" else "normal")
    ax.set_ylabel("CpG-Eliminating sites on array")
    ax.set_ylim(0, array_bottom.max() * 1.16)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(format_count))
    style_axis(ax)


def draw_cumulative_curves(ax, chip):
    """Right half of panels B-E: cumulative CpG loss as the AF cut-off is lowered
    from 1 towards 0. y is the number of CpG sites whose most frequent disrupting
    variant has AF >= x, i.e. the common sites accumulate first, and the x-axis is
    drawn high-to-low so the curve still rises left to right. The right end
    (all frequencies included) matches the bar total on the left."""
    table_chip = table_cumul[table_cumul["Chip"] == chip]
    n_site_chip = dict_chip_to_n_site[chip]

    # Guides at the two frequency boundaries of the stacked bars on the left.
    for af_cut in [0.01, 0.05]:
        ax.axvline(af_cut, color=color_guide, linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)

    list_end = list()
    y_max = 0.0
    for pop in list_pops:
        table_pop = table_chip[table_chip["Superpopulation"] == pop].sort_values("AF_Cutoff")
        array_x = table_pop["AF_Cutoff"].to_numpy()
        # The stored column counts sites with AF <= cut-off; the complement against
        # the all-frequency total gives the common-first count with AF > cut-off.
        array_cumul_le = table_pop["N_CpG_Site_Cumulative"].to_numpy().astype(float)
        array_y = array_cumul_le[-1] - array_cumul_le
        is_korea = pop == "KOREAN"
        ax.plot(array_x, array_y,
                color=dict_pop_to_color[pop],
                linewidth=2.6 if is_korea else 2.0,
                solid_capstyle="round",
                zorder=5 if is_korea else 4)
        # Curves now end at the LOW-frequency end of the axis, which the inverted
        # x-limits put on the right, so the label anchor is the first grid point.
        list_end.append((dict_pop_to_label[pop], array_y[0], dict_pop_to_color[pop]))
        y_max = max(y_max, array_y[0])

    ax.set_xscale("log")
    ax.set_xlim(1.0, 3e-5)          # high frequency on the left, rare on the right
    ax.set_ylim(0, y_max * 1.10)
    ax.set_xlabel("Allele frequency cut-off (high $\\rightarrow$ low)")
    ax.set_ylabel("Cumulative CpG sites lost")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(format_count))
    style_axis(ax, grid_axis='both')

    for af_cut, name in zip([0.01, 0.05], ["0.01", "0.05"]):
        ax.text(af_cut, y_max * 1.12, name, ha="center", va="bottom",
                fontsize=8.5, color=color_ink)

    place_end_labels(ax, list_end)
    ax.text(0.02, 0.97, f"{n_site_chip:,} CpG sites on array",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=8.5, color=color_ink)


def draw_ageclock_panel(ax):
    """Panel E: share of each aging clock's CpG markers lost, per cohort."""
    n_pop = len(list_pops)
    width_bar = 0.8 / n_pop
    array_group = np.arange(len(list_clocks))

    # Bars run horizontally so the panel stands as a tall column beside A-D, and
    # so the cohort names read upright rather than rotated.
    list_tick_pos = list()
    list_tick_label = list()
    for index_pop, pop in enumerate(list_pops):
        array_offset = array_group - 0.4 + width_bar * (index_pop + 0.5)
        table_pop = table_clock[table_clock["Superpopulation"] == pop].set_index("AgeClock")
        array_common = np.array([table_pop.loc[clock, "Pct_Common"] for clock in list_clocks])
        array_rare = np.array([table_pop.loc[clock, "Pct_Rare"] for clock in list_clocks])
        # Frequency is carried by colour; the cohort is named beside each bar.
        ax.barh(array_offset, array_common, height=width_bar * 0.82,
                color=color_af_common, edgecolor=color_surface, linewidth=0.6, zorder=3)
        ax.barh(array_offset, array_rare, left=array_common, height=width_bar * 0.82,
                color=color_af_rare, edgecolor=color_surface, linewidth=0.6, zorder=3)
        for y_value, total in zip(array_offset, array_common + array_rare):
            ax.text(total + 0.45, y_value, f"{total:.1f}", ha="left", va="center",
                    fontsize=9.5, color=color_ink)
        list_tick_pos.extend(array_offset.tolist())
        list_tick_label.extend([dict_pop_to_label[pop]] * len(list_clocks))

    # Only the cohort names sit on the axis; the clock name is placed inside the
    # plot area, right-aligned in the empty space beyond the bars, so the panel
    # needs no wide outer margin.
    ax.set_yticks(list_tick_pos, minor=True)
    ax.set_yticklabels(list_tick_label, minor=True, fontsize=9.5, color=color_ink)
    ax.tick_params(axis='y', which='minor', length=2, width=0.7, pad=2)
    ax.set_yticks([])

    ax.set_xlabel("Clock CpG markers lost (%)")
    # First clock at the top. The strip above the first clock title (-0.80 to
    # -0.5) holds the colour key, clear of the Horvath title.
    ax.set_ylim(len(list_clocks) - 0.45, -0.80)
    ax.set_xlim(0, 34)
    ax.set_xticks([0, 5, 10, 15, 20, 25, 30])
    style_axis(ax, grid_axis='x')

    # The clock name sits inside the plot area, in the blank strip just above each
    # group of bars, so the panel needs no outer margin for it at any width.
    legend_af = ax.legend(
        handles=[Patch(facecolor=color_af_common, edgecolor=color_surface,
                       label="Common variant (AF > 0.01)"),
                 Patch(facecolor=color_af_rare, edgecolor=color_surface,
                       label="Rare variant (AF $\\leq$ 0.01)")],
        loc="upper right", bbox_to_anchor=(1.0, 1.0),
        frameon=False, handlelength=1.8, handleheight=1.1, labelspacing=0.4,
        borderaxespad=0.1, borderpad=0.2, prop={'size': size_key_text_panel},
    )
    ax.add_artist(legend_af)

    for index_clock, clock in enumerate(list_clocks):
        ax.text(0.3, index_clock - 0.43,
                f"{clock}  ({dict_clock_to_n_marker[clock]:,} CpG markers)",
                ha="left", va="bottom",
                fontsize=9.5, fontweight="bold", color=color_ink, zorder=6)


def draw_genotype_beta_panel(ax):
    """Panel A: per-CpG mean beta against the dosage of the disrupting allele."""
    list_data = [array_beta_by_gt[gt] for gt in (0, 1, 2)]
    array_x = np.arange(3) + 1

    parts = ax.violinplot(list_data, positions=array_x, widths=0.80,
                          showextrema=False, showmedians=False)
    for index_gt, body in enumerate(parts["bodies"]):
        body.set_facecolor(list_color_genotype[index_gt])
        body.set_edgecolor(color_surface)
        body.set_linewidth(0.8)
        body.set_alpha(1.0)
        body.set_zorder(3)

    boxes = ax.boxplot(list_data, positions=array_x, widths=0.14,
                       showfliers=False, patch_artist=True, zorder=4,
                       medianprops=dict(color=color_ink, linewidth=1.4),
                       boxprops=dict(facecolor=color_surface, edgecolor=color_ink, linewidth=0.9),
                       whiskerprops=dict(color=color_ink, linewidth=0.9),
                       capprops=dict(color=color_ink, linewidth=0.9))
    del boxes

    for x_value, array_beta in zip(array_x, list_data):
        ax.text(x_value, 1.035, f"Mean = {array_beta.mean():.3f}",
                ha="center", va="bottom", fontsize=9.5, color=color_ink)

    ax.set_xticks(array_x)
    ax.set_xticklabels(["0", "1", "2"])
    ax.set_xlim(0.45, 3.55)
    ax.set_ylim(0, 1.10)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Genotype of CpG-Eliminating variant")
    ax.set_ylabel("Mean beta value per CpG site")
    ax.text(0.985, 0.72,
            f"{len(list_data[0]):,} CpG sites\n($\\geq$ {n_min_genotype} samples in every genotype)",
            transform=ax.transAxes, ha="right", va="top", multialignment="right",
            fontsize=9, color=color_ink)
    style_axis(ax)


#%%
# --- Compose ---------------------------------------------------------------
# 1.6 in taller than the original 18.3 in, all of it spent on the band that
# holds the enlarged colour keys under panels B-E; the panels keep their size.
fig = plt.figure(figsize=(12.0, 19.9))
# Two grids so the band between the A-D block and panel E can hold the shared
# legends at a fixed height, independent of the subplot spacing.
# A-D occupy a narrower left column; E stands as a tall panel to their right.
# Column widths split 0.4 (A-D) / 0.6 (E) of the plotting area; inside A-D the
# cumulative curve gets 65% of the bar chart's width.
# Row 0 of the left column is the genotype panel, spanning both of its cells;
# rows 1-4 are the four arrays.
grid_chip = fig.add_gridspec(
    5, 2, width_ratios=[1.0, 0.78], height_ratios=[0.82, 1.0, 1.0, 1.0, 1.0],
    left=0.085, right=0.535, top=0.960, bottom=0.160,
    hspace=0.50, wspace=0.34,
)
ax_geno = fig.add_subplot(grid_chip[0, :])
draw_genotype_beta_panel(ax_geno)

box_geno = ax_geno.get_position()

dict_panel_to_axes = dict()
for index_chip, chip in enumerate(list_chips):
    ax_bar = fig.add_subplot(grid_chip[index_chip + 1, 0])
    ax_cum = fig.add_subplot(grid_chip[index_chip + 1, 1])
    draw_af_class_bars(ax_bar, chip)
    draw_cumulative_curves(ax_cum, chip)
    dict_panel_to_axes[chip] = (ax_bar, ax_cum)

# Panel F spans exactly the left column: its top edge lines up with panel A's and
# its bottom with panel E's, taken from the drawn axes rather than hard-coded.
box_msa = dict_panel_to_axes["MSA"][0].get_position()
grid_clock = fig.add_gridspec(
    1, 1,
    left=0.685, right=0.985,
    top=box_geno.y1, bottom=box_msa.y0,
)
ax_clock = fig.add_subplot(grid_clock[0, 0])
draw_ageclock_panel(ax_clock)

fig.text(box_geno.x0 - 0.056, box_geno.y1 + 0.012, "A",
         fontsize=15, fontweight="bold", ha="left", va="bottom", color=color_ink)
fig.text(box_geno.x0 + 0.5 * box_geno.width, box_geno.y1 + 0.012,
         "Methylation loss with CpG-Eliminating allele dosage",
         fontsize=11.5, fontweight="bold", ha="center", va="bottom", color=color_ink)

# Panel letter + array name, centred over each pair of subplots.
for index_chip, chip in enumerate(list_chips):
    ax_bar, ax_cum = dict_panel_to_axes[chip]
    box_left = ax_bar.get_position()
    box_right = ax_cum.get_position()
    y_title = box_left.y1 + 0.012
    fig.text(box_left.x0 - 0.056, y_title, "BCDE"[index_chip],
             fontsize=15, fontweight="bold", ha="left", va="bottom", color=color_ink)
    fig.text((box_left.x0 + box_right.x1) / 2, y_title, dict_chip_to_title[chip],
             fontsize=11.5, fontweight="bold", ha="center", va="bottom", color=color_ink)

box_clock = ax_clock.get_position()
fig.text(box_clock.x0 - 0.056, box_clock.y1 + 0.022, "F",
         fontsize=15, fontweight="bold", ha="left", va="bottom", color=color_ink)
fig.text(box_clock.x0 + 0.5 * box_clock.width, box_clock.y1 + 0.016,
         "Epigenetic-clock CpG markers lost",
         fontsize=11.5, fontweight="bold", ha="center", va="bottom", color=color_ink)

# Shared legends, stated once in the band under the B-E block: the AF classes
# key the stacked bars, the cohort hues key the cumulative curves. At this key
# size the two do not fit side by side, so they are stacked and centred on the
# figure. Their heights come from the rendered text, not from guesses: the band
# starts just below the lowest tick label / axis label of the panels above.
def figure_box(artist, renderer):
    return artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())


def style_key_title(legend):
    legend.get_title().set_fontsize(size_key_title)
    legend.get_title().set_fontweight("bold")
    legend.get_title().set_color(color_ink)


renderer = fig.canvas.get_renderer()
list_bottom_axes = list(dict_panel_to_axes["MSA"]) + [ax_clock]
y_key_top = min(ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted()).y0
                for ax in list_bottom_axes) - 0.008

legend_class = fig.legend(
    handles=[Patch(facecolor=dict_af_class_to_color[name], edgecolor=color_surface,
                   label=dict_af_class_to_label[name])
             for name in list_af_class],
    title="Allele frequency",
    loc="upper center", bbox_to_anchor=(0.5, y_key_top),
    frameon=False, ncol=2, handlelength=2.2, handleheight=1.3, columnspacing=2.0,
    labelspacing=0.45, borderpad=0.2, prop={'size': size_key_text},
    alignment="center",
)
style_key_title(legend_class)

legend_pop = fig.legend(
    handles=[plt.Line2D([0], [0], color=dict_pop_to_color[pop], linewidth=5.0,
                        solid_capstyle="butt", label=dict_pop_to_label[pop])
             for pop in list_pops],
    title="Population",
    loc="upper center", bbox_to_anchor=(0.5, figure_box(legend_class, renderer).y0 - 0.006),
    frameon=False, ncol=len(list_pops), handlelength=2.6, columnspacing=2.0,
    borderpad=0.2, prop={'size': size_key_text},
    alignment="center",
)
style_key_title(legend_pop)

box_pop = figure_box(legend_pop, renderer)
if box_pop.y0 < 0.0 or box_pop.x0 < 0.0 or box_pop.x1 > 1.0:
    print(f"WARNING: colour keys overflow the figure ({box_pop}); enlarge figsize or shrink size_key_text")

fig.savefig(path_out_png, dpi=300)
fig.savefig(path_out_pdf)
print(f"Wrote {path_out_png}")
print(f"Wrote {path_out_pdf}")

#%%
