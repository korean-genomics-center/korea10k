#%%
# ---------------------------------------------------------------------------
# Figure: how far the methylation beta value moves away from the homozygous
#         reference genotype at CpG-eliminating variants, split by the allele
#         frequency of the variant.
#
#   Row 1 : delta1 = mean_gt1 - mean_gt0   (heterozygous vs hom. reference)
#   Row 2 : delta2 = mean_gt2 - mean_gt0   (hom. alternate vs hom. reference)
#   Cols  : four allele-frequency groups
#
# Input : Data/Methylation_Change/per_cpg_beta_by_genotype.tsv, one row per CpG
#         position with the per-genotype sample count and mean beta.
#
# The table carries no allele frequency, so it is derived from the genotype
# counts of the same samples:  AF = (n_gt1 + 2 n_gt2) / (2 (n_gt0+n_gt1+n_gt2)).
# That keeps the frequency and the beta values on exactly the same cohort.
#
# A position enters a row only when both genotypes of that contrast have at
# least n_min_sample carriers, so every delta is a difference of two reasonably
# estimated means rather than of one or two individuals.
# ---------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt

path_genotype_beta = ("/BiO/Research/Korea10KGenome/Analysis/Revision/Draw_Figure/Data/"
                      "Methylation_Change/per_cpg_beta_by_genotype.tsv")
dir_figure_out = "/BiO/Research/Korea10KGenome/Analysis/Revision/Draw_Figure/Figures"
os.makedirs(dir_figure_out, exist_ok=True)

path_out_png = os.path.join(dir_figure_out, "Figure.Beta_Change_by_AF.png")
path_out_pdf = os.path.join(dir_figure_out, "Figure.Beta_Change_by_AF.pdf")
path_out_summary = os.path.join(dir_figure_out, "Figure.Beta_Change_by_AF.Summary.tsv")

n_min_sample = 10

# Allele-frequency groups. Upper bounds only; the first group is everything below
# the first bound.
list_af_group = [
    ("AF $\\leq$ 0.01", 0.0, 0.01),
    ("0.01 < AF $\\leq$ 0.05", 0.01, 0.05),
    ("0.05 < AF $\\leq$ 0.20", 0.05, 0.20),
    ("AF > 0.20", 0.20, 1.01),
]
# Blue ordinal ramp of the reference palette, light to dark with rising
# frequency (steps 250 / 350 / 500 / 650).
list_color_af = ["#86b6ef", "#5598e7", "#256abf", "#104281"]

list_contrast = [
    ("delta1", "mean_gt1", "n_gt1", "GT = 1"),
    ("delta2", "mean_gt2", "n_gt2", "GT = 2"),
]

color_surface = "#ffffff"
color_ink = "#0b0b0b"
color_grid = "#e1e0d9"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10.5,
    "axes.labelweight": "bold",
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 9,
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
# --- Load and derive ---------------------------------------------------------
table = pd.read_csv(path_genotype_beta, sep='\t', na_values=["NA"])

array_n_total = (table["n_gt0"] + table["n_gt1"] + table["n_gt2"]).to_numpy()
array_ac = (table["n_gt1"] + 2 * table["n_gt2"]).to_numpy()
table["AF"] = array_ac / (2 * array_n_total)
table["delta1"] = table["mean_gt1"] - table["mean_gt0"]
table["delta2"] = table["mean_gt2"] - table["mean_gt0"]

print(f"{len(table):,} CpG positions, median cohort size {int(np.median(array_n_total)):,}")


def select_contrast(column_n):
    """Positions where both genotypes of the contrast are well sampled."""
    return table[(table["n_gt0"] >= n_min_sample) & (table[column_n] >= n_min_sample)]


#%%
# --- Summary table -----------------------------------------------------------
list_rows = list()
for key_delta, _column_mean, column_n, _label in list_contrast:
    table_use = select_contrast(column_n)
    for name_group, af_low, af_high in list_af_group:
        flag = (table_use["AF"] > af_low) & (table_use["AF"] <= af_high)
        array_delta = table_use.loc[flag, key_delta].to_numpy()
        list_rows.append({
            "Contrast": key_delta,
            "AF_Group": name_group.replace("$\\leq$", "<=").replace("$-$", "-"),
            "N_CpG": len(array_delta),
            "Mean_Delta": round(float(array_delta.mean()), 4) if len(array_delta) else np.nan,
            "Median_Delta": round(float(np.median(array_delta)), 4) if len(array_delta) else np.nan,
            "Pct_Below_Minus_0.2": (round(100.0 * float((array_delta <= -0.2).mean()), 2)
                                    if len(array_delta) else np.nan),
        })

table_summary = pd.DataFrame(list_rows)
table_summary.to_csv(path_out_summary, sep='\t', index=False)
print(table_summary.to_string(index=False))
print(f"\nWrote {path_out_summary}")


#%%
# --- Draw --------------------------------------------------------------------
array_bin = np.linspace(-1.0, 0.5, 76)

fig, array_ax = plt.subplots(
    len(list_contrast), len(list_af_group),
    figsize=(15.0, 7.0), sharex=True, sharey='row',
)

for index_row, (key_delta, _column_mean, column_n, label_row) in enumerate(list_contrast):
    table_use = select_contrast(column_n)
    for index_col, (name_group, af_low, af_high) in enumerate(list_af_group):
        ax = array_ax[index_row, index_col]
        flag = (table_use["AF"] > af_low) & (table_use["AF"] <= af_high)
        array_delta = table_use.loc[flag, key_delta].to_numpy()

        ax.axvline(0.0, color=color_ink, linewidth=0.9, zorder=2)
        if len(array_delta):
            ax.hist(array_delta, bins=array_bin, color=list_color_af[index_col],
                    edgecolor=color_surface, linewidth=0.3, zorder=3)
            median_delta = float(np.median(array_delta))
            ax.axvline(median_delta, color=color_ink, linewidth=1.4,
                       linestyle=(0, (4, 2)), zorder=4)
            ax.text(0.97, 0.95,
                    f"n = {len(array_delta):,}\nmedian = {median_delta:.3f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=9,
                    multialignment="right", color=color_ink)
        else:
            # Homozygous alternates only exist once the variant is common, so the
            # low-frequency cells of the second row are empty by construction.
            ax.text(0.5, 0.5,
                    f"no position with\n$\\geq$ {n_min_sample} carriers",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9.5, color=color_ink)

        ax.set_xlim(-1.0, 0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis='y', color=color_grid, linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(length=3, width=0.7)

        if index_row == 0:
            ax.set_title(name_group, fontsize=11, fontweight="bold", pad=8)
        if index_row == len(list_contrast) - 1:
            ax.set_xlabel("Change in mean beta")
        if index_col == 0:
            ax.set_ylabel("CpG positions")

for index_row in range(len(list_contrast)):
    y_max = max(ax.get_ylim()[1] for ax in array_ax[index_row])
    for ax in array_ax[index_row]:
        ax.set_ylim(0, y_max * 1.10)

# Row labels on the left, outside the axes.
for index_row, (_key, _mean, _n, label_row) in enumerate(list_contrast):
    box = array_ax[index_row, 0].get_position()
    fig.text(0.028, box.y0 + 0.5 * box.height, label_row,
             rotation=90, ha="center", va="center",
             fontsize=11.5, fontweight="bold", color=color_ink)
    fig.text(0.006, box.y1 + 0.010, "AB"[index_row],
             fontsize=15, fontweight="bold", ha="left", va="bottom", color=color_ink)

fig.subplots_adjust(left=0.078, right=0.985, top=0.905, bottom=0.105, wspace=0.16, hspace=0.22)
fig.savefig(path_out_png, dpi=300)
fig.savefig(path_out_pdf)
print(f"Wrote {path_out_png}")
print(f"Wrote {path_out_pdf}")

# %%
