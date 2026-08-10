# %% 
import os
import pickle

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# -------------------- Helper --------------------
def load_pickle(path_pkl):
    with open(path_pkl, "rb") as f:
        return pickle.load(f)

# -------------------- Plotting settings --------------------
plt.rcParams.update({
    "font.size": 15,
    "axes.linewidth": 1.5,
    "axes.labelpad": 6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "figure.dpi": 300
})

display_label_map = {
    "Singleton": "Singleton",
    "Doubleton": "Doubleton",
    "Extremely Rare": "Ultra-Rare",
    "Very Rare": "Very Rare",
    "Rare": "Rare",
    "Common": "Common",
    "Very Common": "Very Common",
    "Total": "Total"
}

# Figure
fig, ax = plt.subplots(figsize=(5, 5))

# Load data
dir_variant_plot = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/plot"
dict_frequency_type_to_mean_line = load_pickle(
    os.path.join(dir_variant_plot, "dict_frequency_type_to_mean_line.pkl")
)
dict_frequency_type_to_saturation_points = load_pickle(
    os.path.join(dir_variant_plot, "dict_frequency_type_to_saturation_points.pkl")
)

# Color palette
dict_frequency_palette = {
    "Total": "#8a8a8aff",
    "Singleton": "#000000ff",
    "Doubleton": "#78004A",
    "Extremely Rare": "#cb181d",
    "Very Rare": "#ff7b00",
    "Rare": "#0B8100",
    "Common": "#0026ff",
    "Very Common": "#3900e3"
}

list_frequency_type = [
    "Total", "Singleton", "Doubleton", "Extremely Rare", 
    "Very Rare", "Rare", "Common", "Very Common"
]

# -------------------- Vertical Boxplot --------------------
# Reverse order for plotting
freq_list_rev = list_frequency_type[::-1]

# Prepare data
data_to_plot = [
    dict_frequency_type_to_saturation_points[freq]
    for freq in freq_list_rev
]

colors = [dict_frequency_palette[freq] for freq in freq_list_rev]

# Create boxplot
box = ax.boxplot(
    data_to_plot,
    vert=True,            # vertical orientation (y-axis = sample count)
    patch_artist=True,
    widths=0.8
)

# Color styling
for patch, col in zip(box["boxes"], colors):
    patch.set_facecolor(col)
    patch.set_edgecolor("black")
    patch.set_linewidth(1.2)

for element in ["whiskers", "caps", "medians"]:
    for item in box[element]:
        item.set_color("black")
        item.set_linewidth(1.2)

for flier, col in zip(box["fliers"], colors):
    flier.set_marker("o")
    flier.set_markerfacecolor(col)
    flier.set_markeredgecolor("white")
    flier.set_markersize(6)

# -------------------- Axes & Labels --------------------
ax.set_xticks(range(1, len(freq_list_rev) + 1))
ax.set_yticks(range(0, 10001, 1000))
ax.set_xticklabels(
    [display_label_map[freq] for freq in freq_list_rev],
    rotation=45,
    ha="right",
    fontsize=plt.rcParams["font.size"]
)

ax.set_ylabel("Sample Count", fontsize=plt.rcParams["font.size"], weight="bold")
ax.set_xlabel("Frequency Category", fontsize=plt.rcParams["font.size"], weight="bold")

ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.8)

from matplotlib.patches import Patch

# Create legend handles
legend_handles = [
    Patch(facecolor=dict_frequency_palette[freq], edgecolor="black", label=display_label_map[freq])
    for freq in freq_list_rev
]

# Add legend
ax.legend(
    handles=legend_handles,
    title="Frequency Class",
    fontsize=plt.rcParams["font.size"]-5,
    title_fontsize=plt.rcParams["font.size"]-4,
    frameon=True,
    loc="lower right",
    bbox_to_anchor=(1, 0.1)
)

sns.despine(ax=ax, top=True, right=True)

plt.tight_layout()
plt.show()
