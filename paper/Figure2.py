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
from korea10k.config import PROJECT_DIR, WORK_DIR


# -------------------- Helper --------------------
def load_pickle(path_pkl):
    with open(path_pkl, "rb") as f:
        return pickle.load(f)

# -------------------- Plotting settings --------------------
plt.rcParams.update({
    "font.size": 18,
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

# -------------------- Figure layout --------------------
fig = plt.figure(figsize=(15, 13))
gs = gridspec.GridSpec(2, 2, 
                       figure=fig,
                       height_ratios=[2.5, 2.5],
                       width_ratios=[1, 1],
                       hspace=0.35, 
                       wspace=0.25)

axA = fig.add_subplot(gs[0:1, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])
for ax in [axA, axB, axC, axD]:
    ax.xaxis.label.set_weight("bold")
    ax.yaxis.label.set_weight("bold")

def widen_axis(ax, factor=1.15, anchor="center"):
    pos = ax.get_position()
    old_width = pos.width
    new_width = old_width * factor

    if anchor == "left":
        new_x0 = pos.x0
    elif anchor == "right":
        new_x0 = pos.x0 - (new_width - old_width)
    elif anchor == "center":
        new_x0 = pos.x0 - (new_width - old_width) / 2
    else:
        raise ValueError("anchor must be 'left', 'center', or 'right'")

    ax.set_position([new_x0, pos.y0, new_width, pos.height])


def heighten_axis(ax, factor=1.15, anchor="bottom"):
    pos = ax.get_position()
    old_height = pos.height
    new_height = old_height * factor

    if anchor == "bottom":
        new_y0 = pos.y0
    elif anchor == "top":
        new_y0 = pos.y0 - (new_height - old_height)
    elif anchor == "center":
        new_y0 = pos.y0 - (new_height - old_height) / 2
    else:
        raise ValueError("anchor must be 'bottom', 'center', or 'top'")

    ax.set_position([pos.x0, new_y0, pos.width, new_height])
    
sns.despine(left=False, bottom=False)

# -------------------- Panel A: Variant count --------------------
dir_variant = f"{WORK_DIR}/genome/variant/Data"
count_freq_reported = load_pickle(os.path.join(dir_variant, "count_freq_reported.pkl"))
count_freq_novel = load_pickle(os.path.join(dir_variant, "count_freq_novel.pkl"))
count_freq_all = load_pickle(os.path.join(dir_variant, "count_freq_all.pkl"))

list_xticks = ["Singleton", "Doubleton", "Extremely Rare", "Very Rare", "Rare", "Common", "Very Common"]
bottom = np.zeros(len(list_xticks))

values_novel = np.array([count_freq_novel[x] for x in list_xticks])
values_reported = np.array([count_freq_reported[x] for x in list_xticks])

axA.bar(list_xticks, values_novel, bottom=bottom, color="firebrick", label="Novel", edgecolor='k', linewidth=1.5, width=0.9, zorder=2)
bottom += values_novel
axA.bar(list_xticks, values_reported, bottom=bottom, color="gray", label="dbSNP", edgecolor='k', linewidth=1.5, width=0.9, zorder=2)

for i, freq in enumerate(list_xticks):
    total = count_freq_all[freq]
    axA.annotate(f"{values_novel[i]/total*100:.2f}", (i, (values_novel[i]/2) if values_novel[i]/total > 0.04 else 0.05e7), color='white', ha="center", va="center", weight="bold", fontsize=14)
    axA.annotate(f"{values_reported[i]/total*100:.2f}", (i, values_novel[i]+ values_reported[i]/2 + (0 if freq != "Common" else 0.05e7)), color='white', ha="center", va="center", weight="bold", fontsize=12)
    axA.annotate(f"{total/1e6:.1f}M", (i, total + 0.01*max(count_freq_all.values())), ha="center", va="bottom", weight="bold", fontsize=14)

axA.set_xticklabels([display_label_map[x] for x in list_xticks],
                    rotation=45, rotation_mode="anchor", ha="right", weight="bold")
axA.set_ylabel("Variant Count", fontsize=plt.rcParams["font.size"], weight="bold")
axA.set_ylim(top=2.7e7)
axA.grid(axis="y", linestyle="--", linewidth=0.8, alpha=0.5)

axA.legend(frameon=False, loc="center right", fontsize=plt.rcParams["font.size"]-2, bbox_to_anchor=(1, 0.5))

total_novel = sum(count_freq_novel.values())
total_reported = sum(count_freq_reported.values())
ax_pie = inset_axes(axA, width="50%", height="50%", loc='upper center')
ax_pie.pie([total_novel, total_reported], colors=["firebrick", "gray"], autopct="%1.1f%%",
           startangle=90, counterclock=False,
           wedgeprops={"linewidth":2, "edgecolor":"k"},
           textprops={"weight":"bold", "fontsize": plt.rcParams["font.size"]-2, "color":"white"})
ax_pie.set_aspect("equal")
ax_pie.text(0.5, 0.92, f"{(total_novel+total_reported)/1e6:.1f}M", ha="center", va="bottom", weight="bold", fontsize=plt.rcParams["font.size"], transform=ax_pie.transAxes)
axA.text(-0.18, 1.0, "A", transform=axA.transAxes, fontsize=plt.rcParams["font.size"]+5, fontweight='bold')
axA.spines['bottom'].set_zorder(10)
for tick in axA.xaxis.get_major_ticks():
    tick.label1.set_zorder(10)
    
# -------------------- Panel B --------------------
dir_variant_plot = f"{WORK_DIR}/genome/shapeit/Resources/Data/plot"
dict_frequency_type_to_mean_line = load_pickle(os.path.join(dir_variant_plot, "dict_frequency_type_to_mean_line.pkl"))
dict_frequency_type_to_saturation_points = load_pickle(os.path.join(dir_variant_plot, "dict_frequency_type_to_saturation_points.pkl"))

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

list_frequency_type = ["Total", "Singleton", "Doubleton", "Extremely Rare", "Very Rare", "Rare", "Common", "Very Common"]

xlim_padd = 500

highlight_cats = ["Extremely Rare"]

max_samples = 0
for frequency_cat in list_frequency_type:
    mean_line = dict_frequency_type_to_mean_line[frequency_cat]["mean"]
    list_samples = list(range(1, len(mean_line)+1))
    
    # Ultra-thin for total curve
    if frequency_cat.lower() in ["total", "all", "all variants"]:
        lw = 2
        zord = 2
    # Highlighted categories
    elif frequency_cat in highlight_cats:
        lw = 8
        zord = 7
    # Normal categories
    else:
        lw = 4
        zord = 5

    axB.plot(list_samples, mean_line, color=dict_frequency_palette[frequency_cat],
         linewidth=lw, label=display_label_map.get(frequency_cat, frequency_cat), alpha=0.7)
    if max(list_samples) > max_samples:
        max_samples = max(list_samples)
axB.set_xlabel("Sample Count", fontsize=plt.rcParams["font.size"], weight="bold")
axB.set_yscale("log")
axB.set_ylabel("Variant Count\n(log-scale)", fontsize=plt.rcParams["font.size"], weight="bold")
axB.grid(axis="both", linestyle="--", linewidth=0.5, alpha=0.5)
axB.legend(frameon=True, fontsize=16, loc="lower right")
axB.set_xlim(right=10001)
axB.text(-0.18, 1.0,"B", transform=axB.transAxes, fontsize=plt.rcParams["font.size"]+5, fontweight='bold')


# -------------------- Panel D --------------------
workdir = f"{WORK_DIR}/genome/shapeit/Results/imputation"
R2_template = os.path.join(workdir, "corr_results/chr{i}.dose.{panel}.concat.corr.maf_10K_cat_added.txt.gz")
def calc_squared_with_minus(val): return (1 if val >= 0 else -1) * val**2

df1 = pd.read_csv(R2_template.format(i=2, panel="1K"), sep="\t")
df2 = pd.read_csv(R2_template.format(i=2, panel="4K"), sep="\t")
df3 = pd.read_csv(R2_template.format(i=2, panel="10K"), sep="\t")
for df in [df1, df2, df3]: df["Rsq"] = df["R"].apply(calc_squared_with_minus)

bin_edges = [0,0.0005,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0]
labels = ["0-0.05","0.05-0.2","0.2-0.5","0.5-1","1-2","2-5","5-10","10-20","20-50","50-100"]
for df in [df1, df2, df3]: df["AF_bin"]=pd.cut(df["ALT_FREQ"], bins=bin_edges, labels=labels)

agg1 = df1.groupby("AF_bin", observed=True)["Rsq"].mean().reset_index()
agg2 = df2.groupby("AF_bin", observed=True)["Rsq"].mean().reset_index()
agg3 = df3.groupby("AF_bin", observed=True)["Rsq"].mean().reset_index()

sns.lineplot(data=agg3, x="AF_bin", y="Rsq", color="crimson", marker="o", markeredgecolor="k", markersize=8, linewidth=4, ax=axC, label="10K", zorder=5)

sns.lineplot(data=agg2, x="AF_bin", y="Rsq", color="limegreen", marker="o", markeredgecolor="k", markersize=8, linewidth=4, ax=axC, label="4K", zorder=4)

sns.lineplot(data=agg1, x="AF_bin", y="Rsq", color="dodgerblue", marker="o", markeredgecolor="k", markersize=8, linewidth=4, ax=axC, label="1K", zorder=3)

axC.set_xlabel("Alt allele frequency (%)", fontsize=plt.rcParams["font.size"], weight="bold")
axC.set_ylabel("Aggregated R²", fontsize=plt.rcParams["font.size"], weight="bold")
axC.set_ylim(0.1, 1.19)
axC.set_xticklabels(labels, rotation=45, ha="right")
axC.legend(frameon=False, loc="upper right", bbox_to_anchor=(1, 0.6))
axC.grid(axis="both", linestyle="--", linewidth=0.5, alpha=0.5)
axC.text(-0.18, 1.0,"C", transform=axC.transAxes, fontsize=plt.rcParams["font.size"]+5, fontweight='bold')

x_positions = range(len(agg1))
cat_ranges = {
    "Low\nfrequency": (0, 3),
    "Common": (3, 6),
    "Very\ncommon": (6, len(labels)-1)
}

for label, (start, end) in cat_ranges.items():
    axC.axvspan(start, end, ymin=0, ymax=0.13, color="white", ec="black", lw=1.0, clip_on=False, zorder=3)

for label, (start, end) in cat_ranges.items():
    mid = (start + end) / 2
    axC.text(mid, 0.065, label, transform=axC.get_xaxis_transform(),
            ha="center", va="center", weight="bold", fontsize=plt.rcParams["font.size"]-3, color="black")

num_variants = [df1.shape[0], df2.shape[0], df3.shape[0]]
num_variants_m = np.array(num_variants) / 1e6

ax_inset = inset_axes(
    axC,
    width="45%",
    height="17%",
    loc="upper left",
    borderpad=2.5
)

bars = ax_inset.barh(
    ["1K", "4K", "10K"],
    num_variants,
    color=["dodgerblue", "limegreen", "crimson"],
    edgecolor="k",
    linewidth=1.2,
    height=1
)

for i, bar in enumerate(bars):
    ax_inset.text(
        bar.get_width() / 2,
        bar.get_y() + bar.get_height() / 2,
        f"{num_variants[i]:,}",
        ha="center", 
        va="center",
        color="white", 
        weight="bold", 
        fontsize=15
    )

ax_inset.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x/1e6:.1f}"))
ax_inset.set_xlim(0, max(num_variants) * 1.05)
ax_inset.xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))

ax_inset.xaxis.tick_top()
ax_inset.xaxis.set_label_position('top')

ax_inset.set_xlabel("Number of imputed variants", fontsize=17, labelpad=13)

ax_inset.text(
    1.15, 1.0, "1e6",
    transform=ax_inset.transAxes,
    ha="right", va="bottom",
    fontsize=13
)

ax_inset.tick_params(axis="y", labelsize=15)
ax_inset.tick_params(axis="x", labelsize=14)


ratio_10K_4K = (agg3["Rsq"] - agg2["Rsq"])/agg2["Rsq"]*100
ratio_10K_1K = (agg3["Rsq"] - agg1["Rsq"])/agg1["Rsq"]*100
x_positions = range(len(agg1))
axD.plot(x_positions, ratio_10K_4K, color="orange", marker="o", markeredgecolor="k", markersize=8, label="10K / 4K")
axD.plot(x_positions, ratio_10K_1K, color="purple", marker="o", markeredgecolor="k", markersize=8, label="10K / 1K")
axD.set_xticks(x_positions)
axD.set_xticklabels(labels, rotation=45, ha="right")
axD.set_ylabel("R² improvement (%)", weight="bold")
axD.set_ylim(-30, 101)
axD.grid(axis="both", linestyle="--", linewidth=0.5, alpha=0.5)
axD.set_xlabel("Alt allele frequency (%)", fontsize=plt.rcParams["font.size"], weight="bold")
axD.legend(frameon=False, fontsize=18, loc="center right", bbox_to_anchor = (1, 0.6))
axD.text(-0.18,1.0,"D", transform=axD.transAxes, fontsize=plt.rcParams["font.size"]+5, fontweight='bold')

x_positions = range(len(agg1))
cat_ranges = {
    "Low\nfrequency": (0, 3),
    "Common": (3, 6),
    "Very\ncommon": (6, len(labels)-1)
}

for label, (start, end) in cat_ranges.items():
    axD.axvspan(start, end, ymin=0, ymax=0.15, color="white", ec="black", lw=1.0, clip_on=False, zorder=3)

# Add text annotations below x-axis
for label, (start, end) in cat_ranges.items():
    mid = (start + end) / 2
    axD.text(mid, 0.075, label, transform=axD.get_xaxis_transform(),
            ha="center", va="center", weight="bold", fontsize=plt.rcParams["font.size"]-3, color="black")

plt.tight_layout()

dir_figure = f"{PROJECT_DIR}/Analysis/Revision/Draw_Figure/Figures"
os.makedirs(dir_figure, exist_ok=True)
fig.savefig(os.path.join(dir_figure, "Figure2.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(dir_figure, "Figure2.pdf"), bbox_inches="tight")

plt.show()

# %%
