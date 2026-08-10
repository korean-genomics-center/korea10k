# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# %%
plt.rcParams["font.size"] = 16
plt.rcParams["axes.linewidth"] = 1.2

workdir = "/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation"
R2_maf_cat_template = os.path.join(
    workdir, "corr_results/chr{i}.dose.{panel}.concat.corr.maf_10K_cat_added.txt.gz"
)

def calc_squared_with_minus(val):
    return (1 if val >= 0 else -1 ) * val ** 2

# Read files
df1 = pd.read_csv(R2_maf_cat_template.format(i=2, panel="1K"), sep="\t")
df1["Rsq"] = df1["R"].apply(calc_squared_with_minus)

df2 = pd.read_csv(R2_maf_cat_template.format(i=2, panel="4K"), sep="\t")
df2["Rsq"] = df2["R"].apply(calc_squared_with_minus)

df3 = pd.read_csv(R2_maf_cat_template.format(i=2, panel="10K"), sep="\t")
df3["Rsq"] = df3["R"].apply(calc_squared_with_minus)

# Keep overlapping variants across panels
# overlap_variants = set(df1["Variant"]).intersection(df2["Variant"], df3["Variant"])
# df1_overlap = df1[df1["Variant"].isin(overlap_variants)]
# df2_overlap = df2[df2["Variant"].isin(overlap_variants)]
# df3_overlap = df3[df3["Variant"].isin(overlap_variants)]

df1_overlap = df1
df2_overlap = df2
df3_overlap = df3


# Define allele frequency bins
bin_edges = [0, 0.0005, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
labels = ["0-0.05", "0.05-0.2", "0.2-0.5", "0.5-1", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100"]

# bin_edges = [0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
# labels = ["0-0.2", "0.2-0.5", "0.5-1", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100"]
df1_overlap["AF_bin"] = pd.cut(df1_overlap["ALT_FREQ"], bins=bin_edges, labels=labels)
df2_overlap["AF_bin"] = pd.cut(df2_overlap["ALT_FREQ"], bins=bin_edges, labels=labels)
df3_overlap["AF_bin"] = pd.cut(df3_overlap["ALT_FREQ"], bins=bin_edges, labels=labels)

# Aggregate mean Rsq per AF bin
agg1 = df1_overlap.groupby("AF_bin", observed=True)["Rsq"].mean().reset_index()
agg2 = df2_overlap.groupby("AF_bin", observed=True)["Rsq"].mean().reset_index()
agg3 = df3_overlap.groupby("AF_bin", observed=True)["Rsq"].mean().reset_index()

# Main Rsq vs AF_bin plot
plt.figure(figsize=(10, 10))
ax = plt.gca()

sns.lineplot(data=agg1, x="AF_bin", y="Rsq", color="dodgerblue", marker="o", markeredgecolor="k",
             markersize=10, linewidth=4, label="1K", ax=ax)
sns.lineplot(data=agg2, x="AF_bin", y="Rsq", color="limegreen", marker="o", markeredgecolor="k",
             markersize=10, linewidth=4, label="4K", ax=ax)
sns.lineplot(data=agg3, x="AF_bin", y="Rsq", color="crimson", marker="o", markeredgecolor="k",
             markersize=10, linewidth=4, label="10K", ax=ax)

ax.set_xlabel("Alt allele frequency (%)", fontsize=plt.rcParams["font.size"]+4)
ax.set_ylabel(r"Aggregated $R^2$", fontsize=plt.rcParams["font.size"]+4)
ax.set_yticks(np.arange(0.25, 0.95, 0.05))
ax.set_ylim(0.30, 0.95)
sns.despine()
ax.grid(axis="both", linestyle="--", linewidth=0.5)
ax.legend(bbox_to_anchor=(1.2, 1), frameon=False)

# Define category boundaries (by index in AF_bin)
x_positions = range(len(agg1))
cat_ranges = {
    "Low frequency": (0, 4),
    "Common": (4, 6),
    "Very common": (6, len(labels)-1)
}

# Add subtle vertical spans for categories
for label, (start, end) in cat_ranges.items():
    ax.axvspan(start, end, ymin=0, ymax=0.07, color="white", ec="black", lw=1.0, clip_on=False, zorder=3)

# Add text annotations below x-axis
for label, (start, end) in cat_ranges.items():
    mid = (start + end) / 2
    ax.text(mid, 0.04, label, transform=ax.get_xaxis_transform(),
            ha="center", va="top", weight="bold", fontsize=13, color="black")

# Set x-tick labels
ax.set_xticks(x_positions)
ax.set_xticklabels(agg1["AF_bin"].tolist())
plt.xticks(rotation=45, rotation_mode="anchor", ha="right")

# Inset subplot: Rsq ratios
ratio_10K_4K = agg3["Rsq"] / agg2["Rsq"]
ratio_4K_1K = agg2["Rsq"] / agg1["Rsq"]
ratio_10K_1K = agg3["Rsq"] / agg1["Rsq"]

w, h = 1, 1
x0 = -0.05
y0 = -0.3

ax_inset = inset_axes(ax, width="35%", height="50%", bbox_to_anchor=(x0, y0, w, h), bbox_transform=ax.transAxes)
ax_inset.plot(x_positions, ratio_10K_4K, color="orange", marker="o", markeredgecolor="k", markersize=6, label="10K / 4K")
ax_inset.plot(x_positions, ratio_4K_1K, color="grey", marker="o", markeredgecolor="k", markersize=6, label="4K / 1K")
ax_inset.plot(x_positions, ratio_10K_1K, color="purple", marker="o", markeredgecolor="k", markersize=6, label="10K / 1K")
ax_inset.set_xticks(x_positions)
ax_inset.set_xticklabels(agg1["AF_bin"], rotation=45, rotation_mode="anchor", ha="right", fontsize=plt.rcParams["font.size"]-4)
ax_inset.set_yticklabels(ax_inset.get_yticklabels(), fontsize=plt.rcParams["font.size"]-4)
ax_inset.set_xlabel("Alt allele frequency (%)", fontsize=plt.rcParams["font.size"])
ax_inset.set_ylabel(r"$R^2$ difference", fontsize=plt.rcParams["font.size"])
ax_inset.grid(axis="both", linestyle="--", linewidth=0.3)
ax_inset.legend(frameon=False, fontsize=plt.rcParams["font.size"]-2)
# sns.despine(ax=ax_inset, top=True, right=True, left=False, bottom=False)
# ax_inset.set_facecolor("#f5f5f5") 

plt.tight_layout()
plt.show()

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import ScalarFormatter

# Example data
dict_imp_variants_cnt = {"1K": df1.shape[0], "4K": df2.shape[0], "10K": df3.shape[0]}
df_imp_variants_cnt = (
    pd.DataFrame.from_dict(dict_imp_variants_cnt, orient="index")
    .reset_index()
    .rename(columns={"index": "project", 0: "count"})
)

plt.figure(figsize=(6, 3))
ax = sns.barplot(
    data=df_imp_variants_cnt[::-1],
    y="project",
    x="count",
    palette=["dodgerblue", "limegreen", "crimson"],
    edgecolor="black",
    linewidth=2
)

# Force x-axis to display as ×1e5 instead of 10⁵
formatter = ScalarFormatter(useMathText=False)
formatter.set_scientific(True)
formatter.set_powerlimits((5, 5))
ax.xaxis.set_major_formatter(formatter)

# Change the offset text (the “1e5” label at the axis)
ax.ticklabel_format(style='sci', axis='x', scilimits=(5, 5))
ax.get_xaxis().get_offset_text().set_fontsize(12)
ax.get_xaxis().get_offset_text().set_text("×e5")

# Annotate each bar with comma-separated counts
for container in ax.containers:
    ax.bar_label(
        container,
        labels=[f"{int(v.get_width()):,}" for v in container],
        label_type="center",
        color="white",
        fontsize=12,
        fontweight="bold"
    )

plt.xlabel("Number of Imputed Variants")
plt.tight_layout()
plt.show()
plt.close()

# %%
for label in labels[:2]:
    plt.hist(df3[df3["AF_bin"] == label]["Rsq"], alpha=0.5, label="1K", bins=10)
    plt.hist(df1[df1["AF_bin"] == label]["Rsq"], alpha=0.5, label="4K", bins=10)
    plt.hist(df2[df2["AF_bin"] == label]["Rsq"], alpha=0.5, label="10K", bins=10)
    # plt.ylim(0, 10000)
    plt.ylabel("Variant Count")
    plt.xlabel(r"Individual $R^2$")
    plt.yscale("log")
    plt.title(f"AF bin: {label}")
    plt.legend()
    plt.show()
    plt.close()
# %%
