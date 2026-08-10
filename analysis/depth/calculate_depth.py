# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec

# %%
path_ori = "/BiO/Research/Korea10KGenome/Resources/MetaData/Sequencing/KOREA10K_DATA_TABLE.xlsx"
path_ori_meta = "/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.VQSR.PASS/Plink_All/Step1_MakePlink.biallelic/chr21.biallelic.psam"
path_new_depth = "/BiO/Access/kyungwhan1998/genome/depthCoverage/KU10K_10243_base_mapped.txt"
path_out = "/BiO/Access/kyungwhan1998/genome/depthCoverage/Korea10K_10239_sequencing_depth.txt"

human_genome_size = 3298912062
list_samples_exclude = ["KU10K-10433", "KU10K-10689", "KU10K-04846", "KU10K-10007"]

df_ori = pd.read_excel(path_ori, sheet_name="Multiomics", engine="openpyxl")
df_ori_meta = pd.read_csv(path_ori_meta, delim_whitespace=True).rename(columns={1:"IID"})
list_samples_all = list(filter(lambda x: "U10K" in x, list(df_ori_meta["IID"])))
list_samples_include = list(set(list_samples_all).difference(list_samples_exclude))
df_new_depth = pd.read_csv(path_new_depth, delim_whitespace=True, header=None)
df_new_depth.columns=["SampleID", "Production_Amount_GB"]
df_new_depth["SampleID"] = df_new_depth["SampleID"].apply(lambda x: str(x).split(".")[0])
df_new_depth = df_new_depth[df_new_depth["SampleID"].isin(list_samples_include)]
df_new_depth["Average_Depth_x"] = df_new_depth["Production_Amount_GB"].apply(lambda x: int(x)/int(human_genome_size))
platforms = ["HiSeq X", "MGI_T7", "NovaSeq"]
df_new_depth["Rounded_Depth"] = df_new_depth["Average_Depth_x"].round().astype(int)

# %%
df_ori_filt = df_ori[["KU10K-ID", "WGS Platform"]].rename(columns={"KU10K-ID": "SampleID"})
df_ori_filt_new_depth = pd.merge(df_new_depth, df_ori_filt, how="inner", on="SampleID")
df_ori_filt_new_depth

# %%
# Count table for bar plots
count_table = df_ori_filt_new_depth.groupby(["Rounded_Depth", "WGS Platform"]).size().unstack(fill_value=0)

# Pie chart data
platforms = ["HiSeq X", "MGI_T7", "NovaSeq"]
platform_colors = {"HiSeq X": "#1f77b4", "MGI_T7": "#ff7f0e", "NovaSeq": "#2ca02c"}
pie_counts = count_table.sum(axis=0).reindex(platforms).tolist()
pie_colors = [platform_colors[p] for p in platforms]

# ----------------------------------------------------------------------
# Create figure with GridSpec
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(10, 5))
gs = GridSpec(1, 2, width_ratios=[2, 1], wspace=0)

# ----------------------------------------------------------------------
# Panel A: Stacked Bar Plot (Average Depth by Platform)
# ----------------------------------------------------------------------
axA = fig.add_subplot(gs[0, 0])
axA_xticks = list(range(0, 101))  # Set x-axis from 0 to 100
bottom = np.zeros(len(axA_xticks))

for platform in platforms:
    platform_counts = count_table.get(platform, pd.Series(dtype=int)).to_dict()
    heights = [platform_counts.get(depth, 0) for depth in axA_xticks]
    axA.bar(
        axA_xticks,
        heights,
        width=1,
        color=platform_colors[platform],
        edgecolor='k',
        bottom=bottom,
        label=platform
    )
    bottom += np.array(heights)

axA.set_xlabel("Average WGS Depth (×)", fontsize=16)
axA.set_ylabel("Sample Count", fontsize=16)
axA.set_xticks(np.arange(0, 101, 10))
axA.set_ylim(top=bottom.max() + 50)
axA.text(-0.12, 1.05, "A", transform=axA.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')
axA.legend(title="WGS Platform", loc="center right", bbox_to_anchor=(0.95, 0.5), frameon=False, fontsize=13, title_fontsize=14)
axA.axvline(x=round(df_ori_filt_new_depth["Average_Depth_x"].mean(), 1), color="firebrick", linestyle="--")
axA.text(round(df_ori_filt_new_depth["Average_Depth_x"].mean(), 1)+1, 990, f"mean depth={round(df_ori_filt_new_depth['Average_Depth_x'].mean(), 1)}x", color="firebrick", weight="bold", fontsize=12)
# ----------------------------------------------------------------------
# Panel B: Pie Chart (Platform Distribution)
# ----------------------------------------------------------------------
axB = fig.add_subplot(gs[0, 1])
wedges, texts, autotexts = axB.pie(
    pie_counts,
    labels=platforms,
    startangle=90,
    colors=pie_colors,
    autopct=lambda p: f"{int(round(p * sum(pie_counts) / 100)):,}\n({p:.1f}%)",
    pctdistance=0.6,
    wedgeprops={"edgecolor": "k", "linewidth": 2},
    textprops={"color": "white", "weight": "bold", "fontsize": 11}
)
axB.text(-0.12, 1.32, "B", transform=axB.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')
sns.despine(top=True, right=True)
plt.show()

