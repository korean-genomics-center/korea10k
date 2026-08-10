# %%
import os

import pandas as pd

# %%
path_ori = "/BiO/Research/Korea10KGenome/Resources/MetaData/Sequencing/KOREA10K_DATA_TABLE.xlsx"
path_add = "/BiO/Research/Korea10KGenome/Resources/Experiment_Sheets/T7_Additional_Sequencing_Experiment_Metadata.txt"
path_ori_meta = "/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.VQSR.PASS/Plink_All/Step1_MakePlink.biallelic/chr21.biallelic.psam"
# path_ori_meta = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.fam"
path_new_depth = "/BiO/Access/kyungwhan1998/genome/depthCoverage/KU10K_10243_base_mapped.txt"
path_out = "/BiO/Access/kyungwhan1998/genome/depthCoverage/Korea10K_10239_sequencing_depth.txt"

human_genome_size = 3298912062
list_samples_exclude = ["KU10K-10433", "KU10K-10689", "KU10K-04846", "KU10K-10007"]

df_ori = pd.read_excel(path_ori, sheet_name="Multiomics", engine="openpyxl")
df_add = pd.read_csv(path_add, delim_whitespace=True)
df_ori_meta = pd.read_csv(path_ori_meta, delim_whitespace=True, header=None).rename(columns={1:"IID"})
df_new_depth = pd.read_csv(path_new_depth, delim_whitespace=True, header=None)
df_new_depth.columns = ["SampleID", "Production Amount (Gbp)"]
df_new_depth["SampleID"] = df_new_depth["SampleID"].apply(lambda x: str(x).replace(".mapping.stats.txt", ""))
df_new_depth_set_ind = df_new_depth.set_index("SampleID")
list_samples_all = list(filter(lambda x: "U10K" in x, list(df_ori_meta["IID"])))
list_samples_include = list(set(list_samples_all).difference(list_samples_exclude))
df_new_depth_filt = df_new_depth_set_ind.loc[list_samples_include].reset_index()
df_new_depth_filt["Average Depth (x)"] = df_new_depth_filt["Production Amount (Gbp)"].apply(lambda x: int(x)/int(human_genome_size))
df_ori_filt = df_ori[df_ori["KU10K-ID"].isin(list_samples_include)].rename(columns={"KU10K-ID": "SampleID"})
df_ori_filt = df_ori_filt.drop(columns=["WGS Depth (avg)"])

df_ori_new_depth_added = pd.merge(df_ori_filt, df_new_depth_filt, how="inner", on="SampleID")
df_ori_new_depth_added.to_csv(path_out, sep="\t", index=False)

# %%
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.displot(df_ori_new_depth_added, x="Average Depth (x)", hue="WGS Platform")
plt.show()
plt.close()

# %%
# 1. Bin the depth into ranges
bins = [0, 10, 20, 30, 40, 50, 100, 200]
labels = ["0–10", "10–20", "20–30", "30–40", "40–50", "50–100", "100+"]
df_ori_new_depth_added["Depth Bin"] = pd.cut(df_ori_new_depth_added["Average Depth (x)"], bins=bins, labels=labels, right=False, include_lowest=True)

# 2. Count samples per platform and depth bin
count_table = df_ori_new_depth_added.groupby(["Depth Bin", "WGS Platform"]).size().unstack(fill_value=0)

# 3. Plot stacked bar chart
ax = count_table.plot(kind="bar", stacked=True, figsize=(4, 5))

# 4. Beautify
plt.title("Sample Count by Sequencing Depth and WGS Platform")
plt.xlabel("Depth Range (× Coverage)")
plt.ylabel("Sample Count")
plt.legend(title="WGS Platform", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()
plt.close()

# %%
plt.rcParams["font.size"] = 20

df_ori_new_depth_added["Rounded Depth"] = df_ori_new_depth_added["Average Depth (x)"].astype(int)

count_table = df_ori_new_depth_added.groupby(["Rounded Depth", "WGS Platform"]).size().unstack(fill_value=0)
count_table = count_table.sort_index()
ax = count_table.plot(kind="bar", stacked=True, figsize=(14, 6), color=sns.color_palette("colorblind"), zorder=3)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=plt.rcParams["font.size"]-7)
plt.xlabel("Sequencing Depth", fontsize=plt.rcParams["font.size"])
plt.ylabel("Sample Count", fontsize=plt.rcParams["font.size"])
plt.legend(title="WGS Platform", bbox_to_anchor=(1.01, 0.6), loc="upper left", frameon=False)
plt.yscale("log")
plt.grid(axis="y", zorder=1)
plt.tight_layout()
plt.show()