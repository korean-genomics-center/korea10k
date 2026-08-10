# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# %%
path_eigenval = "/BiO/Access/kyungwhan1998/genome/pca/Resources/Data/ku10k_1KGP/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.LD_pruning_200kb_0.5.pruned.pca.eigenval"
eigenval = pd.read_csv(path_eigenval, sep="\t", header=None)
eigenval = eigenval.rename(columns={0:"VarianceExplained"})
eigenval["PropVarianceExplained"] = round(eigenval["VarianceExplained"]/sum(eigenval["VarianceExplained"])*100, 1)
eigenval

# %%
path_eigenvec = "/BiO/Access/kyungwhan1998/genome/pca/Resources/Data/ku10k_1KGP/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.LD_pruning_200kb_0.5.pruned.pca.eigenvec"
eigenvec = pd.read_csv(path_eigenvec, sep="\t")
eigenvec = eigenvec.rename(columns={"#FID":"SampleID"})

# %%
list_samples_10K = list(filter(lambda x: "10K" in x, eigenvec["SampleID"].to_list()))
list_pop_10K = ["KOR"]*len(list_samples_10K)
list_superpop_10K = ["EAS"]*len(list_samples_10K)
dict_sample_info_10K = {"SampleID": list_samples_10K, "Population": list_pop_10K, "Superpopulation": list_superpop_10K}
df_sample_info_10K = pd.DataFrame(dict_sample_info_10K)
df_sample_info_10K

# %%
path_sample_info_1KGP = "/BiO/Research/Korea10KGenome/Resources/External_Genome_Data/1KGP/1KGP_30x_GRCh38/20130606_g1k_3202_samples_ped_population.txt"
df_sample_info_1KGP = pd.read_csv(path_sample_info_1KGP, delim_whitespace=True)[["SampleID", "Population", "Superpopulation"]]

# %%
df_sample_info = pd.concat([df_sample_info_10K, df_sample_info_1KGP], axis=0)
df_sample_info

# %%
eigenvec_sample_info_merged = pd.merge(eigenvec, df_sample_info, how="inner", on="SampleID")

# %%
list_samples_exclude = ["KU10K-05465", "KU10K-08832"]
# eigenvec_sample_info_merged_filt = eigenvec_sample_info_merged[~eigenvec_sample_info_merged["SampleID"].isin(list_samples_exclude)]
plt.rcParams["font.size"] = 20
plt.figure(figsize=(15, 15))
sns.scatterplot(
    data=eigenvec_sample_info_merged,
    x="PC1",
    y="PC2",
    hue="Superpopulation",
    palette="tab10",
    alpha=0.6,
    s=100,
    edgecolor="black"
)

eigenvec_sample_info_merged_set_idx = eigenvec_sample_info_merged.set_index("SampleID")

offset = 0.00003
for i, sample  in enumerate(list_samples_exclude):
    v_shift = -offset if i % 2 == 0 else offset
    plt.text(
        eigenvec_sample_info_merged_set_idx.loc[sample, "PC1"] + offset,
        eigenvec_sample_info_merged_set_idx.loc[sample, "PC2"] + v_shift,
        eigenvec_sample_info_merged_set_idx.loc[sample, "IID"],
        fontsize=14,
        fontweight="bold",
        color="black",
        ha="left"
    )

plt.xlabel(f"PC1", fontsize=22, fontweight="bold")
plt.ylabel(f"PC2", fontsize=22, fontweight="bold")
plt.legend(frameon=False, bbox_to_anchor=(1.3, 1))
plt.tight_layout()

plt.grid(axis="both", linestyle="--")
plt.show()
plt.close()

# %%
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Merge PCA + population info
eigenvec_sample_info_merged = pd.merge(eigenvec, df_sample_info, how="inner", on="SampleID")

# Define base color palettes for each superpopulation
superpop_base_colors = {
    "EAS": "YlOrBr",
    "EUR": "Blues",
    "AFR": "Purples",
    "SAS": "Reds",
    "AMR": "Greens"
}

# Desired order for legend stacking
superpop_order = ["AFR", "AMR", "EAS", "EUR", "SAS"]

# Build Population → Superpopulation mapping
pop_to_superpop = dict(zip(df_sample_info["Population"], df_sample_info["Superpopulation"]))

# Create population color dictionary clustered by superpopulation
population_colors = {}
for sp in superpop_order:
    pops_in_sp = sorted([p for p in df_sample_info["Population"].unique() if pop_to_superpop.get(p) == sp])
    if len(pops_in_sp) == 0:
        continue
    cmap = sns.color_palette(superpop_base_colors[sp], len(pops_in_sp))
    for pop, color in zip(pops_in_sp, cmap):
        population_colors[pop] = color

# -------------------------------------------------------------------
# Main PCA plot
# -------------------------------------------------------------------
plt.rcParams["font.size"] = 18
fig, ax = plt.subplots(figsize=(12, 12))

sns.scatterplot(
    data=eigenvec_sample_info_merged,
    x="PC1",
    y="PC2",
    hue="Population",
    palette=population_colors,
    s=100,
    alpha=0.7,
    edgecolor="black",
    linewidth=0.4,
    legend=False  # we'll handle legends manually
)

# Highlight Koreans
subset_kor = eigenvec_sample_info_merged[eigenvec_sample_info_merged["Population"] == "KOR"]
ax.scatter(
    subset_kor["PC1"],
    subset_kor["PC2"],
    s=70,
    c="firebrick",
    edgecolor="black",
    linewidth=1,
    label="KOR (Korea10K)"
)

ax.set_xlabel("PC1", fontsize=22, fontweight="bold")
ax.set_ylabel("PC2", fontsize=22, fontweight="bold")
ax.grid(axis="both", linestyle="--", alpha=0.5)

# -------------------------------------------------------------------
# Create separate legends for each superpopulation
# -------------------------------------------------------------------
legend_ypos = 1.05
legend_gap = 0.225

for i, sp in enumerate(superpop_order):
    pops_in_sp = sorted([p for p in df_sample_info["Population"].unique() if pop_to_superpop.get(p) == sp])
    if not pops_in_sp:
        continue

    handles = [
        plt.Line2D(
            [0], [0],
            marker='o', color='w',
            label=pop,
            markerfacecolor=population_colors[pop],
            markeredgecolor="black",
            markersize=10
        )
        for pop in pops_in_sp
    ]

    # Add the legend to the *figure*, not the axis*
    fig.legend(
        handles,
        pops_in_sp,
        title=f"{sp}",
        title_fontsize=18,
        fontsize=14,
        frameon=False,
        loc='upper left',
        bbox_to_anchor=(1.0, legend_ypos - i * legend_gap)
    )
plt.xlim(0.02, 0.021)
plt.ylim(0.02, 0.023)
plt.tight_layout()
plt.show()
plt.close()


# %%
import matplotlib.pyplot as plt
import seaborn as sns

# List of samples to exclude
list_samples_exclude = ["KU10K-05465", "KU10K-08832"]

# Filter EAS superpopulation
eas_samples = eigenvec_sample_info_merged[
    eigenvec_sample_info_merged["Superpopulation"] == "EAS"
]
eas_samples = eas_samples[~eas_samples["SampleID"].isin(list_samples_exclude)]

# Separate Koreans and other EAS
korean_samples = eas_samples[eas_samples["Population"] == "KOR"]  # "KOR" label
other_eas_samples = eas_samples[eas_samples["Population"] != "KOR"]

# Plot settings
plt.rcParams["font.size"] = 20
plt.figure(figsize=(10, 8))

# Plot other EAS first (grey, lower alpha, lower zorder)
sns.scatterplot(
    data=other_eas_samples,
    x="PC1",
    y="PC2",
    hue="Population",
    palette="Spectral",
    alpha=0.5,
    s=100,
    edgecolor="black",
    zorder=1
)

# Plot Koreans on top (red, higher alpha, top zorder)
sns.scatterplot(
    data=korean_samples,
    x="PC1",
    y="PC2",
    color="firebrick",
    alpha=0.8,
    s=50,
    edgecolor="black",
    zorder=2,
    label="Korean"
)

korean_samples_pc1_sort = korean_samples.sort_values(by=["PC1"])
korean_samples_pc1_filt_out = korean_samples_pc1_sort[korean_samples_pc1_sort["PC1"] > -0.003]
list_samples_qc_fail = korean_samples_pc1_filt_out["SampleID"].to_list()
list_samples_exclude += list_samples_qc_fail

offset = 0.00003
for i, (idx, row) in enumerate(korean_samples_pc1_filt_out.iterrows()):
    v_shift = -offset if i % 2 == 0 else offset
    plt.text(
        row["PC1"] + offset,
        row["PC2"] + v_shift,
        row["SampleID"],
        fontsize=14,
        fontweight="bold",
        color="black",
        ha="left"
    )

# Axis labels
plt.xlabel(f"PC1", fontsize=22, fontweight="bold")
plt.ylabel(f"PC2", fontsize=22, fontweight="bold")
plt.axvline(x=-0.003, color="k", linestyle="--")
# Legend
plt.legend(frameon=False, fontsize=16)

plt.grid(axis="both", linestyle="--")
plt.tight_layout()
plt.show()
plt.close()

# %%
file_samples_exclude = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/remove_samples.list"
with open(file_samples_exclude, mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in list_samples_exclude:
        fw.write(sample + "\t" + sample + "\n")

# %%
eas_samples = eigenvec_sample_info_merged[
    eigenvec_sample_info_merged["Superpopulation"] == "EAS"
]
eas_samples = eas_samples[~eas_samples["SampleID"].isin(list_samples_exclude)]

file_samples_include = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/include_eas_samples.list"
with open(file_samples_include, mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in eas_samples["SampleID"]:
        fw.write("\t".join([sample, sample]) + "\n")

# %%
eas_samples = eigenvec_sample_info_merged[
    eigenvec_sample_info_merged["Superpopulation"] == "EAS"
]
eas_samples = eas_samples[~eas_samples["SampleID"].isin(list_samples_exclude)]

file_samples_include = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/include_eas_samples.list"
with open(file_samples_include, mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in eas_samples["SampleID"]:
        fw.write("\t".join([sample, sample]) + "\n")
        
# %%
kor_samples = eas_samples[eas_samples["Population"] == "KOR"]

file_samples_include = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/Koreans.list"
with open(file_samples_include, mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in kor_samples["SampleID"]:
        fw.write("\t".join([sample, sample]) + "\n")
