# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# %%
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"

# %%
path_eigenval = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.200kb_0.5.pca.eigenval")
eigenval = pd.read_csv(path_eigenval, sep="\t", header=None)
eigenval = eigenval.rename(columns={0:"VarianceExplained"})
eigenval["PropVarianceExplained"] = round(eigenval["VarianceExplained"]/sum(eigenval["VarianceExplained"])*100, 1)
eigenval

# %%
path_eigenvec = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.200kb_0.5.pca.eigenvec")
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
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ---- Paths and lists ----
clust_path = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/cluster.cluster2"

list_samples_annotate = [
    "U10K-00014",
    "KU10K-01224",
    "KU10K-05697",
    "KU10K-05470",
    "U10K-00017",
    "KU10K-02463",
    "KU10K-02029",
    "KU10K-00922",
]

# ---- Load cluster info ----
df_clust = pd.read_csv(
    clust_path,
    delim_whitespace=True,
    names=["FID", "IID", "Cluster"]
)

# ---- Merge PCA and cluster ----
df_merged = eigenvec_sample_info_merged.merge(
    df_clust, left_on="SampleID", right_on="IID", how="left"
)

# ---- Plot ----
plt.figure(figsize=(7, 7))

# # Base scatter (all samples)
# sns.kdeplot(
#     data=df_merged,
#     x="PC1",
#     y="PC2",
#     cmap="viridis",
#     fill=True
# )

sns.scatterplot(
    data=df_merged,
    x="PC1",
    y="PC2",
    s=50,
    linewidth=0.5,
    alpha=0.2,
    edgecolor="black",
    color="white",
    legend=False
)

# Highlighted samples (in red)
df_highlight = df_merged[df_merged["SampleID"].isin(list_samples_annotate)]
sns.scatterplot(
    data=df_highlight,
    x="PC1",
    y="PC2",
    s=20,
    color="red",
    edgecolor="black",
    linewidth=0.6,
    legend=False
)

# ---- Annotate sample names ----
for _, row in df_highlight.iterrows():
    plt.text(
        row["PC1"],
        row["PC2"],
        row["SampleID"],
        color="k",
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="bottom",
        clip_on=True,
    )

# ---- Labels and title ----
plt.title("PCA Colored by Cluster Assignment", fontsize=18, fontweight="bold")
plt.xlabel("PC1", fontsize=16, fontweight="bold")
plt.ylabel("PC2", fontsize=16, fontweight="bold")
plt.tight_layout()
plt.show()



# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ============================================================
# PARAMETERS
# ============================================================

max_pc = 5  # number of PCs to include
palette = "Spectral"
alpha = 0.5
point_size = 30
edgecolor = "black"

# ============================================================
# PREPARE DATA
# ============================================================

# eigenvec_sample_info_merged should have columns: ["SampleID", "Population", "PC1", "PC2", ..., "PC10"]
pcs = [f"PC{i}" for i in range(1, max_pc + 1)]


# Get explained variance (%)
prop_var = eigenval["PropVarianceExplained"].tolist()
if len(prop_var) < max_pc:
    prop_var += [0] * (max_pc - len(prop_var))

# ============================================================
# CREATE TRIANGULAR SUBPLOTS (PAIRWISE PC SCATTER)
# ============================================================

plt.rcParams["font.size"] = 12

fig, axes = plt.subplots(
    nrows=max_pc - 1,
    ncols=max_pc - 1,
    figsize=(12, 10),
    sharex=False,
    sharey=False
)

# Adjust layout spacing
plt.subplots_adjust(wspace=0.2, hspace=0.2)

# ============================================================
# PLOT LOWER TRIANGLE ONLY
# ============================================================

for i in range(1, max_pc):
    for j in range(0, max_pc - 1):
        ax = axes[i - 1, j]

        if j >= i:
            # upper triangle → turn off axis
            ax.axis("off")
            continue

        # Scatterplot for (PCj vs PCi)
        sns.scatterplot(
            data=df_merged,
            x=pcs[j],
            y=pcs[i],
            hue="Cluster",
            palette=palette,
            alpha=alpha,
            s=point_size,
            edgecolor=edgecolor,
            linewidth=0.5,
            ax=ax,
            legend=False
        )

        # Axis labeling (only leftmost & bottom)
        if j == 0:
            ax.set_ylabel(f"{pcs[i]}", fontsize=12, fontweight="bold")
        else:
            ax.set_ylabel("")

        if i == max_pc - 1:
            ax.set_xlabel(f"{pcs[j]}", fontsize=12, fontweight="bold")
        else:
            ax.set_xlabel("")

plt.tight_layout(rect=[0, 0, 0.95, 0.95])
plt.show()
plt.close()

# %%
