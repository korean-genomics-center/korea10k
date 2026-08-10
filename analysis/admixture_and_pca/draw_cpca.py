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
path_eigenval = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.include_eas_samples_only.postmerge_QC_filtered.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.200kb_0.5.pca.eigenval")

eigenval = pd.read_csv(path_eigenval, sep="\t", header=None)
eigenval = eigenval.rename(columns={0:"VarianceExplained"})
eigenval["PropVarianceExplained"] = round(eigenval["VarianceExplained"]/sum(eigenval["VarianceExplained"])*100, 1)
eigenval

# %%
path_eigenvec = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.include_eas_samples_only.postmerge_QC_filtered.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.200kb_0.5.pca.eigenvec")

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
# ============================================================
# PARAMETERS
# ============================================================

max_pc = 10  # number of PCs to include
palette = "Set2"
alpha = 0.5
point_size = 80
edgecolor = "black"

# ============================================================
# PREPARE DATA
# ============================================================

# eigenvec_sample_info_merged should have columns: ["SampleID", "Population", "PC1", "PC2", ..., "PC10"]
pcs = [f"PC{i}" for i in range(1, max_pc + 1)]

# Optionally exclude outliers:
list_samples_exclude = ["KU10K-03207", "KU10K-03220", "KU10K-03341", "KU10K-03238"]
df = eigenvec_sample_info_merged[~eigenvec_sample_info_merged["SampleID"].isin(list_samples_exclude)].copy()

# Get explained variance (%)
prop_var = eigenval["PropVarianceExplained"].tolist()
if len(prop_var) < max_pc:
    prop_var += [0] * (max_pc - len(prop_var))

# ============================================================
# CREATE TRIANGULAR SUBPLOTS (PAIRWISE PC SCATTER)
# ============================================================

sns.set_style("whitegrid")
plt.rcParams["font.size"] = 12

fig, axes = plt.subplots(
    nrows=max_pc - 1,
    ncols=max_pc - 1,
    figsize=(20, 20),
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
            data=df,
            x=pcs[j],
            y=pcs[i],
            hue="Population",
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
            ax.set_ylabel(f"{pcs[i]} ({prop_var[i]:.2f}%)", fontsize=12, fontweight="bold")
        else:
            ax.set_ylabel("")

        if i == max_pc - 1:
            ax.set_xlabel(f"{pcs[j]} ({prop_var[j]:.2f}%)", fontsize=12, fontweight="bold")
        else:
            ax.set_xlabel("")

# ============================================================
# ADD LEGEND OUTSIDE
# ============================================================

# Create a single legend for populations
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, title="Population", bbox_to_anchor=(1.05, 0.5), loc="center left", frameon=False)

plt.tight_layout(rect=[0, 0, 0.95, 0.95])
plt.show()
plt.close()
