# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# %%
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"

# %%
path_eigenval = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.preprocssed.prune_200kb_0.5.pca.eigenval")

eigenval = pd.read_csv(path_eigenval, sep="\t", header=None)
eigenval = eigenval.rename(columns={0:"VarianceExplained"})
eigenval["PropVarianceExplained"] = round(eigenval["VarianceExplained"]/sum(eigenval["VarianceExplained"])*100, 1)
eigenval

# %%
path_eigenvec = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.preprocssed.prune_200kb_0.5.pca.eigenvec")

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
plt.rcParams["font.size"] = 20
plt.figure(figsize=(12, 12))
sns.scatterplot(
    data=eigenvec_sample_info_merged,
    x="PC1",
    y="PC2",
    hue="Population",
    palette="Set2",
    alpha=0.5,
    s=200,
    edgecolor="black"
)

plt.xlabel(f"PC1", fontsize=22, fontweight="bold")
plt.ylabel(f"PC2", fontsize=22, fontweight="bold")
plt.legend(frameon=False, bbox_to_anchor=(1.05, 1))
plt.tight_layout()

plt.show()
plt.close()

# %%
plt.rcParams["font.size"] = 20
plt.figure(figsize=(12, 12))
sns.scatterplot(
    data=eigenvec_sample_info_merged,
    x="PC1",
    y="PC3",
    hue="Population",
    palette="Set2",
    alpha=0.5,
    s=200,
    edgecolor="black"
)

plt.xlabel(f"PC1", fontsize=22, fontweight="bold")
plt.ylabel(f"PC2", fontsize=22, fontweight="bold")
plt.legend(frameon=False, bbox_to_anchor=(1.05, 1))
plt.tight_layout()

plt.show()
plt.close()

# %%
import matplotlib.pyplot as plt
import seaborn as sns

# Filter EAS superpopulation
eas_samples = eigenvec_sample_info_merged[
    eigenvec_sample_info_merged["Superpopulation"] == "EAS"
]
# Separate Koreans and other EAS
korean_samples = eas_samples[eas_samples["Population"] == "KOR"]  # "KOR" label
other_eas_samples = eas_samples[eas_samples["Population"] != "KOR"]

# %%
plt.rcParams["font.size"] = 20
plt.figure(figsize=(12, 12))

# ============================================================
# Plot background samples (other EAS)
# ============================================================
sns.scatterplot(
    data=other_eas_samples,
    x="PC1",
    y="PC2",
    hue="Population",
    palette="Set2",
    alpha=0.5,
    s=100,
    edgecolor=None,
    zorder=1
)

# ============================================================
# Plot Koreans
# ============================================================
sns.scatterplot(
    data=korean_samples,
    x="PC1",
    y="PC2",
    color="red",
    alpha=0.8,
    s=100,
    edgecolor="black",
    linewidth=2,
    zorder=2,
    label="Korean"
)

# ============================================================
# Annotate alternating left/right
# ============================================================
list_samples_annotate = korean_samples[
    np.logical_or(korean_samples["PC2"] < -0.05, korean_samples["PC2"] > 0.03)
]["SampleID"].to_list()

annot_df = korean_samples[korean_samples["SampleID"].isin(list_samples_annotate)]

# Alternate left/right by row index
for i, (_, row) in enumerate(annot_df.iterrows()):
    if i % 2 == 0:
        # Right side label
        x_offset = 0.002
        ha = "left"
    else:
        # Left side label
        x_offset = -0.002
        ha = "right"

    plt.text(
        row["PC1"] + x_offset,
        row["PC2"],
        row["SampleID"],
        fontsize=14,
        fontweight="bold",
        color="black",
        ha=ha,
        va="center",
        zorder=3
    )

# ============================================================
# Axis labels and legend
# ============================================================
plt.xlabel(f"PC1", fontsize=22, fontweight="bold")
plt.ylabel(f"PC2", fontsize=22, fontweight="bold")
plt.legend(frameon=False, fontsize=16)

plt.tight_layout()
plt.show()
plt.close()


# %%
list_samples_korean = set(korean_samples["SampleID"]) - set(list_samples_annotate)

# %%
with open(os.path.join(workdir, "nonKoreans.list"), mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in list_samples_annotate:
        fw.write("\t".join([sample, sample]) + "\n")

with open(os.path.join(workdir, "remove_samples.list"), mode="r") as fr:
    list_samples_already_removed = list(map(lambda x: x.split("\t")[0], fr.readlines()[1:]))

list_samples_already_removed_plus_annot = list_samples_already_removed + list_samples_annotate
with open(os.path.join(workdir, "removePlusnonKoreans.list"), mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in list_samples_already_removed_plus_annot:
        fw.write(sample + "\t" + sample + "\n")

# %%
with open(os.path.join(workdir, "Koreans.list"), mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in list_samples_korean:
        fw.write("\t".join([sample, sample]) + "\n")

# %%
size=500
list_koreans_random_sampled = list(np.random.choice(list(list_samples_korean), size=size, replace=False))
with open(os.path.join(workdir, f"Koreans_random_sampled_{size}.list"), mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in list_koreans_random_sampled:
        fw.write("\t".join([sample, sample]) + "\n")

list_other_eas = other_eas_samples["SampleID"].to_list()
list_koreans_random_sampled_plus_1kgp_eas = list_koreans_random_sampled + list_other_eas
with open(os.path.join(workdir, f"Koreans_random_sampled_{size}+1KGP_EAS.list"), mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in list_koreans_random_sampled_plus_1kgp_eas:
        fw.write("\t".join([sample, sample]) + "\n")

# %%

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# ----------------------------
# Your setup
# ----------------------------
list_sample_annotate = [
    "U10K-00014", "KU10K-01224", "KU10K-05697",
    "KU10K-05470", "U10K-00017", "KU10K-02463",
    "KU10K-02029", "KU10K-00922"
]

plt.rcParams["font.size"] = 15

# Split Koreans vs non-Koreans
df_korean = eigenvec_sample_info_merged[eigenvec_sample_info_merged["Population"] == "KOR"]
df_non_korean = eigenvec_sample_info_merged[eigenvec_sample_info_merged["Population"] != "KOR"]

# Non-Korean palette
unique_pops_non_korean = df_non_korean["Population"].unique()
palette_non_korean = dict(zip(unique_pops_non_korean, sns.color_palette("Set2", n_colors=len(unique_pops_non_korean))))

# Annotated subset
df_annotated = eigenvec_sample_info_merged[eigenvec_sample_info_merged["SampleID"].isin(list_sample_annotate)]

# ----------------------------
# Compute genetic centroid and most central Korean sample (PC1~PC20)
# ----------------------------
pcs = [f"PC{i}" for i in range(1, 21)]
korean_coords = df_korean[pcs]

# Centroid in 20D
centroid = korean_coords.mean(axis=0)

# Euclidean distance in 20D
distances = np.linalg.norm(korean_coords.values - centroid.values, axis=1)
most_central_idx = distances.argmin()
most_central_sample = df_korean.iloc[most_central_idx]
most_central_id = most_central_sample["SampleID"]

# ----------------------------
# Scatter plot
# ----------------------------
plt.figure(figsize=(14, 12))

# Koreans (grey)
sns.scatterplot(
    data=df_korean,
    x="PC1",
    y="PC2",
    color="white",
    edgecolor="grey",
    s=50,
    alpha=0.5,
    legend=False
)

# Non-Koreans (colored by population)
sns.scatterplot(
    data=df_non_korean,
    x="PC1",
    y="PC2",
    hue="Population",
    palette=palette_non_korean,
    edgecolor="black",
    s=50,
    alpha=0.6,
    legend=False
)

# Annotated samples (red)
sns.scatterplot(
    data=df_annotated,
    x="PC1",
    y="PC2",
    color="firebrick",
    edgecolor="black",
    s=150,
    legend=False
)

# Genetic centroid (gold star)
plt.scatter(
    centroid["PC1"], centroid["PC2"],
    color="gold", s=1000, marker="*", edgecolor="black", linewidth=1.5, label="Centroid"
)

# Most central Korean individual (limegreen cross)
plt.scatter(
    most_central_sample["PC1"], most_central_sample["PC2"],
    color="limegreen", s=500, marker="o", edgecolor="black", linewidth=1.5, label="Most Central KOR"
)

# Annotate KOGIC samples
for _, row in df_annotated.iterrows():
    plt.text(
        row["PC1"],
        row["PC2"],
        row["SampleID"],
        fontsize=15,
        fontweight="bold",
        color="k"
    )


plt.text(
    most_central_sample["PC1"], most_central_sample["PC2"], most_central_id,
    fontsize=15, fontweight="bold", color="k"
)

# ----------------------------
# Manual legend
# ----------------------------
legend_elements = [
    mpatches.Patch(facecolor="lightgrey", edgecolor="black", label="KOR")
]
for pop, color in palette_non_korean.items():
    legend_elements.append(mpatches.Patch(facecolor=color, edgecolor="black", label=pop))
legend_elements.append(mpatches.Patch(facecolor="firebrick", edgecolor="black", label="KOGIC"))
legend_elements.append(mpatches.Patch(facecolor="gold", edgecolor="black", label="Centroid"))
legend_elements.append(mpatches.Patch(facecolor="limegreen", edgecolor="black", label="Most Central KOR"))

plt.legend(handles=legend_elements, fontsize=15, bbox_to_anchor=(1.05, 0.7), frameon=False)

plt.xlabel("PC1", fontweight="bold")
plt.ylabel("PC2", fontweight="bold")
plt.tight_layout()
plt.show()
plt.close()


# %%
