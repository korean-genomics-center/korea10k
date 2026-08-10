# %%
import os
import re

import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

# %%
plt.rcParams.update({
    "font.size": 18,
    "axes.linewidth": 1.2,
    "axes.labelpad": 6,
    "xtick.direction": "out",
    "ytick.direction": "out"
})

plt.rcParams.update({
    "font.size": 18,
    "axes.linewidth": 1.2,
    "axes.labelpad": 6,
    "xtick.direction": "out",
    "ytick.direction": "out"
})

height_pca = 100
hspaces = [28,4,28,28,28]
height_adm = 22
height_prop = 40

fig = plt.figure(figsize=(16, 24))
row = sum(hspaces) + height_pca + height_adm*2 + height_prop*3
col = 36
gsfig = gridspec.GridSpec(
    row, col, 
    left=0, right=1, bottom=0,
    top=1, wspace=1, hspace=1)


gs_axA = gsfig[0:height_pca, 0:15]
axA = fig.add_subplot(gs_axA)

gs_axB = gsfig[0:height_pca, 21:36]
axB = fig.add_subplot(gs_axB)

gs_axC = gsfig[height_pca+sum(hspaces[:1]):height_pca+sum(hspaces[:1])+height_adm, 0:36]
gs_axC1 = gsfig[height_pca+height_adm+sum(hspaces[:2]):height_pca+height_adm*2+sum(hspaces[:2]), 0:36]
axC = fig.add_subplot(gs_axC)
axC1 = fig.add_subplot(gs_axC1)

gs_axD = gsfig[height_pca+height_adm*2+sum(hspaces[:3]):height_pca+height_adm*2+height_prop+sum(hspaces[:3]), 0:36]
axD = fig.add_subplot(gs_axD)

gs_axE = gsfig[height_pca+height_adm*2+height_prop+sum(hspaces[:4]):height_pca+height_adm*2+height_prop*2+sum(hspaces[:4]), 0:36]
axE = fig.add_subplot(gs_axE)

gs_axF = gsfig[height_pca+height_adm*2+height_prop*2+sum(hspaces[:5]):height_pca+height_adm*2+height_prop*3+sum(hspaces[:5]), 0:36]
axF = fig.add_subplot(gs_axF)

for ax in [axC1]:
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0 + 0.007, pos.width, pos.height])
    
def widen_axis(ax, factor=1.15):
    pos = ax.get_position()
    new_width = pos.width * factor
    ax.set_position([pos.x0, pos.y0, new_width, pos.height])

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

# widen_axis(axA, factor=1.05)
# heighten_axis(axC, factor=2.2, anchor="top")
# heighten_axis(axC1, factor=2.2, anchor="bottom")

workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
path_eigenval = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.preprocssed.prune_200kb_0.5.pca.eigenval")
path_eigenvec = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.preprocssed.prune_200kb_0.5.pca.eigenvec")
path_sample_info_1KGP = "/BiO/Research/Korea10KGenome/Resources/External_Genome_Data/1KGP/1KGP_30x_GRCh38/20130606_g1k_3202_samples_ped_population.txt"

eigenval = pd.read_csv(path_eigenval, sep="\t", header=None)
eigenval = eigenval.rename(columns={0:"VarianceExplained"})
eigenval["PropVarianceExplained"] = round(eigenval["VarianceExplained"]/sum(eigenval["VarianceExplained"])*100, 1)

eigenvec = pd.read_csv(path_eigenvec, sep="\t")
eigenvec = eigenvec.rename(columns={"#FID":"SampleID"})

list_samples_10K = list(filter(lambda x: "10K" in x, eigenvec["SampleID"].to_list()))
list_pop_10K = ["KOR"]*len(list_samples_10K)
list_superpop_10K = ["EAS"]*len(list_samples_10K)
dict_sample_info_10K = {"SampleID": list_samples_10K, "Population": list_pop_10K, "Superpopulation": list_superpop_10K}
df_sample_info_10K = pd.DataFrame(dict_sample_info_10K)
df_sample_info_1KGP = pd.read_csv(path_sample_info_1KGP, delim_whitespace=True)[["SampleID", "Population", "Superpopulation"]]
df_sample_info = pd.concat([df_sample_info_10K, df_sample_info_1KGP], axis=0)

eigenvec_sample_info_merged = pd.merge(eigenvec, df_sample_info, how="inner", on="SampleID")

eigenvec_sample_info_merged_kor_only = eigenvec_sample_info_merged[eigenvec_sample_info_merged["Population"] == "KOR"]
sns.scatterplot(
    data=eigenvec_sample_info_merged,
    x="PC1",
    y="PC2",
    color="grey", 
    alpha=0.5,
    s=20,
    edgecolor="black",
    legend=False,
    ax=axA
)

# Compute cluster centroids for superpopulation
centroids = eigenvec_sample_info_merged.groupby("Superpopulation")[["PC1", "PC2"]].median()

# Annotate each superpopulation near its cluster with offsets
for superpop, row in centroids.iterrows():
    if superpop == "AMR":
        axA.text(
            row["PC1"]-0.01,
            row["PC2"]+0.001,
            superpop,
            fontsize=plt.rcParams["font.size"]+2,
            fontweight="bold",
            color="k",
            ha="center",
            va="center"
        )
    elif superpop == "AFR":
        axA.text(
            row["PC1"]-0.003,
            row["PC2"]+0.005,
            superpop,
            fontsize=plt.rcParams["font.size"]+2,
            fontweight="bold",
            color="k",
            ha="center",
            va="center"
        )
    
    elif superpop == "EAS":
        axA.text(
            row["PC1"]+0.004,
            row["PC2"]+0.004,
            superpop,
            fontsize=plt.rcParams["font.size"]+2,
            fontweight="bold",
            color="k",
            ha="center",
            va="center"
        )
    else:
        axA.text(
            row["PC1"]-0.003,
            row["PC2"]+0.002,
            superpop,
            fontsize=plt.rcParams["font.size"]+2,
            fontweight="bold",
            color="k",
            ha="center",
            va="center"
        )

def plot_confidence_ellipse(x, y, ax, n_std=1.0, facecolor='none', **kwargs):
    cov = np.cov(x, y)
    mean_x = np.mean(x)
    mean_y = np.mean(y)

    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]

    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)
    
    ellipse = patches.Ellipse(
        (mean_x, mean_y),
        width,
        height,
        angle=theta,
        facecolor=facecolor,
        **kwargs
    )
    ax.add_patch(ellipse)

color_dict = {
    "EAS": "red",
    "EUR": "blue",
    "AFR": "green",
    "AMR": "purple",
    "SAS": "k"
}

superpops = eigenvec_sample_info_merged["Superpopulation"].unique()
for sp in superpops:
    if sp == "SAS":
        sp_data = eigenvec_sample_info_merged[eigenvec_sample_info_merged["Superpopulation"]==sp]
        plot_confidence_ellipse(
            sp_data["PC1"], sp_data["PC2"],
            ax=axA,
            n_std=2.0,
            edgecolor=color_dict.get(sp, "grey"),
            linestyle="-",
            alpha=1,
            linewidth=5
        )

sns.scatterplot(
    data=eigenvec_sample_info_merged_kor_only,
    x="PC1",
    y="PC2",
    color="firebrick",
    alpha=0.8,
    s=20,
    edgecolor="black",
    legend=False,
    ax=axA
)

centroid_KOR = eigenvec_sample_info_merged_kor_only.groupby("Population")[["PC1", "PC2"]].median()

for pop, row in centroid_KOR.iterrows():
    axA.text(
        row["PC1"]+0.001,
        row["PC2"]-0.004,
        pop,
        fontsize=plt.rcParams["font.size"]+8,
        fontweight="bold",
        color="firebrick",
        ha="center",
        va="center"
    )

axA.set_xlabel(f"PC1", fontsize=plt.rcParams["font.size"]+7, weight="bold")
axA.set_ylabel(f"PC2", fontsize=plt.rcParams["font.size"]+7, weight="bold")

axA.text(
    -0.22, 1.02,
    "A",
    transform=axA.transAxes,
    fontsize=plt.rcParams["font.size"]+12,
    fontweight="bold",
    va="top",
    ha="left"
)

workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
path_eigenval = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.200kb_0.5.pca.eigenval")
path_eigenvec = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.200kb_0.5.pca.eigenvec")
path_sample_info_1KGP = "/BiO/Research/Korea10KGenome/Resources/External_Genome_Data/1KGP/1KGP_30x_GRCh38/20130606_g1k_3202_samples_ped_population.txt"

eigenval = pd.read_csv(path_eigenval, sep="\t", header=None)
eigenval = eigenval.rename(columns={0:"VarianceExplained"})
eigenval["PropVarianceExplained"] = round(eigenval["VarianceExplained"]/sum(eigenval["VarianceExplained"])*100, 1)

eigenvec = pd.read_csv(path_eigenvec, sep="\t")
eigenvec = eigenvec.rename(columns={"#FID":"SampleID"})

list_samples_10K = list(filter(lambda x: "10K" in x, eigenvec["SampleID"].to_list()))
list_pop_10K = ["KOR"]*len(list_samples_10K)
list_superpop_10K = ["EAS"]*len(list_samples_10K)
dict_sample_info_10K = {"SampleID": list_samples_10K, "Population": list_pop_10K, "Superpopulation": list_superpop_10K}
df_sample_info_10K = pd.DataFrame(dict_sample_info_10K)
df_sample_info_1KGP = pd.read_csv(path_sample_info_1KGP, delim_whitespace=True)[["SampleID", "Population", "Superpopulation"]]
df_sample_info = pd.concat([df_sample_info_10K, df_sample_info_1KGP], axis=0)

eigenvec_sample_info_merged = pd.merge(eigenvec, df_sample_info, how="inner", on="SampleID")
eigenvec_sample_info_merged_kor_only = eigenvec_sample_info_merged[eigenvec_sample_info_merged["Population"] == "KOR"]

sns.scatterplot(
    data=eigenvec_sample_info_merged,
    x="PC1",
    y="PC2",
    color="grey", 
    alpha=0.5,
    s=20,
    edgecolor="black",
    legend=False,
    ax=axB
)

label_offsets = {
    "JPT": (0.01, 0),
    "CHB": (0, -0.01),
    "CHS": (0, -0.01),
    "CDX": (0, -0.01),
    "KHV": (0, -0.01)
}

centroids = eigenvec_sample_info_merged.groupby("Population")[["PC1", "PC2"]].median()

for ind, row in centroids.iterrows():
    pop = ind
    if pop != "KOR":
        dx, dy = label_offsets.get(pop, (0, 0))
        axB.text(
            row["PC1"] + dx,
            row["PC2"] + dy,
            pop,
            fontsize=plt.rcParams["font.size"]+2,
            fontweight="bold",
            color="black",
            ha="center",
            va="center",
            zorder=5
        )

sns.scatterplot(
    data=eigenvec_sample_info_merged_kor_only,
    x="PC1",
    y="PC2",
    color="firebrick", 
    alpha=0.5,
    s=20,
    edgecolor="black",
    legend=False,
    ax=axB
)

centroid_KOR = eigenvec_sample_info_merged_kor_only.groupby("Population")[["PC1", "PC2"]].median()

for pop, row in centroid_KOR.iterrows():
    axB.text(
        row["PC1"]+0.01,
        row["PC2"]-0.02,
        pop,
        fontsize=plt.rcParams["font.size"]+8,
        fontweight="bold",
        color="firebrick",
        ha="center",
        va="center"
    )

samples_highlight = ["KU10K-00922", "KU10K-02029"]
labels_highlight = ["KOREF1", "KOREF2"]
highlight_rows = eigenvec_sample_info_merged.loc[
    eigenvec_sample_info_merged["SampleID"].isin(samples_highlight)
]

for highlight_row_pc1, highlight_row_pc2, highlight_label in zip(highlight_rows["PC1"], highlight_rows["PC2"], labels_highlight):
    axB.scatter(
        highlight_row_pc1,
        highlight_row_pc2,
        color="red",
        edgecolor="black",
        s=80,
        zorder=10,
        label=str(highlight_label)
    )

    axB.text(
        highlight_row_pc1+0.007,
        highlight_row_pc2-0.003,
        str(highlight_label),
        color="k",
        fontsize=plt.rcParams["font.size"],
        fontweight="bold",
        ha="left",
        va="bottom",
        zorder=11
    )

    if highlight_label in ["KOREF1", "KOREF2"]:
        axB.plot(
            [highlight_row_pc1, highlight_row_pc1 + 0.006],
            [highlight_row_pc2, highlight_row_pc2],
            color="black",
            lw=2,
            zorder=9
        )

axB.set_xlabel(f"PC1", fontsize=plt.rcParams["font.size"]+7, weight="bold")
axB.set_ylabel(f"PC2", fontsize=plt.rcParams["font.size"]+7, weight="bold")

axB.text(
    -0.22, 1.02,
    "B",
    transform=axB.transAxes,
    fontsize=plt.rcParams["font.size"]+12,
    fontweight="bold",
    va="top",
    ha="left"
)

K = 7
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Results"
Qpop_sorted_file = os.path.join(workdir, f"admixture_plot_input_K{K}.txt")
pop_order_file = os.path.join(workdir, "pop_order.txt")

# -------------------------------
# READ DATA
# -------------------------------
Qpop_sorted = (
    pd.read_csv(Qpop_sorted_file, sep="\t", index_col=0)
    .reset_index()
    .rename(columns={"index": "SampleID"})
)
pop_order = [line.strip() for line in open(pop_order_file) if line.strip() and not line.startswith("#")]
Qpop_sorted["Population"] = pd.Categorical(Qpop_sorted["Population"], categories=pop_order, ordered=True)

# -------------------------------
# SORT BY MAX CLUSTER
# -------------------------------
cluster_cols = [f"Cluster{i+1}" for i in range(K)]
Qpop_sorted["MaxCluster"] = Qpop_sorted[cluster_cols].idxmax(axis=1)
Qpop_sorted["MaxValue"] = Qpop_sorted[cluster_cols].max(axis=1)

sorted_frames = []
for pop in pop_order:
    sub = Qpop_sorted[Qpop_sorted["Population"] == pop].copy()
    if sub.empty:
        continue
    sub = sub.sort_values(by=["MaxCluster", "MaxValue"], ascending=[True, False])[::-1]
    sorted_frames.append(sub)
Qpop_sorted = pd.concat(sorted_frames).reset_index(drop=True)

# -------------------------------
# COLORS
# -------------------------------
cluster_colors = sns.color_palette([
    "#984EA3", "#0047A0", "#FF5E00", "#EBB112", "#FFFF33", "#F781BF", "#01A4BA"
])

x = np.arange(len(Qpop_sorted))
bottom = np.zeros(len(Qpop_sorted))

for i, cluster in enumerate(cluster_cols):
    axC.bar(
        x,
        Qpop_sorted[cluster],
        bottom=bottom,
        width=1.0,
        color=cluster_colors[i],
        edgecolor="none",
        label=cluster,
    )
    bottom += Qpop_sorted[cluster].values

# -------------------------------
# POPULATION BOUNDARIES
# -------------------------------
pop_boundaries, pop_labels, pop_labels_names = [], [], []
start = 0
for pop in pop_order:
    subset = Qpop_sorted[Qpop_sorted["Population"] == pop]
    if subset.empty:
        continue
    end = start + len(subset)
    pop_boundaries.append(end)
    pop_labels.append((start + end) / 2)
    pop_labels_names.append(pop)
    start = end

# AXIS STYLE
for axis in ['top', 'bottom', 'left', 'right']:
    axC.spines[axis].set_linewidth(2)
for b in pop_boundaries[:-1]:
    axC.axvline(b, color="black", lw=2)
axC.set_xticks(pop_labels)
axC.set_xticklabels("", fontsize=0)
axC.set_xlabel("")
axC.set_yticks([])
axC.set_xmargin(0)
axC.set_ymargin(0)

# K LABEL
axC.text(-0.01, 0.5, f"K={K}", transform=axC.transAxes,
         fontsize=plt.rcParams["font.size"] + 7, fontweight="bold",
         rotation=0, va="center", ha="right")

# -------------------------------
# SUPERPOPULATION DICTIONARY
# -------------------------------
superpop_dict = {
    "AFR": ["YRI", "LWK", "GWD", "MSL", "ESN", "ACB", "ASW"],
    "EUR": ["CEU", "GBR", "FIN", "IBS", "TSI"],
    "EAS": ["CHB", "CHS", "JPT", "KHV", "CDX", "KOR"],
    "SAS": ["PJL", "BEB", "GIH", "ITU", "STU"],
    "AMR": ["CLM", "MXL", "PUR", "PEL"]
}

# Reverse mapping: population -> superpopulation
pop_to_superpop = {pop: sp for sp, pop_list in superpop_dict.items() for pop in pop_list}

# Convert ordered population list → superpopulation list safely
superpops = [pop_to_superpop.get(p, "OTHER") for p in pop_labels_names]

# -------------------------------
# SUPERPOPULATION RANGES
# -------------------------------
superpop_ranges = []
sp_current = superpops[0]
start = 0
for i in range(1, len(superpops)):
    if superpops[i] != sp_current:
        superpop_ranges.append((sp_current, start, i))
        sp_current = superpops[i]
        start = i
superpop_ranges.append((sp_current, start, len(superpops)))

# -------------------------------
# PLOT SUPERPOPULATION LINES & LABELS
# -------------------------------
y_line = 1.1
y_text = 1.15

# compute population start indices
pop_start_idx = [0] + pop_boundaries[:-1]

for sp, start_idx, end_idx in superpop_ranges:
    # first sample of the first population in superpop
    first_sample_idx = pop_start_idx[start_idx]
    # last sample of the last population in superpop
    last_sample_idx = pop_boundaries[end_idx - 1] - 1

    # convert to x coordinates
    x_start = first_sample_idx
    x_end = last_sample_idx + 1  # +1 to include the last sample

    # draw horizontal line for superpopulation
    axC.plot([x_start+20, x_end-20], [y_line, y_line],
             transform=axC.get_xaxis_transform(),
             color='black', lw=3, clip_on=False)

    # place superpopulation text centered
    axC.text((x_start + x_end) / 2, y_text, sp,
             transform=axC.get_xaxis_transform(),
             ha='center', va='bottom',
             fontsize=plt.rcParams['font.size'] + 4,
             fontweight='bold')
    
K = 9
Qpop_sorted_file = os.path.join(workdir, f"admixture_plot_input_K{K}.txt")

Qpop_sorted = (
    pd.read_csv(Qpop_sorted_file, sep="\t", index_col=0)
    .reset_index()
    .rename(columns={"index": "SampleID"})
)

Qpop_sorted["Population"] = pd.Categorical(Qpop_sorted["Population"], categories=pop_order, ordered=True)

cluster_cols = [f"Cluster{i+1}" for i in range(K)]
Qpop_sorted["MaxCluster"] = Qpop_sorted[cluster_cols].idxmax(axis=1)
Qpop_sorted["MaxValue"] = Qpop_sorted[cluster_cols].max(axis=1)

sorted_frames = []
for pop in pop_order:
    sub = Qpop_sorted[Qpop_sorted["Population"] == pop].copy()
    if sub.empty:
        continue
    sub = sub.sort_values(
        by=["MaxCluster", "MaxValue"], ascending=[True, False]
    )[::-1]
    sorted_frames.append(sub)
Qpop_sorted = pd.concat(sorted_frames).reset_index(drop=True)

cluster_colors = sns.color_palette([
    "#984EA3", "#0047A0", "#EBB112", "#01A4BA", "#FFFF33",
    "#4DAF4A", "#CD2E3A", "#F781BF", "#FF5E00"
])

x = np.arange(len(Qpop_sorted))
bottom = np.zeros(len(Qpop_sorted))

for i, cluster in enumerate(cluster_cols):
    axC1.bar(
        x,
        Qpop_sorted[cluster],
        bottom=bottom,
        width=1.0,
        color=cluster_colors[i],
        edgecolor="none",
        label=cluster,
    )
    bottom += Qpop_sorted[cluster].values

pop_boundaries, pop_labels, pop_labels_names = [], [], []
start = 0
for pop in pop_order:
    subset = Qpop_sorted[Qpop_sorted["Population"] == pop]
    if subset.empty:
        continue
    end = start + len(subset)
    pop_boundaries.append(end)
    pop_labels.append((start + end) / 2)
    pop_labels_names.append(pop)
    start = end

for axis in ['top', 'bottom', 'left', 'right']:
    axC1.spines[axis].set_linewidth(2)

for b in pop_boundaries[:-1]:
    axC1.axvline(b, color="black", lw=2)

axC1.set_xticks(pop_labels)
axC1.set_xticklabels(
    pop_labels_names,
    rotation=45, rotation_mode="anchor", ha="right",
    fontsize=plt.rcParams["font.size"] + 2,
)
axC1.set_xlabel("Individuals (grouped by population)", fontsize=plt.rcParams["font.size"] + 5, weight="bold")
axC1.set_yticks([])
axC1.set_xmargin(0)
axC1.set_ymargin(0)

axC1.text(-0.01, 0.5, "K=9", transform=axC1.transAxes,
          fontsize=plt.rcParams["font.size"] + 7, fontweight="bold",
          rotation=0, va="center", ha="right")

axC.text(
    -0.09, 1.05,
    "C",
    transform=axC.transAxes,
    fontsize=plt.rcParams["font.size"] + 10,
    fontweight="bold",
    va="top",
    ha="left"
)
# %%

# path_Y_10K = "/BiO/Access/kyungwhan1998/genome/yhaplo/output/Korea10K_Male_Only/haplogroups.chrY.biallelic.hg19.liftover.haploid.fixed.male.filtered.txt"
path_Y_10K = "/BiO/Access/kyungwhan1998/genome/ychrom/Korea10K/haplogroups.korea10K.jointcall.removeDup.biallelic.chrY.liftover.txt"
path_Y_1KGP = "/BiO/Access/kyungwhan1998/genome/yhaplo/output/1000Genomes/haplogroups.1000Y.all.txt"

path_plink_fam = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.preprocssed.prune_200kb_0.5.pruned.fam"
df_Y_10K = pd.read_csv(path_Y_10K, delim_whitespace=True, header=None)
df_Y_10K.columns = ["SampleID", "Haplogroup", "Macro_haplogroup", "Simplified_haplogroup"]
df_Y_10K_excl_Ahaplo = df_Y_10K[~df_Y_10K["Macro_haplogroup"].str.startswith("A")]
df_Y_1KGP = pd.read_csv(path_Y_1KGP, delim_whitespace=True, header=None)
df_Y_1KGP.columns = ["SampleID", "Haplogroup", "Macro_haplogroup", "Simplified_haplogroup"]
df_Y = pd.concat([df_Y_10K_excl_Ahaplo, df_Y_1KGP], axis=0)
df_Y.columns = ["SampleID", "Haplogroup", "Macro_haplogroup", "Simplified_haplogroup"]

# Simplify haplogroups
def simplify_haplogroup(h):
    if str(h).startswith("O"):
        match = re.match(r"([A-Z]+\d*)", h)
    else:
        match = re.match(r"([A-Z]+)", h)
    return match.group(1) if match else h

df_Y["Major_Haplogroup"] = df_Y["Simplified_haplogroup"].apply(simplify_haplogroup)

df_Y_sample_info = pd.merge(df_Y, df_sample_info, on="SampleID", how="inner")

with open(path_plink_fam, "r") as fr:
    list_samples_filter_in = [line.split()[0] for line in fr.readlines()]

df_Y_sample_info_filtered = df_Y_sample_info[
    df_Y_sample_info["SampleID"].isin(list_samples_filter_in)
]

pop_freq = (
    df_Y_sample_info_filtered.groupby(["Population", "Major_Haplogroup"])
    .size()
    .reset_index(name="Count")
)
pop_freq["Group"] = pop_freq["Population"]

super_freq = (
    df_Y_sample_info_filtered.groupby(["Superpopulation", "Major_Haplogroup"])
    .size()
    .reset_index(name="Count")
)
super_freq["Group"] = super_freq["Superpopulation"]

freq_df = pd.concat([pop_freq, super_freq], axis=0)
freq_df["Proportions"] = freq_df.groupby("Group")["Count"].transform(lambda x: x / x.sum())

freq_df[freq_df["Population"]=="KOR"].sort_values(by="Count", ascending=False)
freq_df.to_excel("/BiO/Access/kyungwhan1998/genome/paper/Supplementary_Table_Y.xlsx")

haplo_order= (
    freq_df.groupby("Major_Haplogroup")["Proportions"]
    .mean()
    .sort_values(ascending=False)
    .index.tolist()
)
all_haplos = freq_df["Major_Haplogroup"].unique().tolist()
haplo_order = haplo_order+ [h for h in all_haplos if h not in haplo_order]

group_order = ["AFR", "EUR", "CHB", "JPT", "KOR"]
pivot_df = freq_df.pivot(index="Group", columns="Major_Haplogroup", values="Proportions").fillna(0)
pivot_df = pivot_df.reindex([g for g in group_order if g in pivot_df.index])

palette = sns.color_palette("tab20", n_colors=len(haplo_order))
haplo_color_dict = dict(zip(haplo_order, palette))

added_labels = set()
for idx, group in enumerate(pivot_df.index):
    row = pivot_df.loc[group]
    row_sorted = row[row > 0].sort_values(ascending=False)

    left = 0
    for haplo, val in row_sorted.items():
        color = haplo_color_dict.get(haplo, "lightgrey") if val > 0.1 else "lightgrey"

        axD.barh(
            group,
            val,
            left=left,
            height=1,
            color=color,
            linewidth=0.5,
            edgecolor="black"
        )

        if val > 0.1:
            text = f"{haplo} ({val*100:.1f}%)"
            axD.text(
                left + val / 2,
                idx,
                text,
                va="center",
                ha="center",
                color="k",
                fontsize=plt.rcParams["font.size"],
                weight="bold",
            )


        left += val

axD.set_xlim(0, 1)
axD.set_xticks(np.arange(0, 1.1, 0.1))
axD.set_xlabel("Proportions", fontsize=plt.rcParams["font.size"]+5, weight="bold")
axD.set_ylabel("")
axD.tick_params(axis="y", labelsize=plt.rcParams["font.size"]+2)
axD.tick_params(axis="x", labelsize=plt.rcParams["font.size"])

for axis in ["top", "bottom", "left", "right"]:
    axD.spines[axis].set_linewidth(2)
for i, label in enumerate(pivot_df.index):
    axD.axhline(i + 0.5, color="black", lw=2)

axD.set_ymargin(0)

axD.text(
    -0.09, 1.05,
    "D",
    transform=axD.transAxes,
    fontsize=plt.rcParams["font.size"] + 10,
    fontweight="bold",
    va="top",
    ha="left"
)

# %%
path_mito_10K = "/BiO/Access/kyungwhan1998/genome/mitochondria/Resources/Korea10K.merged_chrM.final.sorted.normed.10239Samples_haplogrep3.txt"
path_mito_1KGP = "/BiO/Access/kyungwhan1998/genome/mitochondria/Resources/1KGP_30X_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrM_filtered.recalibrated_variants.haplogrep3.txt"
path_plink_fam = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.preprocssed.prune_200kb_0.5.pruned.fam"

df_mito_10K = pd.read_csv(path_mito_10K, delim_whitespace=True)
df_mito_1KGP = pd.read_csv(path_mito_1KGP, delim_whitespace=True)
df_mito = pd.concat([df_mito_10K, df_mito_1KGP], axis=0)

# Simplify haplogroups to major letter
def simplify_haplogroup(h):
    match = re.match(r"([A-Z]+)", h)
    return match.group(1) if match else h

df_mito["Major_Haplogroup"] = df_mito["Haplogroup"].apply(simplify_haplogroup)

# Merge with sample info (assuming df_sample_info exists)
df_mito_sample_info = pd.merge(df_mito, df_sample_info, on="SampleID", how="inner")

# Filter samples based on PLINK fam
with open(path_plink_fam, "r") as fr:
    list_samples_filter_in = [line.split()[0] for line in fr.readlines()]

df_mito_sample_info_filtered = df_mito_sample_info[
    df_mito_sample_info["SampleID"].isin(list_samples_filter_in)
]

pop_freq = (
    df_mito_sample_info_filtered.groupby(["Population", "Major_Haplogroup"])
    .size()
    .reset_index(name="Count")
)
pop_freq["Group"] = pop_freq["Population"]

super_freq = (
    df_mito_sample_info_filtered.groupby(["Superpopulation", "Major_Haplogroup"])
    .size()
    .reset_index(name="Count")
)
super_freq["Group"] = super_freq["Superpopulation"]

freq_df = pd.concat([pop_freq, super_freq], axis=0)
freq_df["Proportions"] = freq_df.groupby("Group")["Count"].transform(lambda x: x / x.sum())

freq_df[freq_df["Population"]=="KOR"].sort_values(by="Count", ascending=False)
freq_df.to_excel("/BiO/Access/kyungwhan1998/genome/paper/Supplementary_Table_MT.xlsx")

haplo_order= (
    freq_df.groupby("Major_Haplogroup")["Proportions"]
    .mean()
    .sort_values(ascending=False)
    .index.tolist()
)

all_haplos = freq_df["Major_Haplogroup"].unique().tolist()
haplo_order = haplo_order+ [h for h in all_haplos if h not in haplo_order]

group_order = ["AFR", "EUR", "CHB", "JPT", "KOR"]

pivot_df = freq_df.pivot(index="Group", columns="Major_Haplogroup", values="Proportions").fillna(0)
pivot_df = pivot_df.reindex([g for g in group_order if g in pivot_df.index])

palette = sns.color_palette("tab20", n_colors=len(haplo_order))
haplo_color_dict = dict(zip(haplo_order, palette))

added_labels = set()
for idx, group in enumerate(pivot_df.index):  # pivot_df = MT haplogroup pivot table
    row = pivot_df.loc[group]
    row_sorted = row[row > 0].sort_values(ascending=False)

    left = 0
    for haplo, val in row_sorted.items():
        color = haplo_color_dict.get(haplo, "lightgrey") if val > 0.1 else "lightgrey"

        axE.barh(
            group,
            val,
            left=left,
            height=1,
            color=color,
            linewidth=0.5,
            edgecolor="black"
        )

        if val > 0.1:
            text = f"{haplo} ({val*100:.1f}%)"
            axE.text(left + val / 2, idx, text, va="center", ha="center",
                     color="k", fontsize=plt.rcParams["font.size"], weight="bold")
        left += val

axE.set_xlim(0, 1)
axE.set_xticks(np.arange(0, 1.1, 0.1))
axE.set_xlabel("Proportions", fontsize=plt.rcParams["font.size"]+5, weight="bold")
axE.set_ylabel("")
axE.tick_params(axis="y", labelsize=plt.rcParams["font.size"]+2)
axE.tick_params(axis="x", labelsize=plt.rcParams["font.size"])

for axis in ["top", "bottom", "left", "right"]:
    axE.spines[axis].set_linewidth(2)
for i, label in enumerate(pivot_df.index):
    axE.axhline(i + 0.5, color="black", lw=2)

axE.set_ymargin(0)

axE.text(
    -0.09, 1.05,
    "E",
    transform=axE.transAxes,
    fontsize=plt.rcParams["font.size"] + 10,
    fontweight="bold",
    va="top",
    ha="left"
)

path_hla_file = "/BiO/Access/kyungwhan1998/genome/hla/Resources/Data/1000GplusKOR_HLA.tsv"
df_hla = pd.read_csv(path_hla_file, sep="\t")
df_hla_filtered = df_hla[
    df_hla["SampleID"].isin(list_samples_filter_in)
]

hla_alleles = pd.concat([
    df_hla_filtered[['SampleID', 'A_1', 'Population', 'Superpopulation']].rename(columns={'A_1':'Allele'}),
    df_hla_filtered[['SampleID', 'A_2', 'Population', 'Superpopulation']].rename(columns={'A_2':'Allele'})
], axis=0)

# Compute allele frequencies
pop_freq = (
    hla_alleles.groupby(['Population', 'Allele'])
    .size()
    .reset_index(name='Count')
)
pop_freq['Group'] = pop_freq['Population']

super_freq = (
    hla_alleles.groupby(['Superpopulation', 'Allele'])
    .size()
    .reset_index(name='Count')
)
super_freq['Group'] = super_freq['Superpopulation']

freq_df = pd.concat([pop_freq, super_freq], axis=0)
freq_df['Proportions'] = freq_df.groupby('Group')['Count'].transform(lambda x: x / x.sum())

freq_df[freq_df["Population"]=="KOR"].sort_values(by="Count", ascending=False)
freq_df.to_excel("/BiO/Access/kyungwhan1998/genome/paper/Supplementary_Table_HLA_A.xlsx")

allele_order= (
    freq_df.groupby("Allele")["Proportions"]
    .mean()
    .sort_values(ascending=False)
    .index.tolist()
)

all_alleles = freq_df['Allele'].unique().tolist()
allele_order = allele_order+ [a for a in all_alleles if a not in allele_order]

# Group display order
group_order = ["AFR", "EUR", "CHB", "JPT", "KOR"]

# Pivot table
pivot_df = freq_df.pivot(index='Group', columns='Allele', values='Proportions').fillna(0)
pivot_df = pivot_df.reindex([g for g in group_order if g in pivot_df.index])

# Color palette
palette = sns.color_palette("tab20", n_colors=len(allele_order))
allele_color_dict = dict(zip(allele_order, palette))

added_labels = set()

for idx, group in enumerate(pivot_df.index):
    row = pivot_df.loc[group]
    row_sorted = row[row > 0].sort_values(ascending=False)

    left = 0
    for allele, val in row_sorted.items():
        # Only color bars with proportions > 0.05
        color = allele_color_dict.get(allele, 'lightgrey') if val > 0.1 else 'lightgrey'

        # Add label only once for legend if colored
        label = allele if (allele not in added_labels and val > 0.1) else None
        if label:
            added_labels.add(allele)

        axF.barh(
            group,
            val,
            left=left,
            height=1,
            color=color,
            linewidth=0.5,
            edgecolor='black',
            label=label
        )

        # Annotate only colored bars
        if val > 0.1:
            axF.text(
                left + val / 2,
                idx,
                f"{val * 100:.1f}%",
                va='center',
                ha='center',
                color='k',
                fontsize=plt.rcParams['font.size'],
                weight='bold',
            )
        left += val

# Formatting
axF.set_xlim(0, 1)
axF.set_xticks(np.arange(0, 1.1, 0.1))
axF.set_xlabel('Proportions', fontsize=plt.rcParams['font.size']+5, weight="bold")
axF.set_ylabel('')
axF.tick_params(axis='y', labelsize=plt.rcParams['font.size']+2)
axF.tick_params(axis='x', labelsize=plt.rcParams['font.size'])

# Axis and separators
for axis in ['top','bottom','left','right']:
    axF.spines[axis].set_linewidth(2)
for i, label in enumerate(pivot_df.index):
    axF.axhline(i + 0.5, color='black', lw=2)

# Legend: top 10 alleles
handles, labels = axF.get_legend_handles_labels()

axF.legend(
    handles,
    labels,
    loc='upper center',
    bbox_to_anchor=(0.5, -0.5),
    ncol=len(handles),
    title='HLA-A Allele',
    title_fontproperties={"weight": "bold", "size": plt.rcParams["font.size"]+2},
    frameon=True
)

axF.set_ymargin(0)

axF.text(
    -0.09, 1.05,
    "F",
    transform=axF.transAxes,
    fontsize=plt.rcParams["font.size"] + 10,
    fontweight="bold",
    va="top",
    ha="left"
)
