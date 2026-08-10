# %%
import glob
import os
import subprocess

import numpy as np
import pandas as pd
import seaborn as sns

# %%
plink = "/BiO/Share/Tool/plink"
window_size = "200kb"
r2 = "0.5"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
os.makedirs(workdir, exist_ok=True)
k = 500
bfile_in = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}")

# %%
cmd = f"{plink} --bfile {bfile_in} --cluster --K {k} --out {os.path.join(workdir, 'plink')}"

# subprocess.run(cmd, shell=True)

# %%
from collections import Counter

import pandas as pd
from matplotlib import pyplot as plt

# %%
path_clust = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/plink.cluster1"
df_clust = pd.read_csv(path_clust, sep="\t", names = ["Cluster_name", "Sample_IDs"])

n_samples_per_cluster = df_clust["Sample_IDs"].apply(lambda val: len(val.split()))

count_n_samples_per_cluster = Counter(n_samples_per_cluster)

counts = sorted(count_n_samples_per_cluster.keys())

plt.figure(figsize=(10, 4))

plt.bar(
    list(range(len(counts))),
    list(map(count_n_samples_per_cluster.__getitem__, counts)),
    color = "gray"
)
plt.xticks(
    list(range(len(counts))),
    counts,
    fontsize = 10,
    rotation = -90,
    ha = "center"
)
plt.ylabel("Number of Clusters", fontsize=15)
plt.xlabel("Cluster Size (N)", fontsize=15)
plt.margins(x=0)
plt.show()
plt.close()

# %%
samples_per_cluster = df_clust["Sample_IDs"].apply(lambda val: list(map(lambda x: x.split("_")[0],val.split())))
samples_per_cluster


# %%
clust_path = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/cluster.cluster2"

df_ = pd.read_csv(
    clust_path,
    delim_whitespace=True,
    names=["FID", "IID", "Cluster"]
)

# %
from scipy.spatial.distance import cdist

path_pca = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.200kb_0.5.pca.eigenvec")
df_pca = pd.read_csv(path_pca, sep="\t")
n_pcs_use = 5
pc_cols = [c for c in df_pca.columns if c.startswith("PC")][:n_pcs_use]
# === PARSE SAMPLE LISTS ===
def parse_samples(sample_str):
    """Extract IID part from 'FID_IID' space-separated strings."""
    return [s.split("_")[0] for s in sample_str.split() if "_" in s]

df_clust["SampleList"] = df_clust["Sample_IDs"].apply(parse_samples)
df_clust_filt = df_clust[df_clust["SampleList"].apply(len) > 2]
df_clust_filt["SampleLen"] = df_clust_filt["SampleList"].apply(lambda x: len(x))

# === COMPUTE CENTRAL SAMPLE PER CLUSTER ===
central_rows = []

for _, row in df_clust_filt.iterrows():
    cluster_name = row["Cluster_name"]
    samples = row["SampleList"]
    sub = df_pca[df_pca["IID"].isin(samples)].copy()

    if len(sub) == 0:
        print(f"No matching PCA samples for {cluster_name}")
        continue

    # Compute centroid (robust to outliers → median can be used)
    centroid = sub[pc_cols].median().values.reshape(1, -1)

    # Distance to centroid
    dists = cdist(sub[pc_cols], centroid)
    sub["dist_to_centroid"] = dists

    # Find most central sample
    central_sample = sub.loc[sub["dist_to_centroid"].idxmin()]
    central_rows.append({
        "Cluster_name": cluster_name,
        "FID": central_sample["#FID"],
        "IID": central_sample["IID"],
        "Dist_to_centroid": central_sample["dist_to_centroid"]
    })

# === OUTPUT ===
df_central = pd.DataFrame(central_rows)

# %%
list_central_samples = df_central["FID"].to_list()
path_central = f"{workdir}/Koreans_central_sampled_500.list"
with open(path_central, mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in list_central_samples:
        fw.write("\t".join([sample, sample]) + "\n")
        
# %%
path_sample_info_1KGP = "/BiO/Research/Korea10KGenome/Resources/External_Genome_Data/1KGP/1KGP_30x_GRCh38/20130606_g1k_3202_samples_ped_population.txt"
df_sample_info_1KGP = pd.read_csv(path_sample_info_1KGP, delim_whitespace=True)[["SampleID", "Population", "Superpopulation"]]
list_samples_1KGP = df_sample_info_1KGP["SampleID"].to_list()

list_central_samples_plus_1KGP = list_central_samples + list_samples_1KGP
path_central_plus_1KGP = f"{workdir}/Koreans_central_sampled_500+1KGP.list"
with open(path_central_plus_1KGP, mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in list_central_samples_plus_1KGP:
        fw.write("\t".join([sample, sample]) + "\n")
# %%
