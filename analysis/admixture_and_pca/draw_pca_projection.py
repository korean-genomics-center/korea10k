# %%
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# %%
size = 500
path_projection = f"/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.include_eas_samples_only.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.LD_pruning_200kb_0.5.pruned.projection_{size}.sscore"
path_ref = f"/BiO/Access/kyungwhan1998/genome/pca/Resources/Data/ku10k_1KGP/Koreans_random_sampled_{size}.list"

set_korean_ref = set(pd.read_csv(path_ref, sep="\t")["#FID"].to_list())

# %%
proj_df = pd.read_csv(path_projection, sep='\t')

# Standardize column names
proj_df = proj_df.rename(columns=lambda x: x.strip())

proj_df["Projection"] = proj_df["#FID"].apply(lambda x: False if x in set_korean_ref else True)        
        
#%%

# Extract PCs
pcs = [f"PC{i}_AVG" for i in range(1, 11)]

# Rename reference PCs to match projected naming convention
proj_df = proj_df.rename(columns={f"PC{i}": f"PC{i}_AVG" for i in range(1, 11)})
proj_df = proj_df.rename(columns={"#FID": "SampleID"})
list_samples_10K = list(filter(lambda x: "10K" in x, proj_df["SampleID"].to_list()))
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
plot_df = pd.merge(proj_df, df_sample_info, how="inner", on="SampleID")

# %%
# ============================================================
# Plot PC1 vs PC2 with projected samples highlighted
# ============================================================

plt.rcParams["font.size"] = 18
plt.figure(figsize=(12, 12))

sns.scatterplot(
    data=plot_df[plot_df["Projection"] == False],
    x="PC1_AVG",
    y="PC2_AVG",
    color="indigo",
    edgecolor="black",
    alpha=1,
    s=50,
    label="Korean Reference (N=500)",
    zorder=3
)

# Plot background populations
sns.scatterplot(
    data=plot_df[plot_df["Projection"] == True],
    x="PC1_AVG",
    y="PC2_AVG",
    hue="Population",
    palette="Set1",
    alpha=0.5,
    s=100,
    edgecolor="black",
    zorder=1
)

plot_df_kor = plot_df[plot_df["Population"] == "KOR"]
plot_df_kor_proj = plot_df_kor[plot_df_kor["Projection"]]
list_top = plot_df_kor_proj.sort_values(by=["PC1_AVG"]).head(10)["SampleID"].to_list()

list_sample_annotate = """KU10K-05465
KU10K-08832
KU10K-05189
KU10K-08830
KU10K-03238
KU10K-03341
KU10K-03220
KU10K-03207
KU10K-01884
KU10K-03285
KU10K-05027
KU10K-09623
KU10K-09782""".split('\n')

annot_df = plot_df[plot_df["SampleID"].isin(list_sample_annotate)]

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
        row["PC1_AVG"],
        row["PC2_AVG"],
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

plt.title("Projected Samples onto Reference PCA Space", fontsize=24, fontweight="bold")
plt.tight_layout()
plt.show()
plt.close()


# %%
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# %%
size = 500
path_projection = f"/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.Koreans_only.preprocessed.prune_200kb_0.5.projection_include_eas_samples_only.sscore"

# %%
proj_df = pd.read_csv(path_projection, sep='\t')

# Standardize column names
proj_df = proj_df.rename(columns=lambda x: x.strip())      
        
#%%

# Extract PCs
pcs = [f"PC{i}_AVG" for i in range(1, 11)]

# Rename reference PCs to match projected naming convention
proj_df = proj_df.rename(columns={f"PC{i}": f"PC{i}_AVG" for i in range(1, 11)})
proj_df = proj_df.rename(columns={"#FID": "SampleID"})
list_samples_10K = list(filter(lambda x: "10K" in x, proj_df["SampleID"].to_list()))
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
plot_df = pd.merge(proj_df, df_sample_info, how="inner", on="SampleID")

# %%
# ============================================================
# Plot PC1 vs PC2 with projected samples highlighted
# ============================================================

plt.rcParams["font.size"] = 18
plt.figure(figsize=(12, 12))

# Plot background populations
sns.scatterplot(
    data=plot_df,
    x="PC1_AVG",
    y="PC10_AVG",
    hue="Population",
    palette="RdBu_r",
    alpha=0.3,
    s=100,
    edgecolor="black",
    zorder=1
)

# ============================================================
# Axis labels and legend
# ============================================================
plt.xlabel(f"PC1", fontsize=22, fontweight="bold")
plt.ylabel(f"PC2", fontsize=22, fontweight="bold")
plt.legend(frameon=False, fontsize=16)

plt.title("Projected Samples onto Reference PCA Space", fontsize=24, fontweight="bold")
plt.tight_layout()
plt.show()
plt.close()