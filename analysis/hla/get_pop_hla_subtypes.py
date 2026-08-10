# %%
import json

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# %%
# -----------------------------
# 1. Load HLA table
# -----------------------------
hla_path = "/BiO/Access/kyungwhan1998/genome/hla/Resources/Data/1000GplusKOR_HLA.tsv"
hla_df = pd.read_csv(hla_path, sep="\t")

# -----------------------------
# 2. Load population coordinates
# -----------------------------
worldmap_coord_pop_json = "/BiO/Access/kyungwhan1998/genome/yhaplo/output/Korea10K/1000GplusKOR.pop.loc.json"
with open(worldmap_coord_pop_json, "r") as f:
    pop_coords = json.load(f)

# -----------------------------
# 3. Compute allele fractions per population
# -----------------------------
def compute_allele_fractions(hla_df, locus):
    """
    Returns a dict: {pop: {"fractions": [...], "coords": (lat, lon), "alleles": [...]}}
    """
    populations = {}
    all_alleles = sorted(set(hla_df[f"{locus}_1"].dropna().tolist() + hla_df[f"{locus}_2"].dropna().tolist()))

    for pop, pop_df in hla_df.groupby("Population"):
        alleles = pop_df[f"{locus}_1"].tolist() + pop_df[f"{locus}_2"].tolist()
        counts = [alleles.count(a) for a in all_alleles]
        total = sum(counts)
        fractions = [c / total if total > 0 else 0 for c in counts]
        populations[pop] = {
            "fractions": fractions,
            "coords": pop_coords.get(pop, (0, 0)),
            "alleles": all_alleles
        }
    return populations

# -----------------------------
# 4. Plot HLA pie chart world map
# -----------------------------
def plot_hla_map(populations, locus):
    cmap = plt.get_cmap("tab20")
    haplotypes = populations[list(populations.keys())[0]]["alleles"]

    # Assign fixed colors based on allele index
    dict_colors = {hap: cmap(i % 20 / 20) for i, hap in enumerate(haplotypes)}

    # Initialize map
    fig = plt.figure(figsize=(18, 9))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines(linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.LAND, facecolor="honeydew")
    ax.add_feature(cfeature.OCEAN, facecolor="whitesmoke")

    # Plot pie charts for each population
    for pop, info in populations.items():
        lat, lon = info["coords"]
        fractions = info["fractions"]

        # Enlarged inset pie
        ax_inset = inset_axes(
            ax,
            width=0.8,
            height=0.8,
            loc="center",
            bbox_to_anchor=(lon, lat),
            bbox_transform=ccrs.PlateCarree()._as_mpl_transform(ax),
            borderpad=0
        )
        ax_inset.pie(fractions, colors=[dict_colors[a] for a in haplotypes])
        ax_inset.set_aspect("equal")
        ax_inset.axis("off")

        # Population label
        ax.text(lon, lat + 7, pop,
                transform=ccrs.PlateCarree(),
                ha="center", va="bottom",
                fontsize=10, fontweight="bold")

    # Legend
    legend_handles = [
        mpatches.Patch(color=dict_colors[hap], label=hap) for hap in haplotypes
    ]
    legend_ax = fig.add_axes([0.05, 0.05, 0.9, 0.08])
    legend_ax.axis("off")

    legend = legend_ax.legend(
        handles=legend_handles,
        title=f"HLA-{locus} Alleles",
        loc="center",
        ncol=10,
        fontsize=12,
        title_fontsize=15,
        frameon=True,
    )

    # Make legend background opaque
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_edgecolor("black")

    # Title
    fig.suptitle(f"Global Distribution of HLA-{locus} Alleles", fontsize=16, y=0.97)

    plt.tight_layout()
    plt.show()

# -----------------------------
# 5. Run for each HLA locus
# -----------------------------
loci = ["A", "B", "C", "DRB1", "DQB1"]
for locus in loci:
    populations = compute_allele_fractions(hla_df, locus)
    plot_hla_map(populations, locus)

# %%
import json

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# -----------------------------
# 1. Load HLA table
# -----------------------------
hla_path = "/BiO/Access/kyungwhan1998/genome/hla/Resources/Data/1000GplusKOR_HLA.tsv"
hla_df = pd.read_csv(hla_path, sep="\t")

# -----------------------------
# 2. Load population coordinates
# -----------------------------
worldmap_coord_pop_json = "/BiO/Access/kyungwhan1998/genome/yhaplo/output/Korea10K/1000GplusKOR.pop.loc.json"
with open(worldmap_coord_pop_json, "r") as f:
    pop_coords = json.load(f)

# -----------------------------
# 3. Define superpopulation mapping
# -----------------------------
superpop_map = {
    "CEU": "EUR", "TSI": "EUR", "FIN": "EUR", "GBR": "EUR", "IBS": "EUR",
    "YRI": "AFR", "LWK": "AFR", "GWD": "AFR", "MSL": "AFR", "ESN": "AFR",
    "CHB": "EAS", "JPT": "EAS", "CHS": "EAS", "CDX": "EAS", "KHV": "EAS",
    "PEL": "AMR", "MXL": "AMR", "CLM": "AMR", "PUR": "AMR",
    "KOR": "EAS", "JPN": "EAS",
    # Add others as needed
}

# Superpopulation coordinates (map placement)
superpop_coords = {
    "EUR": (54, 15),
    "AFR": (0, 20),
    "EAS": (35, 105),
    "AMR": (10, -80),
    "SAS": (23, 80),
    "OCE": (-20, 135),
}

# -----------------------------
# 4. Compute allele fractions per superpopulation (all alleles)
# -----------------------------
def compute_allele_fractions_superpop_all(hla_df, locus):
    superpopulations = {}
    all_alleles_global = []

    for superpop, group_df in hla_df.groupby(hla_df["Population"].map(superpop_map)):
        if superpop is None:
            continue

        # All alleles in this superpopulation
        alleles = group_df[f"{locus}_1"].dropna().tolist() + group_df[f"{locus}_2"].dropna().tolist()
        if not alleles:
            continue

        # Count all alleles
        allele_counts = pd.Series(alleles).value_counts()
        allele_list = allele_counts.index.tolist()
        counts = allele_counts.values.tolist()

        # Fractions
        total = sum(counts)
        fractions = [c / total for c in counts]

        superpopulations[superpop] = {
            "fractions": fractions,
            "coords": superpop_coords.get(superpop, (0, 0)),
            "alleles": allele_list
        }
        all_alleles_global.extend(allele_list)

    unique_alleles = list(pd.unique(all_alleles_global))
    return superpopulations, unique_alleles

# -----------------------------
# 5. Plot HLA pie chart world map
# -----------------------------
def plot_hla_map(superpopulations, locus, global_alleles):
    cmap = plt.get_cmap("tab20")
    dict_colors = {allele: cmap(i % 20 / 20) for i, allele in enumerate(global_alleles)}

    fig = plt.figure(figsize=(18, 9))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines(linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.LAND, facecolor="honeydew")
    ax.add_feature(cfeature.OCEAN, facecolor="whitesmoke")

    for sp, info in superpopulations.items():
        lat, lon = info["coords"]
        fractions = info["fractions"]
        alleles_local = info["alleles"]

        ax_inset = inset_axes(
            ax,
            width=1.5,
            height=1.5,
            loc="center",
            bbox_to_anchor=(lon, lat),
            bbox_transform=ccrs.PlateCarree()._as_mpl_transform(ax),
            borderpad=0
        )
        ax_inset.pie(fractions, colors=[dict_colors[a] for a in alleles_local])
        ax_inset.set_aspect("equal")
        ax_inset.axis("off")

        # Title above pie
        ax_inset.set_title(sp, fontsize=16, fontweight="bold", pad=5)

    # Legend
    legend_handles = [mpatches.Patch(color=dict_colors[hap], label=hap) for hap in global_alleles]
    legend_ax = fig.add_axes([0.05, 0.05, 0.9, 0.08])
    legend_ax.axis("off")
    legend = legend_ax.legend(
        handles=legend_handles,
        title=f"HLA-{locus} Alleles (All types per Superpopulation)",
        loc="center",
        ncol=10,
        fontsize=10,
        title_fontsize=12,
        frameon=True,
    )
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_edgecolor("black")

    fig.suptitle(f"Global Distribution of HLA-{locus} Alleles (Superpopulation)", fontsize=16, y=0.97)
    plt.tight_layout()
    plt.show()

# -----------------------------
# 6. Run for each HLA locus
# -----------------------------
loci = ["A", "B", "C", "DRB1", "DQB1"]
for locus in loci:
    superpopulations, global_alleles = compute_allele_fractions_superpop_all(hla_df, locus)
    plot_hla_map(superpopulations, locus, global_alleles)

# %%
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

plt.rcParams["font.size"] = 20

# -----------------------------
# 1. Load HLA data
# -----------------------------
hla_path = "/BiO/Access/kyungwhan1998/genome/hla/Resources/Data/1000GplusKOR_HLA.tsv"
hla_df = pd.read_csv(hla_path, sep="\t")

# -----------------------------
# 2. Superpopulation mapping and color coding
# -----------------------------
superpop_map = hla_df.set_index('Population')['Superpopulation'].to_dict()
superpop_colors = {
    "EUR": "royalblue",
    "AFR": "darkorange",
    "EAS": "forestgreen",
    "AMR": "crimson"}

# -----------------------------
# 3. Compute allele frequencies per population for top N alleles
# -----------------------------
def compute_population_top_freqs(hla_df, locus, top_n=10):
    populations = hla_df['Population'].unique()

    # Find top N alleles globally for this locus
    all_alleles = hla_df[f"{locus}_1"].dropna().tolist() + hla_df[f"{locus}_2"].dropna().tolist()
    allele_counts = pd.Series(all_alleles).value_counts()
    top_alleles = allele_counts.index[:top_n].tolist()

    # Frequency matrix: rows = populations, cols = top alleles
    freq_matrix = pd.DataFrame(0, index=populations, columns=top_alleles, dtype=float)

    for pop in populations:
        pop_df = hla_df[hla_df['Population'] == pop]
        alleles = pop_df[f"{locus}_1"].dropna().tolist() + pop_df[f"{locus}_2"].dropna().tolist()
        total = len(alleles)
        for allele in top_alleles:
            count = alleles.count(allele)
            freq_matrix.loc[pop, allele] = count / total  # relative frequency

    return freq_matrix

# -----------------------------
# 4. PCA and plotting per locus with superpopulation legend
# -----------------------------
loci = ["A", "B", "C", "DRB1", "DQB1"]
top_n = 10  # top alleles

for locus in loci:
    freq_matrix = compute_population_top_freqs(hla_df, locus, top_n=top_n)

    # Standardize
    scaler = StandardScaler()
    freq_scaled = scaler.fit_transform(freq_matrix)

    # PCA
    pca = PCA(n_components=10)
    pca_result = pca.fit_transform(freq_scaled)

    # DataFrame with all PCs
    pca_df = pd.DataFrame(pca_result, columns=[f"PC{i+1}" for i in range(10)], index=freq_matrix.index)

    # Plot PC1 vs PC2
    plt.figure(figsize=(8, 7))
    for pop in pca_df.index:
        color = superpop_colors.get(superpop_map.get(pop, ""), "gray")
        plt.scatter(pca_df.loc[pop, "PC1"], pca_df.loc[pop, "PC2"], s=80, color=color, edgecolor="k")
        if pop != "KOR":
            plt.text(pca_df.loc[pop, "PC1"] + 0.02, pca_df.loc[pop, "PC2"] + 0.02,
                    pop, fontsize=14, weight="bold", color="k")
        else:
            plt.text(pca_df.loc[pop, "PC1"] + 0.02, pca_df.loc[pop, "PC2"] + 0.02,
                    pop, fontsize=14, weight="bold", color="firebrick")

    # Legend for superpopulations
    legend_handles = [mpatches.Patch(color=color, label=sp) for sp, color in superpop_colors.items()]
    plt.legend(handles=legend_handles, loc="upper right", frameon=True, fontsize=14)

    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.xticks(np.arange(-5, 6, 1))
    plt.yticks(np.arange(-5, 6, 1))
    plt.title(f"HLA-{locus}: Top {top_n} Alleles")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# %%
