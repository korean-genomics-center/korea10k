# %%
import os

import pandas as pd

K = 9

workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Results"
qfile = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.500_random_Koreans_plus_1KGP.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.500kb_0.2.{K}.Q")
indfile = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.500_random_Koreans_plus_1KGP.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.500kb_0.2.fam")
ind2popfile = os.path.join(workdir, "ind2pop.txt")
pop_order_file = os.path.join(workdir, "pop_order.txt")

# Load sample IDs
samples = pd.read_csv(indfile, sep=r"\s+", header=None)[1].astype(str).tolist()

# Load Q matrix
Q = pd.read_csv(qfile, sep=r"\s+", header=None)
Q.columns = [f"Cluster{i+1}" for i in range(Q.shape[1])]
Q.index = samples

# Load population info
ind2pop = pd.read_csv(ind2popfile, sep="\t", header=None).rename(columns={0: "Population"})
ind2pop.index = samples

# Combine Q matrix with population labels
Qpop = pd.concat([Q, ind2pop], axis=1)

# Population order
pop_order = [line.strip() for line in open(pop_order_file) if line.strip() and not line.startswith("#")]
pop_order = [p for p in pop_order if p in Qpop["Population"].unique()]
Qpop["Population"] = pd.Categorical(Qpop["Population"], categories=pop_order, ordered=True)
Qpop = Qpop.sort_values(["Population"])

# Sort individuals within each population by max cluster proportion
sorted_Qpop = []
for pop in pop_order:
    subset = Qpop[Qpop["Population"] == pop].copy()
    if subset.empty:
        continue
    subset["max_cluster"] = subset[[f"Cluster{i+1}" for i in range(K)]].idxmax(axis=1)
    subset["max_value"] = subset[[f"Cluster{i+1}" for i in range(K)]].max(axis=1)
    
    # Reverse: sort descending by max_cluster, ascending by max_value
    subset = subset.sort_values(["max_cluster", "max_value"], ascending=[True, True])
    sorted_Qpop.append(subset.drop(columns=["max_cluster", "max_value"]))

Qpop_sorted = pd.concat(sorted_Qpop, axis=0)

# Ensure all cluster columns are float
for col in Qpop_sorted.columns:
    try:
        Qpop_sorted[col] = Qpop_sorted[col].astype(float)
    except:
        print(f"column '{col}' does not contain float values!")

# Save to file without scientific notation
path_output = os.path.join(workdir, f"admixture_plot_input_K{K}.txt")
Qpop_sorted.to_csv(path_output, sep="\t", float_format="%.6f", index=True)