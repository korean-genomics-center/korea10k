# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# %%
K=7
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Results"
Qpop_sorted_file = os.path.join(workdir , f"admixture_plot_input_K{str(K)}.txt")
pop_order_file = os.path.join(workdir, "pop_order.txt")

Qpop_sorted = pd.read_csv(Qpop_sorted_file, sep="\t")
pop_order = [line.strip() for line in open(pop_order_file) if line.strip() and not line.startswith("#")]


plt.rcParams["font.size"] = 14
fig, ax = plt.subplots(figsize=(20, 1))
cluster_colors = sns.color_palette("Set2", K)
x = np.arange(len(Qpop_sorted))
bottom = np.zeros(len(Qpop_sorted))
for i, cluster in enumerate([f"Cluster{i+1}" for i in range(K)]):
    ax.bar(
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

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(5)

for b in pop_boundaries[:-1]:
    ax.axvline(b, color="black", lw=5)
ax.set_xticks(pop_labels)
ax.set_xticklabels(pop_labels_names, fontsize=plt.rcParams["font.size"]+4, fontweight="bold")
ax.set_ylabel("Ancestry proportion", fontsize=plt.rcParams["font.size"]+10)
ax.set_xlabel("Individuals (grouped by population)", fontsize=plt.rcParams["font.size"]+10)
ax.set_xmargin(0)
ax.set_ymargin(0)
# %%
