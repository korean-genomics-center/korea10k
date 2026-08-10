# %%
import glob
import math
import os

# %%
dir_admixture = "/BiO/Access/kyungwhan1998/genome/admixture/Results"
pattern_log = "log*.out"
list_log_files = glob.glob(os.path.join(dir_admixture, pattern_log), recursive=True)

# %%
dict_k_cv_error = dict()
for log_file in list_log_files:
    with open(log_file, mode="r") as fr:
        for line in fr:
            if line.startswith("CV error"):
                record = line.rstrip("\n").split(":")
                K = record[0].split("=")[-1].replace(")", "")
                error = float(record[-1].strip())
                dict_k_cv_error[int(K)] =  error

dict_k_cv_error_sorted = dict(sorted(dict_k_cv_error.items(), key=lambda x: x[0]))

# %%
import matplotlib.pyplot as plt

plt.rcParams["font.size"] = 14
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))

error = list(map(lambda x: x * 100, dict_k_cv_error_sorted.values()))[:-2]
k = list(dict_k_cv_error_sorted.keys())[:-2]

# ---------------------------
# Left panel (A)
# ---------------------------
axes[0].scatter(k, error, color="grey", zorder=2)
axes[0].set_xticks(list(range(1, len(list_log_files) + 1, 1)))
axes[0].set_xlabel("K (number of ancestral population)", fontsize=20)
axes[0].set_xlim(0, 18)
axes[0].set_ylabel("Cross-validation Error", fontsize=20)
axes[0].grid(axis="both")

# Add panel label "A"
axes[0].text(
    -0.15, 1.02, "A", transform=axes[0].transAxes,
    fontsize=24, fontweight="bold", va="bottom", ha="left"
)

# ---------------------------
# Right panel (B)
# ---------------------------
axes[1].scatter(k, error, color="grey", zorder=2)
axes[1].set_xticks(list(range(1, len(list_log_files) + 1, 1)))
axes[1].set_xlabel("K (number of ancestral population)", fontsize=20)
axes[1].grid(axis="both")
axes[1].set_xlim(3, 15)
axes[1].set_ylim(25.97, 26.07)

# Add panel label "B"
axes[1].text(
    -0.15, 1.02, "B", transform=axes[1].transAxes,
    fontsize=24, fontweight="bold", va="bottom", ha="left"
)

# Circle K=7 and K=9 points
highlight_k = [7, 9]
for hk in highlight_k:
    if hk in k:
        idx = k.index(hk)
        axes[1].scatter(
            k[idx], error[idx],
            facecolors='none', edgecolors='firebrick',
            s=300, linewidths=3, zorder=3
        )

plt.tight_layout()
plt.show()
plt.close()

# %%
