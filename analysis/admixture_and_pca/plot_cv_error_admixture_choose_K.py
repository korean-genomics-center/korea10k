# %%
import glob
import os

# %%
dir_admixture = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ukbb"
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
plt.scatter(dict_k_cv_error_sorted.keys(), dict_k_cv_error_sorted.values(), color="grey", zorder=2)
plt.xticks(list(range(1, len(list_log_files)+1, 1)))
plt.xlabel("K (number of ancestral population)", fontsize=15)
plt.ylabel("Cross-validation Error", fontsize=15)
plt.grid(axis="both")
plt.show()
plt.close()

# %%
