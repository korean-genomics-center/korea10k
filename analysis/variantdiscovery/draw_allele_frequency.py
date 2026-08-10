# %%
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# %%
def load_pickle(path_pkl):
    with open(file=path_pkl, mode='rb') as f:
        file_load=pickle.load(f)
        
    return file_load
        
# %%
dir_variant = "/BiO/Access/kyungwhan1998/genome/variant/Data"
path_pkl_count_freq_reported = os.path.join(dir_variant, "count_freq_reported.pkl")
path_pkl_count_freq_novel = os.path.join(dir_variant, "count_freq_novel.pkl")
path_pkl_count_freq_all = os.path.join(dir_variant, "count_freq_all.pkl")
path_pkl_dict_vartype_to_color = os.path.join(dir_variant, "dict_vartype_to_count_novel.pkl")

#%%
count_freq_reported = load_pickle(path_pkl_count_freq_reported)
count_freq_novel = load_pickle(path_pkl_count_freq_novel)
count_freq_all = load_pickle(path_pkl_count_freq_all)
dict_vartype_to_count = load_pickle(path_pkl_dict_vartype_to_color)

# %%
# list_xticks = ["Singleton", "Doubleton", "Extremely Rare", "Very Rare", "Rare", "Common", "Very Common"]
# list_xtickslabel = ["Singleton", "Doubleton", "Extremely\nRare", "Very\nRare", "Rare", "Common", "Very\nCommon"]

# plt.rcParams["font.size"] = 15
# plt.rcParams["font.family"] = "Arial"

# # --- layout ---
# fig, ax1 = plt.subplots(figsize=(12, 6), facecolor='w')

# bottom = np.zeros(len(list_xticks))

# # --- Novel variants ---
# values_novel = np.array(list(map(count_freq_novel.__getitem__, list_xticks)))
# ax1.bar(
#     list_xticks, values_novel, bottom=bottom,
#     color="firebrick", label="Novel",
#     edgecolor='k', linewidth=3.5
# )

# for ind, freq in enumerate(list_xticks):
#     cnt_freq = count_freq_novel[freq]
#     cnt_freq_tot = count_freq_all[freq]
#     yloc = max(cnt_freq / 2, 0.02 * max(count_freq_all.values()))
#     ax1.annotate(f"{cnt_freq / cnt_freq_tot * 100:.2f}%", (ind, yloc),
#                  color='w', weight="bold", ha="center", va="center", fontsize=plt.rcParams["font.size"] - 0.5)
# bottom += values_novel

# # --- Reported variants ---
# values_reported = np.array(list(map(count_freq_reported.__getitem__, list_xticks)))
# ax1.bar(
#     list_xticks, values_reported, bottom=bottom,
#     color="gray", label="dbSNP (v157)",
#     edgecolor='k', linewidth=3.5
# )
# for ind, freq in enumerate(list_xticks):
#     cnt_freq = count_freq_reported[freq]
#     cnt_freq_tot = count_freq_all[freq]
#     yloc = cnt_freq / 2 + count_freq_novel[freq] + 0.02 * max(count_freq_all.values())
#     ax1.annotate(f"{cnt_freq / cnt_freq_tot * 100:.2f}%", (ind, yloc),
#                  color='w', weight="bold", ha="center", va="center", fontsize=plt.rcParams["font.size"] - 0.5)

# for i in range(len(list_xticks)):
#     cnt_total = count_freq_all[list_xticks[i]]
#     ax1.annotate(f"{cnt_total / 1e6:.1f}M",
#                  (i, cnt_total + 0.005 * max(count_freq_all.values())),
#                  ha="center", va="bottom", weight="bold")

# # --- Axes settings ---
# ax1.set_ylim(0, 3e7)
# ax1.set_xlim(-0.55, len(list_xticks)-0.55)
# ax1.set_xticks(range(len(list_xticks)))
# ax1.set_xticklabels(list_xtickslabel, weight="bold")
# ax1.tick_params(left=True)
# ax1.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.6, zorder=1)

# sns.despine(ax=ax1, left=False, bottom=False)
# for spine in ax1.spines.values():
#     spine.set_linewidth(2.5)
# ax1.tick_params(width=2.5)
# ax1.set_ylabel("Variant Count", weight="bold")

# total_novel = sum(count_freq_novel.values())
# total_reported = sum(count_freq_reported.values())
# ax_pie = inset_axes(ax1,
#                     width="50%", height="50%",
#                     loc='upper center',
#                     bbox_to_anchor=(0, 0, 1, 1),
#                     bbox_transform=ax1.transAxes,
#                     borderpad=0)

# sizes = [total_novel, total_reported]
# colors = ["firebrick", "gray"]

# ax_pie.pie(sizes, 
#            colors=colors, 
#            autopct="%1.1f%%", 
#            explode=(0, 0.2), 
#            pctdistance=0.7,
#            startangle=90,
#            counterclock=False,
#            wedgeprops={"linewidth":3, "edgecolor":"k"},
#            textprops={"weight":"bold", "fontsize":14, "color": "white", "va": "center", "ha": "center"})

# ax_pie.set_aspect("equal")

# total_all = total_novel + total_reported

# ax_pie.text(0.5, 0.92, f"{total_all/1e6:.1f}M", 
#             ha="center", va="bottom", weight="bold", fontsize=16, color="black",
#             transform=ax_pie.transAxes)

# # --- Legend ---
# handles, labels = ax1.get_legend_handles_labels()
# ax1.legend(handles[::-1], labels[::-1], frameon=False, prop={"weight": "bold"})

# plt.show()
# plt.close()


# %%
list_xticks = ["Singleton", "Doubleton", "Extremely Rare", "Very Rare", "Rare", "Common", "Very Common"]
list_xtickslabel = ["Singleton", "Doubleton", "Extremely\nRare", "Very\nRare", "Rare", "Common", "Very\nCommon"]

list_vartype = ["SNP", "DEL", "INS"]

dict_vartype_to_color = {
    "SNP":sns.color_palette("colorblind")[3],
    "DEL":sns.color_palette("colorblind")[2],
    "INS":sns.color_palette("colorblind")[0]
}

plt.rcParams["font.size"] = 17
plt.rcParams["font.family"] = "Arial"
plt.figure(figsize = (11, 7), facecolor = 'w')

bottom = np.array([0] * len(list_xticks))


for vartype in list_vartype:
    values_novel = np.array(list(map(dict_vartype_to_count[vartype].__getitem__, list_xticks)))
    plt.bar(
        list_xticks, values_novel, bottom = bottom,
        color = dict_vartype_to_color[vartype], label = vartype,
        edgecolor = 'k', linewidth = 1.5
    )

    bottom += values_novel

for i in range(len(list_xticks)):
    cnt_total = count_freq_novel[list_xticks[i]]
    plt.annotate(
        f"{cnt_total/1e3:,.1f}K",
        (i, cnt_total + 0.02e7),
        ha = "center", 
        va = "bottom",
        weight = "bold"
    )

plt.ylabel("Variant Count", weight = "bold")
plt.ylim(bottom = -0.05e7, top = 0.93e7)
plt.xlim(left = -0.5, right = len(list_xticks)-1 + 0.5)
plt.xticks(range(len(list_xtickslabel)), list_xtickslabel, weight = "bold")
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], loc = "upper right", bbox_to_anchor = (0.95, 0.8), prop = {"weight":"bold"})
plt.show()


#%%

list_xticks = ["Singleton", "Doubleton", "Extremely Rare", "Very Rare", "Rare", "Common", "Very Common"]
list_xtickslabel = ["Singleton", "Doubleton", "Extremely\nRare", "Very\nRare", "Rare", "Common", "Very\nCommon"]

list_vartype = ["SNP", "DEL", "INS"]

dict_vartype_to_color = {
    "SNP":sns.color_palette("colorblind")[3],
    "DEL":sns.color_palette("colorblind")[2],
    "INS":sns.color_palette("colorblind")[0]
}

plt.rcParams["font.size"] = 17
plt.rcParams["font.family"] = "Arial"
plt.figure(figsize = (11, 7), facecolor = 'w')

bottom = np.array([0] * len(list_xticks)).astype(np.float64)

for vartype in list_vartype:
    values_novel = np.array(list(map(dict_vartype_to_count[vartype].__getitem__, list_xticks))) / np.array(list(map(count_freq_novel.__getitem__, list_xticks)))
    
    plt.bar(
        list_xticks, values_novel, bottom = bottom,
        color = dict_vartype_to_color[vartype], label = vartype,
        edgecolor = 'k', linewidth = 1.5
    )

    for ind, freq in enumerate(list_xticks):
        cnt_freq = dict_vartype_to_count[vartype][freq]
        cnt_freq_tot = count_freq_novel[freq]
        
        yloc = values_novel[ind] / 2 + bottom[ind]
        plt.annotate(f"{cnt_freq / cnt_freq_tot * 100:.2f}%", (ind, yloc), color = 'w', weight="bold", ha = "center", va = "center", fontsize = plt.rcParams["font.size"] - 0.5)

    bottom += values_novel

plt.ylabel("Ratio of Variant", weight = "bold")
plt.ylim(bottom = -0.05, top = 1.05)
plt.xlim(left = -0.5, right = len(list_xticks)-1 + 0.5)
plt.xticks(range(len(list_xtickslabel)), list_xtickslabel, weight = "bold")
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], loc = "upper left", bbox_to_anchor = (1.01, 0.8), prop = {"weight":"bold"})
plt.show()

# %%
list_xticks = ["Singleton", "Doubleton", "Extremely Rare", "Very Rare", "Rare", "Common", "Very Common"]
list_xtickslabel = ["Singleton", "Doubleton", "Extremely\nRare", "Very\nRare", "Rare", "Common", "Very\nCommon"]
list_vartype = ["SNP", "DEL", "INS"]
dict_vartype_to_color = {
    "SNP": sns.color_palette("colorblind")[3],
    "DEL": sns.color_palette("colorblind")[2],
    "INS": sns.color_palette("colorblind")[0]
}

plt.rcParams["font.size"] = 17
plt.rcParams["font.family"] = "Arial"

# Create figure with 2 subplots side by side
fig, axes = plt.subplots(1, 2, figsize=(18, 7), facecolor='w')

# ----- Left plot: Variant counts -----
ax = axes[0]
bottom = np.zeros(len(list_xticks))

for vartype in list_vartype:
    values = np.array([dict_vartype_to_count[vartype][x] for x in list_xticks])
    ax.bar(list_xticks, values, bottom=bottom, color=dict_vartype_to_color[vartype],
           label=vartype, edgecolor='k', linewidth=1.5)
    bottom += values

for i, xtick in enumerate(list_xticks):
    ax.annotate(f"{count_freq_novel[xtick]/1e3:,.1f}K",
                (i, count_freq_novel[xtick] + 0.02e7),
                ha="center", va="bottom", weight="bold", fontsize=15)

ax.set_ylabel("Variant Count", weight="bold")
ax.set_ylim(-0.05e7, 0.93e7)
ax.set_xticks(range(len(list_xtickslabel)))
ax.set_xticklabels(list_xtickslabel, weight="bold", rotation=45, rotation_mode="anchor", ha="right")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc="upper right", bbox_to_anchor=(0.95, 0.8), prop={"weight":"bold"}, frameon=False)
ax.text(-0.1, 1.02, "A", transform=ax.transAxes, fontsize=25, weight="bold")

# ----- Right plot: Variant ratio -----
ax = axes[1]
bottom = np.zeros(len(list_xticks))

for vartype in list_vartype:
    values = np.array([dict_vartype_to_count[vartype][x] / count_freq_novel[x] for x in list_xticks])
    ax.bar(list_xticks, values, bottom=bottom, color=dict_vartype_to_color[vartype],
           label=vartype, edgecolor='k', linewidth=1.5)
    for ind, xtick in enumerate(list_xticks):
        yloc = values[ind] / 2 + bottom[ind]
        ax.annotate(f"{values[ind]*100:.2f}%", (ind, yloc), color='w', weight="bold",
                    ha="center", va="center", fontsize=plt.rcParams["font.size"] - 4)
    bottom += values

ax.set_ylabel("Ratio of Variant", weight="bold")
ax.set_ylim(-0.05, 1.05)
ax.set_xticks(range(len(list_xtickslabel)))
ax.set_xticklabels(list_xtickslabel, weight="bold", rotation=45, rotation_mode="anchor", ha="right")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc="upper left", bbox_to_anchor=(1.01, 0.6), prop={"weight":"bold"}, frameon=False)
ax.text(-0.15, 1.02, "B", transform=ax.transAxes, fontsize=25, weight="bold")

plt.tight_layout()
plt.show()