#%%
import os
import pickle

import numpy as np
from matplotlib import pyplot as plt


#%%
def load_pickle(path_pkl):
    with open(file=path_pkl, mode='rb') as f:
        file_load=pickle.load(f)
        
    return file_load
        
# %%
dir_variant = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/plot"
path_pkl_dict_frequency_type_to_mean_line = os.path.join(dir_variant, "dict_frequency_type_to_mean_line.pkl")
path_pkl_dict_frequency_type_to_saturation_points = os.path.join(dir_variant, "dict_frequency_type_to_saturation_points.pkl")

#%%
dict_frequency_type_to_mean_line = load_pickle(path_pkl_dict_frequency_type_to_mean_line)
dict_frequency_type_to_saturation_points = load_pickle(path_pkl_dict_frequency_type_to_saturation_points)

# %%
list_frequency_type = ["Total", "Singleton", "Doubleton", "Extremely Rare", "Very Rare", "Rare", "Common", "Very Common"]

dict_frequency_palette = {
    "Total": "#8a8a8aff",
    "Singleton": "#000000ff",
    "Doubleton": "#78004A",
    "Extremely Rare": "#cb181d",
    "Very Rare": "#ff7b00",
    "Rare": "#006E81",
    "Common": "#0088f7",
    "Very Common": "#1203ae"
}

#%%
plt.rcParams["font.size"] = 14

xlim_padd = 500

fig, ax = plt.subplots(2, 1, figsize = (10, 12), gridspec_kw={"height_ratios":[5, 2]}, sharex=True)
plt.subplots_adjust(hspace = 0.05)

highlight_cats = ["Doubleton", "Extremely Rare", "Very Rare"]

max_samples = 0
for frequency_cat in list_frequency_type:
    mean_line = dict_frequency_type_to_mean_line[frequency_cat]["mean"]
    list_samples = list(range(1, len(mean_line)+1))
    
    # Ultra-thin for total curve
    if frequency_cat.lower() in ["total", "all", "all variants"]:
        lw = 1.5
        zord = 2
        alpha = 0.8
    # Highlighted categories
    elif frequency_cat in highlight_cats:
        lw = 4.5
        zord = 7
        alpha = 1.0
    # Normal categories
    else:
        lw = 3
        zord = 5
        alpha = 1.0

    # plt.fill_between(list_samples, mean_line-ci, mean_line+ci, color = dict_frequency_palette[frequency_cat], alpha = 0.4, zorder = 3)
    ax[0].plot(list_samples, mean_line, color = dict_frequency_palette[frequency_cat], linewidth=lw, zorder = 5, label = frequency_cat)
    if max(list_samples) > max_samples:
        max_samples = max(list_samples)

leg = ax[0].legend(
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    title="AF Category",
    title_fontproperties={"weight": "bold", "size": plt.rcParams["font.size"]+2},
    frameon=False
)

for text in leg.get_texts():
    freq_label = text.get_text()
    if freq_label in dict_frequency_palette:
        text.set_color(dict_frequency_palette[freq_label])

ax[0].set_ylabel("Variant Count", fontsize=plt.rcParams["font.size"]+4)
ax[0].set_xlim(0-xlim_padd, max_samples+xlim_padd)
ax[0].set_xticks(np.linspace(0, 10000, 6))
ax[0].set_xticks(np.linspace(0, 10000, 21), minor = True)
ax[0].tick_params(axis = 'x', length = 5, width = 1)
ax[0].tick_params(axis = 'x', which = "minor", length = 2, width = 1)
ax[0].grid(axis="x", which="minor", linewidth=0.2)
ax[0].grid(axis="x", which="major", linewidth=0.5)

for ind, frequency_cat in enumerate(list_frequency_type[::-1]):
    col = dict_frequency_palette[frequency_cat]
    box = ax[1].boxplot(
        dict_frequency_type_to_saturation_points[frequency_cat],
        widths=0.8,
        labels=[frequency_cat],
        positions=[ind],
        flierprops=dict(marker='o', markeredgecolor="w", markeredgewidth=0.3, markerfacecolor=col, markersize=6, linestyle='none', linewidth=2),
        boxprops=dict(linewidth=2.5, color="k"),
        medianprops=dict(linewidth=2.0, color="k"),
        whiskerprops=dict(linewidth=1.8, color="k"),
        capprops=dict(linewidth=1.8, color="k"),
        vert=False,
        patch_artist=True
    )
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(box[element], color=dict_frequency_palette[frequency_cat])
    for patch in box['boxes']:
        patch.set(facecolor=col, linewidth=2, alpha=0.5)

ax[1].set_xlim(0 - xlim_padd, max_samples + xlim_padd)
ax[1].yaxis.tick_right()
ax[1].set_xlabel("Sample Count (with All Variants Discovered)", fontsize=plt.rcParams["font.size"]+4)
ax[1].set_xticks(np.linspace(0, 10000, 6))
ax[1].set_xticks(np.linspace(0, 10000, 21), minor=True)
ax[1].tick_params(axis='x', length=5, width=1)
ax[1].tick_params(axis='x', which="minor", length=2, width=1)
ax[1].grid(axis="x", which="minor", linewidth=0.2)
ax[1].grid(axis="x", which="major", linewidth=0.5)

for tick_label, frequency_cat in zip(ax[1].get_yticklabels(), list_frequency_type[::-1]):
    tick_label.set_color(dict_frequency_palette[frequency_cat])




