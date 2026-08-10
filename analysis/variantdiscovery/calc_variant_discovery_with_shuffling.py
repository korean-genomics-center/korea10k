# %%
import os
import pickle

import numpy as np
import pandas as pd
from scipy import stats

# %%
list_shuffles = list(range(0, 100))
list_frequency_type = ["Total", "Singleton", "Doubleton", "Ultra Rare", "Very Rare", "Rare", "Common", "Very Common"]

dict_frequency_palette = {
    "Total": "#8a8a8aff",
    "Singleton": "#000000ff",
    "Doubleton": "#78004A",
    "Ultra Rare": "#cb181d",
    "Very Rare": "#ff7b00",
    "Rare": "#006E81",
    "Common": "#0088f7",
    "Very Common": "#1203ae"
} 

path_summary_format = "/BiO/Research/Korea10KGenome/Results/VariantDiscoveryStatistics/Summarized_Discovery_Index.remove_nonKorean/Summarized_Discovered_Index.korea10k_samples.RemoveABHetOutlier2STD.mind_0.1.het_3std.kinship_3rd.remove_nonKorean.shuffle_{ind}.txt"

# %%
dict_frequency_type_to_dict_cumsum = dict()
for ind in list_shuffles:
    table_summary_single = pd.read_csv(path_summary_format.format(ind = ind), sep = '\t')
    
    for frequency_type in list_frequency_type:
        if dict_frequency_type_to_dict_cumsum.get(frequency_type) == None:
            dict_frequency_type_to_dict_cumsum[frequency_type] = list()

        table_summary_single_freq = table_summary_single[table_summary_single["Variant_Set"] == frequency_type]
        dict_frequency_type_to_dict_cumsum[frequency_type].append(
            dict(zip(table_summary_single_freq["Sample_Count"], table_summary_single_freq["Variant_Discovery_Count"]))
        )

# %%
def calculate_confidence_interval(data, confidence = 0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t.ppf((1+confidence) / 2., n-1)
    return h
# %%
def get_average_line_from_shuffling(list_dict_cumsum):
    max_samplecount = 0
    for dict_cumsum in list_dict_cumsum:
        if max(dict_cumsum.keys()) > max_samplecount:
            max_samplecount = max(dict_cumsum.keys())

    array_cumsum_by_sample = None
    for dict_cumsum in list_dict_cumsum:
        max_samplecount_this_shuffle = max(dict_cumsum.keys())
        cumsum_singleshuffle = list(map(lambda samplecount: dict_cumsum[min(samplecount, max_samplecount_this_shuffle)], range(1, max_samplecount+1)))
        if isinstance(array_cumsum_by_sample, type(None)):
            array_cumsum_by_sample = np.array([cumsum_singleshuffle])
        else:
            array_cumsum_by_sample = np.append(array_cumsum_by_sample, np.array([cumsum_singleshuffle]), axis = 0)
    
    average_line = np.mean(array_cumsum_by_sample, axis = 0)
    std = np.std(array_cumsum_by_sample, axis = 0)
    ci = np.apply_along_axis(calculate_confidence_interval, axis = 0, arr = array_cumsum_by_sample)
    return average_line, std, ci

def get_saturation_point_from_shuffling(list_dict_cumsum):
    list_saturation_count = list()
    for dict_cumsum in list_dict_cumsum:
        list_saturation_count.append(max(dict_cumsum.keys()))
    return list_saturation_count

# %%
dict_frequency_type_to_mean_line = dict()
dict_frequency_type_to_saturation_points = dict()

for frequency_cat in list_frequency_type:
    mean_line, std, ci = get_average_line_from_shuffling(dict_frequency_type_to_dict_cumsum[frequency_cat])
    dict_frequency_type_to_mean_line[frequency_cat] = {
        "mean":mean_line,
        "std": std,
        "ci":ci
    }
    dict_frequency_type_to_saturation_points[frequency_cat] = get_saturation_point_from_shuffling(dict_frequency_type_to_dict_cumsum[frequency_cat])

# %%
def dump_pickle(file_dump, dict_dump):
    with open(file_dump, 'wb') as fw:
        pickle.dump(dict_dump, fw, protocol=pickle.HIGHEST_PROTOCOL)


dir_freqtype = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/plot"
file_frequency_type_to_mean_line = os.path.join(dir_freqtype, "dict_frequency_type_to_mean_line.pkl")
file_frequency_type_to_saturation_points = os.path.join(dir_freqtype, "dict_frequency_type_to_saturation_points.pkl")

dump_pickle(file_frequency_type_to_mean_line, dict_frequency_type_to_mean_line)
dump_pickle(file_frequency_type_to_saturation_points, dict_frequency_type_to_saturation_points)


# %%
for frequency_cat in list_frequency_type:
    print(frequency_cat, np.min(dict_frequency_type_to_saturation_points[frequency_cat]), np.median(dict_frequency_type_to_saturation_points[frequency_cat]), np.max(dict_frequency_type_to_saturation_points[frequency_cat]))

#%%
file_saturation_points_frequenct_type = os.path.join(dir_freqtype, "saturation_points_frequency_type.txt")
with open(file_saturation_points_frequenct_type, mode="w") as fw:
    header = "\t".join(["AF Category", "Number of variants for saturation", "95%CI"])+"\n"
    fw.write(header)
    for frequency_cat in list_frequency_type:
        ci = calculate_confidence_interval(data = list(dict_frequency_type_to_saturation_points[frequency_cat]))
        mean = np.mean(dict_frequency_type_to_saturation_points[frequency_cat])
        content = "\t".join([str(frequency_cat), str(mean), str(f"{mean-ci:.2f}"), str(f"{mean+ci:.2f}")])+"\n"
        fw.write(content)

# %%
