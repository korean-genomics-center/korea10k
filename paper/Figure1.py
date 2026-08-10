# %%
from collections import Counter

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import FancyArrowPatch, Patch, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# %%
path_cohort = "/BiO/Access/kyungwhan1998/genome/depthCoverage/10243sample_list.xlsx"
df_cohort = pd.read_excel(path_cohort)

# %% [QC]
path_excel = "/BiO/Research/Korea10KGenome/Resources/MetaData/Sequencing/KOREA10K_DATA_TABLE.xlsx"
df_excel = pd.read_excel(path_excel)
dict_rd_id_conv = dict(zip(df_excel["KU10K-ID"], df_excel["RD_ID"]))
dict_rd_id_conv = {k: str(int(v)) for k, v in dict_rd_id_conv.items() if str(v) != "nan"}
dict_kpgp_id_conv = dict(zip(df_excel["KU10K-ID"], df_excel["KPGP_ID"]))
dict_kpgp_id_conv = {k: str(v) for k, v in dict_kpgp_id_conv.items() if str(v) != "nan"}
dict_kpgp_id_conv_filt = {k: str(v) for k, v in dict_kpgp_id_conv.items() if not "U10K" in str(k)}

dict_4k_id_conv = dict()
dict_4k_id_conv.update(dict_rd_id_conv)
dict_4k_id_conv.update(dict_kpgp_id_conv_filt)

# %%
df_cohort["SampleID"] = df_cohort["ID"].apply(lambda x: dict_4k_id_conv.get(x, x))

path_vcf_4k = "/BiO/Store/KOGIC/Jellyfish/KOGIC-KU10K-Genome-2019-01/Results/JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller/Korea4K.4157Samples.VQSR.Filtered.Related.Rare.Diabetes.NonKorean.Inculdes.KOREFs.Filtered.AvgAB.1_0/chr14.recal.vcf"
list_sample_4k = list()
with open(path_vcf_4k, mode="r") as fr:
    for line in fr:
        if str(line).startswith("#CHROM"):
            record = line.rstrip("\n").split("\t")
            samples = record[9:]
            list_sample_4k.extend(samples)
            break

# %%
df_cohort_set_ind = df_cohort.set_index("SampleID")
df_cohort_4k = df_cohort_set_ind.loc[list_sample_4k]
 
# %%
list_samples_exclude = ["KU10K-10433", "KU10K-10689", "KU10K-04846", "KU10K-10007"]
dict_cohort = dict(zip(df_cohort["ID"], df_cohort["모집군"]))

dict_cohort_filt = dict()
for sample, cohort in dict_cohort.items():
    if sample not in list_samples_exclude:
        dict_cohort_filt[sample] = cohort

dict_cnt = dict(Counter(dict_cohort_filt.values()))

total_cnt = 0
diagnosed_cnt = 0
for category, count in dict_cnt.items():
    if category != "일반":
        diagnosed_cnt += count
    total_cnt += count

dict_cnt["질환군"] = diagnosed_cnt
dict_cnt["총합"] = total_cnt

df_pheno = pd.DataFrame.from_dict(dict_cnt, orient="index").reset_index()

df_pheno.columns = ["Category", "Count"]

# %% 
import json

path_translate_dict = "/BiO/Access/kyungwhan1998/genome/depthCoverage/convert_disease_category_kor_to_eng.json"
with open(path_translate_dict, mode="rb") as fr:
    translate_dict = json.load(fr)
df_pheno["Category_EN"] = df_pheno["Category"].map(translate_dict)

# %% Split summary and disease-specific data
summary_df = df_pheno[df_pheno["Category_EN"].isin(["Healthy", "Diagnosed"])]

disease_df = df_pheno.loc[1:15].copy().sort_values("Count", ascending=False)

# %% Assign major disease categories
def assign_group(cat):
    if cat in ["Myocardial infarction", "Angina"]:
        return "Cardiovascular disorder"
    elif cat in ["Depression", "Anxiety Disorder", "Suicide attempt", "Suicidal ideation", "Sleep disorder"]:
        return "Mental & mood disorder"
    elif "cancer" in cat.lower():
        return "Cancer"
    elif "diabetic" in cat.lower():
        return "Metabolic disorder"
    elif cat in ["Rare disorder", "Congenital hearing loss"]:
        return "Others"
    else:
        return "Others"

disease_df["Major_Group"] = disease_df["Category_EN"].apply(assign_group)

# %% Group summary
summary_groups = (
    disease_df.groupby("Major_Group")["Count"].sum()
    .reindex(["Cardiovascular disorder", "Mental & mood disorder", "Cancer", "Metabolic disorder", "Others"])
    .dropna()
    .reset_index()
)

total_count = summary_groups["Count"].sum()
summary_groups["Percent"] = summary_groups["Count"] / total_count * 100

# %% Color palette
group_colors = {
    "Cardiovascular disorder": "#B10101",
    "Mental & mood disorder": "#008FDC",
    "Cancer": "#0F9200",
    "Metabolic disorder": "#F57600",
    "Others": "#8D0088"
}
disease_df["Color"] = disease_df["Major_Group"].map(group_colors)

# %% 
plt.rcParams["font.size"] = 18
fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.5], width_ratios=[1,1], hspace=0.1, wspace=0)

# ----------------------------------------------------------------------
# Panel A: Recruitment Flowchart (upper left)
# ----------------------------------------------------------------------
axA = fig.add_subplot(gs[0, 0])
axA.set_xlim(0, 10)
axA.set_ylim(0, 5)
axA.axis("off")
axA.text(-0.9, 1.05, "A", transform=axA.transAxes,
         fontsize=plt.rcParams["font.size"]+10, fontweight='bold', va='bottom', ha='left')

# ----------------------------------------------------------------------
# Panel B: Healthy vs Diagnosed Pie (upper right)
# ----------------------------------------------------------------------
axB = fig.add_subplot(gs[0, 1])
wedges, texts, autotexts = axB.pie(
    summary_df["Count"], 
    labels=summary_df["Category_EN"], 
    autopct=lambda p: f"{int(round(p * summary_df['Count'].sum() / 100)):,}\n({p:.1f}%)", 
    startangle=90, 
    pctdistance=0.5, 
    colors=["grey", "firebrick"], 
    wedgeprops={"linewidth": 3, 
                "edgecolor": "k"}, 
    textprops={"weight": "bold", 
               "fontsize": plt.rcParams["font.size"], 
               "color": "white"} ) 

axB.set_title( f'Korea10K: {int(df_pheno.loc[df_pheno["Category_EN"] == "Total", "Count"].values[0]):,} samples', fontsize=plt.rcParams["font.size"] + 6, fontweight="bold", ha="center")

legend_handles = [Patch(facecolor=color, edgecolor='k', label=label)
                  for label, color in zip(summary_df["Category_EN"], ["grey", "firebrick"])]


axB.legend(handles=legend_handles,bbox_to_anchor=(1.0, 0.5), loc="center left", fontsize=plt.rcParams["font.size"], frameon=False)
           
axB.text(-0.6, 1.05, "B", transform=axB.transAxes,
         fontsize=plt.rcParams["font.size"]+10, fontweight='bold', va='bottom', ha='left')

# ----------------------------------------------------------------------
# Panel C: Disease Subgroup Bar Chart (bottom row full width)
# ----------------------------------------------------------------------
axC = fig.add_subplot(gs[1, :])
bars = axC.barh(
    disease_df["Category_EN"],
    disease_df["Count"],
    color=disease_df["Color"],
    edgecolor="k"
)
axC.set_xlabel("Number of Samples", fontsize=plt.rcParams["font.size"]+5)
axC.set_title(f"Disease Subgroups ({sum(disease_df['Count']):,} samples)", fontsize=plt.rcParams["font.size"]+5, fontweight="bold")
axC.invert_yaxis()

for i, v in enumerate(disease_df["Count"]):
    axC.text(v + 40, i, f"{v:,}", va="center", fontsize=plt.rcParams["font.size"]+3)

ax_inset = inset_axes(axC, width="35%", height="35%", loc="lower right", borderpad=2)
ax_inset.barh(
    summary_groups["Major_Group"],
    summary_groups["Count"],
    color=[group_colors[g] for g in summary_groups["Major_Group"]],
    edgecolor="k"
)
ax_inset.set_xticks(np.arange(0, 3200, 1000))
ax_inset.invert_yaxis()
ax_inset.set_title("Major Disease Groups", fontsize=plt.rcParams["font.size"]+4, fontweight="bold")

for i, (count, pct) in enumerate(zip(summary_groups["Count"], summary_groups["Percent"])):
    ax_inset.text(
        count + summary_groups["Count"].max() * 0.02,  # dynamic offset
        i,
        f"{count:,} ({pct:.1f}%)",                   # e.g. 2,340 (18.5%)
        va="center",
        fontsize=plt.rcParams["font.size"]
    )

sns.despine(ax=ax_inset, top=True, right=True)

axC.text(-0.45, 1.05, "C", transform=axC.transAxes,
         fontsize=plt.rcParams["font.size"]+10, fontweight='bold', va='bottom', ha='left')

axC.set_ymargin(0)
sns.despine(ax=axC, top=True, right=True)
plt.tight_layout()
plt.show()

# %%