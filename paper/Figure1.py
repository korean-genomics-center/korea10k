# %%
import os
from collections import Counter

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import FancyArrowPatch, Patch, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from korea10k.config import PROJECT_DIR, STORE_DIR, WORK_DIR

from Figure1A import X_BOX_LEFT, draw_panel_a

# %%
path_cohort = f"{WORK_DIR}/genome/depthCoverage/10243sample_list.xlsx"
df_cohort = pd.read_excel(path_cohort)

# %% [QC]
path_excel = f"{PROJECT_DIR}/Resources/MetaData/Sequencing/KOREA10K_DATA_TABLE.xlsx"
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

path_vcf_4k = f"{STORE_DIR}/Jellyfish/KOGIC-KU10K-Genome-2019-01/Results/JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller/Korea4K.4157Samples.VQSR.Filtered.Related.Rare.Diabetes.NonKorean.Inculdes.KOREFs.Filtered.AvgAB.1_0/chr14.recal.vcf"
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
# reindex (not .loc) so VCF samples missing from the cohort table do not abort the run;
# df_cohort_4k is not consumed by any panel below.
df_cohort_4k = df_cohort_set_ind.reindex(list_sample_4k)
 
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

path_translate_dict = f"{WORK_DIR}/genome/depthCoverage/convert_disease_category_kor_to_eng.json"
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
axA.axis("off")
# content is drawn after tight_layout(), once axA has its final size -- see below

# ----------------------------------------------------------------------
# Panel B: Healthy vs Diagnosed Pie (upper right)
# ----------------------------------------------------------------------
axB = fig.add_subplot(gs[0, 1])
wedges, texts, autotexts = axB.pie(
    summary_df["Count"],
    # slice labels are drawn in white and spill over Panel A; the legend below covers them
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
           
text_panel_B = axB.text(-0.6, 1.05, "B", transform=axB.transAxes,
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

text_panel_C = axC.text(-0.45, 1.05, "C", transform=axC.transAxes,
         fontsize=plt.rcParams["font.size"]+10, fontweight='bold', va='bottom', ha='left')

axC.set_ymargin(0)
sns.despine(ax=axC, top=True, right=True)
plt.tight_layout()

# ----------------------------------------------------------------------
# Panel A content (drawn last)
# ----------------------------------------------------------------------
# tight_layout sizes the left column around panel B's external legend, which leaves
# axA at ~30% of the row width. Position it explicitly instead, against the space
# panels B and C actually occupy once rendered:
#   left   flush with the start of panel C's y tick labels
#   right  MARGIN clear of everything panel B draws (its "B" letter reaches furthest)
#   bottom MARGIN clear of everything panel C draws (its "C" letter and title)
# The flowchart scales its font to whatever width it ends up with.
MARGIN_PANEL_A = 0.025  # figure fraction

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
to_fig = fig.transFigure.inverted()

pos_axB = axB.get_position()
bbox_axB = axB.get_tightbbox(renderer).transformed(to_fig)
bbox_axC = axC.get_tightbbox(renderer).transformed(to_fig)

# panel C's visual left edge is its longest tick label, not its spine, so measure it
# directly -- its tight bbox would instead report the far-flung "C" letter
x_ticklabels = min(t.get_window_extent(renderer).x0 for t in axC.get_yticklabels())
x_target = x_ticklabels / (fig.get_figwidth() * fig.dpi)

inset = X_BOX_LEFT / 10.0  # boxes sit this fraction of the axes width in from its edge
x_right = bbox_axB.x0 - MARGIN_PANEL_A
x_left = (x_target - inset * x_right) / (1 - inset)
y_bottom = bbox_axC.y1 + MARGIN_PANEL_A
y_top = pos_axB.y1

axA.set_position([x_left, y_bottom, x_right - x_left, y_top - y_bottom])
draw_panel_a(axA)

# the "A" letter is placed in figure coordinates rather than relative to axA, so it can
# share its left edge with "C" and its top edge with "B" (all three are ha=left/va=bottom)
x_letter = text_panel_C.get_transform().transform(text_panel_C.get_position())[0]
y_letter = text_panel_B.get_transform().transform(text_panel_B.get_position())[1]
x_letter, y_letter = fig.transFigure.inverted().transform((x_letter, y_letter))

fig.text(x_letter, y_letter, "A",
         fontsize=plt.rcParams["font.size"]+10, fontweight="bold", va="bottom", ha="left")

# %% Save
dir_figure = f"{PROJECT_DIR}/Analysis/Revision/Draw_Figure_ver20260929/Figures"
os.makedirs(dir_figure, exist_ok=True)

fig.savefig(os.path.join(dir_figure, "Figure1.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(dir_figure, "Figure1.pdf"), bbox_inches="tight")

plt.show()

# %%