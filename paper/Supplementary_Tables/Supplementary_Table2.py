# %%
import os

import numpy as np
import pandas as pd

# %%
path_excel = "/BiO/Access/kyungwhan1998/genome/Korea10KGenome/Resources/MetaData/LifestyleQuestionnaire_2023/KU10K_LQ_Integrated_Ver2.03.xlsx"

df_lq = pd.read_excel(path_excel, engine="openpyxl")
list_samples_lq = list(df_lq["KU10K-id"])
list_samples_lq

# %%
import os

import pandas as pd

# %%
path_ori = "/BiO/Research/Korea10KGenome/Resources/MetaData/Sequencing/KOREA10K_DATA_TABLE.xlsx"
path_add = "/BiO/Research/Korea10KGenome/Resources/Experiment_Sheets/T7_Additional_Sequencing_Experiment_Metadata.txt"
path_ori_meta = "/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.VQSR.PASS/Plink_All/Step1_MakePlink.biallelic/chr21.biallelic.psam"

list_samples_exclude = ["KU10K-10433", "KU10K-10689", "KU10K-04846", "KU10K-10007"]

df_ori = pd.read_excel(path_ori, sheet_name="Multiomics", engine="openpyxl")
df_add = pd.read_csv(path_add, delim_whitespace=True)
df_ori_meta = pd.read_csv(path_ori_meta, delim_whitespace=True, header=None).rename(columns={1:"IID"})
list_samples_all = list(filter(lambda x: "U10K" in x, list(df_ori_meta["IID"])))
list_samples_include = list(set(list_samples_all).difference(list_samples_exclude))

# %%
list_questions = list(df_lq.columns)
list_questions

# %%
path_codebook = "/BiO/Access/kyungwhan1998/genome/Korea10KGenome/Resources/MetaData/LifestyleQuestionnaire_2023/KU10K_LQ_Integrated_Codebook_Ver2.0.xlsx"
df_codebook = pd.read_excel(path_codebook, engine="openpyxl")
df_codebook_filt = df_codebook[["Category_Eng", "Question_Tag", "Question_Text_Eng"]]
df_codebook_filt_dropdup = df_codebook_filt.drop_duplicates(subset=["Question_Tag"])
df_questions = df_codebook_filt_dropdup.reset_index(drop=True)
df_questions.to_excel("/BiO/Access/kyungwhan1998/genome/paper/Supplementary_Table4.xlsx")
# %%
