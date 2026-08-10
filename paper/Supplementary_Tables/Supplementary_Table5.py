# %%
import os

import pandas as pd

# %%
path_ori = "/BiO/Research/Korea10KGenome/Resources/MetaData/Sequencing/KOREA10K_DATA_TABLE.xlsx"
path_ori_meta = "/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.VQSR.PASS/Plink_All/Step1_MakePlink.biallelic/chr21.biallelic.psam"
path_new_depth = "/BiO/Access/kyungwhan1998/genome/depthCoverage/KU10K_10243_base_mapped.txt"
path_out = "/BiO/Access/kyungwhan1998/genome/depthCoverage/Korea10K_10239_sequencing_depth.txt"

human_genome_size = 3298912062
list_samples_exclude = ["KU10K-10433", "KU10K-10689", "KU10K-04846", "KU10K-10007"]

df_ori = pd.read_excel(path_ori, sheet_name="Multiomics", engine="openpyxl")
df_ori_meta = pd.read_csv(path_ori_meta, delim_whitespace=True).rename(columns={1:"IID"})
list_samples_all = list(filter(lambda x: "U10K" in x, list(df_ori_meta["IID"])))
list_samples_include = list(set(list_samples_all).difference(list_samples_exclude))
df_new_depth = pd.read_csv(path_new_depth, delim_whitespace=True, header=None)
df_new_depth.columns=["SampleID", "Production_Amount_GB"]
df_new_depth["SampleID"] = df_new_depth["SampleID"].apply(lambda x: str(x).split(".")[0])
df_new_depth = df_new_depth[df_new_depth["SampleID"].isin(list_samples_include)]
df_new_depth["Average_Depth_x"] = df_new_depth["Production_Amount_GB"].apply(lambda x: int(x)/int(human_genome_size))
platforms = ["HiSeq X", "MGI_T7", "NovaSeq"]
df_new_depth["Rounded_Depth"] = df_new_depth["Average_Depth_x"].round().astype(int)

# %%
df_ori_filt = df_ori[["KU10K-ID", "WGS Platform"]].rename(columns={"KU10K-ID": "SampleID"})
df_ori_filt_new_depth = pd.merge(df_new_depth, df_ori_filt, how="inner", on="SampleID")
df_ori_filt_new_depth

# %%
path_rename = "/BiO/Access/kyungwhan1998/genome/paper/Korea10K.external.id.random.seed.434770919382672563240420251116.list"
df_rename = pd.read_csv(path_rename, delim_whitespace=True, header=None)
dict_rename = dict(zip(df_rename[0], df_rename[1]))

df_ori_filt_new_depth["SampleID"] = df_ori_filt_new_depth["SampleID"].apply(lambda x: dict_rename.get(str(x), "NA"))

df_ori_filt_new_depth.to_excel("/BiO/Access/kyungwhan1998/genome/paper/Supplementary_Table6.xlsx", index=False)