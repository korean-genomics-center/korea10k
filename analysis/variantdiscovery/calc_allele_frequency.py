#%%
import os
import pickle

import pandas as pd

#%%
table_maf = pd.read_csv("/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.RemoveABHetOutlier2STD.VQSR.PASS/Plink_Final/chromosome_merged.qc_filtered.kinship_filtered.nonKorean_filtered/Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.nonKorean.mac1.QC.acount", delim_whitespace = True)

list_chrname = list(map(lambda val: f"chr{val}", range(1, 23)))
table_varid_conv = pd.concat(
    list(map(lambda chrname: pd.read_csv(f"/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.RemoveABHetOutlier2STD.VQSR.PASS/Plink_Final/chromosome_merged.qc_filtered.kinship_filtered.nonKorean_filtered/dbSNP_Annotate/VarID_Conversion/ConversionTable.PlinkID_to_rsID.{chrname}.tsv", sep = '\t'), list_chrname))
)

#%%

table_maf["ALT_FREQS"] = table_maf["ALT_CTS"] / table_maf["OBS_CT"]

def classify_variant_by_vaf(vaf):
    
    if vaf <= 0.0005:
        return "Ultra Rare"
    if vaf <= 0.001:
        return "Very Rare"
    elif vaf <= 0.01:
        return "Rare"
    elif vaf <= 0.05:
        return "Common"
    else:
        return "Very Common"

table_maf_singleton = table_maf[table_maf["ALT_CTS"] == 1]
table_maf_doubleton = table_maf[table_maf["ALT_CTS"] == 2]
table_maf_over_doubleton = table_maf[table_maf["ALT_CTS"] > 2]

table_maf_singleton["Freq_Annotate"] = "Singleton"
table_maf_doubleton["Freq_Annotate"] = "Doubleton"
table_maf_over_doubleton["Freq_Annotate"] = table_maf_over_doubleton["ALT_FREQS"].apply(classify_variant_by_vaf)

table_maf_freq_annotated = pd.concat([table_maf_singleton, table_maf_doubleton, table_maf_over_doubleton])
#%%
dict_varid_to_rsid = dict(zip(table_varid_conv["VarID"], table_varid_conv["rsID"]))

table_maf_freq_annotated["rsID_Annotate"] = table_maf_freq_annotated["ID"].apply(dict_varid_to_rsid.__getitem__)
#%%
from collections import Counter

table_maf_freq_annotated_reported = table_maf_freq_annotated[table_maf_freq_annotated["rsID_Annotate"] != '.']
table_maf_freq_annotated_novel = table_maf_freq_annotated[table_maf_freq_annotated["rsID_Annotate"] == '.']

count_freq_reported = Counter(table_maf_freq_annotated_reported["Freq_Annotate"])
count_freq_novel = Counter(table_maf_freq_annotated_novel["Freq_Annotate"])
count_freq_all = Counter(table_maf_freq_annotated["Freq_Annotate"])

# %%
def get_vartype_from_refalt(ref, alt):
    if len(ref) == 1:
        if len(alt) == 1:
            return "SNP"
        else:
            return "INS"
    else:
        return "DEL"


table_maf_freq_annotated["VarType"] = list(map(lambda ref, alt: get_vartype_from_refalt(ref, alt), table_maf_freq_annotated["REF"], table_maf_freq_annotated["ALT"]))

table_maf_freq_annotated_reported["VarType"] = list(map(lambda ref, alt: get_vartype_from_refalt(ref, alt), table_maf_freq_annotated_reported["REF"], table_maf_freq_annotated_reported["ALT"]))

table_maf_freq_annotated_novel["VarType"] = list(map(lambda ref, alt: get_vartype_from_refalt(ref, alt), table_maf_freq_annotated_novel["REF"], table_maf_freq_annotated_novel["ALT"]))

# %%
dict_vartype_to_count_novel = dict()
for vartype in ["SNP", "INS", "DEL"]:
    dict_vartype_to_count_novel[vartype] = Counter(table_maf_freq_annotated_novel[table_maf_freq_annotated_novel["VarType"] == vartype]["Freq_Annotate"])


# %%
dir_vartype = "/BiO/Access/kyungwhan1998/genome/variant"
file_vartype = os.path.join(dir_vartype, 'dict_vartype_to_count_novel.pkl')

def dump_pickle(file_dump, dict_dump):
    with open(file_dump, 'wb') as fw:
        pickle.dump(dict_dump, fw, protocol=pickle.HIGHEST_PROTOCOL)

with open(file_vartype, 'wb') as fw:
    pickle.dump(dict_vartype_to_count_novel, fw, protocol=pickle.HIGHEST_PROTOCOL)