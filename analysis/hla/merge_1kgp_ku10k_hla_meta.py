# %%
import os

import numpy as np
import pandas as pd
from korea10k.config import PROJECT_DIR, WORK_DIR

# %%
kgp_hla_clean_path = f"{WORK_DIR}/genome/hla/Resources/Data/1KGP_1113samples_HLA.tsv"
kgp_meta_path = f"{PROJECT_DIR}/Resources/External_Genome_Data/1KGP/1KGP_30x_GRCh38/20130606_g1k_3202_samples_ped_population.txt"
kgp_hla_clean = pd.read_csv(kgp_hla_clean_path, sep="\t")

## Select the first HLA as representative
kgp_hla_clean_select_first = kgp_hla_clean.copy()
list_hla_clean_cols = list(kgp_hla_clean_select_first.columns)[1:]
for col in list_hla_clean_cols:
    kgp_hla_clean_select_first[col] = kgp_hla_clean_select_first[col].apply(lambda x: x.split("/")[0])

kgp_meta = pd.read_csv(kgp_meta_path, delim_whitespace=True)[["SampleID", "Population", "Superpopulation"]]
kgp_hla_clean_meta_added = pd.merge(kgp_hla_clean_select_first, kgp_meta, how="inner", on="SampleID")

# %%
ku10k_hla_clean = f"{WORK_DIR}/genome/hla/Resources/Data/Korea10K_HLA.tsv"
ku10k_hla_clean = pd.read_csv(ku10k_hla_clean, sep="\t").rename(columns={"sample_id": "SampleID"})
ku10k_hla_clean_meta_added = ku10k_hla_clean.copy()
ku10k_hla_clean_meta_added["Population"] = "KOR"
ku10k_hla_clean_meta_added["Superpopulation"] = "EAS"

# %%
list_col_kgp_hla_clean_meta_added = list(kgp_hla_clean_meta_added.columns)

ku10k_kgp_hla_clean_meta_added = pd.concat([ku10k_hla_clean_meta_added[list_col_kgp_hla_clean_meta_added], kgp_hla_clean_meta_added], axis=0)
ku10k_kgp_hla_clean_meta_added.to_csv(f"{WORK_DIR}/genome/hla/Resources/Data/1000GplusKOR_HLA.tsv", sep="\t", index=False)
# %%
