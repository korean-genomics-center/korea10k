# %%
import os

import numpy as np
import pandas as pd

# %%
path_1kgp = "/BiO/Research/Korea10KGenome/Resources/External_Genome_Data/1KGP/1KGP_30x_GRCh38/20130606_g1k_3202_samples_ped_population.txt"
df_1kgp = pd.read_csv(path_1kgp, delim_whitespace=True)
df_1kgp_EAS = df_1kgp[df_1kgp["Superpopulation"] == "EAS"]
list_samples_1kgp_eas = list(df_1kgp_EAS["SampleID"])

# %%
path_kor = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/Koreans.list"
df_kor = pd.read_csv(path_kor, sep="\t")
list_samples_kor = list(df_kor["IID"])

# %%
path_kgp_1kgp_random_eas = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/9000Koreans+1KGPEAS_samples.list"
list_samples_ku10k_1kgp_random_eas = list_samples_kor + list_samples_1kgp_eas
with open(path_kgp_1kgp_random_eas, mode="w") as fw:
    fw.write("\t".join(["#FID", "IID"]) + "\n")
    for sample in list_samples_ku10k_1kgp_random_eas:
        fw.write("\t".join([sample, sample]) + "\n")
