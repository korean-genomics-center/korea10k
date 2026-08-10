# %%
import os
import subprocess

import pandas as pd

# %%
tool_plink2 = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
bfile_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
list_sample_remove = f"{workdir}/removePlusnonKoreans.list"
bfile_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
threads = 50
cmd = f"{tool_plink2} --bfile {bfile_in} --remove {list_sample_remove} --make-bed --out {bfile_out} --threads {threads}"

# subprocess.run(cmd, shell=True)

# %%
tool_plink2 = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
bfile_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
list_sample_keep = f"{workdir}/Koreans_include_eas_samples.list"
bfile_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.include_eas_samples_only.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
threads = 50
cmd = f"{tool_plink2} --bfile {bfile_in} --keep {list_sample_keep} --make-bed --out {bfile_out} --threads {threads}"

# subprocess.run(cmd, shell=True)

# %%
tool_plink2 = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
bfile_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
size = 100
list_sample_keep = f"{workdir}/Koreans_random_sampled_{size}+1KGP_EAS.list"
bfile_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_EAS.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
threads = 50
cmd = f"{tool_plink2} --bfile {bfile_in} --keep {list_sample_keep} --make-bed --out {bfile_out} --threads {threads}"

# subprocess.run(cmd, shell=True)

# %%
tool_plink2 = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
bfile_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
size = 500
# list_sample_keep = f"{workdir}/Koreans.list"
list_sample_keep = f"{workdir}/Koreans_central_sampled_{size}+1KGP.list"
bfile_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_1KGP.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
threads = 50
cmd = f"{tool_plink2} --bfile {bfile_in} --keep {list_sample_keep} --make-bed --out {bfile_out} --threads {threads}"

# subprocess.run(cmd, shell=True)

# %%
tool_plink2 = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
bfile_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
list_sample_keep = f"{workdir}/9011Koreans+50JPT+50CHB_samples.list"
bfile_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9011Koreans+50JPT+50CHB.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
threads = 50
cmd = f"{tool_plink2} --bfile {bfile_in} --keep {list_sample_keep} --make-bed --out {bfile_out} --threads {threads}"

# subprocess.run(cmd, shell=True)

# %%
tool_plink2 = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
bfile_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
list_sample_keep = f"{workdir}/9000Koreans+1KGPEAS_samples.list"
bfile_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
threads = 50
cmd = f"{tool_plink2} --bfile {bfile_in} --keep {list_sample_keep} --make-bed --out {bfile_out} --threads {threads}"

# subprocess.run(cmd, shell=True)

# %%
tool_plink2 = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k"
bfile_in = f"{workdir}/Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd"
list_sample_keep = f"{workdir}/9000Korean_samples.list"
bfile_out = f"{workdir}/9000Korean_samples.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd"
threads = 50
cmd = f"{tool_plink2} --bfile {bfile_in} --keep {list_sample_keep} --make-bed --out {bfile_out} --threads {threads}"

subprocess.run(cmd, shell=True)
