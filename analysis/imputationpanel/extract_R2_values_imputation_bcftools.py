# %%
import os
import subprocess

from joblib import Parallel, delayed


# %%
def extract_R2_values(workdir, bcftools, i, panel):
    vcf = f"{workdir}/chr{i}.recal.selected.reheadered.hg19.liftoverInfiniumOmni2-5-8v1-5_A1_chip_sites.hg38.liftover.imputed.{panel}.vcf.bgz"
    cmd = f"{bcftools} query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ER2\n' {vcf} | grep -v '\.$' > {workdir}/chr{i}_imputation_performance_{panel}.txt"
    subprocess.run(cmd, shell=True)

# %%
bcftools = "/BiO/Access/kyungwhan1998/miniconda3/envs/shapeit/bin/bcftools"
workdir = "/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation"
panel = "1K"
with Parallel(n_jobs = 23) as parallel:
    parallel(delayed(extract_R2_values)(workdir, bcftools, i, panel) for i in range(1, 23, 1))
