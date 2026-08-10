# %%
import os
import subprocess

from joblib import Parallel, delayed
from korea10k.config import PROJECT_DIR, STORE_DIR, WORK_DIR

# %%
dir_vcf = f"{STORE_DIR}/Jellyfish/KOGIC-KU10K-Genome-2019-01/Results/Korea4K.PhasedVCF.OnlyBiallelic"
dir_fa = f"{PROJECT_DIR}/Resources/Reference/chromosome"
script = f"{WORK_DIR}/genome/shapeit/Resources/Scripts/fix_vcf_missing_contigs.py"
input_vcf = os.path.join(dir_vcf, "chr{i}.recal.forPhasing.nonOverlap.phased.corrected.vcf")
fai = os.path.join(dir_fa, "hg38.fa.fai")
output_vcf = os.path.join(dir_vcf, "chr{i}.recal.forPhasing.nonOverlap.phased.corrected.headerFixed.vcf")

def run_fix_vcf(script, input_vcf, fai, output_vcf, i):
    cmd = f"python {script} {input_vcf.format(i=i)} {fai} {output_vcf.format(i=i)}"

    subprocess.run(cmd, shell=True)

with Parallel(n_jobs=23) as parallel:
    parallel(delayed(run_fix_vcf)(script, input_vcf, fai, output_vcf, i) for i in range(1, 23, 1))

