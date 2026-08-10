# %%
import os
import subprocess

# %%
workdir = "/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation"
script = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Scripts/split_chr_by_chunks_bcftools.py"
input_vcf = os.path.join(workdir, "chr{i}.omniChipOnly.imputed.4K.dose.vcf.gz")
region_file = "/BiO/Share/Tool/shapeit5/resources/chunks/b38/20cM/chunks_chr{i}.header.incl.txt"
outdir = os.path.join(workdir)
threads = 20

def run(script, input_vcf, region_file, outdir, threads):
    cmd = f"python {script} {input_vcf} {region_file} {outdir} {threads}"
    
    subprocess.run(cmd, shell=True)


for i in range(1, 23, 1):
    run(script, input_vcf.format(i=i), region_file.format(i=i), outdir, threads) 