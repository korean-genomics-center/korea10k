#!/usr/bin/env python3
import glob
import os
import subprocess

import pandas as pd

# ====================== CONFIGURATION ======================
workdir = "/BiO/Access/kyungwhan1998/genome/shapeit"

corr_script = os.path.join(
    workdir,
    "Resources/Scripts/calc_correlation_imputation_genotype.py"
)

imputed_vcf_template = os.path.join(
    workdir,
    "Results/imputation/chr{chrom}.omniChipOnly.imputed.10K.dose_chunk_{chunk}.vcf.gz"
)
genotyped_vcf_template = os.path.join(
    workdir,
    "Resources/Data/grf_vcf/VCF_MatchedOnly/GenotypeGVCFs_chr_{chrom}.haplotypecaller.chrname_fixed.biallelic.vcf.gz"
)
chunk_file_template = os.path.join(
    workdir,
    "Resources/Data/chunks/b38/20cM/chunks_chr{chrom}.header.incl.txt"
)
sample_file = os.path.join(
    workdir,
    "Resources/Data/grf_vcf/VCF_ChipOnly_MatchedOnly/samples.list"
)

outdir = os.path.join(workdir, "Results/imputation/corr_results")
shell_dir = os.path.join(outdir, "sh")
log_dir = os.path.join(outdir, "log")
tmp_dir = os.path.join(outdir, "tmp")

os.makedirs(outdir, exist_ok=True)
os.makedirs(shell_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)
os.makedirs(tmp_dir, exist_ok=True)

hostname_list = ["Client1", "Client2", "Client3", "Client4", "Client5"]
nthreads = 1
# ===========================================================


def get_qsub_script(cmd, jobname, stderr_path, stdout_path, nthreads=1, hostname=[], hold_jids=[]):
    script = "#!/bin/bash\n"
    script += f"#$ -N {jobname}\n"
    script += "#$ -S /bin/bash\n"
    script += "#$ -V\n"
    if hold_jids:
        script += f"#$ -hold_jid {','.join(hold_jids)}\n"
    if hostname:
        script += f"#$ -l hostname={'|'.join(hostname)}\n"
    script += f"#$ -pe smp {nthreads}\n"
    script += f"#$ -e {stderr_path}\n"
    script += f"#$ -o {stdout_path}\n\n"
    script += "echo \"=================================================\"\n"
    script += "echo \"Started on       : $(date)\"\n"
    script += "echo \"Running on node  : $(hostname)\"\n"
    script += "echo \"Current job ID   : $JOB_ID\"\n"
    script += "echo \"Current job name : $JOB_NAME\"\n"
    script += "echo \"=================================================\"\n\n"
    script += cmd + "\n\n"
    script += "echo \"=================================================\"\n"
    script += "echo \"Finished on      : $(date)\"\n"
    script += "echo \"=================================================\"\n"
    script += "echo \"Done.\"\n"
    return script


def submit_qsub_job(cmd, jobname):
    path_shell = os.path.join(shell_dir, f"{jobname}.sh")
    stderr_path = os.path.join(log_dir, f"{jobname}.stderr")
    stdout_path = os.path.join(log_dir, f"{jobname}.stdout")
    script_content = get_qsub_script(cmd, jobname, stderr_path, stdout_path, nthreads, hostname_list)

    with open(path_shell, "w") as fw:
        fw.write(script_content)
    result = subprocess.run(f"qsub {path_shell}", shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] qsub failed: {result.stderr}")
    else:
        print(f"[INFO] Job submitted: {result.stdout.strip()}")


# -------------------- Main Loop -------------------- #
for chrom in range(1, 23):  # Adjust chromosome range as needed
    chunk_file = chunk_file_template.format(chrom=chrom)
    if not os.path.exists(chunk_file):
        print(f"[WARNING] Chunk file not found: {chunk_file}")
        continue

    df = pd.read_csv(chunk_file, sep="\t")
    genotyped_vcf = genotyped_vcf_template.format(chrom=chrom)

    for _, row in df.iterrows():
        chunk_idx = int(row["chunk"])
        imputed_vcf = imputed_vcf_template.format(chrom=chrom, chunk=chunk_idx)
        if not os.path.exists(imputed_vcf):
            print(f"[WARNING] Missing imputed VCF for chunk {chunk_idx}: {imputed_vcf}")
            continue

        region_str = str(row["input_region"])
        chrom_str, coords = region_str.split(":")
        start, end = map(int, coords.split("-"))

        outfile = os.path.join(outdir, f"chr{chrom}.dose_chunk_{chunk_idx}.10K.corr.txt.gz")

        cmd = (
            f"python {corr_script} "
            f"-i {imputed_vcf} "
            f"-g {genotyped_vcf} "
            f"-s {sample_file} "
            f"-c {chrom} "
            f"-b {start} "
            f"-e {end} "
            f"-o {outfile}"
        )

        jobname = f"chr{chrom}_chunk{chunk_idx}_corr"
        submit_qsub_job(cmd, jobname)
