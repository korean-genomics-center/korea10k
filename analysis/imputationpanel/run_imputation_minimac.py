import os
import subprocess

minimac3 = "/BiO/Share/Tool/Minimac3Executable/bin/Minimac3"
tabix = "/BiO/Access/kyungwhan1998/miniconda3/envs/shapeit/bin/tabix"

# ------------ Utility functions ------------

def get_context_for_qsub_submission(cmd, jobname, stderr_path, stdout_path, nthreads, hostname=[], hold_jids=[]):
    script_header = "#!/bin/bash\n"
    script_header += f"#$ -N {jobname}\n"
    script_header += f"#$ -S /bin/bash\n"
    script_header += f"#$ -V\n"
    if hold_jids:
        script_header += f"#$ -hold_jid {','.join(hold_jids)}\n"
    if hostname:
        script_header += f"#$ -l hostname={'|'.join(hostname)}\n"
    script_header += f"#$ -e {stderr_path}\n"
    script_header += f"#$ -o {stdout_path}\n"
    script_header += f"#$ -pe smp {nthreads}\n\n"

    script_header += "echo \"=================================================\"\n"
    script_header += "echo \"Started on       : $(date)\"\n"
    script_header += "echo \"Running on node  : $(hostname)\"\n"
    script_header += "echo \"Current job ID   : $JOB_ID\"\n"
    script_header += "echo \"Current job name : $JOB_NAME\"\n"
    script_header += "echo \"=================================================\"\n\n"

    script_footer = "\n\n"
    script_footer += "echo \"=================================================\"\n"
    script_footer += "echo \"Finished on      : $(date)\"\n"
    script_footer += "echo \"=================================================\"\n"
    script_footer += "echo \"Done.\"\n"

    return script_header + cmd + script_footer


def write_shell_script_and_run(path_shell, shell_context):
    print(f"[INFO] Writing shell script to: {path_shell}")
    with open(path_shell, 'w') as fw:
        fw.write(shell_context)
    print(f"[INFO] Shell script written. Submitting job...")
    result = subprocess.run(f"qsub {path_shell}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"[ERROR] qsub failed: {result.stderr}")
    else:
        print(f"[INFO] qsub submitted: {result.stdout.strip()}")


def run(minimac3,
        tabix,
        reference_panel_m3vcf, 
        target_vcf, 
        imputed_vcf,
        dir_tmp, 
        dir_shell, 
        dir_log, 
        jobname, 
        n_threads, 
        hostname=[], 
        jobs_wait=[]):


    cmd = f"{minimac3} --refHaps {reference_panel_m3vcf} \
            --haps {target_vcf} \
            --prefix {imputed_vcf} \
            --cpus {n_threads} \
            {tabix} -p vcf {imputed_vcf}"

    os.makedirs(os.path.dirname(imputed_vcf), exist_ok=True)
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)

    path_shell = os.path.join(dir_shell, f"{jobname}.sh")
    stderr_path = os.path.join(dir_log, f"{jobname}.stderr")
    stdout_path = os.path.join(dir_log, f"{jobname}.stdout")

    shell_script = get_context_for_qsub_submission(
        cmd, jobname, stderr_path, stdout_path, n_threads, hostname, jobs_wait
    )

    write_shell_script_and_run(path_shell, shell_script)

# ------------ Main ------------
workdir = "/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation/run"
if __name__ == "__main__":
    for i in range(1, 23, 1):
        # reference_panel_m3vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/reference_panel/1K/chr{i}.imputed.chrname_fixed.m3vcf.gz"
        # reference_panel_m3vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/reference_panel/4K/chr{i}.recal.forPhasing.nonOverlap.phased.corrected.headerFixed.chrname_fixed.m3vcf.gz"
        reference_panel_m3vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/reference_panel/10K/chr{i}.phased.concat.m3vcf.gz"
        
        target_vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/grf_vcf/VCF_ChipOnly_MatchedOnly/chr{i}.omniChipOnly.biallelic.vcf.gz"

        # imputed_vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation/chr{i}.omniChipOnly.imputed.1K"
        # imputed_vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation/chr{i}.omniChipOnly.imputed.4K"
        imputed_vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation/chr{i}.omniChipOnly.imputed.10K"
        
        jobname = f"{os.path.basename(imputed_vcf)}-minimac3-IMPUTE"
        print(f"[INFO] Preparing job: {jobname}")
        
        run(minimac3,
            tabix,
            reference_panel_m3vcf, 
            target_vcf, 
            imputed_vcf,
            os.path.join(workdir, "tmp"),
            os.path.join(workdir, "sh"),
            os.path.join(workdir, "log"),
            jobname,
            n_threads=5,
            hostname=["Client1", "Client2", "Client3", "Client4", "Client5"],
            jobs_wait=[]
        )