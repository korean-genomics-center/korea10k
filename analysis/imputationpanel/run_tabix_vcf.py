# %%
import os
import subprocess

# %%
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


def run(tabix,
        output_vcf,
        dir_tmp, dir_shell, 
        dir_log, 
        jobname, 
        n_threads, 
        hostname=[], 
        jobs_wait=[]):

    cmd = f"{tabix} -p vcf {output_vcf}"

    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)
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

if __name__ == "__main__":
    dir_vcf = "/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation"
    workdir = "/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation/run"
    os.makedirs(workdir, exist_ok=True)
    list_file_vcfs = sorted(list(filter(lambda x: str(x).endswith("10K.dose.vcf.gz"), os.listdir(dir_vcf))))
    for file_vcf in list_file_vcfs:
        chrnum = file_vcf.split(".")[0]
        output_vcf = os.path.join(dir_vcf, f"{chrnum}.omniChipOnly.imputed.10K.dose.vcf.gz")
        jobname = f"{file_vcf}-TABIX-VCF"
        print(f"[INFO] Preparing job: {jobname}")
        
        run(
            tabix,
            output_vcf,
            os.path.join(workdir, "tmp"),
            os.path.join(workdir, "sh"),
            os.path.join(workdir, "log"),
            jobname,
            n_threads=1,
            hostname=["Client1", "Client2", "Client3", "Client4", "Client5"],
            jobs_wait=[]
        )


