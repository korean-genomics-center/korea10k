import os
import subprocess

gatk = "/BiO/Research/Korea10KGenome/Resources/Tools/gatk/gatk-4.6.1.0/gatk"

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


def run(gatk, 
        input_vcf, 
        output_vcf, 
        reject_vcf, 
        reference_fa, 
        reference_chain, 
        dir_tmp, dir_shell, 
        dir_log, 
        jobname, 
        n_threads, 
        hostname=[], 
        jobs_wait=[]):


    cmd = f'{gatk} LiftoverVcf  \
            --java-options " -Xmx10g" \
            --CHAIN {reference_chain} \
            --INPUT {input_vcf}  \
            --OUTPUT {output_vcf} \
            --REFERENCE_SEQUENCE {reference_fa} \
            --REJECT {reject_vcf} \
            --VERBOSITY DEBUG'

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
    # reference_fa = "/BiO/Share/References/UCSC/hg19.fa"
    # reference_chain = "/BiO/Share/References/UCSC/hg38ToHg19.over.chain.gz"
    # input_dir = "/BiO/Research/KOREF_PersonalMultiomicsReference/Results/KPGP_WGS"
    # output_dir = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/kpgp_vcf"
    reference_fa = "/BiO/Share/References/Broad/hg38.fa"
    reference_chain = "/BiO/Share/References/UCSC/hg19ToHg38.over.chain.gz"
    input_dir = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/kpgp_vcf"
    output_dir = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/kpgp_vcf"   
    workdir = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/kpgp_vcf/run"
    os.makedirs(workdir, exist_ok=True)
    list_file_vcfs = sorted(list(filter(lambda x: str(x).endswith(".bgz"), os.listdir(input_dir))))
    list_file_vcfs_filt = list(filter(lambda x: "Infinium" in x, list_file_vcfs))
    for file_vcf in list_file_vcfs_filt:
        basename_vcf = file_vcf.split(".")[0]
        input_vcf = os.path.join(input_dir, file_vcf)
        output_vcf = os.path.join(output_dir, file_vcf.replace(".vcf", ".hg38.liftover.vcf"))
        reject_vcf = os.path.join(output_dir, file_vcf.replace(".vcf", ".hg38.rejected.vcf"))
        jobname = f"{file_vcf}-GATK_Liftover"
        print(f"[INFO] Preparing job: {jobname}")
        
        run(
            gatk, 
            input_vcf, 
            output_vcf, 
            reject_vcf, 
            reference_fa, 
            reference_chain, 
            os.path.join(workdir, "tmp"),
            os.path.join(workdir, "sh"),
            os.path.join(workdir, "log"),
            jobname,
            n_threads=20,
            hostname=["Client1", "Client2", "Client3", "Client4"],
            jobs_wait=[]
        )