# %%
import os
import subprocess
from korea10k.config import STORE_DIR, TOOL_DIR, WORK_DIR

# %%
# chrnum = 22
# phase_common = f"{STORE_DIR}/RNASeq/Tools/shapeit5/phase_common/bin/phase_common"
# input_vcf = f"{WORK_DIR}/genome/shapeit/chr{chrnum}.asterisk_removed.vcf.gz"
# maf = 0.001
# gmap = f"{WORK_DIR}/genome/shapeit/maps/b38/chr{chrnum}.b38.gmap.gz"
# output_vcf = f"{WORK_DIR}/genome/shapeit/chr{chrnum}.asterisk_removed.chunk1.vcf.gz"
# num_threads = 20

# %%
phase_common = f"{STORE_DIR}/RNASeq/Tools/shapeit5/phase_common/bin/phase_common"
workdir = f"{TOOL_DIR}/shapeit5/test"
input_vcf = os.path.join(workdir,"array/target.unrelated.bcf")
region = 1
gmap = os.path.join(workdir, "info/chr1.gmap.gz")
output_vcf = os.path.join(workdir, "tmp/target.phased.bcf")
num_threads = 8               

# %%
def run_cmd(phase_common, input_vcf, region, gmap, output_vcf, num_threads):
    cmd = f"{phase_common} \
    --input {input_vcf} \
    --region {region} \
    --map {gmap} \
    --output {output_vcf} \
    --thread {num_threads}"
    
    subprocess.run(cmd, shell=True)

# %%
run_cmd(phase_common, input_vcf, region, gmap, output_vcf, num_threads)
# %%
