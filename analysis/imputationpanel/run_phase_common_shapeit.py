# %%
import os
import subprocess

# %%
chrnum = 22
phase_common = "/BiO/Share/Tool/shapeit5/static_bins/phase_common_static"
input_vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/chr{chrnum}.biallelic.vcf.gz"
maf = 0.001
gmap = f"/BiO/Access/kyungwhan1998/genome/shapeit/maps/b38/chr{chrnum}.b38.gmap.gz"
output_vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/chr{chrnum}.biallelic.scaffold.bcf"
output_format = "bcf"
num_threads = 40

# %%
def run_cmd(phase_common, input_vcf, maf, region, gmap, output_vcf, num_threads):
    cmd = f"{phase_common} \
    --input {input_vcf} \
    --filter-maf {maf} \
    --region {region} \
    --map {gmap} \
    --output {output_vcf} \
    --output-format {output_format} \
    --thread {num_threads}"
    
    subprocess.run(cmd, shell=True)

# %%
run_cmd(phase_common, input_vcf, maf, f"chr{chrnum}", gmap, output_vcf, num_threads)

