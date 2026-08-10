# %%
import glob
import os
import subprocess

from joblib import Parallel, delayed

# %%
plink2 = "/BiO/Share/Tool/plink2"
fasta = "/BiO/Research/Korea10KGenome/Resources/Reference/chromosome/hg38.fa"
outfmt = "vcf-4.2"
vcf_dosage = "HDS-force"
dir_in = "/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.RemoveABHetOutlier2STD.VQSR.PASS/Plink_Final/chromosome_merged.qc_filtered.kinship_filtered"
list_pgen = list(map(lambda x: ".".join(x.split(".")[:-1]), glob.glob(f"{dir_in}/*.pgen")))
list_vcf = list(map(lambda x: x.replace("pgen", outfmt.split("-")[0]) + f".{vcf_dosage}", list_pgen))
os.makedirs(os.path.dirname(list_vcf[0]), exist_ok=True)

# %%
def run_cmd(plink2, file_pgen, file_vcf, outfmt, vcf_dosage):
    cmd = f"{plink2} --pfile {file_pgen} --export {outfmt} vcf-dosage={vcf_dosage} --out {file_vcf}"
    subprocess.run(cmd, shell=True)
    print("job done!")
# %%
with Parallel(n_jobs=len(list_vcf)) as parallel:
    list_res = parallel(delayed(run_cmd)(plink2, file_pgen, file_vcf, outfmt, vcf_dosage) for file_pgen, file_vcf in zip(list_pgen, list_vcf))

# %%
for res in list_res:
    print(res)