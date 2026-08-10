# %%
import glob
import os
import subprocess

from joblib import Parallel, delayed

# %%
plink2 = "/BiO/Share/Tool/plink2"
# indir = "/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.RemoveABHetOutlier2STD.VQSR.PASS/Plink_Final/chromosome_merged.qc_filtered.kinship_filtered"
# outdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/ku10k"
indir = "/BiO/Research/Korea10KGenome/Results/External_Genome_Data/1KGP/1KGP_30x_GRCh38/Plink_All/Step6_Kinship_Filtering"
outdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/1kgp"
list_pgen = list(map(lambda x: ".".join(x.split(".")[:-1]), glob.glob(f"{indir}/*.pgen")))
list_outprefix = list(map(lambda x: os.path.join(outdir, os.path.basename(x)), list_pgen))
os.makedirs(outdir, exist_ok=True)

# %%
def run_cmd(plink2, inprefix, outprefix):
    cmd = f"{plink2} --pfile {inprefix} --make-bed --out {outprefix}"
    subprocess.run(cmd, shell=True)
    print("job done!")
    
# %%
with Parallel(n_jobs=len(list_outprefix)) as parallel:
    list_res = parallel(delayed(run_cmd)(plink2, inprefix, outprefix) for inprefix, outprefix in zip(list_pgen, list_outprefix))

# %%
for res in list_res:
    print(res)