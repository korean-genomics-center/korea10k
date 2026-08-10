# %%
import os
import subprocess

from joblib import Parallel, delayed

tool_admixture = "/BiO/Share/Tool/admixture_linux-1.3.0/admixture"
path_input = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.LD_pruning_500kb_0.2.pruned.bed"

os.chdir(os.path.dirname(path_input))
list_k = list(range(4, 15, 1))

def run_cmd(tool_admixture, path_input, k):
    cmd = f"{tool_admixture} --cv {path_input} {k} | tee log{k}.out"
    subprocess.run(cmd, shell=True)
    print(f"running K={k}...")

for k in list_k:
    run_cmd(tool_admixture, path_input, k)