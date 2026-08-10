# %%
import os
import subprocess
from korea10k.config import PROJECT_DIR, TOOL_DIR

# %%
workdir = f"{PROJECT_DIR}/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.RemoveABHetOutlier2STD.VQSR.PASS/Plink_Final/chromosome_merged.qc_filtered.kinship_filtered"
plink2 = f"{TOOL_DIR}/plink2"
pfile_in = f"{workdir}/Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd"
af_out = f"{workdir}/Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd_freq"
threads = 50

cmd = f"{plink2} --pfile {pfile_in} --freq --out {af_out} --threads {threads}"
subprocess.run(cmd, shell=True)