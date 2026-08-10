# %%
import subprocess

bcftools = "/BiO/Research/Korea10KGenome/Resources/Tools/bcftools-1.20/bcftools"
path_input = "/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.RemoveABHetOutlier2STD.VQSR.PASS/Plink_Final/chromosome_merged.qc_filtered.kinship_filtered/Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.bcf"

cmd = f"{bcftools} index {path_input}"

subprocess.run(cmd, shell=True)