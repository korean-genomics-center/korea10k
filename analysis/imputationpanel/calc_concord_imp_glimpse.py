# %%
import os
import subprocess

# %%
path_glimpse = "/BiO/Share/Tool/GLIMPSE2_concordance_static"
chr_num = "chr9"
Korea10K_vcf = f"/BiO/Access/kyungwhan1998/genome/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.RemoveABHetOutlier2STD.VQSR.PASS/Plink_Final/chromosome_merged.qc_filtered.kinship_filtered/Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.bcf"
original_vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/kpgp_vcf/{chr_num}.recal.selected.reheadered.vcf.bgz"
imputed_vcf = f"/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation/{chr_num}.recal.selected.reheadered.hg19.liftoverInfiniumOmni2-5-8v1-5_A1_chip_sites.hg38.liftover.imputed.10K.vcf.bgz"
concordance_input = f"/BiO/Access/kyungwhan1998/genome/shapeit/Results/concordance/{chr_num}.10K.input.txt"
os.makedirs(os.path.dirname(concordance_input), exist_ok=True)
concordance_output = f"/BiO/Access/kyungwhan1998/genome/shapeit/Results/concordance/{chr_num}.10K.output.txt"
cmd1 = f'echo "{chr_num} {Korea10K_vcf} {original_vcf} {imputed_vcf}" > {concordance_input}'
cmd2 = f"{path_glimpse} --gt-val \
        --ac-bins 1 5 10 20 30\
        --threads 10 \
        --input  {concordance_input} \
        --output {concordance_output}"

cmd = f"{cmd1} && {cmd2}"
subprocess.run(cmd, shell=True)