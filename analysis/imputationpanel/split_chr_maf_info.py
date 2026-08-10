# %%
import os

from joblib import Parallel, delayed

# %%
path_freq = "/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.RemoveABHetOutlier2STD.VQSR.PASS/Plink_Final/chromosome_merged.qc_filtered.kinship_filtered/Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd_freq.afreq"
outdir = "/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation/maf_info"
os.makedirs(outdir, exist_ok=True)
outfreq = os.path.join(outdir, "chr{i}_maf_10K.txt")
    
def split_chr_maf_info(outfreq, i):
    with open(outfreq.format(i=i), mode="w") as fw:
        with open(path_freq, mode="r") as fr:
            header = fr.readline()
            fw.write(header)
            for line in fr:
                record = line.rstrip("\n").split("\t")
                chrnum = record[0]
                if str(chrnum) == str(i):
                    fw.write(line)

with Parallel(n_jobs=23) as parallel:
    parallel(delayed(split_chr_maf_info)(outfreq, i) for i in range(1, 23, 1))
