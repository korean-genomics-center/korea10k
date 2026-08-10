# %%
import glob
import os
import subprocess

plink2 = "/BiO/Share/Tool/plink2"
# indir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
# bfile_in = os.path.join(indir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
# bfile_out = os.path.join(indir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
# bfile_in = f"{indir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
# bfile_out = f"{indir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
# size = 500
# bfile_in = os.path.join(indir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_1KGP.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
# bfile_out = os.path.join(indir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_1KGP.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")

indir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k"
bfile_in = os.path.join(indir, "9000Korean_samples.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd")
bfile_out = os.path.join(indir, f"9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd")


cmd = f"{plink2} --bfile {bfile_in} --geno 0.01 --maf 0.01 --hwe 1e-6 --make-bed --out {bfile_out}"

subprocess.run(cmd, shell=True)