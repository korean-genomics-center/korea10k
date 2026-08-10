# %%
import subprocess

plink = "/BiO/Share/Tool/plink2"
# workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k"

# bfile_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
# bfile_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"

bfile_in = f"{workdir}/9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd"
window_size = "200kb"
r2 = "0.5"

# prune_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}.prune.in"
# prune_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}.prune.in"
prune_in = f"{workdir}/9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.prune.{window_size}_{r2}.prune.in"
pca_cnt = "20"

# pca_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}.pca"
# pca_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}.pca"
pca_out = f"{workdir}/9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.prune.{window_size}_{r2}.pca"
threads = 100

cmd = f"{plink} --bfile {bfile_in} --extract {prune_in} --pca {pca_cnt} --out {pca_out} --threads {threads}"

subprocess.run(cmd, shell=True)
