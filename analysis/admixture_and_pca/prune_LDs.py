# %%
import glob
import os
import subprocess

# %%
plink2 = "/BiO/Share/Tool/plink2"
window_size = "200kb"
r2 = "0.5"
# window_size = "500kb"
# r2 = "0.2"
# workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k"
os.makedirs(workdir, exist_ok=True)

# %%
def run_ld_pruning(plink2, inprefix, outprefix, window_size, r2, threads):
    cmd = f"{plink2} --bfile {inprefix} --indep-pairwise {window_size} {r2} --out {outprefix} --threads {threads}"
    print(cmd)
    subprocess.run(cmd, shell=True)
 
# %%
def run_extract_ld_pruned(plink2, inprefix, prune_in, outprefix, threads):
    cmd = f"{plink2} --bfile {inprefix} --extract {prune_in} --make-bed --out {outprefix} --threads {threads}"
    print(cmd)
    subprocess.run(cmd, shell=True)
 
# %%
# size = 500
# inprefix = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
# outprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.LD_pruning_{window_size}_{r2}")
# inprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_EAS.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
# outprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_EAS.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}")
# inprefix = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
# outprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}")
# inprefix = os.path.join(workdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
# outprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.preprocssed.prune_{window_size}_{r2}")

inprefix = os.path.join(workdir, "9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd")
outprefix = os.path.join(workdir, f"9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.prune.{window_size}_{r2}")

threads = 50
run_ld_pruning(plink2, inprefix, outprefix, window_size, r2, threads)

# %%
# inprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
# prune_in = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.LD_pruning_{window_size}_{r2}.prune.in")
# outprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.LD_pruning_{window_size}_{r2}.pruned")
# inprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_EAS.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
# prune_in = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_EAS.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}.prune.in")
# outprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_EAS.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}.pruned")
# inprefix = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
# prune_in = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}.prune.in")
# outprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.9000Koreans+1KGPEAS_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}")
# inprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
# prune_in = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.preprocssed.prune_{window_size}_{r2}.prune.in")
# outprefix = os.path.join(workdir, f"Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.excluded_outlier_plus_nonKorean_samples.postmerge_QC_filtered.preprocssed.prune_{window_size}_{r2}.pruned")


inprefix = os.path.join(workdir, "9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd")
prune_in = os.path.join(workdir, f"9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.prune.{window_size}_{r2}.prune.in")
outprefix = os.path.join(workdir, f"9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.prune.{window_size}_{r2}.pruned")

threads = 50
run_extract_ld_pruned(plink2, inprefix, prune_in, outprefix, threads)