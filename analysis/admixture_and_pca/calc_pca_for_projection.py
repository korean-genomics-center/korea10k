# %%
import subprocess

plink = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
size = 500
bfile_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_EAS.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
window_size = "200kb"
r2 = "0.5"
prune_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_EAS.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}.prune.in"
pca_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_EAS.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.{window_size}_{r2}.pca"
threads = 50
def run_pca_ref(plink, bfile_in, prune_in, pca_out, threads):
        cmd = f"{plink} --bfile {bfile_in} \
        --freq counts \
        --extract {prune_in} \
        --pca allele-wts vcols=chrom,ref,alt \
        --out {pca_out} \
        --threads {threads}"

        subprocess.run(cmd, shell=True)

# run_pca_ref(plink, bfile_in, prune_in, pca_out, threads)
# %%
plink = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
bfile_in  = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.include_eas_samples_only.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.LD_pruning_200kb_0.5.pruned"
pca_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.{size}_random_Koreans_plus_EAS.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.prune.200kb_0.5.pca"
projection_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.include_eas_samples_only.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.LD_pruning_200kb_0.5.pruned.projection_{size}"
threads = 50
size = 500

def project_pca_test(plink, bfile_in, pca_in, projection_out, threads):
        cmd = f"{plink} --bfile {bfile_in} \
        --read-freq {pca_in}.acount\
        --score {pca_in}.eigenvec.allele 2 5 header-read no-mean-imputation \
                variance-standardize \
        --score-col-nums 6-15 \
        --out {projection_out} \
        --threads {threads}"
        
        subprocess.run(cmd, shell=True)

# project_pca_test(plink, bfile_in, pca_in, projection_out, threads)


# %%
import subprocess

plink = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
size = 500
bfile_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd"
window_size = "200kb"
r2 = "0.5"
prune_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.preprocssed.prune_{window_size}_{r2}.prune.in"
pca_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.preprocssed.prune_{window_size}_{r2}.pca"
threads = 50
def run_pca_ref(plink, bfile_in, prune_in, pca_out, threads):
        cmd = f"{plink} --bfile {bfile_in} \
        --freq counts \
        --extract {prune_in} \
        --pca allele-wts vcols=chrom,ref,alt \
        --out {pca_out} \
        --threads {threads}"

        subprocess.run(cmd, shell=True)

# run_pca_ref(plink, bfile_in, prune_in, pca_out, threads)
# %%
plink = "/BiO/Share/Tool/plink2"
workdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
bfile_in  = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.include_eas_samples_only.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd.LD_pruning_{window_size}_{r2}.pruned"
pca_in = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.exclude_outlier_samples.Koreans_only.postmerge_QC_filtered.preprocssed.prune_{window_size}_{r2}.pca"
projection_out = f"{workdir}/Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.Koreans_only.preprocessed.prune_200kb_0.5.projection_include_eas_samples_only"
threads = 50

def project_pca_test(plink, bfile_in, pca_in, projection_out, threads):
        cmd = f"{plink} --bfile {bfile_in} \
        --read-freq {pca_in}.acount\
        --score {pca_in}.eigenvec.allele 2 5 header-read no-mean-imputation \
                variance-standardize \
        --score-col-nums 6-15 \
        --out {projection_out} \
        --threads {threads}"
        
        subprocess.run(cmd, shell=True)

project_pca_test(plink, bfile_in, pca_in, projection_out, threads)